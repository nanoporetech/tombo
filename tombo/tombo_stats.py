import sys, os

import re
import h5py
import Queue
import warnings
import pkg_resources

import numpy as np
np.seterr(all='raise')
import multiprocessing as mp

from time import sleep
from scipy import stats
from collections import defaultdict
from itertools import repeat, product
from scipy.spatial.distance import pdist
from scipy.cluster.hierarchy import single, leaves_list

import tombo_helper as th

from c_helper import c_mean_std, c_calc_llh_ratio

VERBOSE = False

PROFILE_SIGNIF = False
PROFILE_EST_KMER = False

DEBUG_EST_ALT = True
DEBUG_EST_BW = 0.05
DEBUG_EST_NUM_KMER_SAVE = 500

DNA_BASES = ['A','C','G','T']
MIN_POSITION_SD = 0.01

#######################################################
##### Pair-wise distance functions for clustering #####
#######################################################

def order_reads(log_r_pvals):
    """
    Compute order of reads based on log p-values
    """
    # get pairwise distances between reads
    # will get some empty slice means warnings, so suppress those
    #   (no seterr for this specific warning)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        r_dists = pdist(log_r_pvals, lambda u, v:
                        np.nanmean(np.sqrt(((u-v)**2))))
    r_dists[np.isnan(r_dists)] = np.nanmax(r_dists) + 1
    # then perform single/min linkage clustering and return the leaf order
    return leaves_list(single(r_dists))

def sliding_window_dist(sig_diffs1, sig_diffs2, slide_span, num_bases):
    """
    Compute distance over the minimum over a sliding window
    """
    return np.sqrt(min(np.sum(np.square(
        sig_diffs1[i1:i1+num_bases] - sig_diffs2[i2:i2+num_bases]))
                       for i1 in range((slide_span * 2) + 1)
                       for i2 in range((slide_span * 2) + 1)))

def euclidian_dist(sig_diffs1, sig_diffs2):
    """
    Compute Euclidean distance
    """
    return np.sqrt(np.sum(np.square(sig_diffs1 - sig_diffs2)))

def get_pairwise_dists(reg_sig_diffs, index_q, dists_q, slide_span=None):
    """
    Compute pairwise distances between a set of signal shifts
    """
    if slide_span > 0:
        num_bases=reg_sig_diffs[0].shape[0] - (slide_span * 2)
    while not index_q.empty():
        try:
            index = index_q.get(block=False)
        except Queue.Empty:
            break

        if slide_span > 0:
            row_dists = np.array(
                [sliding_window_dist(
                    reg_sig_diffs[index], reg_sig_diffs[j],
                    slide_span, num_bases)
                 for j in range(0,index+1)] +
                [0 for _ in range(index+1, len(reg_sig_diffs))])
        else:
            row_dists = np.array(
                [euclidian_dist(reg_sig_diffs[index], reg_sig_diffs[j])
                 for j in range(0,index+1)] +
                [0 for _ in range(index+1, len(reg_sig_diffs))])
        dists_q.put((index, row_dists))

    return


#################################################
##### Significant region selection function #####
#################################################

def get_most_signif_regions(
        all_stats, num_bases, num_regions, qval_thresh=None, unique_pos=True,
        fraction_order=False):
    """
    Select the most significant genome locations based on some criteria
    """
    # applied threshold for scores on each chromosome, so now
    # we include all here
    if qval_thresh is not None:
        num_regions = np.argmax(np.logical_and(
            np.logical_not(np.isnan(all_stats['mt_stat'])),
            all_stats['mt_stat'] > qval_thresh))
        if num_regions == 0:
            sys.stderr.write(
                '*' * 60 + '\nERROR: No regions identified q-value ' +
                'below thresh. Minumum q-value: {:.2g}\n'.format(
                    all_stats[0]['mt_stat']) + '*' * 60 + '\n')
            sys.exit()

    if fraction_order:
        all_stats.sort(order='frac')
    else:
        all_stats.sort(order='stat')
    plot_intervals = []
    used_intervals = defaultdict(set)
    for i, stat in enumerate(all_stats):
        int_start = max(0, stat['pos'] - int(num_bases / 2.0))
        interval_poss = range(int_start, int_start + num_bases)
        if not unique_pos or \
           stat['pos'] not in used_intervals[(stat['chrm'], stat['strand'])]:
            used_intervals[(stat['chrm'], stat['strand'])].update(interval_poss)
            plot_intervals.append(th.intervalData(
                '{:03d}'.format(i), stat['chrm'], int_start,
                int_start + num_bases, stat['strand'],
                'Frac Standard:{1:.2g} Frac. Alternate:{2:.2g}'.format(
                    stat['stat'], stat['frac'], stat['alt_frac'],
                    stat['cov'])))
            if len(plot_intervals) >= num_regions: break

    return plot_intervals


#######################################
##### K-mer model estimation code #####
#######################################

def write_tombo_model(kmer_ref, kmer_ref_fn, central_pos,
                      alt_base=None, alt_name=None):
    """
    Write a tombo model file
    """
    with h5py.File(kmer_ref_fn, 'w') as kmer_ref_fp:
        kmer_ref_fp.create_dataset('model', data=kmer_ref, compression="gzip")
        kmer_ref_fp.attrs['central_pos'] = central_pos
        if alt_base is None:
            kmer_ref_fp.attrs['model_name'] = 'standard'
        else:
            kmer_ref_fp.attrs['model_name'] = alt_name
            kmer_ref_fp.attrs['alt_base'] = alt_base

    return

def parse_tombo_model(kmer_ref_fn):
    """
    Parse a tombo model file
    """
    try:
        with h5py.File(kmer_ref_fn, 'r') as kmer_ref_fp:
            kmer_ref_raw = kmer_ref_fp['model'].value
            central_pos = kmer_ref_fp.attrs['central_pos']
            model_name = kmer_ref_fp.attrs['model_name']
            try:
                alt_base = kmer_ref_fp.attrs['alt_base']
            except:
                alt_base = None
    except:
        sys.stderr.write(
            '********* ERROR *********\n\tInvalid tombo kmer model ' +
            'file provided: ' + str(kmer_ref_fn) + '\n')
        sys.exit()

    kmer_ref = dict((kmer, (kmer_mean, kmer_std))
                    for kmer, kmer_mean, kmer_std in kmer_ref_raw)

    return kmer_ref, central_pos, alt_base, model_name

def parse_tombo_models(alt_fns, ref_upstrm_bases, ref_kmer_width):
    """
    Parse several alternative tombo model files
    """
    alt_refs = []
    for alt_model_fn in alt_fns:
        alt_ref, alt_upstrm_bases, alt_base, alt_name = parse_tombo_model(alt_model_fn)
        if (ref_upstrm_bases != alt_upstrm_bases or
            ref_kmer_width != len(next(alt_ref.iterkeys()))):
            sys.stderr.write(
                '********* WARNING *********\n\tStandard and ' + alt_model_fn +
                ' alternative base models must be estimated using the same ' +
                'k-mer positions.\n')
            continue
        if alt_base is None:
            sys.stderr.write(
                '********* WARNING *********\n\tAlternative model ' + alt_model_fn +
                ' appears to be a standard model and will not be processed.\n')
            continue

        alt_refs.append((alt_name, alt_ref, alt_base))

    return alt_refs

def get_default_standard_ref(raw_read_coverage):
    if th.is_rna(raw_read_coverage):
        standard_ref_fn = th.STANDARD_MODELS['RNA']
    else:
        standard_ref_fn = th.STANDARD_MODELS['DNA']
    # get full filename path with setuptools
    standard_ref_fn = pkg_resources.resource_filename(
        'tombo', 'tombo_models/' + standard_ref_fn)

    return standard_ref_fn

def get_default_standard_ref_from_files(fast5_fns):
    if th.is_rna_from_files(fast5_fns):
        standard_ref_fn = th.STANDARD_MODELS['RNA']
    else:
        standard_ref_fn = th.STANDARD_MODELS['DNA']
    # get full filename path with setuptools
    standard_ref_fn = pkg_resources.resource_filename(
        'tombo', 'tombo_models/' + standard_ref_fn)

    return standard_ref_fn

def get_default_alt_ref(alt_name, raw_read_coverage):
    if th.is_rna(raw_read_coverage):
        samp_type = 'RNA'
        try:
            alt_model_fn = th.ALTERNATE_MODELS['RNA_' + alt_name]
        except KeyError:
            alt_model_fn = None
    else:
        samp_type = 'DNA'
        try:
            alt_model_fn = th.ALTERNATE_MODELS['DNA_' + alt_name]
        except KeyError:
            alt_model_fn = None
    if alt_model_fn is not None:
        # get full filename path with setuptools
        alt_model_fn = pkg_resources.resource_filename(
            'tombo', 'tombo_models/' + alt_model_fn)
    if alt_model_fn is None or not os.path.isfile(alt_model_fn):
        sys.stderr.write(
            '******** WARNING *********\n\tTombo default model for ' +
            alt_name + ' in ' + samp_type + ' does not exists.\n')
        return

    return alt_model_fn

def load_alt_refs(alt_names, raw_read_coverage,
                  ref_upstrm_bases, ref_kmer_width):
    """
    Load several default alternative tombo models
    """
    alt_fns = []
    for alt_name in alt_names:
        alt_model_fn = get_default_alt_ref(alt_name, raw_read_coverage)
        if alt_model_fn is None:
            continue
        alt_fns.append(alt_model_fn)

    return parse_tombo_models(alt_fns, ref_upstrm_bases, ref_kmer_width)

def calc_med_sd(vals):
    """
    Helper function to compute median and standard deviation from a numpy array
    """
    return np.median(vals), np.std(vals)

def get_region_kmer_levels(
        reg_reads, cov_thresh, upstrm_bases, dnstrm_bases,
        cs_cov_thresh, est_mean, region_size, reg_start, strand):
    """
    Compute mean or median and standard deviation for each k-mer
    """
    if cs_cov_thresh is not None:
        # sample reads until requested mean depth of coverage is achieved
        cs_num_bases_thresh = region_size * cs_cov_thresh
        np.random.shuffle(reg_reads)
        cumm_num_bases = np.cumsum([r_data.end - r_data.start
                                    for r_data in reg_reads])
        try:
            cs_num_reads = next((i for i, v in enumerate(cumm_num_bases)
                                 if v >= cs_num_bases_thresh))
            reg_reads = reg_reads[:cs_num_reads]
        except StopIteration:
            # if threshold is not met use all reads from region
            pass
    base_events = th.get_reads_events(reg_reads, strand == '-')
    if len(base_events) == 0:
        return

    # get intervals where coverage is high enough for model estimation
    reg_cov = np.array([len(base_events[pos]) if pos in base_events else 0
                        for pos in xrange(reg_start, reg_start + region_size)])
    cov_intervals = np.where(np.diff(np.concatenate(
        [[False], reg_cov > cov_thresh])))[0]
    if reg_cov[-1] > cov_thresh:
        # if valid coverage goes until end of coverage, add the last
        # interval endpoint
        cov_intervals = np.concatenate([cov_intervals, [region_size]])
    if cov_intervals.shape[0] <= 1:
        return
    cov_intervals = cov_intervals.reshape(-1,2)

    kmer_width = upstrm_bases + dnstrm_bases + 1
    reg_kmer_levels = dict(
        (''.join(kmer),[])
        for kmer in product(DNA_BASES, repeat=kmer_width))
    # upstream and downstream changes the sequence selection
    # depending on the strand
    bb, ab = (upstrm_bases, dnstrm_bases) if strand == '+' else \
             (dnstrm_bases, upstrm_bases)
    for int_start, int_end in cov_intervals:
        int_seq = th.get_seq_from_reads(
            reg_start + int_start, reg_start + int_end, reg_reads)
        if strand == '-':
            int_seq = th.rev_comp(int_seq)
        int_len = int_end - int_start
        for pos in range(int_len - kmer_width + 1):
            pos_kmer = int_seq[pos:pos+kmer_width] \
                       if strand == '+' else \
                          int_seq[int_len-kmer_width-pos:int_len-pos]
            try:
                if est_mean:
                    reg_kmer_levels[pos_kmer].append(c_mean_std(
                        base_events[reg_start+pos+int_start+bb]))
                else:
                    reg_kmer_levels[pos_kmer].append(calc_med_sd(
                        base_events[reg_start+pos+int_start+bb]))
            except KeyError:
                continue

    return reg_kmer_levels

def _est_kmer_model_worker(
        region_q, kmer_level_q, raw_read_coverage, cov_thresh,
        upstrm_bases, dnstrm_bases, cs_cov_thresh, est_mean, region_size):
    while not region_q.empty():
        try:
            chrm, strand, reg_start = region_q.get(block=False)
        except Queue.Empty:
            break

        reg_reads = [r_data for r_data in raw_read_coverage[(chrm, strand)]
                     if not (r_data.start >= reg_start + region_size or
                             r_data.end <= reg_start)]
        if len(reg_reads) == 0:
            if VERBOSE:
                sys.stderr.write('.')
                sys.stderr.flush()
            continue

        reg_kmer_levels = get_region_kmer_levels(
            reg_reads, cov_thresh, upstrm_bases, dnstrm_bases,
            cs_cov_thresh, est_mean, region_size, reg_start, strand)
        if reg_kmer_levels is not None:
            kmer_level_q.put(reg_kmer_levels)
        if VERBOSE:
            sys.stderr.write('.')
            sys.stderr.flush()

    return

if PROFILE_EST_KMER:
    _est_kmer_model_wrapper = _est_kmer_model_worker
    def _est_kmer_model_worker(*args):
        import cProfile
        cProfile.runctx('_est_kmer_model_wrapper(*args)', globals(), locals(),
                        filename='est_kmer_model.prof')
        return

def estimate_kmer_model(
        f5_dirs, corrected_group, basecall_subgroups,
        kmer_ref_fn, cov_thresh, upstrm_bases, dnstrm_bases, min_kmer_obs,
        kmer_specific_sd, cs_cov_thresh, est_mean, region_size, num_processes):
    """
    Estimate a standard tombo k-mer model
    """
    raw_read_coverage = th.parse_fast5s(
        f5_dirs, corrected_group, basecall_subgroups)
    chrm_sizes = th.get_chrm_sizes(raw_read_coverage)

    manager = mp.Manager()
    region_q = manager.Queue()
    kmer_level_q = manager.Queue()
    num_regions = 0
    for chrm, chrm_len in chrm_sizes.items():
        plus_covered = (chrm, '+') in raw_read_coverage
        minus_covered = (chrm, '-') in raw_read_coverage
        for reg_start in range(0, chrm_len, region_size):
            if plus_covered:
                region_q.put((chrm, '+', reg_start))
                num_regions +=1
            if minus_covered:
                region_q.put((chrm, '-', reg_start))
                num_regions +=1

    if VERBOSE: sys.stderr.write(
            'Extracting average kmer levels across ' + str(num_regions) +
            ' regions. (Will print a dot or each batch completed)\n')
    est_args = (
        region_q, kmer_level_q, raw_read_coverage, cov_thresh,
        upstrm_bases, dnstrm_bases,
        cs_cov_thresh, est_mean, region_size)
    est_ps = []
    for p_id in xrange(num_processes):
        p = mp.Process(target=_est_kmer_model_worker, args=est_args)
        p.start()
        est_ps.append(p)

    all_reg_kmer_levels = []
    while any(p.is_alive() for p in est_ps):
        try:
            reg_kmer_levels = kmer_level_q.get(block=False)
            all_reg_kmer_levels.append(reg_kmer_levels)
        except Queue.Empty:
            sleep(1)
            continue
    while not kmer_level_q.empty():
        reg_kmer_levels = kmer_level_q.get(block=False)
        all_reg_kmer_levels.append(reg_kmer_levels)
    if VERBOSE: sys.stderr.write('\n')

    if len(all_reg_kmer_levels) == 0:
        sys.stderr.write('********** ERROR *********\n\tNo genomic positions ' +
                         'contain --minimum-test-reads.\n')
        sys.exit()

    if VERBOSE: sys.stderr.write('Tabulating k-mer model statistics.\n')
    all_kmer_mean_sds = []
    kmer_width = upstrm_bases + dnstrm_bases + 1
    if DEBUG_EST_ALT:
        from sklearn.neighbors import KernelDensity
        kmer_dens = []
        save_x = np.linspace(-5, 5, DEBUG_EST_NUM_KMER_SAVE)[:, np.newaxis]
    for kmer in product(DNA_BASES, repeat=kmer_width):
        kmer = ''.join(kmer)
        kmer_levels = np.concatenate([
            reg_kmer_levels[kmer] for reg_kmer_levels in all_reg_kmer_levels
            if len(reg_kmer_levels[kmer]) > 0])
        if kmer_levels.shape[0] < min_kmer_obs:
            min_obs = min(
                sum(len(reg_levs[''.join(kmer)])
                    for reg_levs in all_reg_kmer_levels)
                for kmer in product(DNA_BASES, repeat=kmer_width))
            sys.stderr.write(
                'ERROR: K-mers represeneted in fewer observations than ' +
                'requested in the provided reads. Consider a shorter ' +
                'k-mer or providing more reads.\n\t' + str(min_obs) +
                ' observations found in least common kmer.\n')
            sys.exit()
        all_kmer_mean_sds.append((kmer, np.median(kmer_levels[:,0]),
                                  np.mean(kmer_levels[:,1])))
        if DEBUG_EST_ALT:
            kmer_kde = KernelDensity(
                kernel='gaussian', bandwidth=DEBUG_EST_BW).fit(
                    kmer_levels[:,0][:,np.newaxis])
            with np.errstate(under='ignore'):
                kmer_dens.append((
                    kmer, np.exp(kmer_kde.score_samples(save_x))))

    if DEBUG_EST_ALT:
        with open('debug_est_standard_ref.density.txt', 'w') as fp:
            fp.write('Kmer\tSignal\tDensity\n')
            fp.write('\n'.join('\t'.join(map(str, (kmer, x, y)))
                               for kmer, dens_i in kmer_dens
                               for x, y in zip(save_x[:,0], dens_i)) + '\n')
    kmer_ref = np.array(
        all_kmer_mean_sds,
        dtype=[('kmer', 'S' + str(kmer_width)),
               ('mean', 'f8'), ('sd', 'f8')])
    if not kmer_specific_sd:
        kmer_ref['sd'] = np.median(kmer_ref['sd'])
    write_tombo_model(kmer_ref, kmer_ref_fn, upstrm_bases)

    return

def parse_kmer_levels(
        all_reads, kmer_width, upstrm_bases, check_min_kmer_batch,
        kmer_obs_thresh, max_kmer_obs, min_kmer_obs_to_est):
    """
    Parse base levels and store grouped by k-mer
    """
    dnstrm_bases = kmer_width - upstrm_bases - 1
    mixed_kmers_levels = dict(
        (''.join(kmer), [])
        for kmer in product(DNA_BASES, repeat=kmer_width))
    # store set of k-mers with enough observations to save on memory footprint
    # while filling more rare k-mers
    completed_kmers = set()
    n_reads = 0
    for r_data in all_reads:
        r_means, r_seq = th.get_read_means_and_seq(r_data)
        if r_means is None: continue
        r_seq = ''.join(r_seq)
        n_reads += 1
        r_kmers = [r_seq[i:i+kmer_width]
                   for i in range(len(r_seq)-kmer_width+1)]
        r_kmer_levels = defaultdict(list)
        for kmer, level in zip(r_kmers, r_means[upstrm_bases:-dnstrm_bases]):
            r_kmer_levels[kmer].append(level)
        # only add observations for k-mers that have not been seen enough times
        # save memory for abundent k-mers
        for kmer in set(r_kmer_levels.keys()).difference(completed_kmers):
            mixed_kmers_levels[kmer].extend(r_kmer_levels[kmer])

        # every check_min_kmer_batch check to see if there are enough
        # observations of each kmer to continue to estimation
        if n_reads % check_min_kmer_batch == 0:
            kmer_levels_counts = sorted([
                (len(kmer_levels), kmer)
                for kmer, kmer_levels in mixed_kmers_levels.iteritems()])
            if kmer_levels_counts[0][0] > kmer_obs_thresh:
                break
            if VERBOSE: sys.stderr.write(
                    '\t' + str(n_reads) + ' reads processed. Current ' +
                    'minimum k-mer observations: ' +
                    str(kmer_levels_counts[0][0]) + ' towards goal of ' +
                    str(kmer_obs_thresh) + '\n')
            # re-compute completed kmers after each batch
            completed_kmers = set()
            for kmer_count, kmer in kmer_levels_counts[::-1]:
                if kmer_count < max_kmer_obs:
                    break
                completed_kmers.add(kmer)

    if min(len(kmer_levels) for kmer_levels in
           mixed_kmers_levels.itervalues()) < kmer_obs_thresh:
        fewest_kmer_obs = min(len(kmer_levels) for kmer_levels in
                              mixed_kmers_levels.itervalues())
        if fewest_kmer_obs < min_kmer_obs_to_est:
            sys.stderr.write(
                '********* ERROR ********\n\tToo few minimal k-mer ' +
                'observations to continue to alternative estimation. ' +
                'Minimal k-mer has ' + str(fewest_kmer_obs) + ' total ' +
                'observations and ' + str(min_kmer_obs_to_est) +
                ' observations per k-mer are required.\n')
            sys.exit()
        sys.stderr.write(
            '********* WARNING ********\n\tToo few minimal k-mer ' +
            'observations. ' +
            'Continuing to estimation using a k-mer with ' +
            str(fewest_kmer_obs) + ' total observations\n')

    return mixed_kmers_levels

def estimate_alt_model(
        f5_dirs, corrected_group, basecall_subgroups,
        standard_ref_fn, alt_base, min_alt_base_frac, kmer_obs_thresh,
        num_sd_thresh, alt_name, check_min_kmer_batch=1000,
        max_kmer_obs=10000, min_kmer_obs_to_est=50):
    """
    Estimate an alternative model from a standard bases only sample with a single,
    known spiked-in alternative base
    """
    from sklearn import mixture

    raw_read_coverage = th.parse_fast5s(
        f5_dirs, corrected_group, basecall_subgroups)
    all_reads = [r_data for cs_reads in raw_read_coverage.values()
                 for r_data in cs_reads]
    np.random.shuffle(all_reads)

    if VERBOSE: sys.stderr.write('Parsing standard model file\n')
    if standard_ref_fn is None:
        standard_ref_fn = get_default_standard_ref(raw_read_coverage)

    standard_ref, upstrm_bases, _, _ = parse_tombo_model(standard_ref_fn)
    kmer_width = len(next(standard_ref.iterkeys()))
    # get mean model SD. most models will be constant, but use mean in case
    model_sd = np.mean(zip(*standard_ref.values())[1])

    if VERBOSE: sys.stderr.write('Parsing k-mer levels from reads\n')
    mixed_kmers_levels = parse_kmer_levels(
        all_reads, kmer_width, upstrm_bases, check_min_kmer_batch,
        kmer_obs_thresh, max_kmer_obs, min_kmer_obs_to_est)

    if VERBOSE: sys.stderr.write(
            'Fitting mixture models to k-mer level distributions\n')
    kmers_mix_fit_data = []
    if DEBUG_EST_ALT:
        from sklearn.neighbors import KernelDensity
        kmer_dens = []
        save_x = np.linspace(-5, 5, DEBUG_EST_NUM_KMER_SAVE)[:, np.newaxis]
    for kmer, norm_levels in mixed_kmers_levels.items():
        # reshape norm levels for sklearn format and conver to numpy array
        norm_levels = np.array(norm_levels).reshape(-1,1)
        clf = mixture.GaussianMixture(n_components=2, covariance_type='tied')
        clf.fit(norm_levels)
        kmers_mix_fit_data.append((kmer, clf.means_[:,0], clf.weights_))
        if DEBUG_EST_ALT:
            kmer_kde = KernelDensity(
                kernel='gaussian', bandwidth=DEBUG_EST_BW).fit(norm_levels)
            with np.errstate(under='ignore'):
                kmer_dens.append((kmer, np.exp(kmer_kde.score_samples(save_x))))
    if VERBOSE: sys.stderr.write('Analyzing fitted mixture distributions\n')
    if DEBUG_EST_ALT:
        with open('debug_est_alt.' + alt_base + '.' + alt_name +
                  '.txt', 'w') as fp:
            fp.write('peaks_diff\tminor_frac\tcontains_alt\tkmer\t' +
                     'sd_width\tmin_frac\n')
            fp.write('\n'.join(
                '\t'.join(map(str, (
                    np.diff(comp_means)[0], min(comp_weights),
                    len(re.findall(alt_base, kmer)) > 0, kmer,
                    model_sd * num_sd_thresh, min_alt_base_frac)))
                for kmer, comp_means, comp_weights in
                kmers_mix_fit_data) + '\n')
        with open('debug_est_alt.' + alt_base + '.' + alt_name +
                  '.density.txt', 'w') as fp:
            fp.write('Kmer\tSignal\tDensity\n')
            fp.write('\n'.join('\t'.join(map(str, (kmer, x, y)))
                               for kmer, dens_i in kmer_dens
                               for x, y in zip(save_x[:,0], dens_i)) + '\n')
    # filter results to get k-mers that contribute to the final alternative model
    alt_model_data = dict(
        (kmer, comp_means) for kmer, comp_means, comp_weights in
        kmers_mix_fit_data
        if np.diff(comp_means)[0] > model_sd * num_sd_thresh and
        min(comp_weights) > min_alt_base_frac and
        len(re.findall(alt_base, kmer)) > 0)

    # determine which mixture component matches is further from the standard
    # model should be the minor component, but definitely should be the
    # one further from the standard base model.
    alt_ref = []
    for kmer, (standard_level, _) in standard_ref.iteritems():
        try:
            comp_levels = alt_model_data[kmer]
            alt_level = comp_levels[np.argmax(np.abs(
                comp_levels - standard_level))]
        except KeyError:
            # if the level did not pass thresholds above use standard base level
            alt_level = standard_level
        alt_ref.append((kmer, alt_level, model_sd))

    alt_ref = np.array(alt_ref, dtype=[('kmer', 'S' + str(kmer_width)),
                                       ('mean', 'f8'), ('sd', 'f8')])

    return alt_ref, upstrm_bases


#########################################
##### Base statistical testing code #####
#########################################

def p_value_to_z_score(pvalue):
    """
    Helper function to convert p-value to z-score
    """
    return -stats.norm.ppf(pvalue)

def z_score_to_p_value(zscore):
    """
    Helper function to convert z-score to p-value
    """
    return stats.norm.cdf(zscore)

def correct_multiple_testing(pvals):
    """
    Use FDR Benjamini-Hochberg multiple testing correction
    """
    pvals = np.asarray(pvals)

    pvals_sortind = np.argsort(pvals)
    pvals_sorted = pvals[pvals_sortind]
    sortrevind = pvals_sortind.argsort()

    nobs = len(pvals)
    ecdffcator = np.arange(1,nobs+1)/float(nobs)
    # ignore underflow values
    with np.errstate(under='ignore'):
        pvals_corrected_raw = pvals_sorted / ecdffcator
    pvals_corrected = np.minimum.accumulate(
        pvals_corrected_raw[::-1])[::-1]
    del pvals_corrected_raw
    pvals_corrected[pvals_corrected>1] = 1

    return pvals_corrected[sortrevind]

def calc_vectorized_fm_pvals(split_pvals, filter_nan=True):
    """
    Compute Fisher's Method p-values in a vectorized fashion
    """
    if filter_nan:
        chi_stats = [np.sum(np.log(base_pvals[~np.isnan(base_pvals)])) * -2
                     for base_pvals in split_pvals]
        chi_shapes = [
            np.sum(~np.isnan(base_pvals)) * 2 for base_pvals in split_pvals]
    else:
        chi_stats = [np.sum(np.log(base_pvals)) * -2
                     for base_pvals in split_pvals]
        chi_shapes = [base_pvals.shape[0] * 2
                     for base_pvals in split_pvals]

    f_pvals = stats.chi2.sf(chi_stats, chi_shapes)
    return f_pvals

def calc_window_fishers_method(pvals, lag):
    """
    Compute Fisher's Method over a moving window across a set of p-values
    """
    assert lag > 0, 'Invalid p-value window provided.'
    width = (lag * 2) + 1
    if pvals.shape[-1] < width:
        raise NotImplementedError, (
            "P-values vector too short for Fisher's Method " +
            "window compuation.")
    with np.errstate(invalid='ignore'):
        pvals = np.maximum(pvals, th.SMALLEST_PVAL)
    log_sums = np.lib.stride_tricks.as_strided(
        np.log(pvals),
        shape=pvals.shape[:-1] + (pvals.shape[-1] - width + 1, width),
        strides=pvals.strides + (pvals.strides[-1],)).sum(-1)
    f_pvals = np.empty(pvals.shape)
    f_pvals[:] = np.NAN
    with np.errstate(invalid='ignore'):
        f_pvals[...,lag:-lag] = stats.chi2.sf(log_sums * -2, width * 2)

    return f_pvals

def calc_window_z_transform(r_means, ref_means, ref_sds, lag):
    """
    Compute Stouffer's Z-transformation across a read
    """
    z_scores = np.abs(r_means - ref_means) / ref_sds
    width = (lag * 2) + 1
    window_z_trans = np.lib.stride_tricks.as_strided(
        z_scores,
        shape=(z_scores.shape[0] - (lag * 2), width),
        strides=z_scores.strides +
        (z_scores.strides[-1],)).sum(1) / np.sqrt(width)
    window_z_trans = np.concatenate([np.array([np.NAN] * lag),
                                     window_z_trans, np.array([np.NAN] * lag)])

    return window_z_trans

def calc_mann_whitney_z_score(samp1, samp2):
    """
    Compute Mann-Whitney z-scores comparing two samples
    """
    s1_len = samp1.shape[0]
    s2_len = samp2.shape[0]
    tot_len = s1_len + s2_len

    all_vals = np.concatenate([samp1, samp2])
    ranks = np.empty(tot_len, int)
    ranks[all_vals.argsort()] = np.arange(1, tot_len + 1)
    s1_ranks_sum = ranks[:s1_len].sum()
    #s2_ranks_sum = ranks[s1_len:].sum()

    u1 = s1_ranks_sum - (s1_len * (s1_len + 1)) / 2
    #u2 = s2_ranks_sum - (s2_len * (s2_len + 1)) / 2

    mu = s1_len * s2_len / 2
    rhou = np.sqrt(s1_len * s2_len * (s1_len + s2_len + 1) / 12)

    z = np.abs(u1 - mu) / rhou

    return z

def add_multiple_testing(all_stats):
    """
    Add multiple testing to a set of statistics
    """
    if len(all_stats) == 0:
        sys.stderr.write(
            '*' * 60 + '\nERROR: No regions contain minimum ' +
            'number of reads.\n' + '*' * 60 + '\n')
        sys.exit()

    # get FDR corrected q-values
    all_stats.sort(order='stat')
    all_stats['mt_stat'] = correct_multiple_testing(all_stats['stat'])

    return all_stats


###################################
##### Model-based re-squiggle #####
###################################

def get_dynamic_prog_params(match_evalue):
    """
    Compute dynamic programming shift parameters from an expected match expected value
    """
    z_shift = stats.halfnorm.expect() + match_evalue
    stay_pen = match_evalue
    return z_shift, stay_pen

def get_begin_nan(arr):
    """
    Find the index of the first NAN value
    """
    tot=0
    for val in iter(arr):
        if not np.isnan(val): break
        tot += 1
    return tot

def get_read_signif_shift_regions(
        z_scores, z_thresh, context_bases, signif_shift_len_thresh=None):
    """
    Identify regions along a read that do not match well with the genomic reference tombo model
    """
    # extend NANs by context_bases to avoid regions extending outside of
    # valid regions over which statistics were computed
    start_nans = get_begin_nan(z_scores) + context_bases
    z_scores[:start_nans] = np.NAN
    end_nans = get_begin_nan(z_scores[::-1]) + context_bases
    if end_nans > 0:
        z_scores[-end_nans:] = np.NAN
    # suppress NAN errors form less than compare
    with np.errstate(invalid='ignore'):
        signif_shift_locs = z_scores > z_thresh
    # find the edges of significantly shifted signal regions
    signif_shift_chngpnts = np.where(np.diff(signif_shift_locs) != 0)[0] + 1

    signif_shift_regs = zip(signif_shift_chngpnts[:-1:2],
                            signif_shift_chngpnts[1::2])
    signif_shift_cntxt_regs = []
    curr_start, curr_end = signif_shift_regs[0]
    for reg_start, reg_end in signif_shift_regs[1:]:
        # if next region overlaps the current region with context
        if reg_start - (context_bases * 2) <= curr_end:
            # extend the current region to cover both regions
            curr_end = reg_end
        else:
            # else check that the region is long enough
            if (signif_shift_len_thresh is None or
                curr_end - curr_start >= signif_shift_len_thresh):
                # and add it to the list of regions to model re-squiggle
                signif_shift_cntxt_regs.append((curr_start - context_bases,
                                                curr_end + context_bases))
            # and set the next region to be the current one
            curr_start, curr_end = reg_start, reg_end

    # add last region
    if (signif_shift_len_thresh is None or
        curr_end - curr_start >= signif_shift_len_thresh):
        signif_shift_cntxt_regs.append((curr_start - context_bases,
                                        curr_end + context_bases))

    return signif_shift_cntxt_regs


################################
##### Base-by-base Testing #####
################################

def calc_llh_ratio(reg_means, reg_ref_means, reg_ref_vars,
                   reg_alt_means, reg_alt_vars):
    """
    Compute log likelihood ratio
    """
    # TODO fix cython implementation
    # compute log likelihood ratio
    # negative value means standard base fits data better
    # positive value means alternative base fits data better
    return (np.sum(np.square(reg_means - reg_alt_means) / reg_alt_vars) +
            np.sum(np.log(reg_alt_vars))) - (
                np.sum(np.square(reg_means - reg_ref_means) / reg_ref_vars) +
                np.sum(np.log(reg_ref_vars)))

def get_reads_ref(ctrl_reg_reads, reverse_strand, reg_start, region_size,
                  min_test_vals, fm_offset):
    """
    Get mean and standard deviation of levels from a sample across the genome
    """
    ctrl_base_events = th.get_reads_events(ctrl_reg_reads, reverse_strand)
    if ctrl_base_events is None:
        raise NotImplementedError
    arr_size = region_size + (fm_offset * 2)
    ctrl_means, ctrl_sds = np.empty(arr_size), np.empty(arr_size)
    ctrl_means[:] = np.NAN
    ctrl_sds[:] = np.NAN
    ctrl_cov = {}
    for pos, pos_events in sorted(ctrl_base_events.iteritems()):
        # if position is past the end of the region return
        if pos - fm_offset >= reg_start + region_size:
            break
        if pos + fm_offset < reg_start:
            continue
        ctrl_cov[pos] = len(pos_events)
        if ctrl_cov[pos] < min_test_vals:
            continue
        ctrl_mean, ctrl_sd = c_mean_std(pos_events)
        ctrl_sd = max(ctrl_sd, MIN_POSITION_SD)
        ctrl_means[pos - reg_start + fm_offset] = ctrl_mean
        ctrl_sds[pos - reg_start + fm_offset] = ctrl_sd

    return ctrl_means, ctrl_sds, ctrl_cov

def get_region_stats(
        chrm, strand, reg_start, reg_reads,
        fm_offset, min_test_vals, region_size, single_read_thresh,
        ctrl_reg_reads, kmer_ref, upstrm_bases, alt_ref, alt_base):
    """
    Compute requested statistics for a specific region of the genome
    """
    if ctrl_reg_reads is None:
        # compute begin and end lag wrt the genome from upstream and downstream
        # which are wrt to the read
        kmer_width = len(next(kmer_ref.iterkeys()))
        dnstrm_bases = kmer_width - upstrm_bases - 1
        begin_lag = upstrm_bases if strand == '+' else dnstrm_bases
        end_lag = dnstrm_bases if strand == '+' else upstrm_bases

    def get_read_comp_stats(r_data, ctrl_means, ctrl_sds):
        def comp_clip_and_flip():
            r_means, _ = th.get_read_means_and_seq(r_data)
            if r_means is None:
                return None, None

            read_start, read_end = r_data.start, r_data.end
            if read_start + fm_offset < reg_start:
                num_start_clip = reg_start - (read_start + fm_offset)
                read_start = reg_start - fm_offset
                if strand == '+':
                    r_means = r_means[num_start_clip:]
                else:
                    r_means = r_means[:-num_start_clip]
            if read_end - fm_offset > reg_start + region_size:
                num_end_clip = (read_end - fm_offset) - (reg_start + region_size)
                read_end = reg_start + region_size + fm_offset
                if strand == '+':
                    r_means = r_means[:-num_end_clip]
                else:
                    r_means = r_means[num_end_clip:]

            # flip means to match genomic positions
            if strand == '-':
                r_means = r_means[::-1]

            return r_means, read_start, read_end

        def get_read_comp_z_score(r_means, read_start, read_end):
            r_z_scores = np.abs(
                r_means - ctrl_means[read_start-reg_start+fm_offset:
                                     read_end-reg_start+fm_offset]) / ctrl_sds[
                                         read_start-reg_start+fm_offset:
                                         read_end-reg_start+fm_offset]

            return r_z_scores

        def get_pvals(r_z_scores):
            # mask out nan z-scores for efficient CDF computation
            r_poss = np.where(np.logical_not(np.isnan(r_z_scores)))[0]
            r_pvals = np.empty(r_z_scores.shape)
            r_pvals[:] = np.NAN
            valid_r_pvals = stats.norm.cdf(-r_z_scores[r_poss]) * 2.0
            r_pvals[r_poss] = valid_r_pvals

            return r_pvals, r_poss

        r_means, read_start, read_end = comp_clip_and_flip()
        r_z_scores = get_read_comp_z_score(r_means, read_start, read_end)

        if np.sum(np.logical_not(np.isnan(r_z_scores))) == 0:
            return None, None
        r_pvals, r_poss = get_pvals(r_z_scores)
        if fm_offset > 0:
            try:
                r_pvals = calc_window_fishers_method(r_pvals, fm_offset)
            except NotImplementedError:
                return None, None
            r_poss = np.where(np.logical_not(np.isnan(r_pvals)))[0]
            r_pvals = r_pvals[r_poss]
        else:
            r_pvals = r_pvals[r_poss]

        r_poss += read_start

        return r_pvals, r_poss

    def clip_and_flip_data(r_data):
        def get_mean_seq():
            r_means, r_seq = th.get_read_means_and_seq(r_data)
            if r_means is None or r_seq is None or len(r_seq) <= kmer_width:
                raise NotImplementedError
            r_seq = ''.join(r_seq)

            read_start, read_end = r_data.start, r_data.end
            # clip read if it extends outside the current genomic region, so
            # stats are only computed within this region
            if read_start + begin_lag + fm_offset < reg_start:
                num_start_clip = reg_start - (read_start + begin_lag + fm_offset)
                read_start = reg_start - begin_lag - fm_offset
                if strand == '+':
                    r_means = r_means[num_start_clip:]
                    r_seq = r_seq[num_start_clip:]
                else:
                    r_means = r_means[:-num_start_clip]
                    r_seq = r_seq[:-num_start_clip]
            if read_end - end_lag - fm_offset > reg_start + region_size:
                num_end_clip = (read_end - end_lag - fm_offset) - (
                    reg_start + region_size)
                read_end = reg_start + region_size + end_lag + fm_offset
                if strand == '+':
                    r_means = r_means[:-num_end_clip]
                    r_seq = r_seq[:-num_end_clip]
                else:
                    r_means = r_means[num_end_clip:]
                    r_seq = r_seq[num_end_clip:]

            # if this read does not cover enough of this region for stat
            # computation raise an error to be handled below
            if len(r_seq) < kmer_width:
                raise NotImplementedError

            return r_means, r_seq, read_start, read_end

        def get_refs(r_seq):
            # get stat lookups from read in native strand then flip if
            # read is on reverse strand
            if strand == '-':
                r_ref_means, r_ref_sds = map(np.array, zip(*[
                    kmer_ref[r_seq[i:i+kmer_width]]
                    for i in range(len(r_seq)-kmer_width+1)][::-1]))
                if alt_ref is not None:
                    r_alt_means, r_alt_sds = map(np.array, zip(*[
                        alt_ref[r_seq[i:i+kmer_width]]
                        for i in range(len(r_seq)-kmer_width+1)][::-1]))
            else:
                r_ref_means, r_ref_sds = map(np.array, zip(*[
                    kmer_ref[r_seq[i:i+kmer_width]]
                    for i in range(len(r_seq)-kmer_width+1)]))
                if alt_ref is not None:
                    r_alt_means, r_alt_sds = map(np.array, zip(*[
                        alt_ref[r_seq[i:i+kmer_width]]
                        for i in range(len(r_seq)-kmer_width+1)]))
            if alt_ref is None:
                return r_ref_means, r_ref_sds, None, None
            return r_ref_means, r_ref_sds, r_alt_means, r_alt_sds

        r_means, r_seq, read_start, read_end = get_mean_seq()
        r_ref_means, r_ref_sds, r_alt_means, r_alt_sds = get_refs(r_seq)

        if strand == '-':
            # reverse means and seq to match genomic order
            r_means = r_means[::-1]
            r_seq = r_seq[::-1]
        # clip means and seq that won't have model info due to only having read
        # sequence
        r_means = r_means[begin_lag:-end_lag]
        r_seq = r_seq[begin_lag:-end_lag]
        read_start += begin_lag
        read_end -= end_lag

        return (r_means, r_seq, r_ref_means, r_ref_sds, read_start,
                read_end, r_alt_means, r_alt_sds)

    def get_read_stats(r_data):
        try:
            (r_means, r_seq, r_ref_means, r_ref_sds,
             read_start, read_end, _, _) = clip_and_flip_data(r_data)
        except NotImplementedError:
            return None, None

        z_scores = np.abs(r_means - r_ref_means) / r_ref_sds
        r_pvals = stats.norm.cdf(-z_scores) * 2.0
        if fm_offset > 0:
            try:
                r_pvals = calc_window_fishers_method(r_pvals, fm_offset)
            except NotImplementedError:
                return None, None

        # ignore errors in max over NAN values if fisher's method was used
        with np.errstate(invalid='ignore'):
            r_pvals = np.maximum(r_pvals, th.SMALLEST_PVAL)

        r_poss = np.arange(read_start, read_end)

        return r_pvals, r_poss

    def get_read_alt_stats(r_data):
        try:
            (r_means, r_seq, r_ref_means, r_ref_sds, read_start, read_end,
             r_alt_means, r_alt_sds) = clip_and_flip_data(r_data)
        except NotImplementedError:
            return None, None
        r_ref_vars = np.square(r_ref_sds)
        r_alt_vars = np.square(r_alt_sds)

        alt_base_poss = []
        log_lh_ratios = []
        # note search space is clipped since all k-mers covering the position
        # of interest must be valid
        for alt_base_pos in re.finditer(alt_base, r_seq[end_lag:-begin_lag]):
            alt_pos = alt_base_pos.start() + end_lag
            alt_base_poss.append(alt_pos + read_start)
            # TODO cython version needs to be checked for bugs
            # probably overflow of typed cython variables
            #pos_lh_ratio = c_calc_llh_ratio(
            pos_lh_ratio = calc_llh_ratio(
                r_means[alt_pos-end_lag:alt_pos+begin_lag],
                r_ref_means[alt_pos-end_lag:alt_pos+begin_lag],
                r_ref_vars[alt_pos-end_lag:alt_pos+begin_lag],
                r_alt_means[alt_pos-end_lag:alt_pos+begin_lag],
                r_alt_vars[alt_pos-end_lag:alt_pos+begin_lag])
            log_lh_ratios.append(pos_lh_ratio)

        return log_lh_ratios, alt_base_poss


    if ctrl_reg_reads is not None:
        try:
            ctrl_means, ctrl_sds, ctrl_cov = get_reads_ref(
                ctrl_reg_reads, strand == '-', reg_start, region_size,
                min_test_vals, fm_offset)
        except NotImplementedError:
            # if there are no events in this window return
            return None

    reg_stats, reg_poss = [], []
    for r_data in reg_reads:
        if ctrl_reg_reads is not None:
            r_stats, r_poss = get_read_comp_stats(r_data, ctrl_means, ctrl_sds)
        elif alt_ref is None:
            r_stats, r_poss = get_read_stats(r_data)
        else:
            r_stats, r_poss = get_read_alt_stats(r_data)
        if r_stats is None: continue
        reg_stats.append(r_stats)
        reg_poss.append(r_poss)

    if len(reg_stats) == 0:
        return None

    reg_stats = np.concatenate(reg_stats)
    reg_poss = np.concatenate(reg_poss)
    # remove nans possibly introduced by fisher's method calculcations
    valid_poss = ~np.isnan(reg_stats)
    reg_poss = reg_poss[valid_poss]
    reg_stats = reg_stats[valid_poss]
    assert reg_poss.shape[0] == reg_stats.shape[0], '\t'.join(map(str, (
        reg_poss.shape[0], reg_stats.shape[0])))

    # get order of all bases from position array
    as_reg_poss = np.argsort(reg_poss)
    # sort all positions from all reads
    reg_poss = reg_poss[as_reg_poss]
    # get unique tested genomic positions across all reads
    us_reg_poss = np.unique(reg_poss)

    # then sort the stats array by genomic position and
    # split into stats by genomic base position
    reg_stats = np.split(
        reg_stats[as_reg_poss],
        np.where(np.concatenate([[0,], np.diff(reg_poss)]) > 0)[0])

    if alt_ref is None:
        reg_combine_read_stats = calc_vectorized_fm_pvals(
            reg_stats, filter_nan=False)
    else:
        with np.errstate(invalid='ignore'):
            reg_combine_read_stats = [
                np.mean(base_stats) for base_stats in reg_stats]

    reg_cov = [base_stats.shape[0] for base_stats in reg_stats]

    reg_frac_standard_base = np.array([
        np.sum(base_stats > single_read_thresh) / float(base_stats.shape[0])
        for base_stats in reg_stats])
    if alt_ref is None:
        reg_frac_alt_base = 1 - reg_frac_standard_base
    else:
        reg_frac_alt_base = np.array([
            np.sum(base_stats < -single_read_thresh) / float(base_stats.shape[0])
            for base_stats in reg_stats])

    if ctrl_reg_reads is None:
        ctrl_cov = repeat(0)
    else:
        ctrl_cov = [ctrl_cov[pos] if pos in ctrl_cov else 0
                    for pos in us_reg_poss]

    reg_stats = np.array(
        [pos_stats for pos_stats in zip(
            reg_combine_read_stats, repeat(np.NAN),
            reg_frac_standard_base, reg_frac_alt_base,
            us_reg_poss, repeat(chrm), repeat(strand), reg_cov, ctrl_cov)
         if pos_stats[-2] >= min_test_vals and
         not np.isnan(pos_stats[0])],
        dtype=[('stat', 'f8'), ('mt_stat', 'f8'),
               ('frac', 'f8'), ('alt_frac', 'f8'),
               ('pos', 'u4'), ('chrm', 'S32'), ('strand', 'S1'),
               ('cov', 'u4'), ('control_cov', 'u4')])
    if reg_stats.shape[0] == 0:
        return None

    return reg_stats

def _test_signif_worker(
        region_q, stats_q, raw_read_coverage, fm_offset, min_test_vals,
        single_read_thresh, region_size,
        ctrl_read_coverage, kmer_ref, upstrm_bases, alt_ref, alt_base):
    ctrl_reg_reads = None
    while not region_q.empty():
        try:
            chrm, strand, reg_start = region_q.get(block=False)
        except Queue.Empty:
            break

        reg_reads = [r_data for r_data in raw_read_coverage[(chrm, strand)]
                     if not (r_data.start >= reg_start + region_size or
                             r_data.end <= reg_start)]
        if len(reg_reads) == 0:
            if VERBOSE:
                sys.stderr.write('.')
                sys.stderr.flush()
            continue

        if ctrl_read_coverage is not None:
            ctrl_reg_reads = [
                r_data for r_data in ctrl_read_coverage[(chrm, strand)]
                if not (r_data.start >= reg_start + region_size or
                        r_data.end <= reg_start)]
        reg_stats = get_region_stats(
            chrm, strand, reg_start, reg_reads,
            fm_offset, min_test_vals, region_size, single_read_thresh,
            ctrl_reg_reads, kmer_ref, upstrm_bases, alt_ref, alt_base)
        if reg_stats is not None:
            stats_q.put(reg_stats)
        if VERBOSE:
            sys.stderr.write('.')
            sys.stderr.flush()

    return

if PROFILE_SIGNIF:
    _test_signif_wrapper = _test_signif_worker
    def _test_signif_worker(*args):
        import cProfile
        cProfile.runctx('_test_signif_wrapper(*args)', globals(), locals(),
                        filename='test_signif.prof')
        return

def test_significance(
        raw_read_coverage, min_test_vals, fm_offset, single_read_thresh,
        region_size, num_processes, ctrl_read_coverage=None,
        kmer_ref=None, upstrm_bases=None, alt_ref=None, alt_base=None):
    """
    Test for significant shifted signal in mutliprocessed batches
    """
    manager = mp.Manager()
    region_q = manager.Queue()
    stats_q = manager.Queue()
    # split chromosomes into separate regions to process independently
    chrm_sizes = th.get_chrm_sizes(raw_read_coverage, ctrl_read_coverage)
    num_regions = 0
    for chrm, chrm_len in chrm_sizes.items():
        plus_covered = (chrm, '+') in raw_read_coverage
        minus_covered = (chrm, '-') in raw_read_coverage
        for reg_start in range(0, chrm_len, region_size):
            if plus_covered:
                region_q.put((chrm, '+', reg_start))
                num_regions += 1
            if minus_covered:
                region_q.put((chrm, '-', reg_start))
                num_regions += 1

    if VERBOSE: sys.stderr.write(
            'Performing significance testing across ' + str(num_regions) +
            ' regions. (Will print a dot for each batch completed)\n')
    test_args = (
        region_q, stats_q, raw_read_coverage, fm_offset, min_test_vals,
        single_read_thresh, region_size,
        ctrl_read_coverage, kmer_ref, upstrm_bases, alt_ref, alt_base)
    test_ps = []
    for p_id in xrange(num_processes):
        p = mp.Process(target=_test_signif_worker, args=test_args)
        p.start()
        test_ps.append(p)

    all_reg_stats = []
    while any(p.is_alive() for p in test_ps):
        try:
            reg_stats = stats_q.get(block=False)
            all_reg_stats.append(reg_stats)
        except Queue.Empty:
            sleep(1)
            continue
    while not stats_q.empty():
        reg_stats = stats_q.get(block=False)
        all_reg_stats.append(reg_stats)
    if VERBOSE: sys.stderr.write('\n')

    if len(all_reg_stats) == 0:
        sys.stderr.write('********** ERROR *********\n\tNo genomic positions ' +
                         'contain --minimum-test-reads.\n')
        sys.exit()
    # put all stats back together
    all_stats = np.concatenate(all_reg_stats)

    if alt_ref is None:
        if VERBOSE: sys.stderr.write('Performing multiple testing correction.\n')
        all_stats = add_multiple_testing(all_stats)

    return all_stats


##########################
##### Statistics I/O #####
##########################

def parse_stats(stats_fn):
    """
    Parse a tombo statistics file
    """
    if stats_fn is None or not os.path.isfile(stats_fn):
            sys.stderr.write(
                '*' * 60 + '\nERROR: No statistics file provided.\n' +
                '*' * 60 + '\n')
            sys.exit()

    try:
        with h5py.File(stats_fn, 'r') as stats_fp:
            all_stats = stats_fp['stats'].value
            try:
                stat_type = stats_fp.attrs['stat_type']
            except:
                # if this is the old stats file assume sample compare
                stat_type = 'sample_compare'
    except:
        sys.stderr.write(
            '*' * 60  + '\nERROR: Attempt to load statistics ' +
            'file failed. May be an old version of statistics ' +
            'file. Try deleting statistics file and ' +
            'recalculating using current tombo version.\n' +
            '*' * 60 + '\n')
        sys.exit()

    return all_stats, stat_type

def write_stats(all_stats, stats_bsnm, stat_type):
    """
    Write a tombo statistics file
    """
    if VERBOSE: sys.stderr.write(
            'Saving signal shift significance testing results.\n')
    if stat_type == 'model_compare':
        # for alternative model testing, write one stats file per
        # alternative model
        for alt_name, alt_stats in all_stats:
            with h5py.File(stats_bsnm + '.' + alt_name +
                           '.tombo.stats', 'w') as stats_fp:
                stats_fp.create_dataset(
                    'stats', data=alt_stats, compression="gzip")
                stats_fp.attrs['stat_type'] = stat_type
    else:
        with h5py.File(stats_bsnm + '.tombo.stats', 'w') as stats_fp:
            stats_fp.create_dataset('stats', data=all_stats, compression="gzip")
            stats_fp.attrs['stat_type'] = stat_type

    return


##########################
##### Main functions #####
##########################

def test_shifts_main(args):
    global VERBOSE
    VERBOSE = not args.quiet
    th.VERBOSE = VERBOSE

    # apply single read threshold defaults
    single_read_thresh = args.single_read_threshold
    if single_read_thresh is None:
        if args.alternate_model_filenames is not None:
            single_read_thresh = 2.0
        else:
            single_read_thresh = 0.01

    raw_read_coverage = th.parse_fast5s(
        args.fast5_basedirs, args.corrected_group, args.basecall_subgroups)
    # if second set of reads is prodived, perform comparison testing
    if args.control_fast5_basedirs is not None:
        stat_type = 'sample_compare'
        if VERBOSE: sys.stderr.write(
                'Performing two-sample comparison significance testing.\n')
        ctrl_read_coverage = th.parse_fast5s(
            args.control_fast5_basedirs, args.corrected_group,
            args.basecall_subgroups)
        all_stats = test_significance(
            raw_read_coverage, args.minimum_test_reads,
            args.fishers_method_context, single_read_thresh,
            args.multiprocess_region_size, args.processes,
            ctrl_read_coverage=ctrl_read_coverage)
    else:
        tb_model_fn = args.tombo_model_filename
        if tb_model_fn is None:
            tb_model_fn = get_default_standard_ref(raw_read_coverage)
        kmer_ref, upstrm_bases, _, _ = parse_tombo_model(tb_model_fn)

        # if no alt model provided perform de novo testing for shifts
        # from a standard model
        if (args.alternate_model_filenames is None and
            args.alternate_bases is None):
            stat_type = 'de_novo'
            if VERBOSE: sys.stderr.write(
                    'Performing de novo model testing against a standard model\n')
            all_stats = test_significance(
                raw_read_coverage, args.minimum_test_reads,
                args.fishers_method_context, single_read_thresh,
                args.multiprocess_region_size, args.processes, kmer_ref=kmer_ref,
                upstrm_bases=upstrm_bases)
        # else perform comparison model testing
        else:
            stat_type = 'model_compare'
            if VERBOSE: sys.stderr.write(
                    'Performing alternative model testing\n')
            kmer_width = len(next(kmer_ref.iterkeys()))
            if args.alternate_model_filenames is not None:
                alt_refs = parse_tombo_models(
                    args.alternate_model_filenames, upstrm_bases, kmer_width)
            else:
                alt_refs = load_alt_refs(
                    args.alternate_bases, raw_read_coverage,
                    upstrm_bases, kmer_width)
            if len(alt_refs) == 0:
                sys.stderr.write(
                    '********* ERROR *********\n\tNo alternative models ' +
                    'successfully loaded\n')
                sys.exit()

            all_stats = []
            for alt_name, alt_ref, alt_base in alt_refs:
                if VERBOSE: sys.stderr.write(
                        'Performing alternative model testing against ' +
                        alt_name + ' model\n')
                all_stats.append((alt_name, test_significance(
                    raw_read_coverage, args.minimum_test_reads, 0,
                    single_read_thresh, args.multiprocess_region_size,
                    args.processes, kmer_ref=kmer_ref,
                    upstrm_bases=upstrm_bases, alt_ref=alt_ref,
                    alt_base=alt_base)))
    # TODO add comparison to processed genome reference determined by
    # deep learning performed on the genomic sequence

    write_stats(all_stats, args.statistics_file_basename, stat_type)

    return

def write_kmer_ref_main(args):
    global VERBOSE
    VERBOSE = not args.quiet
    th.VERBOSE = VERBOSE

    if min(args.upstream_bases, args.downstream_bases) == 0:
        sys.stderr.write(
            '********** ERROR *********\n\tContext upstream and downstream ' +
            'must be greater than 0 for model estimation.\n')
        sys.exit()

    estimate_kmer_model(
        args.fast5_basedirs, args.corrected_group, args.basecall_subgroups,
        args.tombo_model_filename, args.minimum_test_reads,
        args.upstream_bases, args.downstream_bases,
        args.minimum_kmer_observations, args.kmer_specific_sd,
        args.coverage_threshold, args.estimate_mean,
        args.multiprocess_region_size, args.processes)

    return

def write_alt_ref_main(args):
    global VERBOSE
    VERBOSE = not args.quiet
    th.VERBOSE = VERBOSE

    min_alt_frac = args.min_alt_base_percentage / 100.0
    alt_ref, upstrm_bases = estimate_alt_model(
        args.fast5_basedirs, args.corrected_group, args.basecall_subgroups,
        args.tombo_model_filename, args.alternate_model_base, min_alt_frac,
        args.minimum_kmer_observations, args.sd_threshold,
        args.alternate_model_name)
    write_tombo_model(alt_ref, args.alternate_model_filename, upstrm_bases,
                      args.alternate_model_base, args.alternate_model_name)

    return


if __name__ == '__main__':
    raise NotImplementedError, (
        'This is a module. See commands with `tombo -h`')
