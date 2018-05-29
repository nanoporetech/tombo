from __future__ import division, unicode_literals, absolute_import

from builtins import int, range, dict, map, zip, object

import os
import io
import re
import sys
import queue
import random
import pkg_resources

# Future warning from cython in h5py
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
import h5py

import numpy as np
np.seterr(all='raise')

from tqdm import tqdm
from time import sleep
from operator import itemgetter
from scipy import stats, optimize
from collections import defaultdict
from scipy.spatial.distance import pdist
from multiprocessing import Process, Queue, Pipe
from numpy.lib.recfunctions import append_fields
from scipy.cluster.hierarchy import single, leaves_list
from itertools import repeat, product, count, combinations

if sys.version_info[0] > 2:
    unicode = str

# import tombo functions
from . import tombo_helper as th

from .c_helper import c_mean_std, c_apply_outlier_thresh, c_new_means, \
    c_calc_llh_ratio, c_calc_llh_ratio_const_var, \
    c_calc_scaled_llh_ratio_const_var
from ._default_parameters import (
    SMALLEST_PVAL, MIN_POSITION_SD, STANDARD_MODELS, ALTERNATE_MODELS,
    MIN_KMER_OBS_TO_EST, ALT_EST_BATCH, MAX_KMER_OBS, NUM_DENS_POINTS,
    LLR_THRESH, SAMP_COMP_THRESH, DE_NOVO_THRESH, KERNEL_DENSITY_RANGE,
    ROC_PLOT_POINTS, NANOPOLISH_CENTRAL_POS, NUM_READS_FOR_SCALE,
    ROBUST_QUANTS, MAX_POINTS_FOR_THEIL_SEN, NUM_READS_TO_ADJUST_MODEL,
    OCLLHR_SCALE, OCLLHR_HEIGHT, OCLLHR_POWER)

VERBOSE = False

_PROFILE_SIGNIF = False
_PROFILE_EST_REF = False
_PROFILE_ALT_EST = False

_DEBUG_EST_STD = False
_DEBUG_EST_BW = 0.05
_DEBUG_EST_NUM_KMER_SAVE = 500

PER_READ_BLOCKS_QUEUE_LIMIT = 5

DNA_BASES = ['A','C','G','T']

HALF_NORM_EXPECTED_VAL = stats.halfnorm.expect()

STANDARD_MODEL_NAME = 'standard'

SAMP_COMP_TXT = 'sample_compare'
DE_NOVO_TXT = 'de_novo'
ALT_MODEL_TXT = 'model_compare'

ALT_MODEL_SEP_CHAR = '_'

NORM_TYPES = ('none', 'pA', 'pA_raw', 'median', 'robust_median',
              'median_const_scale')

# options specifying testing methods
# assume constant SD in model to save on computation
CONST_SD_MODEL = True


#############################################
##### Pair-wise Distance and Clustering #####
#############################################

def order_reads(log_r_pvals):
    """
    Compute order of reads based on log p-values
    """
    if log_r_pvals.shape[0] == 1:
        return [0]
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
        except queue.Empty:
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


##################################
###### Signal Normalization ######
##################################

def get_valid_cpts(norm_signal, running_stat_width, num_events):
    """
    DEPRECATED. Hook still included in re-squiggle, but commented out.

    Get valid changepoints given largest differences in neighboring
    moving windows

    Note that this method is completely vectorized, but allows segments
    as small as 2 observations. This should be okay R9+, but is problematic
    for <=R7 and RNA
    """
    raw_cumsum = np.cumsum(np.concatenate([[0], norm_signal[:-1]]))
    # get difference between all neighboring running_stat_width regions
    running_diffs = np.abs(
        (2 * raw_cumsum[running_stat_width:-running_stat_width]) -
        raw_cumsum[:-2*running_stat_width] -
        raw_cumsum[2*running_stat_width:])
    not_peaks = np.logical_not(np.logical_and(
        running_diffs > np.concatenate([[0], running_diffs[:-1]]),
        running_diffs > np.concatenate([running_diffs[1:], [0]])))
    running_diffs[not_peaks] = 0
    valid_cpts = np.argsort(
        running_diffs)[::-1][:num_events].astype(np.int64) + running_stat_width

    return valid_cpts

def calc_kmer_fitted_shift_scale(
        prev_shift, prev_scale, r_event_means, r_model_means,
        r_model_inv_vars=None, method='theil_sen'):
    """
    Compute fitted shift and scale parameters based on read sequence
    """
    if method == 'robust':
        def read_lad_objective(x):
            return np.sum(np.abs(((r_event_means - x[0]) / x[1]) -
                                 r_model_means))

        shift_corr_factor, scale_corr_factor = optimize.minimize(
            read_lad_objective, np.array([0,1]), method='nelder-mead',
            options={'xtol': 1e-8}).x
    elif method == 'theil_sen':
        n_points = r_model_means.shape[0]
        # potentially sample points for long reads (>1kb)
        if r_model_means.shape[0] > MAX_POINTS_FOR_THEIL_SEN:
            n_points = MAX_POINTS_FOR_THEIL_SEN
            samp_ind = np.random.choice(
                r_model_means.shape[0], n_points, replace=False)
            r_model_means = r_model_means[samp_ind]
            r_event_means = r_event_means[samp_ind]
        # compute Theil-Sen slope estimator
        # despite computing each diff twice this vectorized solution is about
        # 10X faster than a list comprehension approach
        delta_event = r_event_means[:, np.newaxis] - r_event_means
        delta_model = r_model_means[:, np.newaxis] - r_model_means
        slopes = delta_model[delta_event > 0] / delta_event[delta_event > 0]
        slopes.sort()
        slope = np.median(slopes)
        inter = np.median(r_model_means - (slope * r_event_means))
        if slope == 0:
            raise NotImplementedError(
                'Read failed sequence-based signal re-scaling parameter ' +
                'estimation.')
        # convert to shift and scale parameters (e.g. (obs - shift) / scale)
        scale_corr_factor = 1 / slope
        shift_corr_factor = -inter / slope
    elif method == 'mom':
        model_mean_var = r_model_means * r_model_inv_vars
        # prep kmer model coefficient matrix for the k-mers from this read
        model_mean_var_sum = model_mean_var.sum()
        coef_mat = np.array((
            (r_model_inv_vars.sum(), model_mean_var_sum),
            (model_mean_var_sum, (model_mean_var * r_model_means).sum())))

        # prep dependent values from this reads true events
        r_event_var = r_event_means * r_model_inv_vars
        r_event_var_mean = r_event_var * r_model_means
        dep_vect = np.array((r_event_var.sum(), r_event_var_mean.sum()))

        shift_corr_factor, scale_corr_factor = np.linalg.solve(
            coef_mat, dep_vect)
    else:
        th._error_message_and_exit(
            'Invalid k-mer fitted normalization parameter method: ' + method +
            '\n\t\tValid methods are "robust" and "mom".')

    # apply shift and scale values fitted from kmer conditional model
    shift = prev_shift + (shift_corr_factor * prev_scale)
    scale = prev_scale * scale_corr_factor

    return shift, scale, shift_corr_factor, scale_corr_factor

def estimate_global_scale(fast5_fns, num_reads=NUM_READS_FOR_SCALE):
    if VERBOSE: th._status_message('Estimating global scale parameter.')
    np.random.shuffle(fast5_fns)
    read_mads = []
    if VERBOSE:
        bar = tqdm(total=num_reads, desc='Total reads processed', smoothing=0)
    for fast5_fn in fast5_fns:
        try:
            with h5py.File(fast5_fn, 'r') as fast5_data:
                all_sig = th.get_raw_read_slot(fast5_data)['Signal'].value
            shift = np.median(all_sig)
            read_mads.append(np.median(np.abs(all_sig - shift)))
            if VERBOSE: bar.update(1)
        except:
            continue
        if len(read_mads) >= num_reads:
            break

    if VERBOSE: bar.close()
    if len(read_mads) == 0:
        th._error_message_and_exit(
            'No reads contain raw signal for ' +
            'global scale parameter estimation.')
    if len(read_mads) < num_reads:
        th._warning_message(
            'Few reads contain raw signal for global scale parameter ' +
            'estimation. Results may not be optimal.')

    return np.mean(read_mads)

def normalize_raw_signal(
        all_raw_signal, read_start_rel_to_raw, read_obs_len,
        norm_type=None, channel_info=None, outlier_thresh=None,
        scale_values=None, event_means=None, model_means=None,
        model_inv_vars=None, const_scale=None):
    """
    Apply scaling and windsorizing parameters to normalize raw signal
    """
    if norm_type not in NORM_TYPES and scale_values is None:
        raise NotImplementedError(
            'Normalization type ' + norm_type + ' is not a valid ' +
            'option and shift or scale parameters were not provided.')

    raw_signal = all_raw_signal[read_start_rel_to_raw:
                                read_start_rel_to_raw + read_obs_len]
    if scale_values is None:
        if norm_type == 'none':
            shift, scale = 0, 1
        elif norm_type in ('pA_raw', 'pA'):
            # correct raw signal as described here:
            #     https://community.nanoporetech.com
            #           /posts/squiggle-plot-for-raw-data
            shift, scale = (
                -1 * channel_info.offset,
                channel_info.digitisation / channel_info.range)
            if norm_type == 'pA':
                # perform k-mer model fitted correction as in
                # nanocorr/nanopolish/albacore(pre-RNN)
                shift, scale, _, _ = calc_kmer_fitted_shift_scale(
                    shift, scale, event_means, model_means, model_inv_vars,
                    method='mom')
        elif norm_type == 'median':
            shift = np.median(raw_signal)
            scale = np.median(np.abs(raw_signal - shift))
        elif norm_type == 'median_const_scale':
            assert const_scale is not None
            shift = np.median(raw_signal)
            scale = const_scale
        elif norm_type == 'robust_median':
            shift = np.mean(np.percentile(raw_signal, ROBUST_QUANTS))
            scale = np.median(np.abs(raw_signal - read_robust_med))
    else:
        shift = scale_values.shift
        scale = scale_values.scale

    norm_signal = (raw_signal - shift) / scale

    # windsorize the raw signal
    lower_lim, upper_lim = None, None
    if outlier_thresh is not None or scale_values is not None:
        if outlier_thresh is not None:
            read_med = np.median(norm_signal)
            read_mad = np.median(np.abs(norm_signal - read_med))
            lower_lim = read_med - (read_mad * outlier_thresh)
            upper_lim = read_med + (read_mad * outlier_thresh)
        else:
            lower_lim = scale_values.lower_lim
            upper_lim = scale_values.upper_lim
        norm_signal = c_apply_outlier_thresh(norm_signal, lower_lim, upper_lim)

    return norm_signal, th.scaleValues(shift, scale, lower_lim, upper_lim)


#############################
##### Tombo Model Class #####
#############################

class TomboModel(object):
    """
    Load, store and access Tombo model attributes and sequence-based expected
    mean and standard deviation levels (median normalization only)
    """
    def center_model(self, shift_corr_factor, scale_corr_factor):
        centered_means = {}
        for kmer, k_mean in self.means.items():
            centered_means[kmer] = (
                k_mean * scale_corr_factor) + shift_corr_factor

        self.means = centered_means

        return

    def make_constant_sd(self):
        med_sd = np.median(list(self.sds.values()))
        self.sds = dict((kmer, med_sd) for kmer in self.sds)
        return

    def write_model(self, ref_fn, alt_base=None, alt_name=None):
        # Explicity use btype string names for py3 compatiblity as well as
        # pickle-ability of numpy arrays for consistency. See discussion here:
        # https://github.com/numpy/numpy/issues/2407
        ref_for_file = np.array(
            [(kmer, self.means[kmer], self.sds[kmer]) for kmer in self.means],
            dtype=[(str('kmer'), 'S' + unicode(self.kmer_width)),
                   (str('mean'), 'f8'), (str('sd'), 'f8')])

        with h5py.File(ref_fn, 'w') as ref_fp:
            ref_fp.create_dataset('model', data=ref_for_file, compression="gzip")
            ref_fp.attrs['central_pos'] = self.central_pos
            if alt_base is None:
                ref_fp.attrs['model_name'] = STANDARD_MODEL_NAME
            else:
                ref_fp.attrs['model_name'] = alt_name
                ref_fp.attrs['alt_base'] = alt_base

        return

    def _parse_tombo_model(self):
        """
        Parse a tombo model file
        """
        try:
            with h5py.File(self.ref_fn, 'r') as ref_fp:
                ref_raw = ref_fp['model'].value
                central_pos = ref_fp.attrs['central_pos']
                model_name = ref_fp.attrs['model_name']

                try:
                    model_name = model_name.decode()
                except (AttributeError, TypeError):
                    pass

                try:
                    alt_base = ref_fp.attrs['alt_base']
                except:
                    alt_base = None
                try:
                    alt_base = alt_base.decode()
                except (AttributeError, TypeError):
                    pass

        except:
            th._error_message_and_exit('Invalid tombo model file provided: '
                                       + unicode(self.ref_fn))

        mean_ref = {}
        sd_ref = {}
        for kmer, kmer_mean, kmer_std in ref_raw:
            kmer = kmer.decode()
            mean_ref[kmer] = kmer_mean
            sd_ref[kmer] = kmer_std

        self.means = mean_ref
        self.sds = sd_ref
        self.central_pos = central_pos
        self.alt_base = alt_base
        self.name = model_name

        return

    def _parse_text_model(self):
        """
        Parse a text model file (such as those from nanopolish)
        """
        try:
            mean_ref, sd_ref = {}, {}
            with io.open(self.ref_fn) as fp:
                for line in fp:
                    if line.startswith('#'): continue
                    try:
                        kmer, kmer_mean, kmer_sd = line.split()[:3]
                        kmer_mean, kmer_sd = map(float, (kmer_mean, kmer_sd))
                    except ValueError:
                        # header or other non-kmer field
                        continue
                    mean_ref[kmer] = kmer_mean
                    sd_ref[kmer] = kmer_sd
        except:
            th._error_message_and_exit('Invalid text pA model file provided: '
                                       + unicode(self.ref_fn))

        self.means = mean_ref
        self.sds = sd_ref
        self.central_pos = NANOPOLISH_CENTRAL_POS
        self.alt_base = None
        self.name = STANDARD_MODEL_NAME

        return

    def _load_std_model(self, kmer_ref, central_pos):
        mean_ref = {}
        sd_ref = {}
        for kmer, kmer_mean, kmer_std in kmer_ref:
            # reference may or may not be stored as a numpy array
            try:
                kmer = kmer.decode()
            except AttributeError:
                pass
            mean_ref[kmer] = kmer_mean
            sd_ref[kmer] = kmer_std

        self.means = mean_ref
        self.sds = sd_ref
        self.central_pos = central_pos
        self.alt_base = None
        self.name = STANDARD_MODEL_NAME

        return

    def add_invvar(self):
        self.inv_var = {}
        for kmer, stdev in self.sds.items():
            self.inv_var[kmer] = 1 / (stdev * stdev)

        return

    def __init__(self, ref_fn, is_text_model=False, kmer_ref=None,
                 central_pos=None, minimal_startup=False):
        if ref_fn is None:
            assert kmer_ref is not None and central_pos is not None
            self._load_std_model(kmer_ref, central_pos)
        else:
            self.ref_fn = th.resolve_path(ref_fn)
            if is_text_model:
                self._parse_text_model()
            else:
                self._parse_tombo_model()

        self.kmer_width = len(next(k for k in self.means))
        self.is_std_model = (self.name == STANDARD_MODEL_NAME and
                             self.alt_base is None)
        self.is_alt_model = not self.is_std_model

        if not minimal_startup:
            self.add_invvar()


############################
##### Model Estimation #####
############################

def parse_tombo_models(alt_fns, std_ref):
    """
    Parse several alternative tombo model files
    """
    alt_refs = {}
    for alt_model_fn in alt_fns:
        alt_ref = TomboModel(alt_model_fn)
        if (std_ref.central_pos != alt_ref.central_pos or
            std_ref.kmer_width != alt_ref.kmer_width):
            th._warning_message(
                'Standard and ' + alt_model_fn + ' alternative base ' +
                'models must be estimated using the same k-mer positions.')
            continue
        if not alt_ref.is_alt_model:
            th._warning_message(
                'Alternative model ' + alt_model_fn + ' appears to be a ' +
                'standard model and will not be processed.')
            continue

        if alt_ref.name in alt_refs:
            th._warning_message(
                alt_ref.name + ' alternative model found in more than one ' +
                'model file. Ignoring: ' + alt_model_fn)
            continue
        alt_refs[alt_ref.name] = alt_ref

    return alt_refs

def get_default_standard_ref(raw_read_coverage, bio_samp_type=None):
    if bio_samp_type is not None:
        standard_ref_fn = STANDARD_MODELS[bio_samp_type]
    elif th.is_rna(raw_read_coverage):
        if VERBOSE: th._status_message(
                'Using default canonical ***** RNA ***** model.')
        standard_ref_fn = STANDARD_MODELS['RNA']
    else:
        if VERBOSE: th._status_message(
                'Using default canonical ***** DNA ***** model.')
        standard_ref_fn = STANDARD_MODELS['DNA']
    # get full filename path with setuptools
    standard_ref_fn = pkg_resources.resource_filename(
        'tombo', 'tombo_models/' + standard_ref_fn)

    return standard_ref_fn, bio_samp_type

def get_default_standard_ref_from_files(fast5_fns, bio_samp_type=None):
    if bio_samp_type is not None:
        standard_ref_fn = STANDARD_MODELS[bio_samp_type]
    elif th.is_rna_from_files(fast5_fns):
        if VERBOSE: th._status_message(
                'Using default canonical ***** RNA ***** model.')
        standard_ref_fn = STANDARD_MODELS['RNA']
        bio_samp_type = 'RNA'
    else:
        if VERBOSE: th._status_message(
                'Using default canonical ***** DNA ***** model.')
        standard_ref_fn = STANDARD_MODELS['DNA']
        bio_samp_type = 'DNA'
    # get full filename path with setuptools
    standard_ref_fn = pkg_resources.resource_filename(
        'tombo', 'tombo_models/' + standard_ref_fn)

    return standard_ref_fn, bio_samp_type

def _print_alt_models():
    alt_model_types = [tuple(mod_name.split(ALT_MODEL_SEP_CHAR))
                       for mod_name in ALTERNATE_MODELS.keys()]
    alt_bio_samps = ['',] + sorted(set(list(zip(*alt_model_types))[0]))
    alt_mods = list(set(list(zip(*alt_model_types))[1]))
    row_format ="{:<10}" * (len(alt_bio_samps)) + '\n'
    sys.stderr.write(row_format.format(*alt_bio_samps))
    for alt_mod in alt_mods:
        has_mod = [alt_mod,]
        for bio_samp in alt_bio_samps[1:]:
            has_mod.append(' X' if (bio_samp, alt_mod) in alt_model_types else '')
        sys.stderr.write(row_format.format(*has_mod))

    return

def get_default_alt_ref(alt_name, raw_read_coverage, bio_samp_type=None):
    if bio_samp_type is not None:
        try:
            alt_model_fn = ALTERNATE_MODELS[
                bio_samp_type + ALT_MODEL_SEP_CHAR + alt_name]
        except KeyError:
            alt_model_fn = None
    elif th.is_rna(raw_read_coverage):
        bio_samp_type = 'RNA'
        try:
            alt_model_fn = ALTERNATE_MODELS['RNA' + ALT_MODEL_SEP_CHAR + alt_name]
        except KeyError:
            alt_model_fn = None
    else:
        bio_samp_type = 'DNA'
        try:
            alt_model_fn = ALTERNATE_MODELS['DNA' + ALT_MODEL_SEP_CHAR + alt_name]
        except KeyError:
            alt_model_fn = None
    if alt_model_fn is not None:
        # get full filename path with setuptools
        alt_model_fn = pkg_resources.resource_filename(
            'tombo', 'tombo_models/' + alt_model_fn)
    if alt_model_fn is None or not os.path.isfile(alt_model_fn):
        th._warning_message(
            'Tombo default model for ' + alt_name + ' in ' +
            bio_samp_type + ' does not exists.')
        return None, None

    return alt_model_fn, bio_samp_type

def load_alt_refs(alt_names, raw_read_coverage, std_ref, bio_samp_type=None):
    """
    Load several default alternative tombo models
    """
    alt_fns = []
    for alt_name in alt_names:
        alt_model_fn, _ = get_default_alt_ref(
            alt_name, raw_read_coverage, bio_samp_type)
        if alt_model_fn is None:
            continue
        alt_fns.append(alt_model_fn)

    return parse_tombo_models(alt_fns, std_ref)

def get_ref_from_seq(seq, std_ref, rev_strand=False, alt_ref=None):
    seq_kmers = [seq[i:i + std_ref.kmer_width]
                 for i in range(len(seq) - std_ref.kmer_width + 1)]
    # get stat lookups from seq on native strand then flip if rev_strand
    if rev_strand:
        seq_kmers = seq_kmers[::-1]

    try:
        ref_means = np.array([std_ref.means[kmer] for kmer in seq_kmers])
        ref_sds = np.array([std_ref.sds[kmer] for kmer in seq_kmers])
    except KeyError:
        th._error_message_and_exit(
            'Invalid sequence encountered from genome sequence.')
    if alt_ref is None:
        alt_means, alt_sds = None, None
    else:
        alt_means = np.array([alt_ref.means[kmer] for kmer in seq_kmers])
        alt_sds = np.array([alt_ref.sds[kmer] for kmer in seq_kmers])

    return ref_means, ref_sds, alt_means, alt_sds

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
    base_events = th.get_reads_events(reg_reads)
    if len(base_events) == 0:
        return

    # get intervals within the region where coverage is high enough
    # for model estimation
    reg_cov = np.array([len(base_events[pos]) if pos in base_events else 0
                        for pos in range(reg_start, reg_start + region_size)])
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
            reg_start + int_start - bb, reg_start + int_end + ab, reg_reads)
        if strand == '-':
            int_seq = th.comp_seq(int_seq)
        int_len = int_end - int_start
        for pos in range(int_len):
            pos_kmer = int_seq[pos:pos+kmer_width]
            if strand == '-':
                pos_kmer = pos_kmer[::-1]
            try:
                if est_mean:
                    reg_kmer_levels[pos_kmer].append(c_mean_std(
                        base_events[reg_start+pos+int_start]))
                else:
                    reg_kmer_levels[pos_kmer].append(calc_med_sd(
                        base_events[reg_start+pos+int_start]))
            except KeyError:
                continue

    return reg_kmer_levels

def _est_kmer_model_worker(
        region_q, kmer_level_q, progress_q, raw_read_coverage, cov_thresh,
        upstrm_bases, dnstrm_bases, cs_cov_thresh, est_mean, region_size):
    while not region_q.empty():
        try:
            chrm, strand, reg_start = region_q.get(block=False)
        except queue.Empty:
            # sometimes throws false empty error with get(block=False)
            if not region_q.empty():
                continue
            break

        reg_reads = [r_data for r_data in raw_read_coverage[(chrm, strand)]
                     if not (r_data.start >= reg_start + region_size or
                             r_data.end <= reg_start)]
        if len(reg_reads) == 0:
            progress_q.put(1)
            continue

        reg_kmer_levels = get_region_kmer_levels(
            reg_reads, cov_thresh, upstrm_bases, dnstrm_bases,
            cs_cov_thresh, est_mean, region_size, reg_start, strand)
        if reg_kmer_levels is not None:
            kmer_level_q.put(reg_kmer_levels)
        progress_q.put(1)

    return

if _PROFILE_EST_REF:
    _est_kmer_model_wrapper = _est_kmer_model_worker
    def _est_kmer_model_worker(*args):
        import cProfile
        cProfile.runctx('_est_kmer_model_wrapper(*args)', globals(), locals(),
                        filename='est_kmer_model.prof')
        return

def extract_kmer_levels(
        raw_read_coverage, region_size, cov_thresh, upstrm_bases, dnstrm_bases,
        cs_cov_thresh, est_mean, num_processes):
    chrm_sizes = th.get_chrm_sizes(raw_read_coverage)

    region_q = Queue()
    kmer_level_q = Queue()
    progress_q = Queue()
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

    est_args = (
        region_q, kmer_level_q, progress_q, raw_read_coverage, cov_thresh,
        upstrm_bases, dnstrm_bases, cs_cov_thresh, est_mean, region_size)
    est_ps = []
    for p_id in range(num_processes):
        p = Process(target=_est_kmer_model_worker, args=est_args)
        p.start()
        est_ps.append(p)

    if VERBOSE:
        th._status_message('Extracting average k-mer levels.')
        bar = tqdm(total=num_regions, smoothing=0)
    all_reg_kmer_levels = []
    while any(p.is_alive() for p in est_ps):
        try:
            reg_kmer_levels = kmer_level_q.get(block=False)
            all_reg_kmer_levels.append(reg_kmer_levels)
        except queue.Empty:
            try:
                iter_proc = progress_q.get(block=False)
                if VERBOSE: bar.update(iter_proc)
            except queue.Empty:
                sleep(1)
                continue
    while not kmer_level_q.empty():
        reg_kmer_levels = kmer_level_q.get(block=False)
        all_reg_kmer_levels.append(reg_kmer_levels)
    if VERBOSE: bar.close()

    if len(all_reg_kmer_levels) == 0:
        th._error_message_and_exit(
            'No genomic positions contain --minimum-test-reads. Consider ' +
            'setting this option to a lower value.')

    return all_reg_kmer_levels

def tabulate_kmer_levels(kmer_width, all_reg_kmer_levels, min_kmer_obs):
    if VERBOSE: th._status_message('Tabulating k-mer model statistics.')
    all_kmer_mean_sds = []
    if _DEBUG_EST_STD:
        kmer_dens = []
        save_x = np.linspace(KERNEL_DENSITY_RANGE[0], KERNEL_DENSITY_RANGE[1],
                             _DEBUG_EST_NUM_KMER_SAVE)
    for kmer in product(DNA_BASES, repeat=kmer_width):
        kmer = ''.join(kmer)
        try:
            kmer_levels = np.concatenate([
                reg_kmer_levels[kmer] for reg_kmer_levels in all_reg_kmer_levels
                if len(reg_kmer_levels[kmer]) > 0])
        except ValueError:
            th._error_message_and_exit(
                'At least one k-mer is not covered at any poitions by ' +
                '--minimum-test-reads.\n\t\tConsider fitting to a smaller ' +
                'k-mer via the --upstream-bases and --downstream-bases, ' +
                'or lowering --minimum-test-reads.\n\t\tNote that this may ' +
                'result in a lower quality model.')
        if kmer_levels.shape[0] < min_kmer_obs:
            min_obs = min(
                sum(len(reg_levs[''.join(kmer)])
                    for reg_levs in all_reg_kmer_levels)
                for kmer in product(DNA_BASES, repeat=kmer_width))
            th._error_message_and_exit(
                'K-mers represeneted in fewer observations than ' +
                'requested in the provided reads. Consider a shorter ' +
                'k-mer or providing more reads.\n\t' + unicode(min_obs) +
                ' observations found in least common kmer.')
        all_kmer_mean_sds.append((kmer, np.median(kmer_levels[:,0]),
                                  np.median(kmer_levels[:,1])))
        if _DEBUG_EST_STD:
            kmer_kde = stats.gaussian_kde(
                kmer_levels[:,0],
                bw_method=_DEBUG_EST_BW / kmer_levels[:,0].std(ddof=1))
            with np.errstate(under='ignore'):
                kmer_dens.append((kmer, kmer_kde.evaluate(save_x)))

    if _DEBUG_EST_STD:
        with io.open('debug_est_standard_ref.density.txt', 'wt') as fp:
            fp.write('Kmer\tSignal\tDensity\n')
            fp.write('\n'.join('\t'.join(map(str, (kmer, x, y)))
                               for kmer, dens_i in kmer_dens
                               for x, y in zip(save_x, dens_i)) + '\n')
    return all_kmer_mean_sds

def center_model_to_median_norm(
        raw_read_coverage, init_ref, max_reads=NUM_READS_TO_ADJUST_MODEL):
    upstrm_bases = init_ref.central_pos
    dnstrm_bases = init_ref.kmer_width - init_ref.central_pos - 1
    def get_read_corr_factors(r_data):
        with h5py.File(r_data.fn, 'r+') as fast5_data:
            all_raw_signal = th.get_raw_read_slot(fast5_data)['Signal'].value
            event_starts, r_seq = th.get_multiple_slots_read_centric(
                fast5_data, ('start', 'base'), r_data.corr_group)

        if r_data.rna:
            all_raw_signal = all_raw_signal[::-1]
        norm_signal, scale_values = normalize_raw_signal(
            all_raw_signal, 0, all_raw_signal.shape[0], 'median', None, None)

        event_starts = event_starts.astype(np.int64)
        rsrtr = r_data.read_start_rel_to_raw + event_starts[upstrm_bases]
        # since the last segment end wasn't extracted the events are already
        # clipped by one, so deal with event boundary clipping here
        if dnstrm_bases > 1:
            event_starts = event_starts[upstrm_bases:-(dnstrm_bases - 1)]
        else:
            assert dnstrm_bases == 1, (
                'Must have at least one upstream and downstream base for ' +
                'a Tombo model.')
            event_starts = event_starts[upstrm_bases:]
        # reset event starts to 0
        event_starts -= event_starts[0]

        norm_signal = norm_signal[rsrtr:rsrtr + event_starts[-1]]

        r_seq = b''.join(r_seq).decode()
        r_ref_means = get_ref_from_seq(r_seq, init_ref)[0]

        (_, _, shift_corr_factor,
         scale_corr_factor) = calc_kmer_fitted_shift_scale(
             scale_values.shift, scale_values.scale,
             c_new_means(norm_signal, event_starts), r_ref_means,
             method='theil_sen')

        return shift_corr_factor, scale_corr_factor


    all_shift_corr_factors, all_scale_corr_factors = [], []
    all_reads = [r_data for cs_reads in raw_read_coverage.values()
                 for r_data in cs_reads]
    random.shuffle(all_reads)
    not_enough_reads = True
    for r_data in all_reads:
        try:
            r_shift_corr_factor, r_scale_corr_factor = get_read_corr_factors(
                r_data)
            all_shift_corr_factors.append(r_shift_corr_factor)
            all_scale_corr_factors.append(r_scale_corr_factor)
            if len(all_scale_corr_factors) >= max_reads:
                not_enough_reads = False
                break
        except:
            continue

    if not_enough_reads:
        if len(all_shift_corr_factors) == 0:
            th._error_message_and_exit(
                'No reads succcessfully processed for sequence-based ' +
                'normalization parameter re-fitting.')
        th._warning_message(
            'Fewer reads succcessfully processed for sequence-based ' +
            'normalization parameter re-fitting than requested.')

    # compute median shift and scale correction factors
    # scale parameter should be taken in log space, but median performs
    # the same computation
    med_shift_corr_factor = np.median(all_shift_corr_factors)
    med_scale_corr_factor = np.median(all_scale_corr_factors)

    th._status_message('Shift and scale adjustments to match model to ' +
                       'median normalization: ' + str(med_shift_corr_factor) +
                       " " + str(med_scale_corr_factor))
    init_ref.center_model(med_shift_corr_factor, med_scale_corr_factor)

    return init_ref

def estimate_kmer_model(
        f5_dirs, corrected_group, basecall_subgroups,
        kmer_ref_fn, cov_thresh, upstrm_bases, dnstrm_bases, min_kmer_obs,
        kmer_specific_sd, cs_cov_thresh, est_mean, region_size, num_processes):
    """
    Estimate a standard tombo k-mer model
    """
    raw_read_coverage = th.parse_fast5s(
        f5_dirs, corrected_group, basecall_subgroups)
    all_reg_kmer_levels = extract_kmer_levels(
        raw_read_coverage, region_size, cov_thresh, upstrm_bases, dnstrm_bases,
        cs_cov_thresh, est_mean, num_processes)

    all_kmer_mean_sds = tabulate_kmer_levels(
        upstrm_bases + dnstrm_bases + 1, all_reg_kmer_levels, min_kmer_obs)

    # adjust model to match median normalization best via Theil-Sen optimizer fit
    # this will increase the accuracy of median normalized re-squiggle results
    # and should reduce the need for (or number of) iterative re-squiggle runs
    init_ref = TomboModel(
        ref_fn=None, kmer_ref=all_kmer_mean_sds, central_pos=upstrm_bases)

    centered_ref = center_model_to_median_norm(raw_read_coverage, init_ref)

    if not kmer_specific_sd:
        centered_ref.make_constant_sd()
    centered_ref.write_model(kmer_ref_fn)

    return


########################################
##### Alternative Model Estimation #####
########################################

def _parse_base_levels_worker(
        reads_q, kmer_level_q, kmer_width, central_pos, completed_kmers):
    dnstrm_bases = kmer_width - central_pos - 1
    proc_kmer_levels = dict(
        (''.join(kmer), [])
        for kmer in product(DNA_BASES, repeat=kmer_width))
    while not reads_q.empty():
        try:
            r_fn, corr_slot = reads_q.get(block=False)
        except queue.Empty:
            # sometimes throws false empty error with get(block=False)
            if not reads_q.empty():
                continue
            break

        with h5py.File(r_fn, 'r') as fast5_data:
            r_means, r_seq = th.get_multiple_slots_read_centric(
                fast5_data, ['norm_mean', 'base'], corr_slot)
        if r_means is None: continue
        r_seq = b''.join(r_seq).decode()
        for kmer, level in zip(
                (r_seq[i:i + kmer_width]
                 for i in range(len(r_seq) - kmer_width + 1)),
                r_means[central_pos:-dnstrm_bases]):
            if kmer in completed_kmers: continue
            proc_kmer_levels[kmer].append(level)

    kmer_level_q.put(proc_kmer_levels)

    return

def get_batch_kmer_levels(
        reads_q, kmer_level_q, all_reads, parse_levels_batch_size,
        std_ref, completed_kmers, num_processes):
    no_more_reads = False
    try:
        for _ in range(parse_levels_batch_size):
            r_data = next(all_reads)
            reads_q.put((r_data.fn, r_data.corr_group))
    except StopIteration:
        no_more_reads = True

    base_lev_args = (reads_q, kmer_level_q, std_ref.kmer_width,
                     std_ref.central_pos, completed_kmers)
    base_lev_ps = []
    for p_id in range(num_processes):
        p = Process(target=_parse_base_levels_worker, args=base_lev_args)
        p.start()
        base_lev_ps.append(p)

    batch_kmer_levels = []
    while any(p.is_alive() for p in base_lev_ps):
        try:
            proc_kmer_levels = kmer_level_q.get(block=False)
            batch_kmer_levels.append(proc_kmer_levels)
        except queue.Empty:
            sleep(1)
            continue
    while not kmer_level_q.empty():
        proc_kmer_levels = kmer_level_q.get(block=False)
        batch_kmer_levels.append(proc_kmer_levels)

    return batch_kmer_levels, no_more_reads

def parse_base_levels(
        all_reads, std_ref, parse_levels_batch_size, kmer_obs_thresh,
        max_kmer_obs, min_kmer_obs_to_est, num_processes):
    """
    Parse base levels and store grouped by k-mer
    """
    reads_q = Queue()
    kmer_level_q = Queue()

    all_kmer_levels = dict(
        (''.join(kmer), [])
        for kmer in product(DNA_BASES, repeat=std_ref.kmer_width))
    # store set of k-mers with enough observations to save on memory footprint
    # while filling more rare k-mers
    completed_kmers = set()
    if VERBOSE:
        all_reads_bar = tqdm(total=len(all_reads), smoothing=0,
                             desc='Number of total reads used', leave=True)
        min_kmer_bar = tqdm(total=kmer_obs_thresh, smoothing=0,
                            desc='K-mer with fewest observations', leave=True)
        curr_min_kmer_count = 0
    all_reads = iter(all_reads)
    n_batches = 0
    while True:
        batch_kmer_levels, no_more_reads = get_batch_kmer_levels(
            reads_q, kmer_level_q, all_reads, parse_levels_batch_size,
            std_ref, completed_kmers, num_processes)

        # only add observations for k-mers that have not been seen enough times
        # save memory for abundent k-mers
        batch_total_kmers = []
        for kmer in set(all_kmer_levels).difference(completed_kmers):
            all_kmer_levels[kmer].extend((
                kmer_level for proc_kmer_levels in batch_kmer_levels
                for kmer_level in proc_kmer_levels[kmer]))
            batch_total_kmers.append(len(all_kmer_levels[kmer]))
            if batch_total_kmers[-1] > max_kmer_obs:
                completed_kmers.add(kmer)
        if VERBOSE:
            if no_more_reads:
                all_reads_bar.update(all_reads_bar.total - all_reads_bar.n)
            else:
                all_reads_bar.update(parse_levels_batch_size)
            new_min_kmer_count = min(batch_total_kmers)
            min_kmer_bar.update(new_min_kmer_count - curr_min_kmer_count)
            curr_min_kmer_count = new_min_kmer_count
            sleep(0.1)
        else:
            curr_min_kmer_count = min(batch_total_kmers)
        if curr_min_kmer_count > kmer_obs_thresh or no_more_reads:
            break

        n_batches += 1

    if VERBOSE:
        all_reads_bar.close()
        min_kmer_bar.close()

    fewest_kmer_obs = min(len(kmer_levels) for kmer_levels in
                          all_kmer_levels.values())
    if fewest_kmer_obs < kmer_obs_thresh:
        if fewest_kmer_obs < min_kmer_obs_to_est:
            th._error_message_and_exit(
                'Too few minimal k-mer observations to continue to ' +
                'alternative estimation. Minimal k-mer has ' +
                unicode(fewest_kmer_obs) + ' total observations and ' +
                unicode(min_kmer_obs_to_est) +
                ' observations per k-mer are required.')
        th._warning_message(
            'Requested minimal k-mer observations not found in all reads. ' +
            'Continuing to estimation using a k-mer with ' +
            unicode(fewest_kmer_obs) + ' total observations')

    return all_kmer_levels

def write_kmer_densities_file(dens_fn, kmer_dens, save_x):
    with io.open(dens_fn, 'wt') as fp:
        fp.write('Kmer\tSignal\tDensity\n')
        fp.write('\n'.join('\t'.join(map(str, (kmer, x, y)))
                           for kmer, dens_i in kmer_dens.items()
                           for x, y in zip(save_x, dens_i)) + '\n')

    return

def parse_kmer_densities_file(dens_fn):
    kmer_dens_raw = defaultdict(list)
    with io.open(dens_fn) as dens_fp:
        # read in header
        dens_fp.readline()
        for line in dens_fp:
            kmer, _, dens_i = line.split()
            kmer_dens_raw[kmer].append(float(dens_i))

    kmer_dens = {}
    first_len = None
    for kmer, dens_i in kmer_dens_raw.items():
        if first_len is None: first_len = len(dens_i)
        if len(dens_i) != first_len:
            th._error_message_and_exit('Density file is valid.')
        kmer_dens[kmer] = np.array(dens_i)

    return kmer_dens

def est_kernel_density(
        raw_read_coverage, std_ref, kmer_obs_thresh,
        density_basename, save_x, kernel_dens_bw, num_processes,
        alt_or_stnd_name='alt', parse_levels_batch_size=ALT_EST_BATCH,
        max_kmer_obs=MAX_KMER_OBS, min_kmer_obs_to_est=MIN_KMER_OBS_TO_EST):
    all_reads = [r_data for cs_reads in raw_read_coverage.values()
                 for r_data in cs_reads]
    np.random.shuffle(all_reads)
    base_levels = parse_base_levels(
        all_reads, std_ref, parse_levels_batch_size, kmer_obs_thresh,
        max_kmer_obs, min_kmer_obs_to_est, num_processes)

    if VERBOSE: th._status_message('Fitting kernel densities for k-mer levels.')
    kmer_dens = {}
    for kmer, norm_levels in base_levels.items():
        norm_levels = np.array(norm_levels)
        kmer_kde = stats.gaussian_kde(
            norm_levels, bw_method=kernel_dens_bw / norm_levels.std(ddof=1))
        with np.errstate(under='ignore'):
            kmer_dens[kmer] = kmer_kde.evaluate(save_x)

    if density_basename is not None:
        write_kmer_densities_file(
            density_basename + '.' + alt_or_stnd_name + '_density.txt',
            kmer_dens, save_x)

    return kmer_dens

def estimate_kmer_densities(
        f5_dirs, control_dirs, corrected_group, basecall_subgroups,
        standard_ref_fn, bio_samp_type, kmer_obs_thresh, density_basename,
        kernel_dens_bw, save_x, num_processes):
    raw_read_coverage = th.parse_fast5s(
        f5_dirs, corrected_group, basecall_subgroups)
    cntrl_read_coverage = th.parse_fast5s(
        control_dirs, corrected_group, basecall_subgroups)

    if VERBOSE: th._status_message('Parsing standard model file.')
    if standard_ref_fn is None:
        standard_ref_fn, bio_samp_type = get_default_standard_ref(
            raw_read_coverage, bio_samp_type)
    std_ref = TomboModel(standard_ref_fn)

    if VERBOSE: th._status_message('Parsing base levels from alternative reads.')
    alt_dens = est_kernel_density(
        raw_read_coverage, std_ref, kmer_obs_thresh, density_basename,
        save_x, kernel_dens_bw, num_processes, 'alternate')
    if VERBOSE: th._status_message('Parsing base levels from standard reads.')
    std_dens = est_kernel_density(
        cntrl_read_coverage, std_ref, kmer_obs_thresh, density_basename,
        save_x, kernel_dens_bw, num_processes, 'control')

    return alt_dens, std_dens, std_ref

def load_kmer_densities(
        alt_dens_fn, std_dens_fn, f5_dirs, corrected_group,
        basecall_subgroups, std_ref_fn, bio_samp_type):
    if VERBOSE: th._status_message('Parsing standard model file.')
    if std_ref_fn is None:
        if f5_dirs is None and bio_samp_type is None:
            th._error_message_and_exit(
                'Must provide a FAST5s directory, a canonical model ' +
                'file or spcify the biological sample type.')
        raw_read_coverage = None
        if f5_dirs is not None:
            raw_read_coverage = th.parse_fast5s(
                f5_dirs, corrected_group, basecall_subgroups)
        std_ref_fn, bio_samp_type = get_default_standard_ref(
            raw_read_coverage, bio_samp_type)
    std_ref = TomboModel(std_ref_fn)

    if VERBOSE: th._status_message('Parsing density files.')
    alt_dens = parse_kmer_densities_file(alt_dens_fn)
    std_dens = parse_kmer_densities_file(std_dens_fn)
    num_dens_points = next(v for v in alt_dens.values()).shape[0]
    if num_dens_points != next(v for v in std_dens.values()).shape[0]:
        th._error_message_and_exit(
            'Alternative and standard density ' +
            'estimates do not correspond.')

    save_x = np.linspace(KERNEL_DENSITY_RANGE[0], KERNEL_DENSITY_RANGE[1],
                         num_dens_points)

    return alt_dens, std_dens, std_ref, save_x

def isolate_alt_density(
        alt_dens, std_dens, alt_base, alt_frac_pctl, std_ref, save_x):
    def calc_mean(dens):
        return np.average(save_x[dens>1e-10], weights=dens[dens>1e-10])


    # estimate density shift from k-mers without the alternate base
    no_alt_std_means, no_alt_mean_diffs = [], []
    for kmer in std_dens:
        if alt_base in kmer: continue
        no_alt_std_means.append(calc_mean(std_dens[kmer]))
        kmer_alt_mean = calc_mean(alt_dens[kmer])
        no_alt_mean_diffs.append(kmer_alt_mean - no_alt_std_means[-1])
    calc_offset = np.poly1d(np.polyfit(no_alt_std_means, no_alt_mean_diffs, 2))
    save_x_unit = save_x[1] - save_x[0]

    shifted_alt_dens = {}
    for kmer, kmer_alt_dens in alt_dens.items():
        est_offset = int(calc_offset(calc_mean(std_dens[kmer])) / save_x_unit)
        # if std density mean is estimated to be greater
        if est_offset < 0:
            # shift alt dens to the right
            shifted_alt_dens[kmer] = np.concatenate([
                [0.0,] * -est_offset, kmer_alt_dens[:est_offset]])
        else:
            # else shift alt dens to the left
            shifted_alt_dens[kmer] = np.concatenate([
                kmer_alt_dens[est_offset:], [0.0,] * est_offset])

    def get_peak_frac(kmer_std_dens, kmer_alt_dens):
        # find highest peak in standard density
        std_peak = np.argmax(kmer_std_dens)

        # find closest peak in alternative density
        alt_local_peaks = np.where(np.concatenate([
            [False,], np.logical_and(
                kmer_alt_dens[1:-1] > kmer_alt_dens[:-2],
                kmer_alt_dens[1:-1] > kmer_alt_dens[2:]) + [False,]]))[0]
        matched_alt_peak = alt_local_peaks[
            np.argmin(abs(alt_local_peaks - std_peak))]
        return kmer_alt_dens[matched_alt_peak] / kmer_std_dens[std_peak]


    # estimate the alternative base incorporation rate
    std_frac = np.percentile([
        get_peak_frac(std_dens[kmer], shifted_alt_dens[kmer])
        for kmer in std_dens if kmer.count(alt_base) == 1], alt_frac_pctl)
    if VERBOSE: th._status_message(
            'Alternative base incorporation rate estimate: ' +
            unicode(1 - std_frac))
    if std_frac >= 1:
        th._warning_message(
            'Alternative base incorporation rate ' +
            'estimate is approximately 0. Consider lowering ' +
            '--alt-fraction-percentile.')

    # get mean model SD. most models will be constant, but use mean in case
    model_sd = np.mean(list(std_ref.sds.values()))
    # subtract off fraction of standard density from alt density
    # to estimate mean of isolated alternative distribution
    alt_ref = []
    for kmer, std_level in std_ref.means.items():
        if kmer.count(alt_base) == 0:
            alt_ref.append((kmer, std_level, model_sd))
            continue
        # assuming random incorporation the prortion of standard base
        # observations at this k-mer is the standard fraction raised
        # to the number of alt_base occurences in the k-mer
        kmer_std_frac = std_frac**kmer.count(alt_base)
        with np.errstate(under='ignore'):
            diff_dens = shifted_alt_dens[kmer] - (
                std_dens[kmer] * kmer_std_frac)
            diff_dens[diff_dens < 0] = 0
            alt_level = np.average(save_x, weights=diff_dens)
        alt_ref.append((kmer, alt_level, model_sd))

    alt_ref = TomboModel(
        ref_fn=None, kmer_ref=alt_ref, central_pos=std_ref.central_pos,
        minimal_startup=True)

    return alt_ref

def estimate_alt_model(
        f5_dirs, control_dirs, corrected_group, basecall_subgroups,
        std_ref_fn, bio_samp_type, alt_base, alt_frac_pctl,
        kmer_obs_thresh, density_basename, kernel_dens_bw, alt_dens_fn,
        std_dens_fn, num_processes, num_dens_points=NUM_DENS_POINTS):
    """
    Estimate an alternative model from a sample with a single,
    known, randomly-incorporated alternative base
    """
    if alt_dens_fn is None or std_dens_fn is None:
        save_x = np.linspace(KERNEL_DENSITY_RANGE[0], KERNEL_DENSITY_RANGE[1],
                             num_dens_points)
        alt_dens, std_dens, std_ref = estimate_kmer_densities(
            f5_dirs, control_dirs, corrected_group, basecall_subgroups,
            std_ref_fn, bio_samp_type, kmer_obs_thresh, density_basename,
            kernel_dens_bw, save_x, num_processes)
    else:
        alt_dens, std_dens, std_ref, save_x = load_kmer_densities(
            alt_dens_fn, std_dens_fn, f5_dirs, corrected_group,
            basecall_subgroups, std_ref_fn, bio_samp_type)

    if VERBOSE: th._status_message('Isolating alternative base distribtuions.')
    # perform alternative density isolation algorithm
    alt_ref = isolate_alt_density(
        alt_dens, std_dens, alt_base, alt_frac_pctl, std_ref, save_x)

    return alt_ref

if _PROFILE_ALT_EST:
    _est_alt_wrapper = estimate_alt_model
    def estimate_alt_model(*args):
        import cProfile
        cProfile.runctx('_est_alt_wrapper(*args)', globals(), locals(),
                        filename='est_alt_model.prof')
        return None


####################################
##### Core Statistical Testing #####
####################################

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
    ecdffcator = np.arange(1, nobs + 1) / nobs
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
        raise NotImplementedError(
            "P-values vector too short for Fisher's Method " +
            "window compuation.")
    with np.errstate(invalid='ignore'):
        pvals = np.maximum(pvals, SMALLEST_PVAL)
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

def compute_accuracy_rates(stat_has_mod, num_plot_points=ROC_PLOT_POINTS):
    """
    Given a list or numpy array of true/false values, function returns
    num_plot_point evenly spaced values along the true positive, false
    positive and presicion arrays
    """
    tp_cumsum = np.cumsum(stat_has_mod)
    tp_rate = tp_cumsum / tp_cumsum[-1]
    fp_cumsum = np.cumsum(np.logical_not(stat_has_mod))
    fp_rate = fp_cumsum / fp_cumsum[-1]

    precision = tp_cumsum / np.arange(1, len(stat_has_mod) + 1, dtype=float)

    # trim to number of requested points
    tp_rate = tp_rate[np.linspace(0, tp_rate.shape[0] - 1,
                                  num_plot_points).astype(np.int64)]
    fp_rate = fp_rate[np.linspace(0, fp_rate.shape[0] - 1,
                                  num_plot_points).astype(np.int64)]
    precision = precision[np.linspace(0, precision.shape[0] - 1,
                                      num_plot_points + 1).astype(np.int64)][1:]

    return tp_rate, fp_rate, precision

def get_motif_stats(motif, stats, genome_index):
    stat_has_mod = []
    for stat_seq in stats.iter_stat_seqs(
            genome_index, motif.mod_pos - 1,
            motif.motif_len - motif.mod_pos):
        if motif.motif_pat.match(stat_seq) is not None:
            stat_has_mod.append(True)
        # don't include sites that aren't at the base of interest
        elif stat_seq[motif.mod_pos - 1] == motif.mod_base:
            stat_has_mod.append(False)

    return compute_accuracy_rates(stat_has_mod)


#########################################
##### Local Model-based Re-squiggle #####
#########################################

def get_dynamic_prog_params(match_evalue):
    """
    Compute dynamic programming shift parameters from an expected match
    expected value
    """
    z_shift = HALF_NORM_EXPECTED_VAL + match_evalue
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
    Identify regions along a read that do not match well with the genomic
    reference tombo model
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
    curr_start, curr_end = next(signif_shift_regs)
    for reg_start, reg_end in signif_shift_regs:
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


##########################
##### Statistics I/O #####
##########################

def write_stats(
        all_reg_stats, stats_bsnm, stat_type, min_test_vals, alt_name=None):
    """
    Write a tombo statistics file
    """
    if VERBOSE: th._status_message(
            'Saving signal shift significance testing results.')
    def convert_reg_stats(reg_stats):
        # get all unique fasta record names to store in HDF5 attributes and
        # encode as integers in the stats numpy table
        chrms_lookup = dict(zip(
            sorted(set(map(itemgetter(2), reg_stats))), count()))

        np_stats = []
        for (reg_frac_standard_base, reg_poss, chrm, strand,
             reg_cov, ctrl_cov, valid_cov) in reg_stats:
            np_stats.append(np.array(
                [pos_stats for pos_stats in zip(
                    reg_frac_standard_base, reg_poss,
                    repeat(chrms_lookup[chrm]), repeat(strand),
                    reg_cov, ctrl_cov, valid_cov)
                 if not np.isnan(pos_stats[0])],
                dtype=[
                    (str('frac'), 'f8'), (str('pos'), 'u4'), (str('chrm'), 'u4'),
                    (str('strand'), 'S1'), (str('cov'), 'u4'),
                    (str('control_cov'), 'u4'), (str('valid_cov'), 'u4')]))

        np_stats = np.concatenate(np_stats)

        np_stats = np_stats[np.greater_equal(
            np_stats['valid_cov'], min_test_vals)]

        return np_stats, chrms_lookup

    def write_stats_data(stats_fp, stats, stat_type, chrms_lookup):
        stats_fp.create_dataset('stats', data=stats, compression="gzip")
        stats_fp.attrs['stat_type'] = stat_type

        chrms_subgrp = stats_fp.create_group('chromosome_values')
        for chrm, chrm_val in chrms_lookup.items():
            chrms_subgrp.attrs[chrm] = chrm_val

        return

    stats_fn = stats_bsnm + '.tombo.stats' if alt_name is None else \
               stats_bsnm + '.' + alt_name + '.tombo.stats'
    all_reg_stats, chrms_lookup = convert_reg_stats(all_reg_stats)
    with h5py.File(stats_fn, 'w') as stats_fp:
        write_stats_data(stats_fp, all_reg_stats, stat_type, chrms_lookup)

    return

class TomboStats(object):
    """
    Parse and retrieve relevant information from a standard (not per-read) Tombo
    statistics file.
    """
    def _parse_stats(self):
        if self.stats_fn is None or not os.path.isfile(self.stats_fn):
            th._error_message_and_exit(
                'Statistics file not provided or provided file does not exist.')

        try:
            with h5py.File(self.stats_fn, 'r') as stats_fp:
                self.stats = stats_fp['stats'].value
                self.stat_type = stats_fp.attrs['stat_type']
                try:
                    self.chrms_lookup = dict(
                        (chrm_val, chrm_name) for chrm_name, chrm_val in
                        stats_fp['chromosome_values'].attrs.items())
                    self.has_chrm_lookup = True
                except:
                    self.has_chrm_lookup = False
                    th._warning_message(
                        'Old version of Tombo used to create statistics ' +
                        'file. Upgrading to the current version suggested ' +
                        'for best results.')
        except:
            th._error_message_and_exit(
                'Attempt to load statistics file failed. May be an old ' +
                'version of statistics file. Try deleting statistics ' +
                'file and re-calculating using current tombo version.')

        return

    def __init__(self, stats_fn):
        """
        Parse a standard Tombo statistics file.
        """
        self.stats_fn = stats_fn
        self.has_damp_frac = False
        self.cov_damp_counts = None
        self.has_stat_dict = False
        self.stat_dict = None
        self._parse_stats()

        return

    def filter_coverage(self, min_reads):
        """
        Filter statistics at locations with less than `min_reads` coverage.
        """
        self.stats = self.stats[np.logical_and(
            self.stats['valid_cov'] >= min_reads,
            np.logical_or(self.stat_type != SAMP_COMP_TXT,
                          self.stats['control_cov'] >= min_reads))]

        return

    def is_empty(self):
        """
        Is the statistics table empty?
        """
        return self.stats.shape[0] == 0

    def calc_damp_fraction(self, cov_damp_counts):
        """
        Compute dampened fraction of unmodified reads using provided
        un-modified and modified pseudo-counts from cov_damp_counts
        """
        self.has_damp_frac = True
        self.cov_damp_counts = cov_damp_counts
        damp_frac = np.empty(self.stats.shape[0])
        damp_frac[:] = np.nan
        non_mod_counts = np.round(self.stats['frac'] * self.stats['valid_cov'])
        # compute dampened fraction of modified reads by adding psuedo-counts
        # to the modified and un-modified counts (equivalent to a beta prior
        # on the fraction estimation as a binomial variable)
        damp_frac = (non_mod_counts + cov_damp_counts[0]) / (
            self.stats['valid_cov'] + sum(cov_damp_counts))
        damp_name = 'damp_frac' if sys.version_info[0] > 2 else b'damp_frac'
        self.stats = append_fields(self.stats, damp_name, damp_frac)

        return

    def order_by_frac(self, cov_damp_counts=None):
        """
        Order statistics table via fraction of unmodified reads

        If cov_damp_counts is provided or has been previously provided
        fractions will be dampened accordingly
        """
        if cov_damp_counts is None and not self.has_damp_frac:
            self.stats.sort(order=str('frac'))
        else:
            self.calc_damp_fraction(cov_damp_counts)
            self.stats.sort(order=str('damp_frac'))

        return

    def order_by_pos(self):
        """
        Order statistics table by chrmomosome, then strand and then position
        """
        self.stats.sort(order=['chrm', 'strand', 'pos'])
        return

    def _get_chrm_name(self, pos_stat):
        if self.has_chrm_lookup:
            return self.chrms_lookup[pos_stat['chrm']]
        return pos_stat['chrm'].decode()

    def iter_stat_seqs(self, genome_index, before_bases, after_bases,
                       include_pos=False):
        """
        Iterate through statistics table in current order returning the genomic
        sequence surrounding each position.

        `include_position` option will yeild (pos_seq, chrm, strand, start, end)
        for each record.
        """
        for pos_stat in self.stats:
            chrm, strand, pos = (self._get_chrm_name(pos_stat),
                                 pos_stat['strand'].decode(),
                                 pos_stat['pos'])
            if strand == '+':
                start, end = pos - before_bases, pos + after_bases + 1
                if start < 0: continue
                stat_seq = genome_index.get_seq(chrm, start, end)
            else:
                start, end = pos - after_bases, pos + before_bases + 1
                if start < 0: continue
                stat_seq = th.rev_comp(genome_index.get_seq(chrm, start, end))
            if include_pos:
                yield stat_seq, chrm, strand, start, end
            else:
                yield stat_seq

        return

    def iter_fracs(self):
        """
        Iterate through statistics table yeilding
        (chrm, strand, pos, frac, damp_frac).
        """
        for pos_stat in self.stats:
            yield (
                self._get_chrm_name(pos_stat), pos_stat['strand'].decode(),
                pos_stat['pos'], pos_stat['frac'],
                pos_stat['damp_frac'] if self.has_damp_frac else None,
                pos_stat['valid_cov'])

        return

    def get_most_signif_regions(self, num_bases, num_regions, unique_pos=True,
                                cov_damp_counts=None):
        """
        Select regions centered on locations with the largest fraction
        of modified bases
        """
        self.order_by_frac(cov_damp_counts)
        selected_regs = []
        used_intervals = defaultdict(set)
        for i, pos_stat in enumerate(self.stats):
            int_start = max(0, pos_stat['pos'] - int(num_bases / 2.0))
            chrm = self._get_chrm_name(pos_stat)
            strand = pos_stat['strand'].decode()
            if (not unique_pos or
                pos_stat['pos'] not in used_intervals[(chrm, strand)]):
                used_intervals[(chrm, strand)].update(
                    range(int_start, int_start + num_bases))
                int_text = 'Est. Frac. Alternate: {0:.2g}  Coverage: {1}'.format(
                    1 - pos_stat[str('damp_frac')], pos_stat[str('valid_cov')]) \
                    if self.has_damp_frac else \
                       'Frac. Alternate: {0:.2g}  Coverage: {1}'.format(
                           1 - pos_stat[str('frac')], pos_stat[str('valid_cov')])
                selected_regs.append(th.intervalData(
                    '{:03d}'.format(i), chrm, int_start,
                    int_start + num_bases, strand, int_text))
                if len(selected_regs) >= num_regions: break

        if len(selected_regs) == 0:
            th._error_message_and_exit(
                'No locations identified. Most likely an empty statistics file.')
        if len(selected_regs) < num_regions:
            th._warning_message(
                'Fewer unique significant locations more than [--num-bases]/2 ' +
                'apart were identified. Continuing with ' +
                str(len(selected_regs)) + ' unique locations.')

        return selected_regs

    def create_stat_dict(self, dict_batch_size=10000):
        """
        Create random access to fraction modified values by position

        Fraction will be dampened if cov_damp_counts was previously provided
        Access dictionary will be stored in the stat_dict slot
        """
        self.has_stat_dict = True
        self.dict_batch_size = dict_batch_size
        s_stats = np.sort(self.stats, order=['chrm', 'strand', 'pos'])
        self.stat_dict = {}
        # split at chromosome/strand switches
        for cs_stats in np.split(s_stats, np.where(np.logical_or(
                s_stats['strand'][:-1] != s_stats['strand'][1:],
                np.diff(s_stats['chrm']) != 0))[0] + 1):
            for batch_stats in np.split(cs_stats, np.where(np.diff(
                    np.floor_divide(
                        cs_stats['pos'], dict_batch_size)) != 0)[0] + 1):
                batch_fracs = 1 - (
                    batch_stats[str('damp_frac')] if self.has_damp_frac else
                    batch_stats[str('frac')])
                self.stat_dict[(
                    self._get_chrm_name(batch_stats[0]),
                    batch_stats[0]['strand'].decode(),
                    np.floor_divide(batch_stats[0]['pos'],
                                    dict_batch_size))] = (
                                        batch_fracs, batch_stats['pos'])

        return

    def get_pos_frac(self, chrm, strand, pos, missing_value=None):
        """
        Obtain statistic value from the requested genomic position
        """
        # TODO: Add a get_reg_pos function and only get the reg values
        # once. Just need to handle edge of batch cases
        if not self.has_stat_dict:
            self.create_stat_dict()
        try:
            reg_fracs, reg_poss = self.stat_dict[(
                chrm, strand, np.floor_divide(pos, self.dict_batch_size))]
            pos_index = np.where(reg_poss == pos)[0]
            if len(pos_index) != 1: raise KeyError
            pos_frac = reg_fracs[pos_index[0]]
        except KeyError:
            pos_frac = missing_value

        return pos_frac


class PerReadStats(object):
    def __init__(self, per_read_stats_fn, stat_type=None, region_size=None):
        """
        Open per-read statistics file. If stat_type and region_size are provided
        the file is opened for writing, else it is opened for random access.

        WARNING: If stat_type and region_size are provided the current file's
        contents will be deleted.
        """
        if stat_type is None or region_size is None:
            # open file for reading
            try:
                self._fp = h5py.File(per_read_stats_fn, 'r')
                self.stat_type = self._fp.attrs['stat_type']
                self.region_size = self._fp.attrs['block_size']
                self.per_read_blocks = self._fp['Statistic_Blocks']
                blocks_index = defaultdict(list)
                for block_name, block_data in self.per_read_blocks.items():
                    blocks_index[
                        (block_data.attrs['chrm'],
                         block_data.attrs['strand'])].append((
                             block_data.attrs['start'],
                             block_data.attrs['start'] + self.region_size,
                             block_name))
                self.blocks_index = dict(blocks_index)
                self.iter_all_cs = iter(list(self.blocks_index))
                self.iter_curr_cs = next(self.iter_all_cs)
                self.iter_curr_cs_blocks = iter(
                    self.blocks_index[self.iter_curr_cs])
            except:
                th._error_message_and_exit(
                    'Non-existent or invalid per-read statistics file provided.')
        else:
            # set class attributes
            self.stat_type = stat_type
            self.region_size = region_size
            self.curr_block_num = 0

            # try to remove file for overwriting old results
            try:
                os.remove(per_read_stats_fn)
            except:
                pass
            # open file for writing
            self._fp = h5py.File(per_read_stats_fn, 'w')

            # save attributes to file and open stats blocks group
            self._fp.attrs['stat_type'] = stat_type
            self._fp.attrs['block_size'] = region_size
            self.per_read_blocks = self._fp.create_group('Statistic_Blocks')

        self.are_pvals = self.stat_type != ALT_MODEL_TXT

        return

    def write_per_read_block(
            self, per_read_block, read_id_lookup, chrm, strand, start):
        """
        Write region statistics block to file.
        """
        try:
            block_data = self.per_read_blocks.create_group(
                'Block_' + unicode(self.curr_block_num))
            self.curr_block_num += 1
        except:
            th._warning_message(
                'Per-read statistics file not opened for writing.')
            return

        block_data.attrs['chrm'] = chrm
        block_data.attrs['strand'] = strand
        block_data.attrs['start'] = start
        block_data.create_dataset(
            'block_stats', data=per_read_block, compression="gzip")
        # add lookup dict for read_id slot stored in table to save space and
        # avoid memory leak due to vlen slots in h5py datasets
        read_id_grp = block_data.create_group('read_ids')
        for read_id, read_id_val in read_id_lookup.items():
            read_id_grp.attrs[read_id] = read_id_val

        self._fp.flush()

        return

    def get_region_per_read_stats(self, interval_data, num_reads=None):
        """
        Get per-read statistics from the specifed interval for a random selection
        of num_reads.
        """
        try:
            cs_blocks = self.blocks_index[(
                interval_data.chrm, interval_data.strand)]
        except KeyError:
            return

        int_block_stats = []
        for block_data in cs_blocks:
            if (interval_data.end < block_data[0] or
                interval_data.start > block_data[1]): continue
            block_stats = self.per_read_blocks[
                block_data[2]]['block_stats'].value
            reg_poss = block_stats['pos']
            reg_read_stats = block_stats['stat']
            # convert read_ids back into strings
            block_read_id_lookup = dict([
                (read_id_val, read_id) for read_id, read_id_val in
                self.per_read_blocks[block_data[2]][
                    'read_ids'].attrs.items()])
            reg_read_ids = [
                block_read_id_lookup[r_id] for r_id in block_stats['read_id']]
            int_block_stats.append(np.array(
                list(zip(reg_poss, reg_read_stats, reg_read_ids)),
                dtype=[(str('pos'), 'u4'), (str('stat'), 'f8'),
                       (str('read_id'), object)]))

        if len(int_block_stats) == 0:
            return

        if len(int_block_stats) == 1:
            int_block_stats = int_block_stats[0]
        else:
            int_block_stats = np.concatenate(int_block_stats)

        all_int_stats = int_block_stats[
            (int_block_stats['pos'] >= interval_data.start) &
            (int_block_stats['pos'] < interval_data.end)]

        read_ids = set(all_int_stats['read_id'])
        if num_reads is not None and num_reads < len(read_ids):
            int_plot_reads = set(random.sample(read_ids, num_reads))
            all_int_stats = all_int_stats[
                np.array([read_id in int_plot_reads
                          for read_id in all_int_stats['read_id']])]

        return all_int_stats

    def __iter__(self):
        """
        Iterator over all statistics blocks, yeilding chrm, strand,
        start, end, block_stats
        """
        self.iter_all_cs = iter(list(self.blocks_index))
        self.iter_curr_cs = next(self.iter_all_cs)
        self.iter_curr_cs_blocks = iter(
            self.blocks_index[self.iter_curr_cs])
        return self

    def __next__(self):
        try:
            next_start, next_end, next_block_name = next(
                self.iter_curr_cs_blocks)
        except StopIteration:
            # move to next chromosome and strand
            # this will raise a second StopIteration
            # when the end of the blocks is hit
            self.iter_curr_cs = next(self.iter_all_cs)
            self.iter_curr_cs_blocks = iter(
                self.blocks_index[self.iter_curr_cs])
            next_start, next_end, next_block_name = next(
                self.iter_curr_cs_blocks)

        chrm, strand = self.iter_curr_cs
        next_block_stats = self.per_read_blocks[
            next_block_name]['block_stats'].value
        return chrm, strand, next_start, next_end, next_block_stats

    def close(self):
        self._fp.close()
        return


################################
##### Base-by-base Testing #####
################################

def apply_per_read_thresh(
        pr_stats_fn, single_read_thresh, min_test_vals, lower_thresh):
    if VERBOSE: th._status_message(
            'Loading and aggregating per-read statistics.')
    all_reg_stats = []
    pr_stats = PerReadStats(pr_stats_fn)
    for chrm, strand, start, end, block_stats in pr_stats:
        block_stats.sort(order=str('pos'))
        reg_base_stats = np.split(
            block_stats['stat'], np.where(np.concatenate(
                [[0,], np.diff(block_stats['pos'])]) > 0)[0])
        reg_poss = np.unique(block_stats['pos'])

        reg_cov = [base_stats.shape[0] for base_stats in reg_base_stats]
        if lower_thresh is not None:
            # filter base statistics that fall between the upper and lower
            # stat threshold for the log likelihood statistic
            reg_base_stats = [
                base_stats[np.logical_or(
                    base_stats <= lower_thresh,
                    base_stats >= single_read_thresh)]
                for base_stats in reg_base_stats]
            valid_cov = [base_stats.shape[0] for base_stats in reg_base_stats]
        elif pr_stats.stat_type == ALT_MODEL_TXT:
            # filter base statistics that fall between the upper and lower
            # stat threshold for the log likelihood statistic
            reg_base_stats = [
                base_stats[np.abs(base_stats) >= single_read_thresh]
                for base_stats in reg_base_stats]
            valid_cov = [base_stats.shape[0] for base_stats in reg_base_stats]
        else:
            valid_cov = reg_cov

        ctrl_cov = repeat(0)

        reg_frac_standard_base = np.array([
            np.greater_equal(
                base_stats, single_read_thresh).sum() / base_stats.shape[0]
            if base_stats.shape[0] > 0 else np.NAN
            for base_stats in reg_base_stats])

        reg_stats = (reg_frac_standard_base, reg_poss, chrm, strand,
                     reg_cov, ctrl_cov, valid_cov)
        all_reg_stats.append(reg_stats)

    if len(all_reg_stats) == 0:
        th._error_message_and_exit(
            'No genomic positions contain --minimum-test-reads.')

    return all_reg_stats, pr_stats.stat_type

def get_reads_ref(ctrl_reg_reads, reg_start, region_size,
                  min_test_vals, fm_offset):
    """
    Get mean and standard deviation of levels from a sample across the genome
    """
    ctrl_base_events = th.get_reads_events(ctrl_reg_reads)
    if ctrl_base_events is None:
        raise NotImplementedError
    arr_size = region_size + (fm_offset * 2)
    ctrl_means, ctrl_sds = np.empty(arr_size), np.empty(arr_size)
    ctrl_means[:] = np.NAN
    ctrl_sds[:] = np.NAN
    ctrl_cov = {}
    for pos, pos_events in sorted(ctrl_base_events.items()):
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

def compute_sample_compare_read_stats(
        r_data, ctrl_means, ctrl_sds, fm_offset, reg_start, region_size):
    """
    Compute signficance statistics using comparison of two sequenceing samples
    method for a single read within a specified genomic region.
    """
    def comp_clip_and_flip():
        with h5py.File(r_data.fn, 'r') as fast5_data:
            r_means = th.get_single_slot_read_centric(
                fast5_data, 'norm_mean', r_data.corr_group)
            read_id = th.get_raw_read_slot(fast5_data).attrs['read_id']
        if r_means is None:
            raise NotImplementedError(
                'Read does not contain re-squiggled level means.')

        read_start, read_end = r_data.start, r_data.end
        if read_start + fm_offset < reg_start:
            num_start_clip = reg_start - (read_start + fm_offset)
            read_start = reg_start - fm_offset
            if r_data.strand == '+':
                r_means = r_means[num_start_clip:]
            else:
                r_means = r_means[:-num_start_clip]
        if read_end - fm_offset > reg_start + region_size:
            num_end_clip = (read_end - fm_offset) - (reg_start + region_size)
            read_end = reg_start + region_size + fm_offset
            if r_data.strand == '+':
                r_means = r_means[:-num_end_clip]
            else:
                r_means = r_means[num_end_clip:]

        # flip means to match genomic positions
        if r_data.strand == '-':
            r_means = r_means[::-1]

        return r_means, read_start, read_end, read_id

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

    r_means, read_start, read_end, read_id = comp_clip_and_flip()
    r_z_scores = get_read_comp_z_score(r_means, read_start, read_end)

    if np.sum(np.logical_not(np.isnan(r_z_scores))) == 0:
        raise NotImplementedError('No valid z-scores in read.')
    r_pvals, r_poss = get_pvals(r_z_scores)
    if fm_offset > 0:
        r_pvals = calc_window_fishers_method(r_pvals, fm_offset)
        r_poss = np.where(np.logical_not(np.isnan(r_pvals)))[0]
        r_pvals = r_pvals[r_poss]
    else:
        r_pvals = r_pvals[r_poss]

    r_poss += read_start

    return r_pvals, r_poss, read_id

def compute_de_novo_read_stats(
        r_data, gnm_begin_lag, gnm_end_lag, fm_offset, reg_start,
        region_size, std_ref):
    """
    Compute signficance statistics using de novo comparison to a canonical model
    method for a single read within a specified genomic region.
    """
    def de_novo_clip_and_flip():
        with h5py.File(r_data.fn, 'r') as fast5_data:
            r_means, r_seq = th.get_multiple_slots_read_centric(
                fast5_data, ['norm_mean', 'base'], r_data.corr_group)
            read_id = th.get_raw_read_slot(fast5_data).attrs['read_id']

        if r_means is None or r_seq is None:
            raise NotImplementedError(
                'Read does not contain valid re-squiggled data.')
        r_seq = b''.join(r_seq).decode()

        read_start, read_end = r_data.start, r_data.end
        # clip read if it extends outside the current genomic region, so
        # stats are only computed within this region
        if read_start + gnm_begin_lag + fm_offset < reg_start:
            num_start_clip = reg_start - (read_start + gnm_begin_lag + fm_offset)
            read_start = reg_start - gnm_begin_lag - fm_offset
            if r_data.strand == '+':
                r_means = r_means[num_start_clip:]
                r_seq = r_seq[num_start_clip:]
            else:
                r_means = r_means[:-num_start_clip]
                r_seq = r_seq[:-num_start_clip]
        if read_end - gnm_end_lag - fm_offset > reg_start + region_size:
            num_end_clip = (read_end - gnm_end_lag - fm_offset) - (
                reg_start + region_size)
            read_end = reg_start + region_size + gnm_end_lag + fm_offset
            if r_data.strand == '+':
                r_means = r_means[:-num_end_clip]
                r_seq = r_seq[:-num_end_clip]
            else:
                r_means = r_means[num_end_clip:]
                r_seq = r_seq[num_end_clip:]

        # if this read does not cover enough of this region for stat
        # computation raise an error to be handled below
        if len(r_seq) < std_ref.kmer_width:
            raise NotImplementedError(
                'Read does not contain information in this region.')

        r_ref_means, r_ref_sds, _, _ = get_ref_from_seq(
            r_seq, std_ref, r_data.strand == '-')

        if r_data.strand == '-':
            # reverse means to match genomic order
            r_means = r_means[::-1]
        # clip means that don't have testing model data
        # note that this is still extended by fm_offset
        r_means = r_means[gnm_begin_lag:-gnm_end_lag]
        read_start += gnm_begin_lag
        read_end -= gnm_end_lag

        return (r_means, r_ref_means, r_ref_sds, read_start,
                read_end, read_id)


    (r_means, r_ref_means, r_ref_sds,
     read_start, read_end, read_id) = de_novo_clip_and_flip()

    z_scores = np.abs(r_means - r_ref_means) / r_ref_sds
    r_pvals = stats.norm.cdf(-z_scores) * 2.0
    if fm_offset > 0:
        r_pvals = calc_window_fishers_method(r_pvals, fm_offset)

    # ignore errors in max over NAN values if fisher's method was used
    with np.errstate(invalid='ignore'):
        r_pvals = np.maximum(r_pvals, SMALLEST_PVAL)

    r_poss = np.arange(read_start, read_end)

    return r_pvals, r_poss, read_id

def calc_llh_ratio(reg_means, reg_ref_means, reg_ref_vars,
                   reg_alt_means, reg_alt_vars):
    """
    Compute log likelihood ratio

    This is about 10X slower than the cython version in tombo.c_helper, but
    has been kept for debugging purposes
    """
    # compute log likelihood ratio
    # positive value means standard base fits data better
    # negative value means alternative base fits data better
    return (np.sum(np.square(reg_means - reg_alt_means) / reg_alt_vars) +
            np.sum(np.log(reg_alt_vars))) - (
                np.sum(np.square(reg_means - reg_ref_means) / reg_ref_vars) +
                np.sum(np.log(reg_ref_vars)))

def compute_alt_model_read_stats(
        r_data, gnm_begin_lag, gnm_end_lag, reg_start, region_size,
        std_ref, alt_ref, use_standard_llhr):
    """
    Compute signficance statistics using comparison of read signal to canonical
    and alternative models method for a single read within a specified genomic
    region.
    """
    motif_width = gnm_begin_lag + gnm_end_lag + 1
    def alt_clip_and_flip():
        with h5py.File(r_data.fn, 'r') as fast5_data:
            r_means, r_seq = th.get_multiple_slots_read_centric(
                fast5_data, ['norm_mean', 'base'], r_data.corr_group)
            read_id = th.get_raw_read_slot(fast5_data).attrs['read_id']

        if r_means is None or r_seq is None:
            raise NotImplementedError(
                'Read does not contain valid re-squiggled data.')
        r_seq = b''.join(r_seq).decode()

        read_start = r_data.start
        # clip read if it extends outside the current genomic region, so
        # stats are only computed within this region
        if read_start + motif_width - 1 < reg_start:
            num_start_clip = reg_start - (read_start + motif_width - 1)
            read_start = reg_start - (motif_width - 1)
            if r_data.strand == '+':
                r_means = r_means[num_start_clip:]
                r_seq = r_seq[num_start_clip:]
            else:
                r_means = r_means[:-num_start_clip]
                r_seq = r_seq[:-num_start_clip]
        if r_data.end - (motif_width - 1) > reg_start + region_size:
            num_end_clip = (r_data.end - (motif_width - 1)) - (
                reg_start + region_size)
            if r_data.strand == '+':
                r_means = r_means[:-num_end_clip]
                r_seq = r_seq[:-num_end_clip]
            else:
                r_means = r_means[num_end_clip:]
                r_seq = r_seq[num_end_clip:]

        # if this read does not cover enough of this region for stat
        # computation raise an error to be handled below
        if len(r_seq) < std_ref.kmer_width:
            raise NotImplementedError(
                'Read does not contain information in this region.')

        r_ref_means, r_ref_sds, r_alt_means, r_alt_sds = get_ref_from_seq(
            r_seq, std_ref, r_data.strand == '-', alt_ref)

        if r_data.strand == '-':
            # reverse means and seq to match genomic order
            r_means = r_means[::-1]
            r_seq = r_seq[::-1]
        # clip means to individual tested positions
        r_means = r_means[gnm_begin_lag:-gnm_end_lag]
        # trim seq to positions with valid llh ratio test results
        # this is shorter than the means and model
        r_seq = r_seq[(motif_width - 1):-(motif_width - 1)]
        read_start += motif_width - 1

        return (r_means, r_seq, r_ref_means, r_ref_sds, read_start,
                r_alt_means, r_alt_sds, read_id)


    (r_means, r_seq, r_ref_means, r_ref_sds, read_start,
     r_alt_means, r_alt_sds, read_id) = alt_clip_and_flip()
    r_ref_vars = np.square(r_ref_sds)
    r_alt_vars = np.square(r_alt_sds)

    alt_base_poss = []
    log_lh_ratios = []
    # note search space is clipped since all k-mers covering the position
    # of interest must be valid
    for alt_base_pos in re.finditer(alt_ref.alt_base, r_seq):
        alt_pos = alt_base_pos.start()
        alt_base_poss.append(alt_pos + read_start)
        if CONST_SD_MODEL:
            const_var = r_ref_vars[alt_pos]
            if use_standard_llhr:
                pos_lh_ratio = c_calc_llh_ratio_const_var(
                    r_means[alt_pos:alt_pos + motif_width],
                    r_ref_means[alt_pos:alt_pos + motif_width],
                    r_alt_means[alt_pos:alt_pos + motif_width],
                    const_var)
            else:
                pos_lh_ratio = c_calc_scaled_llh_ratio_const_var(
                    r_means[alt_pos:alt_pos + motif_width],
                    r_ref_means[alt_pos:alt_pos + motif_width],
                    r_alt_means[alt_pos:alt_pos + motif_width],
                    const_var, OCLLHR_SCALE, OCLLHR_HEIGHT, OCLLHR_POWER)
        else:
            if use_standard_llhr:
                pos_lh_ratio = c_calc_llh_ratio(
                    r_means[alt_pos:alt_pos + motif_width],
                    r_ref_means[alt_pos:alt_pos + motif_width],
                    r_ref_vars[alt_pos:alt_pos + motif_width],
                    r_alt_means[alt_pos:alt_pos + motif_width],
                    r_alt_vars[alt_pos:alt_pos + motif_width])
            else:
                raise NotImplementedError(
                    'Variable SD scaled likelihood ratio not implemented.')
        log_lh_ratios.append(pos_lh_ratio)

    return np.array(log_lh_ratios), np.array(alt_base_poss), read_id

def compute_read_stats(
        chrm, strand, reg_start, reg_reads, fm_offset, min_test_vals,
        region_size, single_read_thresh, lower_thresh, ctrl_reg_reads, std_ref,
        alt_ref, use_standard_llhr, per_read_q, stat_type):
    if stat_type == SAMP_COMP_TXT:
        ctrl_means, ctrl_sds, ctrl_cov = get_reads_ref(
            ctrl_reg_reads, reg_start, region_size,
            min_test_vals, fm_offset)
    else:
        ctrl_cov = None
        # compute begin and end lag wrt the genome from upstream and downstream
        # which are wrt to the read
        dnstrm_bases = std_ref.kmer_width - std_ref.central_pos - 1
        gnm_begin_lag = std_ref.central_pos if strand == '+' else dnstrm_bases
        gnm_end_lag = dnstrm_bases if strand == '+' else std_ref.central_pos

    reg_read_stats, reg_poss, reg_ids = [], [], []
    for r_data in reg_reads:
        try:
            if stat_type == SAMP_COMP_TXT:
                r_stats, r_poss, read_id = compute_sample_compare_read_stats(
                    r_data, ctrl_means, ctrl_sds, fm_offset, reg_start,
                    region_size)
            elif stat_type == DE_NOVO_TXT:
                r_stats, r_poss, read_id = compute_de_novo_read_stats(
                    r_data, gnm_begin_lag, gnm_end_lag, fm_offset,
                    reg_start, region_size, std_ref)
            else:
                r_stats, r_poss, read_id = compute_alt_model_read_stats(
                    r_data, gnm_begin_lag, gnm_end_lag, reg_start, region_size,
                    std_ref, alt_ref, use_standard_llhr)
        except NotImplementedError:
            continue
        if r_stats is None: continue
        reg_read_stats.append(r_stats)
        reg_poss.append(r_poss)
        reg_ids.append(read_id)

    if len(reg_read_stats) == 0:
        raise NotImplementedError

    if per_read_q is not None:
        # compile read_ids vector for per-read output
        reg_ids = [(r_id, r_poss.shape[0])
                   for r_id, r_poss, in zip(reg_ids, reg_poss)]

    reg_read_stats = np.concatenate(reg_read_stats)
    reg_poss = np.concatenate(reg_poss)
    # remove nans possibly introduced by fisher's method calculcations
    valid_poss = ~np.isnan(reg_read_stats)
    reg_poss = reg_poss[valid_poss]
    reg_read_stats = reg_read_stats[valid_poss]
    assert reg_poss.shape[0] == reg_read_stats.shape[0], '\t'.join(map(str, (
        reg_poss.shape[0], reg_read_stats.shape[0])))

    if per_read_q is not None:
        valid_reg_ids = [
            rep_r_id for rep_r_id, is_valid in zip(
                [rep_r_id for r_id, r_len in reg_ids
                 for rep_r_id in repeat(r_id, r_len)], valid_poss) if is_valid]
        read_id_lookup = dict((
            (read_id, read_id_val)
            for read_id_val, read_id in enumerate(set(valid_reg_ids))))
        conv_reg_ids = np.array([
            read_id_lookup[r_id] for r_id in valid_reg_ids])
        assert conv_reg_ids.shape[0] == reg_poss.shape[0]
        per_read_block = np.array(
            list(zip(reg_poss, reg_read_stats, conv_reg_ids)),
            dtype=[(str('pos'), 'u4'), (str('stat'), 'f8'),
                   (str('read_id'), 'u4')])
        per_read_q.put((
            per_read_block, read_id_lookup, chrm, strand, reg_start))

    # get order of all bases from position array
    as_reg_poss = np.argsort(reg_poss)
    # sort all positions from all reads
    reg_poss = reg_poss[as_reg_poss]
    # get unique tested genomic positions across all reads
    us_reg_poss = np.unique(reg_poss)

    # then sort the stats array by genomic position and
    # split into stats by genomic base position
    reg_base_stats = np.split(
        reg_read_stats[as_reg_poss],
        np.where(np.concatenate([[0,], np.diff(reg_poss)]) > 0)[0])

    reg_cov = [base_stats.shape[0] for base_stats in reg_base_stats]

    if lower_thresh is not None:
        # filter base statistics that fall between the upper and lower
        # stat threshold for the log likelihood statistic
        reg_base_stats = [
            base_stats[np.logical_or(base_stats <= lower_thresh,
                                     base_stats >= single_read_thresh)]
            for base_stats in reg_base_stats]
        valid_cov = [base_stats.shape[0] for base_stats in reg_base_stats]
    elif stat_type == ALT_MODEL_TXT:
        # filter base statistics that fall between the upper and lower
        # stat threshold for the log likelihood statistic
        reg_base_stats = [base_stats[np.abs(base_stats) >= single_read_thresh]
                          for base_stats in reg_base_stats]
        valid_cov = [base_stats.shape[0] for base_stats in reg_base_stats]
    else:
        valid_cov = reg_cov

    if stat_type == SAMP_COMP_TXT:
        ctrl_cov = [ctrl_cov[pos] if pos in ctrl_cov else 0
                    for pos in reg_poss]
    else:
        # convert to list since python2 repeat objects can't be pickled
        ctrl_cov = list(repeat(0, reg_poss.shape[0]))

    return reg_base_stats, us_reg_poss, reg_cov, ctrl_cov, valid_cov

def get_region_stats(
        chrm, strand, reg_start, reg_reads, fm_offset, min_test_vals,
        region_size, single_read_thresh, lower_thresh, ctrl_reg_reads, std_ref,
        alt_ref, use_standard_llhr, per_read_q, stat_type):
    """
    Compute requested statistics for a specific region of the genome
    """
    try:
        (reg_base_stats, reg_poss,
         reg_cov, ctrl_cov, valid_cov) = compute_read_stats(
             chrm, strand, reg_start, reg_reads, fm_offset, min_test_vals,
             region_size, single_read_thresh, lower_thresh, ctrl_reg_reads,
             std_ref, alt_ref, use_standard_llhr, per_read_q, stat_type)
    except NotImplementedError:
        return None

    if reg_poss.shape[0] == 0:
        return None

    reg_frac_standard_base = np.array([
        np.greater_equal(
            base_stats, single_read_thresh).sum() / base_stats.shape[0]
        if base_stats.shape[0] > 0 else np.NAN
        for base_stats in reg_base_stats])

    reg_stats = (reg_frac_standard_base, reg_poss, chrm, strand,
                 reg_cov, ctrl_cov, valid_cov)

    return reg_stats

def _test_signif_worker(
        region_q, stats_q, progress_q, per_read_q, raw_read_coverage, fm_offset,
        min_test_vals, single_read_thresh, lower_thresh, region_size,
        ctrl_read_coverage, std_ref, alt_ref, use_standard_llhr, stat_type):
    ctrl_reg_reads = None
    while not region_q.empty():
        try:
            chrm, strand, reg_start = region_q.get(block=False)
        except queue.Empty:
            # sometimes throws false empty error with get(block=False)
            if not region_q.empty():
                continue
            break

        reg_reads = [r_data for r_data in raw_read_coverage[(chrm, strand)]
                     if not (r_data.start >= reg_start + region_size or
                             r_data.end <= reg_start)]
        if len(reg_reads) == 0:
            progress_q.put(1)
            continue

        if ctrl_read_coverage is not None:
            ctrl_reg_reads = [
                r_data for r_data in ctrl_read_coverage[(chrm, strand)]
                if not (r_data.start >= reg_start + region_size or
                        r_data.end <= reg_start)]
        reg_stats = get_region_stats(
            chrm, strand, reg_start, reg_reads, fm_offset, min_test_vals,
            region_size, single_read_thresh, lower_thresh, ctrl_reg_reads,
            std_ref, alt_ref, use_standard_llhr, per_read_q, stat_type)
        if reg_stats is not None:
            stats_q.put(reg_stats)
        progress_q.put(1)

    return

if _PROFILE_SIGNIF:
    _test_signif_wrapper = _test_signif_worker
    def _test_signif_worker(*args):
        import cProfile
        cProfile.runctx('_test_signif_wrapper(*args)', globals(), locals(),
                        filename='test_signif.prof')
        return


##############################################
########## Testing Multi-processing ##########
##############################################

def _get_progress_queue(progress_q, prog_conn, num_regions):
    th._status_message(
        'Performing modified base detection across genomic regions.')
    bar = tqdm(total=num_regions, smoothing=0)

    tot_num_rec_proc = 0
    while True:
        try:
            iter_val = progress_q.get(block=False)
            tot_num_rec_proc += iter_val
            bar.update(iter_val)
        except queue.Empty:
            if prog_conn.poll():
                break
            sleep(0.1)
            continue

    bar.close()
    prog_conn.send(tot_num_rec_proc)

    return

def _get_stats_queue(stats_q, stats_conn, min_test_reads, stats_file_bn,
                     alt_name, stat_type):
    # TODO convert to a TomboStats class that writes each batch to file as
    # they are received
    all_reg_stats = []
    while True:
        try:
            reg_stats = stats_q.get(block=False)
            all_reg_stats.append(reg_stats)
        except queue.Empty:
            if stats_conn.poll():
                sleep(0.1)
                break
            sleep(0.1)
            continue

    # Clear leftover values from queues
    while not stats_q.empty():
        reg_stats = stats_q.get(block=False)
        all_reg_stats.append(reg_stats)

    if len(all_reg_stats) == 0:
        th._error_message_and_exit(
            'No genomic positions contain --minimum-test-reads.')

    write_stats(all_reg_stats, stats_file_bn,
                stat_type, min_test_reads, alt_name)
    stats_conn.send(True)

    return

def _get_per_read_queue(
        per_read_q, per_read_conn, per_read_fn, stat_type, region_size):
    per_read_stats = PerReadStats(per_read_fn, stat_type, region_size)

    while True:
        try:
            per_read_block = per_read_q.get(block=False)
            per_read_stats.write_per_read_block(*per_read_block)
            del per_read_block
        except queue.Empty:
            if per_read_conn.poll():
                sleep(0.1)
                break
            sleep(0.1)
            continue

    # Clear leftover values from queues
    while not per_read_q.empty():
        per_read_block = per_read_q.get(block=False)
        per_read_stats.write_per_read_block(*per_read_block)
        del per_read_block
    per_read_stats.close()

    # indicate that the process has closed
    per_read_conn.send(True)

    return

def test_significance(
        raw_read_coverage, min_test_vals, fm_offset, single_read_thresh,
        lower_thresh, region_size, num_processes, per_read_bn, stat_type,
        min_test_reads, stats_file_bn,
        ctrl_read_coverage=None, std_ref=None, alt_ref=None,
        use_standard_llhr=False, alt_name=None):
    """
    Test for significant shifted signal in mutliprocessed batches
    """
    region_q = Queue()
    stats_q = Queue()
    progress_q = Queue()
    per_read_q = Queue(PER_READ_BLOCKS_QUEUE_LIMIT) \
                 if per_read_bn else None
    # split chromosomes into separate regions to process independently
    chrm_sizes = th.get_chrm_sizes(raw_read_coverage, ctrl_read_coverage)
    num_regions = 0
    for chrm, chrm_len in chrm_sizes.items():
        # only process regions covered by both samples if control
        # reads are provided
        plus_covered = (
            (chrm, '+') in raw_read_coverage and
            (ctrl_read_coverage is None or (chrm, '+') in ctrl_read_coverage))
        minus_covered = (
            (chrm, '-') in raw_read_coverage and
            (ctrl_read_coverage is None or (chrm, '-') in ctrl_read_coverage))
        for reg_start in range(0, chrm_len, region_size):
            if plus_covered:
                region_q.put((chrm, '+', reg_start))
                num_regions += 1
            if minus_covered:
                region_q.put((chrm, '-', reg_start))
                num_regions += 1

    test_args = (
        region_q, stats_q, progress_q, per_read_q, raw_read_coverage, fm_offset,
        min_test_vals, single_read_thresh, lower_thresh, region_size,
        ctrl_read_coverage, std_ref, alt_ref, use_standard_llhr, stat_type)
    test_ps = []
    for p_id in range(num_processes):
        p = Process(target=_test_signif_worker, args=test_args)
        p.start()
        test_ps.append(p)

    # start queue getter processes
    if VERBOSE:
        main_prog_conn, prog_conn = Pipe()
        prog_p = Process(target=_get_progress_queue,
                            args=(progress_q, prog_conn, num_regions))
        prog_p.daemon = True
        prog_p.start()

    # main region stats queue getter
    main_stats_conn, stats_conn = Pipe()
    stats_p = Process(target=_get_stats_queue, args=(
        stats_q, stats_conn, min_test_reads, stats_file_bn, alt_name, stat_type))
    stats_p.daemon = True
    stats_p.start()

    # per-read stats queue getter
    if per_read_bn is not None:
        if stat_type == ALT_MODEL_TXT:
            per_read_fn = per_read_bn + '.' + alt_name + '.tombo.per_read_stats'
        else:
            per_read_fn = per_read_bn + '.tombo.per_read_stats'
        main_per_read_conn, per_read_conn = Pipe()
        per_read_p = Process(
            target=_get_per_read_queue,
            args=(per_read_q, per_read_conn, per_read_fn, stat_type, region_size))
        per_read_p.daemon = True
        per_read_p.start()

    # wait for test processes to finish
    for test_p in test_ps:
        test_p.join()

    # in a very unlikely case the progress queue could die while the
    # main process remains active and thus we would have a deadlock here
    if VERBOSE and prog_p.is_alive():
        # send signal to getter queue to finish and return results
        main_prog_conn.send(True)
        # returns total number of processed reads if that is needed
        main_prog_conn.recv()

    if per_read_bn is not None:
        main_per_read_conn.send(True)
        main_per_read_conn.recv()

    main_stats_conn.send(True)
    main_stats_conn.recv()

    return


##########################
##### Main Functions #####
##########################

def _test_shifts_de_novo_main(
        args, lower_thresh, single_read_thresh, bio_samp_type, raw_read_coverage):
    tb_model_fn = args.tombo_model_filename
    if bio_samp_type is None:
        bio_samp_type = 'RNA' if th.is_rna(raw_read_coverage) else 'DNA'
    if tb_model_fn is None:
        tb_model_fn, bio_samp_type = get_default_standard_ref(
            raw_read_coverage, bio_samp_type)
    std_ref = TomboModel(tb_model_fn)

    stat_type = DE_NOVO_TXT
    lower_thresh, single_read_thresh = (
        (lower_thresh, single_read_thresh) if single_read_thresh
        is not None else DE_NOVO_THRESH[bio_samp_type])
    if VERBOSE: th._status_message(
            'Performing de novo model testing against canonical model.')
    test_significance(
        raw_read_coverage, args.minimum_test_reads,
        args.fishers_method_context, single_read_thresh, lower_thresh,
        args.multiprocess_region_size, args.processes,
        args.per_read_statistics_basename, stat_type,
        args.minimum_test_reads, args.statistics_file_basename,
        std_ref=std_ref)

    return

def _test_shifts_alt_main(
        args, lower_thresh, single_read_thresh, bio_samp_type, raw_read_coverage):
    tb_model_fn = args.tombo_model_filename
    if bio_samp_type is None:
        bio_samp_type = 'RNA' if th.is_rna(raw_read_coverage) else 'DNA'
    if tb_model_fn is None:
        tb_model_fn, bio_samp_type = get_default_standard_ref(
            raw_read_coverage, bio_samp_type)
    std_ref = TomboModel(tb_model_fn)

    stat_type = ALT_MODEL_TXT
    lower_thresh, single_read_thresh = (
        (lower_thresh, single_read_thresh) if single_read_thresh
        is not None else LLR_THRESH[bio_samp_type])
    if VERBOSE: th._status_message('Performing alternative model testing.')
    if args.alternate_model_filenames is not None:
        alt_refs = parse_tombo_models(
            args.alternate_model_filenames, std_ref)
    else:
        alt_refs = load_alt_refs(
            args.alternate_bases, raw_read_coverage,
            std_ref, bio_samp_type)
    if len(alt_refs) == 0:
        th._error_message_and_exit('No alternative models successfully loaded.')

    for alt_name, alt_ref in alt_refs.items():
        if VERBOSE: th._status_message(
                'Performing alternative model testing against ' +
                alt_name + ' model.')
        test_significance(
            raw_read_coverage, args.minimum_test_reads, 0,
            single_read_thresh, lower_thresh,
            args.multiprocess_region_size, args.processes,
            args.per_read_statistics_basename, stat_type,
            args.minimum_test_reads, args.statistics_file_basename,
            std_ref=std_ref, alt_ref=alt_ref,
            use_standard_llhr=args.standard_log_likelihood_ratio,
            alt_name=alt_name)

    return

def _test_shifts_samp_comp_main(
        args, lower_thresh, single_read_thresh, bio_samp_type, raw_read_coverage):
    stat_type = SAMP_COMP_TXT
    if single_read_thresh is None:
        if bio_samp_type is None:
            bio_samp_type = 'RNA' if th.is_rna(raw_read_coverage) else 'DNA'
        lower_thresh, single_read_thresh = SAMP_COMP_THRESH[bio_samp_type]
    if VERBOSE: th._status_message(
            'Performing two-sample comparison significance testing.')
    ctrl_read_coverage = th.parse_fast5s(
        args.control_fast5_basedirs, args.corrected_group,
        args.basecall_subgroups)
    test_significance(
        raw_read_coverage, args.minimum_test_reads,
        args.fishers_method_context, single_read_thresh, lower_thresh,
        args.multiprocess_region_size, args.processes,
        args.per_read_statistics_basename, stat_type,
        args.minimum_test_reads, args.statistics_file_basename,
        ctrl_read_coverage=ctrl_read_coverage)

    return

def _test_shifts_main(args):
    global VERBOSE
    VERBOSE = not args.quiet
    th.VERBOSE = VERBOSE

    # Check extra requirements for alternative model testing
    if args.action_command == 'alternative_model':
        if args.print_available_models:
            _print_alt_models()
            sys.exit()
        if args.fast5_basedirs is None or args.statistics_file_basename is None:
            th._error_message_and_exit(
                'Must provide both a set of FAST5 read files ' +
                '(--fast5-basedirs) and an output file basename ' +
                '(--statistics-file-basename).')
        if (args.alternate_model_filenames is None and
            args.alternate_bases is None):
            th._error_message_and_exit(
                'Must provide an alterntive model against which to test.\n\t' +
                'Run with --print-available-models option to see possible ' +
                'values for the --alternate-bases option.')

    if args.single_read_threshold is None:
        lower_thresh = None
        single_read_thresh = None
    elif len(args.single_read_threshold) == 1:
        single_read_thresh = args.single_read_threshold[0]
        lower_thresh = None
    else:
        if len(args.single_read_threshold) > 2:
            th._warning_message(
                'Only 1 or 2 values may be passed as single-read ' +
                'thresholds. Only using the first 2 options provided.')
        lower_thresh = args.single_read_threshold[0]
        single_read_thresh = args.single_read_threshold[1]

    try:
        # sample compare does not have bio_sample_type in the namespace
        bio_samp_type = args.bio_sample_type
    except AttributeError:
        bio_samp_type = None

    raw_read_coverage = th.parse_fast5s(
        args.fast5_basedirs, args.corrected_group, args.basecall_subgroups)

    if args.action_command == 'de_novo':
        _test_shifts_de_novo_main(
            args, lower_thresh, single_read_thresh, bio_samp_type,
            raw_read_coverage)
    elif args.action_command == 'alternative_model':
        _test_shifts_alt_main(
            args, lower_thresh, single_read_thresh, bio_samp_type,
            raw_read_coverage)
    elif args.action_command == 'sample_compare':
        _test_shifts_samp_comp_main(
            args, lower_thresh, single_read_thresh, bio_samp_type,
            raw_read_coverage)
    else:
        th._error_message_and_exit('Invalid Tombo detect_modifications command.')

    return

def _aggregate_per_read_main(args):
    global VERBOSE
    VERBOSE = not args.quiet
    th.VERBOSE = VERBOSE

    if len(args.single_read_threshold) == 1:
        lower_thresh = None
        single_read_thresh = args.single_read_threshold[0]
    else:
        if len(args.single_read_threshold) > 2:
            th._warning_message(
                'Only 1 or 2 values may be passed as single-read ' +
                'thresholds. Only using the first 2 options provided.')
        lower_thresh = args.single_read_threshold[0]
        single_read_thresh = args.single_read_threshold[1]

    all_reg_stats, stat_type = apply_per_read_thresh(
        args.per_read_statistics_filename, single_read_thresh,
        args.minimum_test_reads, lower_thresh)

    write_stats(all_reg_stats, args.statistics_file_basename, stat_type,
                args.minimum_test_reads)

    return

def _est_ref_main(args):
    global VERBOSE
    VERBOSE = not args.quiet
    th.VERBOSE = VERBOSE

    if min(args.upstream_bases, args.downstream_bases) == 0:
        th._error_message_and_exit(
            'Context upstream and downstream must be greater ' +
            'than 0 for model estimation.')

    estimate_kmer_model(
        args.fast5_basedirs, args.corrected_group, args.basecall_subgroups,
        args.tombo_model_filename, args.minimum_test_reads,
        args.upstream_bases, args.downstream_bases,
        args.minimum_kmer_observations, args.kmer_specific_sd,
        args.coverage_threshold, args.estimate_mean,
        args.multiprocess_region_size, args.processes)

    return

def _est_alt_ref_main(args):
    global VERBOSE
    VERBOSE = not args.quiet
    th.VERBOSE = VERBOSE

    alt_ref = estimate_alt_model(
        args.fast5_basedirs, args.control_fast5_basedirs,
        args.corrected_group, args.basecall_subgroups,
        args.tombo_model_filename, args.bio_sample_type,
        args.alternate_model_base, args.alt_fraction_percentile,
        args.minimum_kmer_observations, args.save_density_basename,
        args.kernel_density_bandwidth, args.alternate_density_filename,
        args.control_density_filename, args.processes)
    # returns None when profiling method
    if alt_ref is None: return
    alt_ref.write_model(args.alternate_model_filename,
                        args.alternate_model_base, args.alternate_model_name)

    return

def _estimate_scale_main(args):
    global VERBOSE
    VERBOSE = not args.quiet
    th.VERBOSE = VERBOSE

    if VERBOSE: th._status_message('Getting files list.')
    try:
        if not os.path.isdir(args.fast5_basedir):
            th._error_message_and_exit(
                'Provided [fast5-basedir] is not a directory.')
        fast5_basedir = (
            args.fast5_basedir if args.fast5_basedir.endswith('/') else
            args.fast5_basedir + '/')
        fast5_fns = th.get_files_list(fast5_basedir)
    except OSError:
        th._error_message_and_exit(
            'Reads base directory, a sub-directory or an old (hidden) ' +
            'index file does not appear to be accessible. Check ' +
            'directory permissions.')
    if len(fast5_fns) < 1:
        th._error_message_and_exit(
            'No files identified in the specified ' +
            'directory or within immediate subdirectories.')

    th._status_message('Global scaling estimate: ' +
                       unicode(estimate_global_scale(fast5_fns)))

    return


if __name__ == '__main__':
    raise NotImplementedError(
        'This is a module. See commands with `tombo -h`')
