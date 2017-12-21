import sys, os

import re
import h5py
import string
import fnmatch

import numpy as np

from glob import glob
from operator import itemgetter
from itertools import izip, repeat, islice
from collections import defaultdict, namedtuple

from _version import TOMBO_VERSION
from c_helper import c_new_mean_stds, c_new_means, c_apply_outlier_thresh

SMALLEST_PVAL = 1e-50

readData = namedtuple('readData', (
    'start', 'end', 'filtered', 'read_start_rel_to_raw',
    'strand', 'fn', 'corr_group', 'rna'))
intervalData = namedtuple('intervalData', (
    'reg_id', 'chrm', 'start', 'end', 'strand', 'reg_text', 'reads', 'seq'))
# set default values for strand, text, reads and seq
intervalData.__new__.__defaults__ = (None, '', None, None)

channelInfo = namedtuple(
    'channelInfo',
    ('offset', 'range', 'digitisation', 'number', 'sampling_rate'))
scaleValues = namedtuple(
    'scaleValues',
    ('shift', 'scale', 'lower_lim', 'upper_lim'))
genomeLoc = namedtuple(
    'genomeLoc', ('Start', 'Strand', 'Chrom'))

NORM_TYPES = ('none', 'pA', 'pA_raw', 'median', 'robust_median',
              'median_const_scale')
STANDARD_MODELS = {'DNA':'tombo.DNA.model',
                   'RNA':'tombo.RNA.200mV.model'}
ALTERNATE_MODELS = {'DNA_5mC':'tombo.DNA.5mC.model',}

# single base conversion for motifs
SINGLE_LETTER_CODE = {
    'A':'A', 'C':'C', 'G':'G', 'T':'T', 'B':'[CGT]',
    'D':'[AGT]', 'H':'[ACT]', 'K':'[GT]', 'M':'[AC]',
    'N':'[ACGT]', 'R':'[AG]', 'S':'[CG]', 'V':'[ACG]',
    'W':'[AT]', 'Y':'[CT]'}

FN_SPACE_FILLER = '|||'
FASTA_NAME_JOINER = ':::'

VERBOSE = False

# got quantiles from analysis of stability after shift-only normalization
robust_quantiles = (46.5, 53.5)


######################################
###### Various Helper Functions ######
######################################

COMP_BASES = string.maketrans('ACGT', 'TGCA')
def comp_seq(seq):
    """
    Complement DNA sequence
    """
    return seq.translate(COMP_BASES)
def rev_comp(seq):
    """
    Reverse complement DNA sequence
    """
    return seq.translate(COMP_BASES)[::-1]

U_TO_T = string.maketrans('U', 'T')
def rev_transcribe(seq):
    """
    Convert U bases to T
    """
    return seq.translate(U_TO_T)

def parse_fasta(fasta_fn):
    """
    Parse a fasta file for sequence extraction (mostly for BAM processing)
    """
    # Tried Biopython index and that opened the fail again for each
    # record access request and was thus far too slow

    # could consider a conditional dependence on pyfaix if on-disk
    # indexing is required for larger genomes
    # testing shows that human genome only takes 3.2 GB with raw parser
    # so raw parsing is probably fine
    fasta_fp = open(fasta_fn)

    fasta_records = {}
    curr_id = None
    curr_seq = ''
    for line in fasta_fp:
        if line.startswith('>'):
            if (curr_id is not None and
                curr_seq is not ''):
                fasta_records[curr_id] = curr_seq
            curr_seq = ''
            curr_id = line.replace(">","").strip().split()[0]
        else:
            curr_seq += line.strip()

    # add last record
    if (curr_id is not None and
        curr_seq is not ''):
        fasta_records[curr_id] = curr_seq

    fasta_fp.close()

    return fasta_records

def get_chrm_sizes(raw_read_coverage, raw_read_coverage2=None):
    """
    Get covered chromosome sizes from a set of reads
    """
    strand_chrm_sizes = defaultdict(list)
    for (chrm, strand), cs_read_cov in \
        raw_read_coverage.iteritems():
        strand_chrm_sizes[chrm].append(max(
            r_data.end for r_data in cs_read_cov))
    if raw_read_coverage2 is not None:
        for (chrm, strand), cs_read_cov in \
            raw_read_coverage2.iteritems():
            strand_chrm_sizes[chrm].append(max(
                r_data.end for r_data in cs_read_cov))

    return dict((chrm, max(strnd_sizes))
                for chrm, strnd_sizes in
                strand_chrm_sizes.iteritems())

def parse_motif(motif):
    """
    Parse a single letter code motif into a pattern for matching
    """
    invalid_chars = re.findall(
        '[^' + ''.join(SINGLE_LETTER_CODE.keys()) + ']',
        motif)
    if len(invalid_chars) > 0:
       sys.stderr.write(
           '********* ERROR *********\n\tInvalid characters in motif: ' +
           ', '.join(invalid_chars) + '\n')
       sys.exit()

    return re.compile(''.join(
        SINGLE_LETTER_CODE[letter] for letter in motif))


#############################################
###### Automatic Sample-Type Detection ######
#############################################

def is_read_rna(fast5_data):
    """
    Determine if a read is RNA or DNA
    """
    # check both experiment type and kit slots for "rna"
    exp_type, exp_kit = None, None
    try:
        exp_type = fast5_data['UniqueGlobalKey/context_tags'].attrs[
            'experiment_type']
        # remove the word internal since it contains rna.
        exp_type = exp_type.replace('internal', '')
    except:
        pass
    try:
        exp_kit = fast5_data['UniqueGlobalKey/context_tags'].attrs[
            'experiment_kit']
        # remove the word internal since it contains rna.
        exp_kit = exp_kit.replace('internal', '')
    except:
        pass

    if exp_type is None and exp_kit is None:
        rna = False
    else:
        rna = (
            (exp_type is not None and re.search('rna', exp_type) is not None) or
            (exp_kit is not None and re.search('rna', exp_kit) is not None))

    return rna

def is_rna(raw_read_coverage, n_reads=10):
    """
    Determine if a set of reads are RNA or DNA from a small sample
    """
    proc_reads = 0
    for cs_reads in raw_read_coverage.itervalues():
        for r_data in cs_reads:
            if not r_data.rna:
                return False
            proc_reads += 1
            if proc_reads >= n_reads:
                break
        if proc_reads >= n_reads:
            break
    return True

def is_rna_from_files(fast5_fns, n_reads=10):
    """
    Determine if a set of files are RNA or DNA from a small sample
    """
    proc_reads = 0
    for fast5_fn in fast5_fns:
        try:
            fast5_data = h5py.File(fast5_fn, 'r')
            if not is_read_rna(fast5_data):
                return False
            proc_reads += 1
        except:
            continue
        if proc_reads >= n_reads:
            break
    return True


#########################################
###### Index File/Filter Functions ######
#########################################

def get_index_fn(fast5s_dir, corr_grp):
    """
    Get the filename for the requested directory and corrected group
    """
    # dirname should come in with a trailing slash
    split_dir = os.path.split(fast5s_dir[:-1])
    return os.path.join(split_dir[0], "." + split_dir[1] +
                        "." + corr_grp + '.tombo.index')

def prep_index_data(
        fast5_fn, genome_loc, read_start_rel_to_raw, segs,
        corrected_group, subgroup, rna, obs_filter):
    """
    Prepare data for storage in the index file
    """
    if obs_filter is None:
        is_filtered = False
    else:
        base_lens = np.diff(segs)
        is_filtered = any(np.percentile(base_lens, pctl) > thresh
                          for pctl, thresh in obs_filter)
    mapped_end = genome_loc.Start + len(segs) - 1

    return ((genome_loc.Chrom, genome_loc.Strand), (
        fast5_fn, genome_loc.Start, mapped_end, read_start_rel_to_raw,
        corrected_group, subgroup, is_filtered, rna))

def write_index_file(all_index_data, index_fn, basedir):
    """
    Write index file
    """
    try:
        import cPickle as pickle
    except:
        import pickle
    index_data = defaultdict(list)
    for chrm_strand, (fn, start, end, rsrtr, c_grp, s_grp,
                      filtered, rna) in all_index_data:
        # clip the basedir off the FAST5 filename in case later functions are
        # called from another relative path
        from_base_fn = fn.replace(basedir, '')
        index_data[chrm_strand].append((
            from_base_fn, start, end, rsrtr, c_grp, s_grp, filtered, rna))

    with open(index_fn, 'w') as index_fp:
        pickle.dump(dict(index_data), index_fp)

    return

def clear_filters(fast5s_dir, corr_grp):
    """
    Clear filters applied to this directories index files
    """
    fast5s_dir = (fast5s_dir if fast5s_dir.endswith('/') else
                  fast5s_dir + '/')
    index_fn = get_index_fn(fast5s_dir, corr_grp)
    try:
        import cPickle as pickle
    except:
        import pickle
    try:
        with open(index_fn, 'rb') as index_fp:
            index_data = pickle.load(index_fp)
    except IOError:
        sys.stderr.write(
            '******** ERRROR *******\n\tFilters can only be applied to runs ' +
            'with a Tombo index file. Re-run resquiggle without the ' +
            '--skip-index option to apply filters.\n')
        sys.exit()
    new_index_data = []
    for chrm_strand, cs_raw_data in index_data.iteritems():
        new_index_data.extend([(chrm_strand, (
            from_base_fn, start, end, rsrtr, corr_grp, s_grp, False, rna))
                               for from_base_fn, start, end, rsrtr, c_grp,
                               s_grp, filtered, rna in cs_raw_data])

    write_index_file(new_index_data, index_fn, fast5s_dir)

    return

def parse_obs_filter(obs_filter):
    """
    Parse observations per base formatted filtering
    """
    if len(obs_filter) < 1:
        return None

    # parse obs_filter
    try:
        obs_filter = [map(int, pctl_nobs.split(':'))
                      for pctl_nobs in obs_filter]
    except:
        raise RuntimeError, 'Invalid format for observation filter'

    if any(pctl < 0 or pctl > 100 for pctl in zip(*obs_filter)[0]):
       sys.stderr.write(
           '********* ERROR ********* Invalid percentile value. ' +
           ' *********\n')
       sys.exit()

    return obs_filter

def filter_reads(fast5s_dir, corr_grp, obs_filter):
    """
    Filter reads based on some observation per base threshold criteria
    """
    def read_is_stuck(fast5_fn, s_grp):
        try:
            fast5_data = h5py.File(fast5_fn, 'r')
            event_data = fast5_data['/Analyses/' + corr_grp + '/' + s_grp +
                                     '/Events'].value
            events_end = event_data[-1]['start'] + event_data[-1]['length']
            base_lens = np.diff(np.concatenate([
                event_data['start'], [events_end,]]))
            return any(np.percentile(base_lens, pctl) > thresh
                       for pctl, thresh in obs_filter)
        except:
            return True

    fast5s_dir = (fast5s_dir if fast5s_dir.endswith('/') else
                  fast5s_dir + '/')
    index_fn = get_index_fn(fast5s_dir, corr_grp)
    try:
        import cPickle as pickle
    except:
        import pickle
    try:
        with open(index_fn, 'rb') as index_fp:
            index_data = pickle.load(index_fp)
    except IOError:
        sys.stderr.write(
            '******** ERRROR *******\n\tFilters can only be applied to runs ' +
            'with a Tombo index file. Re-run resquiggle without the ' +
            '--skip-index option to apply filters.\n')
        sys.exit()
    filt_index_data = []
    num_reads, num_filt_reads = 0, 0
    for chrm_strand, cs_raw_data in index_data.iteritems():
        cs_filt_reads = [
            (chrm_strand, (
                from_base_fn, start, end, rsrtr, corr_grp, s_grp,
                read_is_stuck(fast5s_dir + '/' + from_base_fn, s_grp), rna))
            for from_base_fn, start, end, rsrtr, c_grp,
            s_grp, filtered, rna in cs_raw_data if not filtered]
        num_reads += len(cs_raw_data)
        num_filt_reads += sum([1 for i_data in cs_filt_reads if i_data[1][6]])
        filt_index_data.extend(cs_filt_reads)

    sys.stderr.write(
        'Filtered ' + str(num_filt_reads) +
        ' reads due to observations per base filter from a ' +
        'total of ' + str(num_reads) + ' reads in ' + fast5s_dir + '.\n')

    write_index_file(filt_index_data, index_fn, fast5s_dir)

    return

def filter_reads_for_coverage(fast5s_dir, corr_grp, frac_to_filter):
    fast5s_dir = (fast5s_dir if fast5s_dir.endswith('/') else
                  fast5s_dir + '/')
    index_fn = get_index_fn(fast5s_dir, corr_grp)
    try:
        import cPickle as pickle
    except:
        import pickle
    try:
        with open(index_fn, 'rb') as index_fp:
            index_data = pickle.load(index_fp)
    except IOError:
        sys.stderr.write(
            '******** ERRROR *******\n\tFilters can only be applied to runs ' +
            'with a Tombo index file. Re-run resquiggle without the ' +
            '--skip-index option to apply filters.\n')
        sys.exit()
    unfilt_data = []
    unfilt_reads_cov = []
    prev_filt_data = []
    for chrm_strand, cs_raw_data in index_data.iteritems():
        max_end = max(end for (_, _, end, _, _, _, _, _) in cs_raw_data)
        cs_coverage = np.zeros(max_end, dtype=np.int_)
        for (from_base_fn, start, end, rsrtr, c_grp,
             s_grp, filtered, rna) in cs_raw_data:
            if filtered: continue
            cs_coverage[start:end] += 1
        # now go through and compute coverage as well
        for (from_base_fn, start, end, rsrtr, c_grp,
             s_grp, filtered, rna) in cs_raw_data:
            if filtered:
                prev_filt_data.append((chrm_strand, (
                    from_base_fn, start, end, rsrtr,
                    c_grp, s_grp, filtered, rna)))
                continue
            # add approximate coverage from middle of read
            # faster than mean over the whole read
            unfilt_reads_cov.append(cs_coverage[start + int((end - start)/2)])
            unfilt_data.append((chrm_strand, (
                from_base_fn, start, end, rsrtr, c_grp, s_grp, filtered, rna)))

    num_reads = len(unfilt_data)
    num_filt_reads = int(frac_to_filter * num_reads)
    sys.stderr.write(
        'Filtered ' + str(num_filt_reads) +
        ' reads due to observations per base filter from a ' +
        'total of ' + str(num_reads) + ' reads in ' + fast5s_dir + '.\n')

    # create probabilities array with coverage values normalized to sum to 1
    unfilt_reads_cov = np.array(unfilt_reads_cov, dtype=np.float)
    unfilt_reads_p = unfilt_reads_cov / unfilt_reads_cov.sum()
    filt_indices = np.random.choice(
        num_reads, size=num_filt_reads, replace=False, p=unfilt_reads_p)
    filt_index_data = [
        (chrm_strand, (from_base_fn, start, end, rsrtr, c_grp, s_grp, True, rna))
        for (chrm_strand, (from_base_fn, start, end, rsrtr, c_grp, s_grp, _, rna))
        in itemgetter(*filt_indices)(unfilt_data)]
    unfilt_index_data = list(itemgetter(*list(set(range(num_reads)).difference(
        filt_indices)))(unfilt_data))

    write_index_file(prev_filt_data + filt_index_data + unfilt_index_data,
                     index_fn, fast5s_dir)

    return


#####################################
###### FAST5 Parsing Functions ######
#####################################

def annotate_with_fastqs(fastq_fns, fast5s_read_ids, fastq_slot):
    if VERBOSE: sys.stderr.write('Annotating FAST5s with sequence from FASTQs.\n')
    for fastq_fn in fastq_fns:
        n_recs = 0
        been_warned_ids = False
        with open(fastq_fn) as fastq_fp:
            while True:
                fastq_rec = list(islice(fastq_fp, 4))
                # if record contains fewer than 4 lines this indicates the
                # EOF, so move to next file
                if len(fastq_rec) != 4: break
                # if sequence identifier line does not start with "@" or quality
                # score line does not start with a "+" the file may be
                # corrupted, so don't process any more records
                if (re.match('@', fastq_rec[0]) is None or
                    re.match('\+', fastq_rec[2]) is None):
                    sys.stderr.write(
                        '********* WARNING ********\n\tSuccessfully parsed ' +
                        str(n_recs) + 'FASTQ records from ' + fastq_fn +
                        ' before encountering an invalid record. The rest of ' +
                        'this file will not be processed.\n')
                    break

                # extract read_id from fastq (which should be the first text
                # after the "@" record delimiter up to the first white space or
                # underscore
                read_id = fastq_rec[0].split()[0].split('_')[0][1:]
                if not been_warned_ids and read_id not in fast5s_read_ids:
                    been_warned_ids = True
                    sys.stderr.write(
                        '********* WARNING ********\n\tSome records from ' +
                        fastq_fn + ' contain read identifiers not found ' +
                        'in any FAST5 files.\n')
                    continue

                with h5py.File(fast5s_read_ids[read_id]) as fast5_data:
                    bc_slot = fast5_data[fastq_slot]
                    bc_slot.create_dataset('Fastq', data=''.join(fastq_rec))

    return

def get_files_list(fast5s_dir):
    """
    Get all fast5 files recursively listed below the directory
    """
    all_fast5s = []
    # walk through directory structure searching for fast5 files
    for root, _, fns in os.walk(fast5s_dir):
        for fn in fnmatch.filter(fns, '*.fast5'):
            all_fast5s.append(os.path.join(root, fn))

    return all_fast5s

def prep_fast5_for_fastq(
        fast5_data, basecall_group, basecall_subgroup, overwrite):
    try:
        read_id = fast5_data[
            '/Raw/Reads/'].values()[0].attrs['read_id']
    except:
        return None

    # if Analyses group doesn't exist yet, create it
    try:
        analyses_grp = fast5_data['/Analyses']
    except:
        analyses_grp = fast5_data.create_group('Analyses')

    # create Fastq slot, unless value exists and --overwrite is not set
    try:
        bc_grp = analyses_grp[basecall_group]
    except:
        bc_grp = analyses_grp.create_group(basecall_group)
        bc_subgrp = bc_grp.create_group(basecall_subgroup)
    else:
        if overwrite:
            del analyses_grp[basecall_group]
            bc_grp = analyses_grp.create_group(basecall_group)
            bc_subgrp = bc_grp.create_group(basecall_subgroup)
        else:
            raise NotImplementedError, (
                basecall_group + ' exists and --overwrite is not set.')

    return read_id

def get_read_ids_and_prep_fastq_slot(
        fast5s_dir, basecall_group, basecall_subgroup, overwrite):
    """
    Extract read id from /Raw group and prep fastq slots for annotation with
    associated FASTQ files.
    """
    if VERBOSE: sys.stderr.write(
            'Preparing reads and extracting read identifiers.\n')
    been_warned_overwrite = False
    been_warned_unique = False
    fast5s_read_ids = {}
    # walk through directory structure searching for fast5 files
    for root, _, fns in os.walk(fast5s_dir):
        for fn in fnmatch.filter(fns, '*.fast5'):
            fast5_fn = os.path.join(root, fn)
            with h5py.File(fast5_fn) as fast5_data:
                try:
                    read_id = prep_fast5_for_fastq(
                        fast5_data, basecall_group, basecall_subgroup, overwrite)
                except NotImplementedError:
                    if VERBOSE and not been_warned_overwrite:
                        been_warned_overwrite = True
                        sys.stderr.write(
                            '********* WARNING ********\n\tBasecalls exsit in ' +
                            basecall_group + ' slot. Set --overwrite option ' +
                            'to overwrite these basecalls in this slot.\n')
                    continue
            if read_id is None:
                continue
            if read_id in fast5s_read_ids:
                # Warn non-unique read_ids in directory
                if VERBOSE and not been_warned_unique:
                    been_warned_unique = True
                    sys.stderr.write(
                        '******** WARNING *********\n\tMultiple FAST5 files ' +
                        'contain the same read identifiers. Ensure that ' +
                        'FAST5 files are from a single run.\n')
                continue

            fast5s_read_ids[read_id] = fast5_fn

    return fast5s_read_ids

def parse_fast5s_wo_index(
        fast5_basedirs, corrected_group, basecall_subgroups, rna):
    """
    Parse re-squiggled reads data from a list of fast5 directories
    """
    files = [fn for fast5_basedir in fast5_basedirs
             for fn in get_files_list(fast5_basedir)]
    raw_read_coverage = defaultdict(list)
    for read_fn in files:
        try:
            read_data = h5py.File(read_fn, 'r')
        except IOError:
            # probably truncated file
            continue
        for basecall_subgroup in basecall_subgroups:
            try:
                corr_data = read_data['/'.join((
                    '/Analyses', corrected_group, basecall_subgroup))]
            except:
                continue

            try:
                align_data = dict(corr_data['Alignment'].attrs.items())
                read_start_rel_to_raw = corr_data['Events'].attrs[
                    'read_start_rel_to_raw']
            except:
                # don't warn here since errored out reads will have get here, but
                # not have alignment and events stored, so just skip these reads
                continue
            raw_read_coverage[(
                align_data['mapped_chrom'],
                align_data['mapped_strand'])].append(
                    readData(
                        align_data['mapped_start'], align_data['mapped_end'],
                        False, read_start_rel_to_raw,
                        align_data['mapped_strand'], read_fn,
                        corrected_group + '/' + basecall_subgroup, rna))

        read_data.close()

    return dict(raw_read_coverage)

def convert_index(index_data, fast5s_dir, corr_grp, new_corr_grp):
    """
    Convert an index and save under a new corrected group. Mostly for
    model_resquiggle
    """
    new_index_data = []
    for (chrm, strand), cs_raw_data in index_data.iteritems():
        for (from_base_fn, start, end, rsrtr, c_grp, s_grp,
             filtered, rna) in cs_raw_data:
            if c_grp != corr_grp: continue
            new_index_data.append(((chrm, strand), (
                from_base_fn, start, end, rsrtr,
                new_corr_grp, s_grp, filtered, rna)))

    fast5s_dir = (fast5s_dir if fast5s_dir.endswith('/') else
                  fast5s_dir + '/')
    new_index_fn = get_index_fn(fast5s_dir, new_corr_grp)
    write_index_file(new_index_data, new_index_fn, fast5s_dir)

    return

def parse_fast5s_w_index(fast5s_dir, corr_grp, subgroups, new_corr_grp):
    """
    Use index file to parse information about a set of reads
    """
    fast5s_dir = (fast5s_dir if fast5s_dir.endswith('/') else
                  fast5s_dir + '/')
    index_fn = get_index_fn(fast5s_dir, corr_grp)
    try:
        import cPickle as pickle
    except:
        import pickle
    with open(index_fn, 'rb') as index_fp:
        index_data = pickle.load(index_fp)
    raw_read_coverage = {}
    for (chrm, strand), cs_raw_data in index_data.iteritems():
        # TODO temporary check that filtered is a boolean value so that old
        # len_percentiles slots will be handled correctly should be removed
        cs_data = [
            readData(start, end, filtered, rsrtr, strand,
                     os.path.join(fast5s_dir, from_base_fn),
                     corr_grp + '/' + s_grp, rna)
            for from_base_fn, start, end, rsrtr, c_grp,
            s_grp, filtered, rna in cs_raw_data
            if c_grp == corr_grp and s_grp in subgroups and
            not (isinstance(filtered, bool) and filtered)]
        raw_read_coverage[(chrm, strand)] = cs_data
    if new_corr_grp is not None:
        # convert corrected group to new corrected group for
        # model re-squiggle
        convert_index(index_data, fast5s_dir, corr_grp, new_corr_grp)

    return raw_read_coverage

def merge_cov(w_index_covs, wo_index_cov):
    """
    Merge coverage from serveral parsed sets of data
    """
    all_covs = w_index_covs + [wo_index_cov,]
    raw_read_coverage = defaultdict(list)
    for chrm_strand in set([cs for d_cov in all_covs
                            for cs in d_cov.keys()]):
        for dir_cov in all_covs:
            if chrm_strand not in dir_cov: continue
            raw_read_coverage[chrm_strand].extend(dir_cov[chrm_strand])

    return dict(raw_read_coverage)

def parse_fast5s(fast5_basedirs, corrected_group, basecall_subgroups,
                 new_corr_grp=None, rna=False):
    """
    Parse data from a list of re-squiggle fast5 directories
    """
    if VERBOSE: sys.stderr.write('Parsing tombo index file(s).\n')
    wo_index_dirs = []
    w_index_covs = []
    warn_index = False
    # determine if index exists for each directory and load appropriately
    for fast5s_dir in fast5_basedirs:
        fast5s_dir = (fast5s_dir if fast5s_dir.endswith('/') else
                      fast5s_dir + '/')
        if os.path.exists(get_index_fn(fast5s_dir, corrected_group)):
            try:
                w_index_covs.append(parse_fast5s_w_index(
                    fast5s_dir, corrected_group, basecall_subgroups,
                    new_corr_grp))
            except:
                sys.stderr.write('WARNING: Failed to parse tombo index ' +
                                 'file for ' + fast5s_dir + ' directory.\n')
                wo_index_dirs.append(fast5s_dir)
        else:
            if not warn_index:
                sys.stderr.write(
                    'WARNING: Index file does not exist for one or more ' +
                    'directories. For optimal performance, either re-run ' +
                    're-squiggle without --skip-index flag or point to ' +
                    'top level fast5 directory of recursive directories.\n')
                warn_index = True
            wo_index_dirs.append(fast5s_dir)
    wo_index_cov = parse_fast5s_wo_index(
        wo_index_dirs, corrected_group, basecall_subgroups, rna)
    raw_read_coverage = merge_cov(w_index_covs, wo_index_cov)

    return raw_read_coverage


###################################
###### Statistical Functions ######
###################################

# Some of these functions should likely be re-factored to tombo_stats

def parse_pore_model(pore_model_fn):
    """
    Parse pore model for pA normalization (Deprecated)
    """
    pore_model = {'mean':{}, 'inv_var':{}}
    with open(pore_model_fn) as fp:
        for line in fp:
            if line.startswith('#'): continue
            try:
                kmer, lev_mean, lev_stdev = line.split()[:3]
                lev_mean, lev_stdev = map(float, (lev_mean, lev_stdev))
            except ValueError:
                # header or other non-kmer field
                continue
            pore_model['mean'][kmer] = lev_mean
            pore_model['inv_var'][kmer] = 1 / (lev_stdev * lev_stdev)

    return pore_model

def calc_kmer_fitted_shift_scale(pore_model, events_means, events_kmers):
    """
    Compute fitted shift and scale parameters for pA normalization
    """
    r_model_means = np.array([pore_model['mean'][kmer]
                              for kmer in events_kmers])
    r_model_inv_vars = np.array([pore_model['inv_var'][kmer]
                                 for kmer in events_kmers])
    model_mean_var = r_model_means * r_model_inv_vars
    # prep kmer model coefficient matrix for the k-mers from this read
    model_mean_var_sum = model_mean_var.sum()
    coef_mat = np.array((
        (r_model_inv_vars.sum(), model_mean_var_sum),
        (model_mean_var_sum, (model_mean_var * r_model_means).sum())))

    # prep dependent values from this reads true events
    r_event_var = events_means * r_model_inv_vars
    r_event_var_mean = r_event_var * r_model_means
    dep_vect = np.array((r_event_var.sum(), r_event_var_mean.sum()))

    shift, scale = np.linalg.solve(coef_mat, dep_vect)

    return shift, scale

def get_valid_cpts(norm_signal, min_base_obs, num_events):
    """
    Get valid changepoints given largest differences in neighboring moving windows

    Note that this method is completely vectorized, but allows segments
    as small as 2 observations. This should be okay R9+, but is problematic
    for <=R7 and RNA
    """
    raw_cumsum = np.cumsum(np.concatenate([[0], norm_signal[:-1]]))
    # get difference between all neighboring min_base_obs regions
    running_diffs = np.abs(
        (2 * raw_cumsum[min_base_obs:-min_base_obs]) -
        raw_cumsum[:-2*min_base_obs] -
        raw_cumsum[2*min_base_obs:])
    not_peaks = np.logical_not(np.logical_and(
        running_diffs > np.concatenate([[0], running_diffs[:-1]]),
        running_diffs > np.concatenate([running_diffs[1:], [0]])))
    running_diffs[not_peaks] = 0
    valid_cpts = np.argsort(
        running_diffs)[::-1][:num_events].astype(np.int32) + min_base_obs

    return valid_cpts

def estimate_global_scale(fast5_fns, num_reads=500):
    sys.stderr.write('Estimating global scale parameter\n')
    np.random.shuffle(fast5_fns)
    read_mads = []
    for fast5_fn in fast5_fns:
        try:
            with h5py.File(fast5_fn, 'r') as fast5_data:
                all_sig = fast5_data['/Raw/Reads'].values()[0]['Signal'].value
            shift = np.median(all_sig)
            read_mads.append(np.median(np.abs(all_sig - shift)))
        except:
            continue
        if len(read_mads) >= num_reads:
            break

    if len(read_mads) == 0:
        sys.stderr.write(
            '******** ERROR *********\n\tNo reads contain raw signal for ' +
            'global scale parameter estimation.\n')
        sys.exit()
    if len(read_mads) < num_reads:
        sys.stderr.write(
            '******** WARNING *********\n\tFew reads contain raw signal for ' +
            'global scale parameter estimation. Results may not be optimal.\n')

    return np.mean(read_mads)

def normalize_raw_signal(
        all_raw_signal, read_start_rel_to_raw, read_obs_len,
        norm_type=None, channel_info=None, outlier_thresh=None,
        shift=None, scale=None, lower_lim=None, upper_lim=None,
        pore_model=None, event_means=None, event_kmers=None,
        const_scale=None):
    """
    Apply scaling and windsorizing parameters to normalize raw signal
    """
    if norm_type not in NORM_TYPES and (shift is None or scale is None):
        raise NotImplementedError, (
            'Normalization type ' + norm_type + ' is not a valid ' +
            'option and shift or scale parameters were not provided.')

    raw_signal = np.array(all_raw_signal[
        read_start_rel_to_raw:
        read_start_rel_to_raw + read_obs_len])
    if shift is None or scale is None:
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
                # nanocorr/nanopolish/ONT
                fit_shift, fit_scale = calc_kmer_fitted_shift_scale(
                    pore_model, event_means, event_kmers)
                # apply shift and scale values fitted from kmer
                # conditional model after raw DAC scaling
                shift = shift + (fit_shift * scale)
                scale = scale * fit_scale
            # print fitted shift and scale for comparisons
            #print 'shift: ' + str(fit_shift) + \
            #  '\tscale: ' + str(fit_scale)
        elif norm_type == 'median':
            shift = np.median(raw_signal)
            scale = np.median(np.abs(raw_signal - shift))
        elif norm_type == 'median_const_scale':
            assert const_scale is not None
            shift = np.median(raw_signal)
            scale = const_scale
        elif norm_type == 'robust_median':
            shift = np.mean(np.percentile(
                raw_signal, robust_quantiles))
            scale = np.median(np.abs(raw_signal - read_robust_med))

    raw_signal = (raw_signal - shift) / scale

    if outlier_thresh is not None or (
            lower_lim is not None and upper_lim is not None):
        if outlier_thresh is not None:
            read_med = np.median(raw_signal)
            read_mad = np.median(np.abs(raw_signal - read_med))
            lower_lim = read_med - (read_mad * outlier_thresh)
            upper_lim = read_med + (read_mad * outlier_thresh)
        raw_signal = c_apply_outlier_thresh(raw_signal, lower_lim, upper_lim)

    return raw_signal, scaleValues(shift, scale, lower_lim, upper_lim)


###########################################
###### Events Table Access Functions ######
###########################################

def get_multiple_slots_read_centric(r_data, slot_names):
    """
    Extract read-centric slot_names from this read's Events table
    """
    try:
        with h5py.File(r_data.fn, 'r') as read_data:
            # note that it's more efficient to try to access the slot
            # and except the error that check if the slot exists first
            r_events = read_data['/'.join((
                '/Analyses', r_data.corr_group, 'Events'))].value
    except:
        # probably truncated file or events don't exist
        return [None,] * len(slot_names)

    return [r_events[slot_name] for slot_name in slot_names]

def get_single_slot_read_centric(r_data, slot_name):
    """
    Extract read-centric slot_name from this read's Events table
    """
    try:
        with h5py.File(r_data.fn, 'r') as read_data:
            # note that it's more efficient to try to access the slot
            # and except the error that check if the slot exists first
            r_events = read_data['/'.join((
                '/Analyses', r_data.corr_group, 'Events'))].value
            r_slot_values = r_events[slot_name]
    except:
        # probably truncated file or events don't exist
        return None

    return r_slot_values

def get_single_slot_genome_centric(r_data, slot_name):
    """
    Extract genome-centric slot_name from this read's Events table
    """
    r_slot_values = get_single_slot_read_centric(r_data, slot_name)
    if r_slot_values is None:
        return None

    if ((r_data.strand == '-' and not r_data.rna) or
        (r_data.strand == '+' and r_data.rna)):
        r_slot_values = r_slot_values[::-1]

    return r_slot_values

def get_mean_slot_genome_centric(cs_reads, chrm_len, slot_name):
    """
    Get the mean over all reads at each covered genomic location for this slots value
    """
    base_sums = np.zeros(chrm_len)
    base_cov = np.zeros(chrm_len, dtype=np.int_)
    for r_data in cs_reads:
        # extract read means data so data across all chrms is not
        # in RAM at one time
        g_slot_values = get_single_slot_genome_centric(r_data, slot_name)
        if g_slot_values is None: continue

        base_sums[r_data.start:
                  r_data.start + len(g_slot_values)] += g_slot_values
        base_cov[r_data.start:r_data.start + len(g_slot_values)] += 1

    return base_sums / base_cov

def get_all_mean_slot_values(raw_read_coverage, chrm_sizes, slot_name):
    """
    Get the mean over all reads at each covered genomic location for this slots value over all covered chromosomes and strands
    """
    # ignore divide by zero errors that occur where there is no
    # coverage. Need to correct nan values after subtracting two sets of
    # coverage so leave as nan for now
    old_err_settings = np.seterr(all='ignore')
    # take the mean over all signal overlapping each base
    all_mean_values = {}
    for chrm, strand in [(c, s) for c in chrm_sizes.keys()
                         for s in ('+', '-')]:
        if (chrm, strand) in raw_read_coverage:
            cs_mean_values = get_mean_slot_genome_centric(
                raw_read_coverage[(chrm, strand)], chrm_sizes[chrm], slot_name)
        else:
            cs_mean_values = np.empty(chrm_sizes[chrm])
            cs_mean_values[:] = np.nan
        all_mean_values[(chrm, strand)] = cs_mean_values
    _ = np.seterr(**old_err_settings)

    return all_mean_values

def get_all_mean_levels(raw_read_coverage, chrm_sizes):
    """
    Helper function to compute genome location mean levels
    """
    return get_all_mean_slot_values(raw_read_coverage, chrm_sizes, 'norm_mean')

def get_all_mean_stdev(raw_read_coverage, chrm_sizes):
    """
    Helper function to compute genome location mean levels
    """
    return get_all_mean_slot_values(raw_read_coverage, chrm_sizes, 'norm_stdev')

def get_all_mean_lengths(raw_read_coverage, chrm_sizes):
    """
    Helper function to compute genome location mean lengths
    """
    return get_all_mean_slot_values(raw_read_coverage, chrm_sizes, 'length')


###########################################
###### Special Data Access Functions ######
###########################################

def get_channel_info(read_data):
    """
    Get channel information for a read
    """
    try:
        fast5_info = read_data['UniqueGlobalKey/channel_id'].attrs
    except:
        raise RuntimeError, ("No channel_id group in HDF5 file. " +
                             "Probably mux scan HDF5 file.")

    channel_info = channelInfo(
        fast5_info['offset'], fast5_info['range'],
        fast5_info['digitisation'], fast5_info['channel_number'],
        fast5_info['sampling_rate'].astype('int_'))

    return channel_info

def get_raw_signal(r_data, int_start, int_end):
    """
    Extract raw signal from where this read overlaps a particular genomic region
    """
    with h5py.File(r_data.fn, 'r') as fast5_data:
        # retrieve shift and scale computed in correction script
        corr_subgrp = fast5_data['/Analyses/' + r_data.corr_group]
        event_data = corr_subgrp['Events'].value
        events_end = event_data[-1]['start'] + event_data[-1]['length']
        segs = np.concatenate([event_data['start'], [events_end,]])

        shift = corr_subgrp.attrs['shift']
        scale = corr_subgrp.attrs['scale']
        lower_lim = corr_subgrp.attrs['lower_lim']
        upper_lim = corr_subgrp.attrs['upper_lim']
        all_sig = fast5_data['/Raw/Reads'].values()[0]['Signal'].value

    rsrtr = r_data.read_start_rel_to_raw
    if r_data.rna:
        # reverse raw signal for RNA
        all_sig = all_sig[::-1]

    if r_data.strand == "-":
        segs = (segs[::-1] * -1) + segs[-1]

    if int_start < r_data.start:
        # handle reads that start in the middle of the interval
        start_offset = r_data.start - int_start
        overlap_seg_data = segs[:int_end - r_data.start + 1]
    else:
        start_offset = 0
        skipped_bases = int_start - r_data.start
        overlap_seg_data = segs[
            skipped_bases:int_end - r_data.start + 1]

    num_reg_obs = overlap_seg_data[-1] - overlap_seg_data[0]
    if r_data.strand == "+":
        reg_start_rel_raw = rsrtr + overlap_seg_data[0]
        r_sig, _ = normalize_raw_signal(
            all_sig, reg_start_rel_raw, num_reg_obs, shift=shift,
            scale=scale, lower_lim=lower_lim, upper_lim=upper_lim)
    else:
        reg_start_rel_raw = rsrtr + segs[-1] - overlap_seg_data[-1]
        r_sig, _ = normalize_raw_signal(
            all_sig, reg_start_rel_raw, num_reg_obs, shift=shift,
            scale=scale, lower_lim=lower_lim, upper_lim=upper_lim)
        r_sig = r_sig[::-1]

    return r_sig, overlap_seg_data, start_offset

def parse_read_correction_data(r_data):
    """
    Parse correction data from an event resquiggled read
    """
    try:
        with h5py.File(r_data.fn, 'r') as fast5_data:
            raw_grp = fast5_data['/Raw/Reads'].values()[0]
            corr_grp = fast5_data['/Analyses/' + r_data.corr_group]
            events_grp = corr_grp['Events']
            event_data = events_grp.value

            read_id = raw_grp.attrs['read_id']
            signal_data = raw_grp['Signal'].value

            raw_offset = events_grp.attrs['read_start_rel_to_raw']
            shift, scale, lower_lim, upper_lim = [
                corr_grp.attrs[attr_name] for attr_name in (
                    'shift', 'scale', 'lower_lim', 'upper_lim')]

            old_segs = corr_grp['Alignment/read_segments'].value
            old_align_vals = corr_grp['Alignment/read_alignment'].value
            new_align_vals = corr_grp['Alignment/genome_alignment'].value
    except:
        return None

    events_end = event_data['start'][-1] + event_data['length'][-1]
    new_segs = np.concatenate([event_data['start'], [events_end,]])

    if r_data.rna:
        signal_data = signal_data[::-1]

    return (read_id, signal_data, raw_offset, shift, scale, lower_lim,
            upper_lim, old_segs, old_align_vals, new_align_vals,
            events_end, new_segs)

def get_all_read_data(r_data):
    """
    Extract most relevant read data from this read
    """
    try:
        with h5py.File(r_data.fn, 'r') as read_data:
            # note that it's more efficient to try to access the slot
            # and except the error that check if the slot exists first
            corr_subgrp = read_data['/Analyses/' + r_data.corr_group]
            algn_subgrp = dict(corr_subgrp['Alignment'].attrs.items())
            event_data = corr_subgrp['Events'].value
            r_attrs = dict(corr_subgrp.attrs.items())
            all_sig = read_data['/Raw/Reads'].values()[0]['Signal'].value
    except:
        # probably truncated file or Events slot doesn't exist
        return None

    if r_data.rna:
        all_sig = all_sig[::-1]
    r_means, r_seq = event_data['norm_mean'], event_data['base']

    events_end = event_data[-1]['start'] + event_data[-1]['length']
    segs = np.concatenate([event_data['start'], [events_end,]]).astype(np.int32)
    r_sig, scale_values = normalize_raw_signal(
        all_sig, r_data.read_start_rel_to_raw, segs[-1] - segs[0],
        shift=r_attrs['shift'], scale=r_attrs['scale'],
        lower_lim=r_attrs['lower_lim'], upper_lim=r_attrs['upper_lim'])

    return (r_means, r_seq, r_sig, segs, scale_values,
            r_attrs['norm_type'], r_attrs['outlier_threshold'],
            genomeLoc(algn_subgrp['mapped_start'], algn_subgrp['mapped_strand'],
                      algn_subgrp['mapped_chrom']))

def get_coverage(raw_read_coverage):
    """
    Get genome coverage for a set of reads
    """
    if VERBOSE: sys.stderr.write('Calculating read coverage.\n')
    read_coverage = {}
    for (chrm, strand), reads_data in raw_read_coverage.items():
        max_end = max(r_data.end for r_data in reads_data)
        chrm_coverage = np.zeros(max_end, dtype=np.int_)
        for r_data in reads_data:
            chrm_coverage[r_data.start:r_data.end] += 1
        read_coverage[(chrm, strand)] = chrm_coverage

    return read_coverage

def get_reads_events(cs_reads, rev_strand):
    """
    Extract read base levels split by genomic position
    """
    # note that this function assumes that all reads come from the same
    # chromosome and strand
    cs_base_means = []
    cs_read_start_ends = []
    for r_data in cs_reads:
        # extract read means data so data across all chrms is not
        # in RAM at one time
        read_means = get_single_slot_genome_centric(r_data, 'norm_mean')
        if read_means is None: continue
        assert read_means.shape[0] == r_data.end - r_data.start, (
            'Read found with mismatching mapping location and ' +
            'signal information.')
        cs_base_means.append(read_means)
        cs_read_start_ends.append((r_data.start, r_data.end))

    if len(cs_base_means) == 0: return None

    chrm_signal = np.concatenate(cs_base_means)
    chrm_pos = np.concatenate([np.arange(start, end)
                               for start, end in cs_read_start_ends])
    # get order of all bases from position array
    as_chrm_pos = np.argsort(chrm_pos)
    # then sort the signal array by genomic position and
    # split into event means by base
    split_poss = np.where(
        np.concatenate([[0,], np.diff(
            chrm_pos[as_chrm_pos])]) > 0)[0]
    cs_base_events = dict(zip(
        np.unique(chrm_pos[as_chrm_pos]),
        np.split(chrm_signal[as_chrm_pos], split_poss)))

    return cs_base_events

def update_seq(r_data, reg_base_data, int_start, int_end):
    """
    Update the sequence for the region based on this read
    """
    r_seq = ''.join(get_single_slot_read_centric(r_data, 'base'))
    if r_seq is None:
        # probably a corrupt file so return that the region is only
        # up to the start of this read so the next valid read will be added
        # to the sequence
        return reg_base_data, max(0, r_data.start - int_start)
    if r_data.strand == '-':
        r_seq = rev_comp(r_seq)

    # if read starts before the interval
    if r_data.start <= int_start:
        r_end_overlap = r_data.end - int_start
        # if read covers the whole interval
        if r_data.end > int_end:
            r_end_clip = r_data.end - int_end
            reg_base_data = r_seq[-r_end_overlap:-r_end_clip]
            return reg_base_data, len(reg_base_data)
        # end of read overlaps beginning of interval
        reg_base_data[:r_end_overlap] = r_seq[-r_end_overlap:]
        return reg_base_data, r_end_overlap
    # read doesn't cover the beginning of region
    if r_data.end > int_end:
        # beginning of read covers to the end of the region
        r_begin_overlap = int_end - r_data.start
        reg_base_data[-r_begin_overlap:] = r_seq[:r_begin_overlap]
        return reg_base_data, len(reg_base_data)
    # first read is completely contained in the interval
    r_len = r_data.end - r_data.start
    r_int_start = r_data.start - int_start
    reg_base_data[r_int_start:r_int_start + r_len] = r_seq
    return reg_base_data, r_int_start + r_len

def get_seq_from_reads(int_start, int_end, reg_reads):
    """
    Extract the forward strand genomic sequence for an interval from a set of reads
    """
    # handle case where no read overlaps whole region
    # let each read contibute its sequence and fill the rest
    # with dashes
    reg_base_data = ['-'] * (int_end - int_start)
    if len(reg_reads) == 0:
        return ''.join(reg_base_data)
    # get region sequence by moving through reads that
    # cover the region, but only extract seqeunce from the
    # (close to) minimal reads
    s_reg_reads = sorted(reg_reads, key=lambda r: (r.start, r.end))
    # begin filling sequence with first (by start pos) read
    reg_base_data, curr_cov_pos = update_seq(
        s_reg_reads.pop(0), reg_base_data, int_start, int_end)
    # if there was only one read return now
    if len(s_reg_reads) == 0 or curr_cov_pos >= int_end - int_start:
        return ''.join(reg_base_data)

    # get next read (by start pos)
    curr_read = s_reg_reads.pop(0)
    for next_read in s_reg_reads:
        # once the next read start passes the region covered thus far
        # add the sequence for the saved curr_read to the reg sequence
        if next_read.start >= curr_cov_pos:
            # add read with curr longest end position to the region seq
            reg_base_data, curr_cov_pos = update_seq(
                curr_read, reg_base_data, int_start, int_end)
            curr_read = next_read
            # if the whole interval is covered return the sequence
            if curr_cov_pos >= int_end - int_start:
                return ''.join(reg_base_data)
            continue
        if next_read.end > curr_read.end:
            curr_read = next_read

    reg_base_data, _ = update_seq(
        curr_read, reg_base_data, int_start, int_end)

    return ''.join(reg_base_data)

def add_reg_seq(all_reg_data):
    """
    Add the region sequence to the region data by extraction from a minimal set of reads
    """
    all_reg_base_data = []
    for reg_data in all_reg_data:
        # add text to each regions data
        all_reg_base_data.append(intervalData(
            reg_data.reg_id, reg_data.chrm, reg_data.start, reg_data.end,
            reg_data.strand, reg_data.reg_text, reg_data.reads,
            get_seq_from_reads(reg_data.start, reg_data.end, reg_data.reads)))

    return all_reg_base_data

def get_region_reads(
        plot_intervals, raw_read_coverage, filter_no_cov=True, add_seq=True):
    """
    Get all reads overlapping a set of intervals
    """
    def get_c_s_data(chrm, strand, start, end):
        # get all reads intersecting the interval
        if (chrm, strand) in raw_read_coverage:
            return [
                r_data for r_data in raw_read_coverage[(chrm, strand)]
                if not (r_data.start >= end or r_data.end <= start)]
        return []


    all_reg_data = []
    for int_i in plot_intervals:
        # get all reads that overlap this interval
        # note that this includes partial overlaps as these contribute
        # to coverage and other statistics so can't really restrict to
        # full coverage as previous versions of code did
        if int_i.strand is None:
            # if strand is None, get data from both strands
            all_reg_data.append(intervalData(
                int_i.reg_id, int_i.chrm, int_i.start, int_i.end, int_i.strand,
                int_i.reg_text,
                get_c_s_data(int_i.chrm, '+', int_i.start, int_i.end) +
                get_c_s_data(int_i.chrm, '-', int_i.start, int_i.end)))
        else:
            all_reg_data.append(intervalData(
                int_i.reg_id, int_i.chrm, int_i.start, int_i.end, int_i.strand,
                int_i.reg_text,
                get_c_s_data(int_i.chrm, int_i.strand, int_i.start, int_i.end)))

    if add_seq:
        all_reg_data = add_reg_seq(all_reg_data)
    if not filter_no_cov:
        return all_reg_data

    # filter out no coverage regions
    all_reg_data = [
        reg_data for reg_data in all_reg_data if len(reg_data.reads) > 0]

    no_cov_regions = [
        (len(reg_data.reads) == 0, str(reg_data.chrm) + ':' + str(reg_data.start))
        for reg_data in all_reg_data]
    if any(no_cov[0] for no_cov in no_cov_regions):
        sys.stderr.write(
            '**** WARNING **** No coverage in regions: ' +
            '; '.join([reg for no_cov, reg in no_cov_regions
                       if no_cov]) + '\n')

    return all_reg_data

def get_region_sequences(
        plot_intervals, raw_read_coverage1, raw_read_coverage2=None):
    """
    Get the sequence for a set of intervals from a set of reads
    """
    all_reg_data = get_region_reads(
        plot_intervals, raw_read_coverage1, filter_no_cov=False, add_seq=False)
    if raw_read_coverage2 is not None:
        all_reg_data2 = get_region_reads(
            plot_intervals, raw_read_coverage2, filter_no_cov=False,
            add_seq=False)
        all_reg_data = [
            intervalData(r1.reg_id, r1.chrm, r1.start, r1.end, r1.strand,
                         r1.reg_text, r1.reads + r2.reads)
            for r1, r2 in zip(all_reg_data, all_reg_data2)]
    all_reg_data = add_reg_seq(all_reg_data)

    return all_reg_data


###################################
###### FAST5 Write Functions ######
###################################

def prep_fast5(fast5_fn, corr_grp, overwrite, in_place, bc_grp=None):
    """
    Prepare a read for re-squiggle processing (This deletes old info for this reads
    """
    # several checks to prepare the FAST5 file for correction before
    # processing to save compute
    if not in_place:
        return ('Not currently implementing new hdf5 file writing',
                fast5_fn)
    # check that the file is writeable before trying to correct
    if not os.access(fast5_fn, os.W_OK):
        return ('FAST5 file is not writable', fast5_fn)

    try:
        # create group to store data
        with h5py.File(fast5_fn, 'r+') as fast5_data:
            try:
                analyses_grp = fast5_data['/Analyses']
            except:
                return 'Analysis group not found at root of FAST5', fast5_fn
            try:
                # check that the requested basecalls group exsists
                if bc_grp is not None:
                    _ = analyses_grp[bc_grp]
            except:
                return 'Basecall group not found at [--basecall-group]', fast5_fn

            try:
                corr_grp_ptr = analyses_grp[corr_grp]
                if not overwrite:
                    return (
                        "Tombo data exsists in [--corrected-group] and " +
                        "[--overwrite] is not set", fast5_fn)
                del analyses_grp[corr_grp]
            except:
                # if the corr_grp isn't there we will write it now, but
                # it's more efficient to try than to check if the slot is there
                pass

            corr_grp = analyses_grp.create_group(corr_grp)
            corr_grp.attrs['tombo_version'] = TOMBO_VERSION
            corr_grp.attrs['basecall_group'] = bc_grp
    except:
        return (
            'Error opening or writing to fast5 file', fast5_fn)

    return

def write_error_status(
        filename, corrected_group, basecall_subgroup, error_text):
    """
    Write error message for a read into the FAST5 file
    """
    fast5_data = h5py.File(filename, 'r+')
    analysis_grp = fast5_data['/Analyses']
    corr_grp = analysis_grp[corrected_group]
    if basecall_subgroup is not None:
        # add subgroup matching subgroup from original basecalls
        corr_subgrp = corr_grp.create_group(basecall_subgroup)
        corr_subgrp.attrs['status'] = error_text
    else:
        corr_grp.attrs['status'] = error_text

    return

def write_new_fast5_group(
        filename, genome_location, read_start_rel_to_raw,
        new_segs, align_seq, norm_signal, scale_values, corrected_group,
        basecall_subgroup, norm_type, outlier_thresh, compute_sd,
        alignVals=None, align_info=None, old_segs=None, rna=False):
    """
    Write new fast5 group with re-squiggle data
    """
    try:
        # compute event data before accessing fast5 file
        if compute_sd:
            norm_means, norm_stds = c_new_mean_stds(norm_signal, new_segs)
        else:
            norm_means = c_new_means(norm_signal, new_segs)
            norm_stds = repeat(np.NAN)

        event_data = np.array(
            zip(norm_means, norm_stds,
                new_segs[:-1], np.diff(new_segs), list(align_seq)),
            dtype=[('norm_mean', '<f8'), ('norm_stdev', '<f8'),
                   ('start', '<u4'), ('length', '<u4'), ('base', 'S1')])

        if alignVals is not None:
            np_read_align = np.chararray(len(alignVals))
            np_read_align[:] = zip(*alignVals)[0]
            np_genome_align = np.chararray(len(alignVals))
            np_genome_align[:] = zip(*alignVals)[1]
    except:
        raise
        raise NotImplementedError, 'Error computing new events.'

    try:
        fast5_data = h5py.File(filename, 'r+')
    except:
        raise NotImplementedError, (
            'Error opening file for new group writing. This should ' +
            'have been caught during the alignment phase. Check that ' +
            'there are no other tombo processes or processes ' +
            'accessing these HDF5 files running simultaneously.')

    try:
        analysis_grp = fast5_data['/Analyses']
        corr_grp = analysis_grp[corrected_group]
        # add subgroup matching subgroup from original basecalls
        corr_subgrp = corr_grp.create_group(basecall_subgroup)
        corr_subgrp.attrs['status'] = 'success'
        corr_subgrp.attrs['rna'] = rna
        corr_subgrp.attrs['shift'] = scale_values.shift
        corr_subgrp.attrs['scale'] = scale_values.scale
        corr_subgrp.attrs['lower_lim'] = scale_values.lower_lim
        corr_subgrp.attrs['upper_lim'] = scale_values.upper_lim
        corr_subgrp.attrs['norm_type'] = norm_type
        corr_subgrp.attrs['outlier_threshold'] = outlier_thresh

        # store alignment statistics
        corr_alignment = corr_subgrp.create_group('Alignment')
        corr_alignment.attrs['mapped_start'] = genome_location.Start
        corr_alignment.attrs['mapped_end'] \
            = genome_location.Start + len(new_segs) - 1
        corr_alignment.attrs['mapped_strand'] = genome_location.Strand
        corr_alignment.attrs['mapped_chrom'] = genome_location.Chrom

        if align_info is not None:
            corr_alignment.attrs['clipped_bases_start'] = align_info.ClipStart
            corr_alignment.attrs['clipped_bases_end'] = align_info.ClipEnd
            corr_alignment.attrs['num_insertions'] = align_info.Insertions
            corr_alignment.attrs['num_deletions'] = align_info.Deletions
            corr_alignment.attrs['num_matches'] = align_info.Matches
            corr_alignment.attrs['num_mismatches'] = align_info.Mismatches

        if alignVals is not None:
            corr_alignment.create_dataset(
                'read_alignment', data=np_read_align, compression="gzip")
            corr_alignment.create_dataset(
                'genome_alignment', data=np_genome_align, compression="gzip")
        if old_segs is not None:
            # store old segmentation in order to plot "correction process"
            corr_alignment.create_dataset(
                'read_segments', data=old_segs, compression="gzip")

        # Add Events to data frame with event means, SDs and lengths
        corr_events = corr_subgrp.create_dataset(
            'Events', data=event_data, compression="gzip")
        corr_events.attrs[
            'read_start_rel_to_raw'] = read_start_rel_to_raw
    except:
        raise NotImplementedError, (
            'Error writing resquiggle information back into fast5 file.')

    try:
        fast5_data.flush()
        fast5_data.close()
    except:
        raise NotImplementedError, (
            'Error closing fast5 file after writing resquiggle ' +
            'information.')

    return


###################################
###### Filter Main Functions ######
###################################

def clear_filters_main(args):
    global VERBOSE
    VERBOSE = not args.quiet

    for fast5s_dir in args.fast5_basedirs:
        clear_filters(fast5s_dir, args.corrected_group)

    return

def filter_stuck_main(args):
    global VERBOSE
    VERBOSE = not args.quiet

    obs_filter = parse_obs_filter(args.obs_per_base_filter)
    for fast5s_dir in args.fast5_basedirs:
        filter_reads(fast5s_dir, args.corrected_group, obs_filter)

    return

def filter_coverage_main(args):
    global VERBOSE
    VERBOSE = not args.quiet

    if not 0 < args.percent_to_filter < 100:
        sys.stderr.write(
            '********** ERROR ********\n\t--percent-to-filter must be between ' +
            '0 and 100.\n')
        sys.exit()

    for fast5s_dir in args.fast5_basedirs:
        filter_reads_for_coverage(
            fast5s_dir, args.corrected_group, args.percent_to_filter / 100.0)

    return


##################################
###### Annotate FAST5s Main ######
##################################

def annotate_reads_with_fastq_main(args):
    global VERBOSE
    VERBOSE = not args.quiet

    fast5s_basedir = (
        args.fast5_basedir if args.fast5_basedir.endswith('/') else
        args.fast5_basedir + '/')
    fast5s_read_ids = get_read_ids_and_prep_fastq_slot(
        fast5s_basedir, args.basecall_group, args.basecall_subgroup,
        args.overwrite)
    fastq_slot = '/'.join(('/Analyses', args.basecall_group,
                           args.basecall_subgroup))
    annotate_with_fastqs(args.fastq_filenames, fast5s_read_ids, fastq_slot)

    return


if __name__ == '__main__':
    raise NotImplementedError, (
        'This is a module. See commands with `tombo -h`')
