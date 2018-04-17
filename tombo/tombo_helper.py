from __future__ import division, unicode_literals, absolute_import

from builtins import int, range, dict, map, zip

import os
import io
import re
import sys
import random
import fnmatch

# Future warning from cython in h5py
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
import h5py

import numpy as np

from glob import glob
from operator import itemgetter
from itertools import repeat, islice
from collections import defaultdict, namedtuple

if sys.version_info[0] > 2:
    unicode = str

# import tombo functions
from ._version import TOMBO_VERSION
from .c_helper import c_new_mean_stds, c_new_means, c_apply_outlier_thresh
from ._default_parameters import ROBUST_QUANTS, NUM_READS_FOR_SCALE

VERBOSE = False


################################
###### Global Namedtuples ######
################################

alignInfo = namedtuple(
    'alignInfo',
    ('ID', 'Subgroup', 'ClipStart', 'ClipEnd',
     'Insertions', 'Deletions', 'Matches', 'Mismatches'))

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

# single base conversion for motifs
SINGLE_LETTER_CODE = {
    'A':'A', 'C':'C', 'G':'G', 'T':'T', 'B':'[CGT]',
    'D':'[AGT]', 'H':'[ACT]', 'K':'[GT]', 'M':'[AC]',
    'N':'[ACGT]', 'R':'[AG]', 'S':'[CG]', 'V':'[ACG]',
    'W':'[AT]', 'Y':'[CT]'}
INVALID_BASES = re.compile('[^ACGT]')


######################################
###### Various Helper Functions ######
######################################

def _warning_message(message):
    sys.stderr.write(
        '*' * 20 + ' WARNING ' + '*' * 20 + '\n\t' +
        message + '\n')
    return

def _error_message_and_exit(message):
    sys.stderr.write(
        '*' * 20 + ' ERROR ' + '*' * 20 + '\n\t' +
        message + '\n')
    sys.exit()
    return

COMP_BASES = dict(zip(map(ord, 'ACGT'), map(ord, 'TGCA')))
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

U_TO_T = {ord('U'):ord('T')}
def rev_transcribe(seq):
    """
    Convert U bases to T
    """
    return seq.translate(U_TO_T)

def get_chrm_sizes(raw_read_coverage, raw_read_coverage2=None):
    """
    Get covered chromosome sizes from a set of reads
    """
    strand_chrm_sizes = defaultdict(list)
    for (chrm, strand), cs_read_cov in \
        raw_read_coverage.items():
        strand_chrm_sizes[chrm].append(max(
            r_data.end for r_data in cs_read_cov))
    if raw_read_coverage2 is not None:
        for (chrm, strand), cs_read_cov in \
            raw_read_coverage2.items():
            strand_chrm_sizes[chrm].append(max(
                r_data.end for r_data in cs_read_cov))

    return dict((chrm, max(strnd_sizes))
                for chrm, strnd_sizes in
                strand_chrm_sizes.items())

def parse_genome_locations(genome_locs, default_strand=None):
    parsed_locs = []
    for chrm_pos_strand in genome_locs:
        # strip off any quotes and return up to the first 3 values
        split_vals = chrm_pos_strand.replace('"', '').replace(
            "'", "").split(':')[:3]
        # default to plus strand if not specified
        if len(split_vals) == 1:
            _error_message_and_exit(
                'Invalid genome location provided: ' + chrm_pos_strand +
                '\n\t\tTry adding quotation marks around specified genome ' +
                'locations (especially for sequence identifiers with ' +
                'special characters).')
        elif len(split_vals) == 2:
            parsed_locs.append((
                split_vals[0], split_vals[1], default_strand))
        else:
            parsed_locs.append(split_vals)

    return parsed_locs

class TomboMotif(object):
    def _parse_motif(self, rev_comp_motif=False):
        """
        Parse a single letter code motif into a pattern for matching
        """
        conv_motif = ''.join(SINGLE_LETTER_CODE[letter]
                             for letter in self.raw_motif)
        if rev_comp_motif:
            # reverse complement and then flip any group brackets
            conv_motif = rev_comp(conv_motif).translate({
                ord('['):']', ord(']'):'['})
        return re.compile(conv_motif)

    def __init__(self, raw_motif, mod_pos=None):
        invalid_chars = re.findall(
            '[^' + ''.join(SINGLE_LETTER_CODE) + ']', raw_motif)
        if len(invalid_chars) > 0:
            _error_message_and_exit(
                'Invalid characters in motif: ' + ', '.join(invalid_chars))

        # basic motif parsing
        self.raw_motif = raw_motif
        self.motif_len = len(raw_motif)
        self.motif_pat = self._parse_motif()
        self.rev_comp_pat = self._parse_motif(True)

        self.is_palindrome = self.motif_pat == self.rev_comp_pat

        # parse modified position from motif if provided
        self.mod_pos = mod_pos
        if mod_pos is None:
            self.mod_base = None
        else:
            self.mod_base = raw_motif[mod_pos - 1]
            if INVALID_BASES.match(self.mod_base):
                _warning_message(
                    'Provided modified position is not a single base, which ' +
                    'is likely an error. Specified modified base is one of: ' +
                    ' '.join(SINGLE_LETTER_CODE[self.mod_base][1:-1]))


def invalid_seq(seq):
    return bool(INVALID_BASES.search(seq))


###########################
###### FASTA Parsing ######
###########################

class Fasta(object):
    """
    Fasta sequence format wrapper class.

    Will load faidx via pyfaidx package if installed, else the fasta will be
    loaded into memory for sequence extraction
    """
    def _load_in_mem(self):
        genome_index = {}
        curr_id = None
        curr_seq = []
        with io.open(self.fasta_fn) as fasta_fp:
            for line in fasta_fp:
                if line.startswith('>'):
                    if (curr_id is not None and
                        curr_seq is not []):
                        genome_index[curr_id] = ''.join(curr_seq)
                    curr_seq = []
                    curr_id = line.replace(">","").split()[0]
                else:
                    curr_seq.append(line.strip())
            # add last record
            if (curr_id is not None and
                curr_seq is not []):
                genome_index[curr_id] = ''.join(curr_seq)

        return genome_index

    def __init__(self, fasta_fn, dry_run=False, force_in_mem=False):
        self.fasta_fn = fasta_fn
        try:
            if force_in_mem: raise ImportError
            import pyfaidx
            self.has_pyfaidx = True
            self.index = pyfaidx.Faidx(fasta_fn)
        except:
            self.has_pyfaidx = False
            if not dry_run:
                self.index = self._load_in_mem()

    def get_seq(self, chrm, start=None, end=None):
        if self.has_pyfaidx:
            if not (start or end):
                return self.index.fetch(
                    chrm, 1, self.index.index[chrm].rlen).seq.upper()
            if start < 0 or end > self.index.index[chrm].rlen:
                raise NotImplementedError(
                    'Encountered invalid genome sequence request.')
            return self.index.fetch(chrm, start + 1, end).seq.upper()
        return self.index[chrm][start:end].upper()

    def iter_chrms(self):
        if self.has_pyfaidx:
            for chrm in self.index.index:
                yield unicode(chrm)
        else:
            for chrm in self.index:
                yield chrm

    def __contains__(self, chrm):
        if self.has_pyfaidx:
            return chrm in self.index.index
        return chrm in self.index


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
        try:
            exp_type = exp_type.decode()
        except (AttributeError, TypeError):
            pass
        # remove the word internal since it contains rna.
        exp_type = exp_type.replace('internal', '')
    except:
        pass
    try:
        exp_kit = fast5_data['UniqueGlobalKey/context_tags'].attrs[
            'experiment_kit']
        try:
            exp_kit = exp_kit.decode()
        except (AttributeError, TypeError):
            pass
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
    for cs_reads in raw_read_coverage.values():
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
            with h5py.File(fast5_fn, 'r') as fast5_data:
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
    # if directory comes with trailing slash, remove for processing
    if fast5s_dir.endswith('/'):
        fast5s_dir = fast5s_dir[:-1]
    split_dir = os.path.split(fast5s_dir)
    return os.path.join(split_dir[0], "." + split_dir[1] +
                        "." + corr_grp + '.tombo.index')

def get_lock_fn(fast5s_dir):
    """
    Get filename for the lock file to indicate that this directory
    is currently being processed. This file should be saved to be deleted later.
    """
    # if directory comes with trailing slash, remove for processing
    if fast5s_dir.endswith('/'):
        fast5s_dir = fast5s_dir[:-1]
    split_dir = os.path.split(fast5s_dir)
    return os.path.join(split_dir[0], "." + split_dir[1] + '.tombo.lock')

def _is_lock_file(lock_fn):
    lock_bn = os.path.split(lock_fn)[1]
    if (len(lock_bn) == 0 or
        not lock_bn.startswith('.') or
        not lock_bn.endswith('.tombo.lock')):
        return False
    return True

def prep_index_data(
        fast5_fn, genome_loc, read_start_rel_to_raw, segs,
        corr_grp, subgroup, rna, obs_filter):
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
        corr_grp, subgroup, is_filtered, rna))

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

    with io.open(index_fn, 'wb') as index_fp:
        # note protocol 2 for py2/3 compatibility
        pickle.dump(dict(index_data), index_fp, protocol=2)

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
        with io.open(index_fn, 'rb') as index_fp:
            index_data = pickle.load(index_fp)
    except IOError:
        _error_message_and_exit(
            'Filters can only be applied to runs ' +
            'with a Tombo index file. Re-run resquiggle without the ' +
            '--skip-index option to apply filters.')
    new_index_data = []
    for chrm_strand, cs_raw_data in index_data.items():
        new_index_data.extend([(chrm_strand, (
            from_base_fn, start, end, rsrtr, corr_grp, s_grp, False, rna))
                               for from_base_fn, start, end, rsrtr, c_grp,
                               s_grp, filtered, rna in cs_raw_data])

    write_index_file(new_index_data, index_fn, fast5s_dir)
    sys.stderr.write('All filters successfully cleared!\n')

    return

def parse_obs_filter(obs_filter):
    """
    Parse observations per base formatted filtering
    """
    if len(obs_filter) < 1:
        return None

    # parse obs_filter
    try:
        obs_filter = [list(map(int, pctl_nobs.split(':')))
                      for pctl_nobs in obs_filter]
    except:
        raise RuntimeError('Invalid format for observation filter')

    if any(pctl < 0 or pctl > 100 for pctl in map(itemgetter(0), obs_filter)):
       _error_message_and_exit('Invalid percentile value.')

    return obs_filter

def filter_reads(fast5s_dir, corr_grp, obs_filter):
    """
    Filter reads based on some observation per base threshold criteria
    """
    def read_is_stuck(fast5_fn, s_grp):
        try:
            with h5py.File(fast5_fn, 'r') as fast5_data:
                base_lens = fast5_data['/Analyses/' + corr_grp + '/' + s_grp +
                                       '/Events']['length']
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
        with io.open(index_fn, 'rb') as index_fp:
            index_data = pickle.load(index_fp)
    except IOError:
        sys.stderr.write(
            '******** ERRROR *******\n\tFilters can only be applied to runs ' +
            'with a Tombo index file. Re-run resquiggle without the ' +
            '--skip-index option to apply filters.\n')
    filt_index_data = []
    num_reads, num_filt_reads = 0, 0
    for chrm_strand, cs_raw_data in index_data.items():
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
        'Filtered ' + unicode(num_filt_reads) +
        ' reads due to observations per base filter from a ' +
        'total of ' + unicode(num_reads) + ' reads in ' + fast5s_dir + '.\n')

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
        with io.open(index_fn, 'rb') as index_fp:
            index_data = pickle.load(index_fp)
    except IOError:
        sys.stderr.write(
            '******** ERRROR *******\n\tFilters can only be applied to runs ' +
            'with a Tombo index file. Re-run resquiggle without the ' +
            '--skip-index option to apply filters.\n')
    unfilt_data = []
    unfilt_reads_cov = []
    prev_filt_data = []
    for chrm_strand, cs_raw_data in index_data.items():
        max_end = max(end for (_, _, end, _, _, _, _, _) in cs_raw_data)
        cs_coverage = np.zeros(max_end, dtype=np.int64)
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
            unfilt_reads_cov.append(cs_coverage[start + ((end - start) // 2)])
            unfilt_data.append((chrm_strand, (
                from_base_fn, start, end, rsrtr, c_grp, s_grp, filtered, rna)))

    num_reads = len(unfilt_data)
    num_filt_reads = int(frac_to_filter * num_reads)
    sys.stderr.write(
        'Filtered ' + unicode(num_filt_reads) +
        ' reads due to observations per base filter from a ' +
        'total of ' + unicode(num_reads) + ' reads in ' + fast5s_dir + '.\n')

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
        with io.open(fastq_fn) as fastq_fp:
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
                    _warning_message(
                        'Successfully parsed ' + unicode(n_recs) +
                        'FASTQ records from ' + fastq_fn + ' before ' +
                        'encountering an invalid record. The rest of ' +
                        'this file will not be processed.')
                    break

                # extract read_id from fastq (which should be the first text
                # after the "@" record delimiter up to the first white space or
                # underscore
                read_id = fastq_rec[0].split()[0].split('_')[0][1:]
                if read_id not in fast5s_read_ids:
                    if not been_warned_ids:
                        been_warned_ids = True
                        _warning_message(
                            'Some records from ' + fastq_fn + ' contain read ' +
                            'identifiers not found in any FAST5 files.')
                    continue

                with h5py.File(fast5s_read_ids[read_id]) as fast5_data:
                    bc_slot = fast5_data[fastq_slot]
                    bc_slot.create_dataset(
                        'Fastq', data=''.join(fastq_rec),
                        dtype=h5py.special_dtype(vlen=unicode))

    return

def reads_contain_basecalls(fast5_fns, bc_grp, num_reads):
    test_fns = random.sample(
        fast5_fns, num_reads) if len(fast5_fns) > num_reads else fast5_fns
    for fast5_fn in test_fns:
        try:
            with h5py.File(fast5_fn, 'r') as fast5_data:
                fast5_data['/Analyses/' + bc_grp]
        except:
            continue
        # if the basecall group is accessible for a single file return true
        return True

    # else if all tested reads did not contain the basecall group return False
    return False

def get_files_list(fast5s_dir):
    """
    Get all fast5 files recursively below this directory
    """
    all_fast5s = []
    # walk through directory structure searching for fast5 files
    for root, _, fns in os.walk(fast5s_dir):
        for fn in fnmatch.filter(fns, '*.fast5'):
            all_fast5s.append(os.path.join(root, fn))

    return all_fast5s

def clear_tombo_locks(lock_fns):
    """
    Clear all lock files
    """
    for lock_fn in lock_fns:
        # safegaurd against incorrect file passed to this function
        if not _is_lock_file(lock_fn):
            continue
        # try to remove lock files so that is some have been corrupted or
        # otherwise cause an error at least all those accessible will be cleared
        try:
            os.remove(lock_fn)
        except:
            pass

    return

def get_files_list_and_lock_dirs(fast5s_dir, ignore_locks):
    """
    Get all fast5 files recursively below this directory and add a Tombo lock
    file to indicate that this directory is currently being re-squiggled
    """
    all_fast5s = []
    lock_fns = []
    # walk through directory structure searching for fast5 files
    for root, _, fns in os.walk(fast5s_dir):
        lock_fn = get_lock_fn(root)
        if not ignore_locks and os.path.exists(lock_fn):
            clear_tombo_locks(lock_fns)
            _error_message_and_exit(
                'This set of reads is currently being processed by another ' +
                'resquiggle command. Multiple resquiggle commands cannot be ' +
                'run concurrently on a set of reads to avoid corrupting ' +
                'read files. If you are sure this set of reads is not being ' +
                'processed by another command (usually caused by previous ' +
                'unexpected exit) set the --ignore-read-locks flag.')
        lock_fns.append(lock_fn)
        try:
            # create empty file indicating this directory is locked
            open(lock_fn, 'w').close()
        except:
            clear_tombo_locks(lock_fns)
            _error_message_and_exit(
                'Could not write tombo lock file. Check that you have write ' +
                'permission within the specified [fast5_basedir].')

        for fn in fnmatch.filter(fns, '*.fast5'):
            all_fast5s.append(os.path.join(root, fn))

    return all_fast5s, lock_fns

def get_raw_read_slot(fast5_data):
    try:
        raw_read_slot = list(fast5_data['/Raw/Reads'].values())[0]
    except:
        raise NotImplementedError(
            'Raw data is not found in /Raw/Reads/Read_[read#]')

    return raw_read_slot

def prep_fast5_for_fastq(
        fast5_data, basecall_group, basecall_subgroup, overwrite):
    try:
        read_id = get_raw_read_slot(fast5_data).attrs['read_id']
        try:
            read_id = read_id.decode()
        except (AttributeError, TypeError):
            pass
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
            raise NotImplementedError(
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
                        _warning_message(
                            'Basecalls exsit in ' + basecall_group + ' slot. ' +
                            'Set --overwrite option to overwrite these ' +
                            'basecalls in this slot.')
                    continue
            if read_id is None:
                continue
            if read_id in fast5s_read_ids:
                # Warn non-unique read_ids in directory
                if VERBOSE and not been_warned_unique:
                    been_warned_unique = True
                    _warning_message(
                        'Multiple FAST5 files contain the same read ' +
                        'identifiers. Ensure that FAST5 files are from ' +
                        'a single run.')
                continue

            fast5s_read_ids[read_id] = fast5_fn

    return fast5s_read_ids

def parse_fast5s_wo_index(fast5_basedirs, corr_grp, bc_subgrps, rna):
    """
    Parse re-squiggled reads data from a list of fast5 directories
    """
    def get_read_data(read_fn, fast5_data, bc_subgrp):
        corr_data = fast5_data['/'.join(('/Analyses', corr_grp, bc_subgrp))]

        align_data = dict(corr_data['Alignment'].attrs.items())
        read_start_rel_to_raw = corr_data['Events'].attrs[
            'read_start_rel_to_raw']
        chrm = align_data['mapped_chrom']
        strand = align_data['mapped_strand']
        try:
            chrm = chrm.decode()
            strand = strand.decode()
        except:
            pass

        return chrm, strand, readData(
            align_data['mapped_start'], align_data['mapped_end'],
            False, read_start_rel_to_raw, strand, read_fn,
            corr_grp + '/' + bc_subgrp, rna)


    files = [fn for fast5_basedir in fast5_basedirs
             for fn in get_files_list(fast5_basedir)]
    raw_read_coverage = defaultdict(list)
    for read_fn in files:
        try:
            with h5py.File(read_fn, 'r') as fast5_data:
                for bc_subgrp in bc_subgrps:
                    chrm, strand, r_data = get_read_data(
                        read_fn, fast5_data, bc_subgrp)
                    raw_read_coverage[(chrm, strand)].append(r_data)
        except:
            # ignore errors and process all reads that don't error
            continue

    return dict(raw_read_coverage)

def convert_index(index_data, fast5s_dir, corr_grp, new_corr_grp):
    """
    Convert an index and save under a new corrected group. Mostly for
    model_resquiggle
    """
    new_index_data = []
    for (chrm, strand), cs_raw_data in index_data.items():
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
    with io.open(index_fn, 'rb') as index_fp:
        index_data = pickle.load(index_fp)
    raw_read_coverage = {}
    for (chrm, strand), cs_raw_data in index_data.items():
        cs_data = [
            readData(start, end, filtered, rsrtr, strand,
                     os.path.join(fast5s_dir, from_base_fn),
                     corr_grp + '/' + s_grp, rna)
            for from_base_fn, start, end, rsrtr, c_grp,
            s_grp, filtered, rna in cs_raw_data
            if c_grp == corr_grp and s_grp in subgroups and not filtered]
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
    for chrm_strand in set([cs for d_cov in all_covs for cs in d_cov]):
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
                raise
                _warning_message(
                    'Failed to parse tombo index file for ' +
                    fast5s_dir + ' directory.')
                wo_index_dirs.append(fast5s_dir)
        else:
            if not warn_index:
                _warning_message(
                    'Tombo index file does not exist for one or more ' +
                    'directories. If --skip-index was not set for ' +
                    're-squiggle command, ensure that the specified ' +
                    'directory is the same as for the re-squiggle command.\n')
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
    with io.open(pore_model_fn) as fp:
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

def get_valid_cpts(norm_signal, running_stat_width, num_events):
    """
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

def estimate_global_scale(fast5_fns, num_reads=NUM_READS_FOR_SCALE):
    sys.stderr.write('Estimating global scale parameter\n')
    np.random.shuffle(fast5_fns)
    read_mads = []
    for fast5_fn in fast5_fns:
        try:
            with h5py.File(fast5_fn, 'r') as fast5_data:
                all_sig = get_raw_read_slot(fast5_data)['Signal'].value
            shift = np.median(all_sig)
            read_mads.append(np.median(np.abs(all_sig - shift)))
        except:
            continue
        if len(read_mads) >= num_reads:
            break

    if len(read_mads) == 0:
        _error_message_and_exit(
            'No reads contain raw signal for ' +
            'global scale parameter estimation.')
    if len(read_mads) < num_reads:
        _warning_message(
            'Few reads contain raw signal for global scale parameter ' +
            'estimation. Results may not be optimal.')

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
        raise NotImplementedError(
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

def get_multiple_slots_read_centric(r_data, slot_names, corr_grp=None):
    """
    Extract read-centric slot_names from this read's Events table
    """
    try:
        do_close = False
        if not isinstance(r_data, h5py.File):
            do_close = True
            corr_grp = r_data.corr_group
            r_data = h5py.File(r_data.fn, 'r')
        event_slot_name = '/'.join(('/Analyses', corr_grp, 'Events'))
        # note that it's more efficient to try to access the slot
        # and except the error that check if the slot exists first
        r_event_data = r_data[event_slot_name].value
        if do_close: r_data.close()
    except:
        # probably truncated file or events don't exist
        return [None,] * len(slot_names)

    return [r_event_data[slot_name] for slot_name in slot_names]

def get_single_slot_read_centric(r_data, slot_name, corr_grp=None):
    """
    Extract read-centric slot_name from this read's Events table
    """
    try:
        # if r_data is an open h5py object then don't open the filename
        do_close = False
        if not isinstance(r_data, h5py.File):
            do_close = True
            corr_grp = r_data.corr_group
            r_data = h5py.File(r_data.fn, 'r')
        # note that it's more efficient to try to access the slot
        # and except the error that check if the slot exists first
        r_slot_values = r_data['/'.join(('/Analyses', corr_grp, 'Events'))][
            slot_name]
        if do_close: r_data.close()
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

    if r_data.strand == '-':
        r_slot_values = r_slot_values[::-1]

    return r_slot_values

def get_mean_slot_genome_centric(cs_reads, chrm_len, slot_name):
    """
    Get the mean over all reads at each covered genomic location for this
    slots value
    """
    base_sums = np.zeros(chrm_len)
    base_cov = np.zeros(chrm_len, dtype=np.int64)
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
    Get the mean over all reads at each covered genomic location for this
    slots value over all covered chromosomes and strands
    """
    # ignore divide by zero errors that occur where there is no
    # coverage. Need to correct nan values after subtracting two sets of
    # coverage so leave as nan for now
    old_err_settings = np.seterr(all='ignore')
    # take the mean over all signal overlapping each base
    all_mean_values = {}
    for chrm, strand in [(c, s) for c in chrm_sizes for s in ('+', '-')]:
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

def get_channel_info(fast5_data):
    """
    Get channel information for a read
    """
    try:
        fast5_info = fast5_data['UniqueGlobalKey/channel_id'].attrs
    except:
        raise NotImplementedError("No channel_id group in HDF5 file. " +
                                  "Probably mux scan HDF5 file.")

    channel_info = channelInfo(
        fast5_info['offset'], fast5_info['range'],
        fast5_info['digitisation'], fast5_info['channel_number'],
        fast5_info['sampling_rate'].astype(np.int64))

    return channel_info

def get_raw_signal(r_data, int_start, int_end):
    """
    Extract raw signal from where this read overlaps a particular genomic region
    """
    with h5py.File(r_data.fn, 'r') as fast5_data:
        # retrieve shift and scale computed in correction script
        corr_subgrp = fast5_data['/Analyses/' + r_data.corr_group]
        event_starts = corr_subgrp['Events']['start']
        events_end = event_starts[-1] + corr_subgrp['Events']['length'][-1]
        segs = np.concatenate([event_starts, [events_end,]])

        shift = corr_subgrp.attrs['shift']
        scale = corr_subgrp.attrs['scale']
        lower_lim = corr_subgrp.attrs['lower_lim']
        upper_lim = corr_subgrp.attrs['upper_lim']
        all_sig = get_raw_read_slot(fast5_data)['Signal'].value

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
            corr_grp = fast5_data['/Analyses/' + r_data.corr_group]
            events_grp = corr_grp['Events']
            event_starts = events_grp['start']
            events_end = event_starts[-1] + events_grp['length'][-1]
            new_segs = np.concatenate([event_starts, [events_end,]])

            raw_grp = get_raw_read_slot(fast5_data)
            read_id = raw_grp.attrs['read_id']
            try:
                read_id = read_id.decode()
            except (AttributeError, TypeError):
                pass
            signal_data = raw_grp['Signal'].value

            raw_offset = events_grp.attrs['read_start_rel_to_raw']
            shift, scale, lower_lim, upper_lim = [
                corr_grp.attrs[attr_name] for attr_name in (
                    'shift', 'scale', 'lower_lim', 'upper_lim')]

            old_segs = corr_grp['Alignment/read_segments'].value
            old_align_vals = list(map(
                lambda x: x.decode(),
                corr_grp['Alignment/read_alignment'].value))
            new_align_vals = list(map(
                lambda x: x.decode(),
                corr_grp['Alignment/genome_alignment'].value))
    except:
        return None

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
        with h5py.File(r_data.fn, 'r') as fast5_data:
            # note that it's more efficient to try to access the slot
            # and except the error that check if the slot exists first
            corr_subgrp = fast5_data['/Analyses/' + r_data.corr_group]
            algn_subgrp = dict(corr_subgrp['Alignment'].attrs.items())
            event_data = corr_subgrp['Events'].value
            r_attrs = dict(corr_subgrp.attrs.items())
            all_sig = get_raw_read_slot(fast5_data)['Signal'].value
    except:
        # probably truncated file or Events slot doesn't exist
        return None

    if r_data.rna:
        all_sig = all_sig[::-1]
    r_means = event_data['norm_mean']
    r_seq = b''.join(event_data['base']).decode()

    events_end = event_data[-1]['start'] + event_data[-1]['length']
    segs = np.concatenate([event_data['start'], [events_end,]]).astype(np.int64)
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
        if len(reads_data) == 0: continue
        max_end = max(r_data.end for r_data in reads_data)
        chrm_coverage = np.zeros(max_end, dtype=np.int64)
        for r_data in reads_data:
            chrm_coverage[r_data.start:r_data.end] += 1
        read_coverage[(chrm, strand)] = chrm_coverage

    return read_coverage

def get_reads_events(cs_reads):
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
    r_seq = b''.join(get_single_slot_read_centric(r_data, 'base')).decode()
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
    Extract the forward strand genomic sequence for an interval from
    a set of reads
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
    Add the region sequence to the region data by extraction from a minimal
    set of reads
    """
    all_reg_base_data = []
    for reg_data in all_reg_data:
        # add text to each regions data
        all_reg_base_data.append(reg_data._replace(seq=get_seq_from_reads(
            reg_data.start, reg_data.end, reg_data.reads)))

    return all_reg_base_data

def get_region_reads(
        plot_intervals, raw_read_coverage, filter_no_cov=True, add_seq=True):
    """
    Get all reads overlapping a set of intervals
    """
    def get_c_s_data(chrm, strand, start, end):
        # get all reads intersecting the interval
        if (chrm, strand) in raw_read_coverage:
            r_data = raw_read_coverage[(chrm, strand)][0]
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
            all_reg_data.append(int_i._replace(
                reads=get_c_s_data(int_i.chrm, '+', int_i.start, int_i.end) +
                get_c_s_data(int_i.chrm, '-', int_i.start, int_i.end)))
        else:
            all_reg_data.append(int_i._replace(
                reads=get_c_s_data(int_i.chrm, int_i.strand,
                                   int_i.start, int_i.end)))

    if add_seq:
        all_reg_data = add_reg_seq(all_reg_data)
    if not filter_no_cov:
        return all_reg_data

    # filter out no coverage regions
    all_reg_data = [
        reg_data for reg_data in all_reg_data if len(reg_data.reads) > 0]

    no_cov_regions = [
        (len(reg_data.reads) == 0, unicode(reg_data.chrm) + ':' +
         unicode(reg_data.start))
        for reg_data in all_reg_data]
    if any(no_cov[0] for no_cov in no_cov_regions):
        _warning_message(
            'No coverage in regions: ' + '; '.join([
                reg for no_cov, reg in no_cov_regions if no_cov]))

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
        all_reg_data = [r1._replace(reads=r1.reads + r2.reads)
                        for r1, r2 in zip(all_reg_data, all_reg_data2)]
    all_reg_data = add_reg_seq(all_reg_data)

    return all_reg_data


###################################
###### FAST5 Write Functions ######
###################################

def prep_fast5(fast5_fn, corr_grp, overwrite, in_place,
               bc_grp=None, return_fp=False):
    """
    Prepare a read for re-squiggle processing (This deletes old re-squiggle
    info for this read)
    """
    def try_close_prep_err(fast5_data, err_str):
        try:
            fast5_data.close()
        except:
            pass
        return (err_str, fast5_fn)

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
        fast5_data = h5py.File(fast5_fn, 'r+')
        try:
            analyses_grp = fast5_data['/Analyses']
        except:
            return try_close_prep_err(
                fast5_data, 'Analysis group not found at root of FAST5')
        try:
            # check that the requested basecalls group exsists
            if bc_grp is not None:
                analyses_grp[bc_grp]
        except:
            return try_close_prep_err(
                fast5_data, 'Basecall group not found at [--basecall-group]')

        try:
            corr_grp_ptr = analyses_grp[corr_grp]
            if not overwrite:
                return try_close_prep_err(
                    fast5_data, "Tombo data exsists in [--corrected-group] " +
                    "and [--overwrite] is not set")
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

    if return_fp:
        return fast5_data

    try:
        fast5_data.close()
    except:
        return 'Error closing fast5 file', fast5_fn

    return

def write_error_status(
        filename, corrected_group, basecall_subgroup, error_text):
    """
    Write error message for a read into the FAST5 file
    """
    with h5py.File(filename, 'r+') as fast5_data:
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
        fast5_data, genome_location, read_start_rel_to_raw,
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

        # had to shift to names formats numpy array specification due to
        # python2 numpy unicode issues. See discussion here:
        # https://github.com/numpy/numpy/issues/2407
        event_data = np.array(
            list(zip(norm_means, norm_stds,
                     new_segs[:-1], np.diff(new_segs), list(align_seq))),
            dtype=[(str('norm_mean'), 'f8'), (str('norm_stdev'), 'f8'),
                   (str('start'), 'u4'), (str('length'), 'u4'),
                   (str('base'), 'S1')])

        if alignVals is not None:
            r_align_vals, g_align_vals = zip(*alignVals)
            np_read_align = np.chararray(len(alignVals))
            np_read_align[:] = r_align_vals
            np_genome_align = np.chararray(len(alignVals))
            np_genome_align[:] = g_align_vals
    except:
        raise NotImplementedError('Error computing new events')

    do_close = False
    if not isinstance(fast5_data, h5py.File):
        try:
            fast5_data = h5py.File(fast5_data, 'r+')
            do_close = True
        except:
            raise NotImplementedError(
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
        raise NotImplementedError(
            'Error writing resquiggle information back into fast5 file.')

    if do_close:
        try:
            fast5_data.close()
        except:
            raise NotImplementedError(
                'Error closing fast5 file after writing resquiggle information.')

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
        _error_message_and_exit(
            '--percent-to-filter must be between 0 and 100.')

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
    raise NotImplementedError(
        'This is a module. See commands with `tombo -h`')
