from __future__ import division, unicode_literals, absolute_import

from builtins import int, range, dict, map, zip

import os
import io
import re
import sys
import queue
import random

# Future warning from cython in h5py
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
import h5py

import numpy as np

from tqdm import tqdm
from glob import glob
from time import sleep
from time import strftime
from operator import itemgetter
from itertools import repeat, islice
from multiprocessing import Process, Queue
from collections import defaultdict, namedtuple

if sys.version_info[0] > 2:
    unicode = str

# import tombo functions
from ._version import TOMBO_VERSION
from .c_helper import c_new_mean_stds, c_new_means
from ._default_parameters import PHRED_BASE

VERBOSE = False

_ITER_QUEUE_LIMIT = 1000
_PROC_UPDATE_INTERVAL = 100

_MAX_FASTQ_QUEUE_SIZE = 10000
_SEQ_SUMMARY_FN_FIELD = 'filename'
_SEQ_SUMMARY_ID_FIELD = 'read_id'

# warning messages for annotate with fastqs over multiple processes,
# requiring passing warning codes to only print warning once.
_WARN_ID_VAL = 'ids'
_WARN_IO_VAL = 'io'
_WARN_MISMATCH_VAL = 'mismatch'
_WARN_OVRWRT_VAL = 'overwrite'
_WARN_UNIQ_VAL = 'uniq'
_WARN_CODES = (_WARN_ID_VAL, _WARN_IO_VAL, _WARN_MISMATCH_VAL, _WARN_OVRWRT_VAL)
_WARN_CODES_PREP = (_WARN_OVRWRT_VAL, _WARN_UNIQ_VAL)


################################
###### Global Namedtuples ######
################################

alignInfo = namedtuple(
    'alignInfo',
    ('ID', 'Subgroup', 'ClipStart', 'ClipEnd',
     'Insertions', 'Deletions', 'Matches', 'Mismatches'))

readData = namedtuple('readData', (
    'start', 'end', 'filtered', 'read_start_rel_to_raw', 'strand', 'fn',
    'corr_group', 'rna', 'sig_match_score', 'mean_q_score'))
# set default values for sig_match_score and q_score
readData.__new__.__defaults__ = (None, None)

intervalData = namedtuple('intervalData', (
    'reg_id', 'chrm', 'start', 'end', 'strand', 'reg_text', 'reads', 'seq'))
""" itervalData - A Tombo namedtuple containing information about a genomic intervar

.. py:attribute:: reg_id

    Region ID - string type

.. py:attribute:: chrm

    Chromosome name - string type

.. py:attribute:: start

    0-based start position - integer type

.. py:attribute:: end

    1-based (or open interval) end position - integer type

.. py:attribute:: strand

    Interval strand ('+', '-' or None). Default: None - string type

.. py:attribute:: reg_test

    Some text describing a region. Used for plot titles. Default: '' - string type

.. py:attribute:: reads

    A list of readData values. Default: None - list type

.. py:attribute:: seq

    The genomic sequence for a region. Default: None - string type

"""
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

def _status_message(message, indent=False):
    pre_str = '\t' if indent else ''
    sys.stderr.write(pre_str + strftime('[%H:%M:%S] ') + message + '\n')
    sys.stderr.flush()
    return

def _warning_message(message):
    sys.stderr.write(
        '*' * 20 + ' WARNING ' + '*' * 20 + '\n\t' +
        message + '\n')
    sys.stderr.flush()
    return

def _error_message_and_exit(message):
    sys.stderr.write(
        '*' * 20 + ' ERROR ' + '*' * 20 + '\n\t' +
        message + '\n')
    sys.exit()
    return

def resolve_path(fn_path):
    """
    Helper function to resolve relative and linked paths that might
    give other packages problems.
    """
    return os.path.realpath(os.path.expanduser(fn_path))

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
    for (chrm, strand), cs_read_cov in raw_read_coverage.items():
        try:
            strand_chrm_sizes[chrm].append(max(
                    r_data.end for r_data in cs_read_cov))
        except ValueError:
            continue
    if raw_read_coverage2 is not None:
        for (chrm, strand), cs_read_cov in raw_read_coverage2.items():
            try:
                strand_chrm_sizes[chrm].append(max(
                        r_data.end for r_data in cs_read_cov))
            except ValueError:
                continue

    chrm_sizes = {}
    for chrm, strnd_sizes in strand_chrm_sizes.items():
        try:
            chrm_sizes[chrm] = max(strnd_sizes)
        except ValueError:
            continue

    return chrm_sizes

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

def parse_genome_regions(all_regs_text):
    parsed_regs = defaultdict(list)
    include_whole_chrms = set()
    for reg_text in all_regs_text:
        try:
            chrm_reg = reg_text.replace('"', '').replace("'", "").split(':')
            if len(chrm_reg) == 1:
                chrm = chrm_reg[0]
                reg_pos = None
                include_whole_chrms.add(chrm)
            elif len(chrm_reg) == 2:
                chrm, reg_pos = chrm_reg
                reg_pos = list(map(lambda x: int(x.replace(',','')),
                                   reg_pos.split('-')))
            else:
                raise NotImplementedError
        except:
            _error_message_and_exit(
                'Invalid [--include-region] format.')

        parsed_regs[chrm].append(reg_pos)

    parsed_regs = dict(parsed_regs)
    for chrm in include_whole_chrms:
        parsed_regs[chrm] = None

    return parsed_regs

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

    def _index_contains_uridines(self, n_chrms=10, n_bases=1000):
        if self.has_pyfaidx:
            # check first N bases of the first M chrms for U characters
            for chrm in islice(self.index.index, n_chrms):
                if re.search('U', self.get_seq(
                        chrm, 1, n_bases, error_end=False)):
                    return True
        else:
            for chrm in islice(self.index.keys(), n_chrms):
                if re.search('U', self.get_seq(
                        chrm, 1, n_bases, error_end=False)):
                    return True
        return False

    def __init__(self, fasta_fn, dry_run=False, force_in_mem=False,
                 assume_dna_base=False):
        self.fasta_fn = resolve_path(fasta_fn)
        self.has_rna_bases = False
        try:
            if force_in_mem: raise ImportError
            import pyfaidx
            self.has_pyfaidx = True
            try:
                self.index = pyfaidx.Faidx(self.fasta_fn)
            except UnicodeDecodeError:
                _error_message_and_exit(
                    'FASTA file does not appear to be formatted correctly.')
        except:
            self.has_pyfaidx = False
            if not dry_run:
                self.index = self._load_in_mem()
        if not dry_run:
            self.has_rna_bases = (assume_dna_base or
                                  self._index_contains_uridines())

    def get_seq(self, chrm, start=None, end=None, error_end=True):
        """
        Extract sequence from a specific genomic region.

        Note if provided, start and end must both be provided or they will
        be ignored.
        """
        if self.has_pyfaidx:
            if not (start or end):
                r_seq = self.index.fetch(
                    chrm, 1, self.index.index[chrm].rlen).seq.upper()
            elif (start < 0 or start > self.index.index[chrm].rlen or (
                    error_end and (
                        end < 0 or end > self.index.index[chrm].rlen))):
                raise NotImplementedError(
                    'Encountered invalid genome sequence request.')
            else:
                r_seq = self.index.fetch(chrm, start + 1, end).seq.upper()
        else:
            if (start is not None and (
                    start < 0 or start > len(self.index[chrm]))) or (
                        error_end and end is not None and
                        (end < 0 or
                         end > len(self.index[chrm]))):
                raise NotImplementedError(
                    'Encountered invalid genome sequence request.')
            r_seq = self.index[chrm][start:end].upper()

        if self.has_rna_bases:
            r_seq = rev_transcribe(r_seq)

        return r_seq

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

def load_index_data(fast5s_dir, corr_grp):
    fast5s_dir = (fast5s_dir if fast5s_dir.endswith('/') else
                  fast5s_dir + '/')
    index_fn = get_index_fn(fast5s_dir, corr_grp)
    try:
        import cPickle as pickle
    except:
        import pickle
    with io.open(index_fn, 'rb') as index_fp:
        raw_index_data = pickle.load(index_fp)

    num_index_vals = len(next(iter(raw_index_data.values()))[0])
    if num_index_vals == 8:
        def convert_r_data(from_base_fn, start, end, rsrtr,
                           c_grp, s_grp, filtered, rna):
            return readData(start, end, filtered, rsrtr, strand,
                            os.path.join(fast5s_dir, from_base_fn),
                            corr_grp + '/' + s_grp, rna)
    elif num_index_vals == 10:
        def convert_r_data(
                from_base_fn, start, end, rsrtr, c_grp, s_grp, filtered, rna,
                sig_match_score, mean_q_score):
            return readData(start, end, filtered, rsrtr, strand,
                            os.path.join(fast5s_dir, from_base_fn),
                            corr_grp + '/' + s_grp, rna,
                            sig_match_score, mean_q_score)
    else:
        raise NotImplementedError('Invalid Tombo index file.')

    raw_read_coverage = {}
    for (chrm, strand), cs_raw_data in raw_index_data.items():
        cs_data = [convert_r_data(*r_data) for r_data in cs_raw_data]
        # don't add chrm/strand if all reads are filtered
        if len(cs_data) > 0:
            raw_read_coverage[(chrm, strand)] = cs_data

    return fast5s_dir, index_fn, raw_read_coverage

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
        corr_grp, subgroup, rna, is_filtered=False, sig_match_score=None,
        mean_q_score=None):
    """
    Prepare data for storage in the index file
    """
    mapped_end = genome_loc.Start + len(segs) - 1

    return ((genome_loc.Chrom, genome_loc.Strand), readData(
        genome_loc.Start, mapped_end, is_filtered, read_start_rel_to_raw,
        genome_loc.Strand, fast5_fn, corr_grp + '/' + subgroup, rna,
        sig_match_score, mean_q_score))

def write_index_file(all_index_data, index_fn, basedir):
    """
    Write index file
    """
    try:
        import cPickle as pickle
    except:
        import pickle
    index_data = defaultdict(list)
    for chrm_strand, rd in all_index_data:
        # clip the basedir off the FAST5 filename in case later functions are
        # called from another relative path
        from_base_fn = rd.fn.replace(basedir, '')
        index_data[chrm_strand].append((
            from_base_fn, rd.start, rd.end, rd.read_start_rel_to_raw,
            rd.corr_group.split('/')[0], rd.corr_group.split('/')[-1],
            rd.filtered, rd.rna, rd.sig_match_score, rd.mean_q_score))

    with io.open(index_fn, 'wb') as index_fp:
        # note protocol 2 for py2/3 compatibility
        pickle.dump(dict(index_data), index_fp, protocol=2)

    return

def clear_filters(fast5s_dir, corr_grp):
    """
    Clear filters applied to this directories index files
    """
    _status_message('Loading index data.')
    try:
        fast5s_dir, index_fn, index_data = load_index_data(fast5s_dir, corr_grp)
    except IOError:
        _error_message_and_exit(
            'Filters can only be applied to runs ' +
            'with a Tombo index file. Re-run resquiggle without the ' +
            '--skip-index option to apply filters.')

    _status_message('Clearing all filters.')
    new_index_data = []
    for chrm_strand, cs_raw_data in index_data.items():
        new_index_data.extend([(chrm_strand, rd._replace(filtered=False))
                               for rd in cs_raw_data])

    write_index_file(new_index_data, index_fn, fast5s_dir)
    _status_message('All filters successfully cleared!')

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

def filter_reads_for_stuck(fast5s_dir, corr_grp, obs_filter):
    """
    Filter reads based on some observation per base threshold criteria
    """
    def read_is_stuck(fast5_fn, s_grp):
        try:
            with h5py.File(fast5_fn, 'r') as fast5_data:
                base_lens = fast5_data['/Analyses/' + s_grp + '/Events']['length']
                return any(np.percentile(base_lens, pctl) > thresh
                           for pctl, thresh in obs_filter)
        except:
            raise
            return True

    _status_message('Loading index data.')
    try:
        fast5s_dir, index_fn, index_data = load_index_data(fast5s_dir, corr_grp)
    except IOError:
        _error_message_and_exit(
            'Filters can only be applied to runs with a Tombo index file. ' +
            'Re-run resquiggle without the --skip-index option to apply ' +
            'filters.')

    _status_message('Filtering stuck reads.')
    filt_index_data = []
    prev_unfilt_reads, num_filt_reads = 0, 0
    for chrm_strand, cs_raw_data in index_data.items():
        prev_unfilt_reads += len(cs_raw_data) - sum([
            rd.filtered for rd in cs_raw_data])
        cs_filt_reads = [(chrm_strand, rd._replace(
            filtered = rd.filtered or read_is_stuck(rd.fn, rd.corr_group)))
                         for rd in cs_raw_data]
        num_filt_reads += sum([i_data[1].filtered for i_data in cs_filt_reads])
        filt_index_data.extend(cs_filt_reads)

    _status_message(
        'Filtered ' + unicode(num_filt_reads) + ' reads due to observations ' +
        'per base filter from a total of ' + unicode(prev_unfilt_reads) +
        ' reads in ' + fast5s_dir + '.')

    write_index_file(filt_index_data, index_fn, fast5s_dir)

    return

def filter_reads_for_coverage(fast5s_dir, corr_grp, frac_to_filter):
    _status_message('Loading index data.')
    try:
        fast5s_dir, index_fn, index_data = load_index_data(fast5s_dir, corr_grp)
    except IOError:
        _error_message_and_exit(
            'Filters can only be applied to runs with a Tombo index file. ' +
            'Re-run resquiggle without the --skip-index option to apply ' +
            'filters.')

    _status_message('Filtering reads to obtain more uniform coverage.')
    unfilt_data = []
    unfilt_reads_cov = []
    prev_filt_data = []
    for chrm_strand, cs_raw_data in index_data.items():
        # compute coverage
        max_end = max(rd.end for rd in cs_raw_data)
        cs_coverage = np.zeros(max_end, dtype=np.int64)
        for rd in cs_raw_data:
            if rd.filtered: continue
            cs_coverage[rd.start:rd.end] += 1
        # assign coverage value to each read
        for rd in cs_raw_data:
            if rd.filtered:
                prev_filt_data.append((chrm_strand, rd))
                continue
            # add approximate coverage from middle of read
            # faster than mean over the whole read
            unfilt_reads_cov.append(cs_coverage[
                rd.start + ((rd.end - rd.start) // 2)])
            unfilt_data.append((chrm_strand, rd))

    num_reads = len(unfilt_data)
    if num_reads == 0:
        _error_message_and_exit(
            'No unfiltered reads present in current Tombo index.')
    num_filt_reads = int(frac_to_filter * num_reads)
    _status_message(
        'Filtering ' + unicode(num_filt_reads) +
        ' reads due even coverage filter from a total of ' +
        unicode(num_reads) + ' reads in ' + fast5s_dir + '.')

    # create probabilities array with coverage values normalized to sum to 1
    unfilt_reads_cov = np.array(unfilt_reads_cov, dtype=np.float)
    unfilt_reads_p = unfilt_reads_cov / unfilt_reads_cov.sum()
    # randomly chose reads to filter
    filt_indices = np.random.choice(
        num_reads, size=num_filt_reads, replace=False, p=unfilt_reads_p)
    filt_index_data = [
        (chrm_strand, rd._replace(filtered=True))
        for chrm_strand, rd in itemgetter(*filt_indices)(unfilt_data)]
    unfilt_index_data = list(itemgetter(*list(set(range(num_reads)).difference(
        filt_indices)))(unfilt_data))

    write_index_file(prev_filt_data + filt_index_data + unfilt_index_data,
                     index_fn, fast5s_dir)

    return

def filter_reads_for_qscore(fast5s_dir, bc_grp, corr_grp, q_score_thresh):
    """
    Filter reads based on mean q-score
    """
    def read_fails_q_score(fast5_fn, s_grp):
        try:
            with h5py.File(fast5_fn, 'r') as fast5_data:
                r_q_scores = fast5_data['/Analyses/' + bc_grp + '/' + s_grp +
                                       '/Fastq'].value.decode().split('\n')[3]
                if sys.version_info[0] > 2:
                    return np.mean(
                        [q_val - PHRED_BASE for q_val in
                         r_q_scores.encode('ASCII')]) < q_score_thresh
                else:
                    return np.mean(
                        [ord(q_val) - PHRED_BASE for q_val in
                         r_q_scores.encode('ASCII')]) < q_score_thresh
        except:
            return True

    _status_message('Loading index data.')
    try:
        fast5s_dir, index_fn, index_data = load_index_data(fast5s_dir, corr_grp)
    except IOError:
        _error_message_and_exit(
            'Filters can only be applied to runs with a Tombo index file. ' +
            'Re-run resquiggle without the --skip-index option to ' +
            'apply filters.')

    _status_message('Filtering reads below a mean q-score cutoff.')
    filt_index_data = []
    num_filt_reads, prev_unfilt_reads = 0, 0
    for chrm_strand, cs_raw_data in index_data.items():
        cs_prev_filt_reads = sum([
            rd.filtered for rd in cs_raw_data])
        prev_unfilt_reads += len(cs_raw_data) - cs_prev_filt_reads
        cs_filt_reads = [
            (chrm_strand, rd._replace(
                # if q_score was previously stored use that else get
                # q-score from fast5
                filtered = rd.filtered or (
                    read_fails_q_score(rd.fn, rd.corr_group.split('/')[-1])
                    if rd.mean_q_score is None else
                    rd.mean_q_score < q_score_thresh)))
            for rd in cs_raw_data]
        num_filt_reads += sum([i_data[1].filtered
                               for i_data in cs_filt_reads]) - cs_prev_filt_reads
        filt_index_data.extend(cs_filt_reads)

    _status_message(
        'Filtered ' + unicode(num_filt_reads) + ' reads due to q-score ' +
        'filter from a total of ' + unicode(prev_unfilt_reads) + ' reads in ' +
        fast5s_dir + '.')

    write_index_file(filt_index_data, index_fn, fast5s_dir)

    return

def filter_reads_for_signal_matching(fast5s_dir, corr_grp, sig_match_thresh):
    """
    Filter reads based on mean half z-score matching to expected levels
    """
    def read_fails_matching_score(fast5_fn, corr_group):
        try:
            with h5py.File(fast5_fn, 'r') as fast5_data:
                return fast5_data['/Analyses/' + corr_group].attrs[
                    'signal_match_score'] > sig_match_thresh
        except:
            return True

    _status_message('Loading index data.')
    try:
        fast5s_dir, index_fn, index_data = load_index_data(fast5s_dir, corr_grp)
    except IOError:
        _error_message_and_exit(
            'Filters can only be applied to runs with a Tombo index file. ' +
            'Re-run resquiggle without the --skip-index option to ' +
            'apply filters.')

    _status_message('Filtering reads above a signal matching score threshold.')
    filt_index_data = []
    num_filt_reads, prev_unfilt_reads = 0, 0
    for chrm_strand, cs_raw_data in index_data.items():
        cs_prev_filt_reads = sum([rd.filtered for rd in cs_raw_data])
        prev_unfilt_reads += len(cs_raw_data) - cs_prev_filt_reads
        cs_filt_reads = [
            (chrm_strand, rd._replace(
                # if sig_match_score was previously stored use that else get
                # sig_match_score from fast5
                filtered = rd.filtered or (
                    read_fails_matching_score(rd.fn, rd.corr_group)
                    if rd.sig_match_score is None else
                    rd.sig_match_score > sig_match_thresh)))
            for rd in cs_raw_data]
        num_filt_reads += sum([i_data[1].filtered for i_data in
                               cs_filt_reads]) - cs_prev_filt_reads
        filt_index_data.extend(cs_filt_reads)

    _status_message(
        'Filtered ' + unicode(num_filt_reads) +
        ' reads due to signal matching filter from a total of ' +
        unicode(prev_unfilt_reads) + ' reads in ' + fast5s_dir + '.')

    write_index_file(filt_index_data, index_fn, fast5s_dir)

    return

def filter_reads_for_genome_pos(fast5s_dir, corr_grp, include_regs):
    """
    Filter reads to include or exclude genomic regions
    """
    def read_not_included(start, end, chrm_include_regs):
        if chrm_include_regs is None:
            return False
        return not any((start >= i_start and end <= i_end)
                       for i_start, i_end in chrm_include_regs)

    _status_message('Loading index data.')
    try:
        fast5s_dir, index_fn, index_data = load_index_data(fast5s_dir, corr_grp)
    except IOError:
        _error_message_and_exit(
            'Filters can only be applied to runs with a Tombo index file. ' +
            'Re-run resquiggle without the --skip-index option to ' +
            'apply filters.')

    _status_message('Filtering reads outside of the specified genomic location.')
    filt_index_data = []
    num_filt_reads, prev_unfilt_reads = 0, 0
    for (chrm, strand), cs_raw_data in index_data.items():
        cs_prev_filt_reads = sum([rd.filtered for rd in cs_raw_data])
        prev_unfilt_reads += len(cs_raw_data) - cs_prev_filt_reads
        do_filter_cs_reads = chrm not in include_regs
        cs_filt_reads = [((chrm, strand), rd._replace(
            filtered = rd.filtered or do_filter_cs_reads or read_not_included(
                rd.start, rd.end, include_regs[chrm])))
                         for rd in cs_raw_data]
        num_filt_reads += sum([i_data[1].filtered for i_data in
                               cs_filt_reads]) - cs_prev_filt_reads
        filt_index_data.extend(cs_filt_reads)

    _status_message(
        'Filtered ' + unicode(num_filt_reads) +
        ' reads due to genomic position filter from a total of ' +
        unicode(prev_unfilt_reads) + ' reads in ' + fast5s_dir + '.')

    write_index_file(filt_index_data, index_fn, fast5s_dir)

    return


#####################################
###### FAST5 Parsing Functions ######
#####################################

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
        for fn in fns:
            if not fn.endswith('.fast5'): continue
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
    ignore_locks_mess = (
        'This set of reads is currently being processed by another ' +
        'resquiggle command. Multiple resquiggle commands cannot be ' +
        'run concurrently on a set of reads to avoid corrupting ' +
        'read files. If you are sure this set of reads is not being ' +
        'processed by another command (usually caused by previous ' +
        'unexpected exit) set the --ignore-read-locks flag.')
    all_fast5s = []
    lock_fns = []
    try:
        # walk through directory structure searching for fast5 files
        for root, _, fns in os.walk(fast5s_dir):
            lock_fn = get_lock_fn(root)
            if not ignore_locks and os.path.exists(lock_fn):
                clear_tombo_locks(lock_fns)
                _error_message_and_exit(ignore_locks_mess)
            lock_fns.append(lock_fn)
            # create empty file indicating this directory is locked
            open(lock_fn, 'w').close()

            for fn in fns:
                if not fn.endswith('.fast5'): continue
                all_fast5s.append(os.path.join(root, fn))
    except:
        clear_tombo_locks(lock_fns)
        _error_message_and_exit(
            'Unexpected error during file enumeration. Check that you have ' +
            'write permission within the specified [fast5_basedir].')

    return all_fast5s, lock_fns

def get_raw_read_slot(fast5_data):
    try:
        raw_read_slot = list(fast5_data['/Raw/Reads'].values())[0]
    except:
        raise NotImplementedError(
            'Raw data is not found in /Raw/Reads/Read_[read#]')

    return raw_read_slot

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
        for rd in cs_raw_data:
            if rd.corr_group.split('/')[0] != corr_grp: continue
            new_index_data.append(((chrm, strand), rd._replace(
                corr_group=new_corr_grp)))

    fast5s_dir = (fast5s_dir if fast5s_dir.endswith('/') else
                  fast5s_dir + '/')
    new_index_fn = get_index_fn(fast5s_dir, new_corr_grp)
    write_index_file(new_index_data, new_index_fn, fast5s_dir)

    return

def parse_fast5s_w_index(fast5s_dir, corr_grp, subgroups, new_corr_grp):
    """
    Use index file to parse information about a set of reads
    """
    try:
        fast5s_dir, index_fn, index_data = load_index_data(fast5s_dir, corr_grp)
    except UnicodeDecodeError:
        _warning_message(
            'Invalid Tombo index file.\n\t\tThis occurs most often when the ' +
            're-squiggle command was completed using a Tombo build against ' +
            'a different python version (2 or 3).')
        raise
    raw_read_coverage = {}
    for (chrm, strand), cs_raw_data in index_data.items():
        cs_data = [
            rd for rd in cs_raw_data
            if rd.corr_group.split('/')[0] == corr_grp and
            rd.corr_group.split('/')[-1] in subgroups and not rd.filtered]
        # don't add chrm/strand if all reads are filtered
        if len(cs_data) > 0:
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
                 new_corr_grp=None, rna=False, sample_name=None):
    """
    Parse data from a list of re-squiggle fast5 directories
    """
    if VERBOSE:
        status_mess = ('Parsing Tombo index file(s).' if sample_name is None else
                       'Parsing ' + sample_name + ' Tombo index file(s).')
        _status_message(status_mess)
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
                _warning_message(
                    'Failed to parse tombo index file for ' + fast5s_dir +
                    ' directory. Creating index from FAST5 files.')
                wo_index_dirs.append(fast5s_dir)
        else:
            if not warn_index:
                _warning_message(
                    'Tombo index file does not exist for one or more ' +
                    'directories.\n\t\tIf --skip-index was not set for ' +
                    're-squiggle command, ensure that the specified ' +
                    'directory is the same as for the re-squiggle command.\n')
                warn_index = True
            wo_index_dirs.append(fast5s_dir)
    wo_index_cov = parse_fast5s_wo_index(
        wo_index_dirs, corrected_group, basecall_subgroups, rna)
    raw_read_coverage = merge_cov(w_index_covs, wo_index_cov)

    return raw_read_coverage


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

def iter_mean_slot_values(raw_read_coverage, chrm_sizes, slot_name,
                          raw_read_coverage2=None):
    """
    Iterate through chromosomes and strands yielding mean slots values over
    all reads at each covered genomic location.

    Generator returns chrmosome, strand, cs_mean_values tuples (3 return values).

    If a second raw_read_coverage object is included another cs_mean_values
    array is yeilded (4 return values)
    """
    # ignore divide by zero errors that occur where there is no
    # coverage. Need to correct nan values after subtracting two sets of
    # coverage so leave as nan for now
    old_err_settings = np.seterr(all='ignore')
    for chrm, strand in [(c, s) for c in sorted(chrm_sizes) for s in ('+', '-')]:
        if raw_read_coverage2 is None:
            if (chrm, strand) not in raw_read_coverage: continue
            cs_mean_values = get_mean_slot_genome_centric(
                raw_read_coverage[(chrm, strand)], chrm_sizes[chrm], slot_name)
            yield chrm, strand, cs_mean_values
        else:
            cs_mean_values, cs_mean_values2 = None, None
            if (chrm, strand) in raw_read_coverage:
                cs_mean_values = get_mean_slot_genome_centric(
                    raw_read_coverage[(chrm, strand)], chrm_sizes[chrm],
                    slot_name)
            if (chrm, strand) in raw_read_coverage2:
                cs_mean_values2 = get_mean_slot_genome_centric(
                    raw_read_coverage2[(chrm, strand)], chrm_sizes[chrm],
                    slot_name)
            if cs_mean_values is None and cs_mean_values2 is None: continue
            yield chrm, strand, cs_mean_values, cs_mean_values2

    _ = np.seterr(**old_err_settings)

    return

def get_largest_signal_differences(
        raw_read_coverage1, raw_read_coverage2, num_regions, num_bases):
    chrm_sizes = get_chrm_sizes(raw_read_coverage1, raw_read_coverage2)

    all_largest_diff_poss = []
    for chrm, strand, cs_sig_means1, cs_sig_means2 in iter_mean_slot_values(
            raw_read_coverage1, chrm_sizes, 'norm_mean', raw_read_coverage2):
        if cs_sig_means1 is None or cs_sig_means2 is None: continue
        chrm_diffs = np.nan_to_num(np.abs(cs_sig_means1 - cs_sig_means2))
        chrm_max_diff_regs = np.argsort(chrm_diffs)[::-1][:num_regions]
        all_largest_diff_poss.extend((
            chrm_diffs[pos], max(pos - int(num_bases / 2.0), 0),
            chrm, strand) for pos in chrm_max_diff_regs)

    return sorted(all_largest_diff_poss, reverse=True)[:num_regions]

def get_signal_differences(raw_read_coverage1, raw_read_coverage2):
    """
    Helper function to compute all signal differences
    """
    chrm_sizes = get_chrm_sizes(raw_read_coverage1, raw_read_coverage2)

    all_diffs =  {}
    for chrm, strand, cs_sig_means1, cs_sig_means2 in iter_mean_slot_values(
            raw_read_coverage1, chrm_sizes, 'norm_mean', raw_read_coverage2):
        if cs_sig_means1 is None or cs_sig_means2 is None: continue
        all_diffs[(chrm, strand)] = np.nan_to_num(cs_sig_means1 - cs_sig_means2)

    return all_diffs


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

        scale_values = scaleValues(
            corr_subgrp.attrs['shift'], corr_subgrp.attrs['scale'],
            corr_subgrp.attrs['lower_lim'], corr_subgrp.attrs['upper_lim'])
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

    # trim and flip raw signal (perform normalization outside of this
    # function in order to avoid circular import between tombo_helper
    # and tombo_stats)
    num_reg_obs = overlap_seg_data[-1] - overlap_seg_data[0]
    if r_data.strand == "+":
        reg_start_rel_raw = rsrtr + overlap_seg_data[0]
        r_sig = all_sig[reg_start_rel_raw:reg_start_rel_raw + num_reg_obs]
    else:
        reg_start_rel_raw = rsrtr + segs[-1] - overlap_seg_data[-1]
        r_sig = all_sig[reg_start_rel_raw:reg_start_rel_raw + num_reg_obs]
        r_sig = r_sig[::-1]

    return r_sig, overlap_seg_data, start_offset, scale_values

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
            scale_values = scaleValues([
                corr_grp.attrs[attr_name] for attr_name in (
                    'shift', 'scale', 'lower_lim', 'upper_lim')])

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

    return (read_id, signal_data, raw_offset, scale_values, old_segs,
            old_align_vals, new_align_vals, events_end, new_segs)

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

    return (r_means, r_seq, all_sig, segs, r_data.read_start_rel_to_raw,
            r_attrs['norm_type'], r_attrs['outlier_threshold'],
            genomeLoc(algn_subgrp['mapped_start'], algn_subgrp['mapped_strand'],
                      algn_subgrp['mapped_chrom']))

def get_coverage(raw_read_coverage):
    """
    Get genome coverage for a set of reads
    """
    if VERBOSE: _status_message('Calculating read coverage.')
    read_coverage = {}
    for (chrm, strand), reads_data in raw_read_coverage.items():
        if len(reads_data) == 0: continue
        max_end = max(r_data.end for r_data in reads_data)
        chrm_coverage = np.zeros(max_end, dtype=np.int64)
        for r_data in reads_data:
            chrm_coverage[r_data.start:r_data.end] += 1
        read_coverage[(chrm, strand)] = chrm_coverage

    return read_coverage

def get_coverage_regions(raw_read_coverage, raw_read_coverage2=None):
    """
    Get genome coverage for a set of reads
    """
    if VERBOSE: _status_message('Calculating read coverage.')
    all_chrm_strands = (
        raw_read_coverage.keys() if raw_read_coverage2 is None else
        set(raw_read_coverage).union(raw_read_coverage2))
    for chrm, strand in sorted(all_chrm_strands):
        if raw_read_coverage2 is None:
            reads_data = raw_read_coverage[(chrm, strand)]
        else:
            reads_data = []
            if (chrm, strand) in raw_read_coverage:
                reads_data += raw_read_coverage[(chrm, strand)]
            if (chrm, strand) in raw_read_coverage2:
                reads_data += raw_read_coverage2[(chrm, strand)]

        if len(reads_data) == 0: continue
        max_end = max(r_data.end for r_data in reads_data)
        cs_cov = np.zeros(max_end, dtype=np.int64)
        for r_data in reads_data:
            cs_cov[r_data.start:r_data.end] += 1

        cs_cov_starts = np.concatenate([
            [0,], np.where(np.diff(cs_cov))[0] + 1,
            [cs_cov.shape[0],]])
        cs_cov = cs_cov[cs_cov_starts[:-1]]
        yield chrm, strand, cs_cov, cs_cov_starts

    return

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
    read_bases = get_single_slot_read_centric(r_data, 'base')
    if read_bases is None:
        _warning_message(
            'Unable to extract data from read. Potentially corrupted file ' +
            'or invalid Tombo index file for this directory.')
        return reg_base_data, max(0, r_data.start - int_start)
    r_seq = b''.join(read_bases).decode()

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
        alignVals=None, align_info=None, old_segs=None, rna=False,
        sig_match_score=None):
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
        if sig_match_score is not None:
            corr_subgrp.attrs['signal_match_score'] = sig_match_score
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
        raise
        raise NotImplementedError(
            'Error writing resquiggle information back into fast5 file.')

    if do_close:
        try:
            fast5_data.close()
        except:
            raise NotImplementedError(
                'Error closing fast5 file after writing resquiggle information.')

    return


####################################
###### Annotate Raw Functions ######
####################################

def _prep_fast5_for_fastq(fast5_data, bc_grp_name, bc_subgrp_name, overwrite):
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
        bc_grp = analyses_grp[bc_grp_name]
        bc_subgrp = analyses_grp[bc_subgrp_name]
    except:
        try:
            bc_grp = analyses_grp.create_group(bc_grp_name)
            bc_subgrp = bc_grp.create_group(bc_subgrp_name)
        except:
            if overwrite:
                del analyses_grp[bc_grp_name]
                bc_grp = analyses_grp.create_group(bc_grp_name)
                bc_subgrp = bc_grp.create_group(bc_subgrp_name)
            else:
                raise NotImplementedError(
                    bc_grp_name + ' exists and --overwrite is not set.')

    return read_id

def _annotate_with_fastqs_worker(
        fastq_rec_q, fast5s_read_ids, fastq_slot, fq_slot_prepped,
        prog_q, warn_q, bc_grp_name, bc_subgrp_name, overwrite):
    been_warned = dict((warn_code, False) for warn_code in _WARN_CODES)
    num_recs_proc = 0
    while True:
        fastq_rec = fastq_rec_q.get()
        if fastq_rec is None:
            break

        # extract read_id from fastq (which should be the first text after
        # the "@" record delimiter up to the first white space or underscore
        read_id = fastq_rec[0].split()[0].split('_')[0][1:]
        if read_id not in fast5s_read_ids:
            if not been_warned[_WARN_ID_VAL]:
                been_warned[_WARN_ID_VAL] = True
                warn_q.put(_WARN_ID_VAL)
            continue

        try:
            with h5py.File(fast5s_read_ids[read_id], 'r+') as fast5_data:
                if not fq_slot_prepped:
                    try:
                        file_parsed_id = _prep_fast5_for_fastq(
                            fast5_data, bc_grp_name, bc_subgrp_name, overwrite)
                    except NotImplementedError:
                        if not been_warned[_WARN_OVRWRT_VAL]:
                            been_warned[_WARN_OVRWRT_VAL] = True
                            warn_q.put(_WARN_OVRWRT_VAL)
                        continue
                    if read_id != file_parsed_id:
                        if not been_warned[_WARN_MISMATCH_VAL]:
                            been_warned[_WARN_MISMATCH_VAL] = True
                            warn_q.put(_WARN_MISMATCH_VAL)
                        continue
                bc_slot = fast5_data[fastq_slot]
                # add sequence to fastq slot
                bc_slot.create_dataset(
                    'Fastq', data=''.join(fastq_rec),
                    dtype=h5py.special_dtype(vlen=unicode))

                # progress q update
                num_recs_proc += 1
                if num_recs_proc % _PROC_UPDATE_INTERVAL == 0:
                    prog_q.put(_PROC_UPDATE_INTERVAL)
        except:
            if not been_warned[_WARN_IO_VAL]:
                been_warned[_WARN_IO_VAL] = True
                warn_q.put(_WARN_IO_VAL)
            continue

    # add last number of records reported from this process
    prog_q.put(num_recs_proc % _PROC_UPDATE_INTERVAL)

    return

def _feed_seq_records_worker(fastq_fns, fastq_rec_q):
    for fastq_fn in fastq_fns:
        n_recs = 0
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
                        ' FASTQ records from ' + fastq_fn + ' before ' +
                        'encountering an invalid record. The rest of ' +
                        'this file will not be processed.')
                    break
                n_recs += 1
                fastq_rec_q.put(fastq_rec)

    return

def _get_ann_queues(prog_q, warn_q, been_warned):
    iter_added = 0
    while True:
        try:
            iter_added += prog_q.get(block=False)
        except queue.Empty:
            break

    while True:
        try:
            warn_val = warn_q.get(block=False)
        except queue.Empty:
            break

        if warn_val == _WARN_ID_VAL:
            if not been_warned[_WARN_ID_VAL]:
                _warning_message(
                    'Some records contain read identifiers not found in ' +
                    'any FAST5 files or sequencing summary files.')
            been_warned[_WARN_ID_VAL] = True
        elif warn_val == _WARN_IO_VAL:
            if not been_warned[_WARN_IO_VAL]:
                _warning_message(
                    'Some read files that could not be accessed.')
            been_warned[_WARN_IO_VAL] = True
        elif warn_val == _WARN_MISMATCH_VAL:
            if not been_warned[_WARN_MISMATCH_VAL]:
                _warning_message(
                    'Read ID found in sequencing summary and FAST5 ' +
                    'file are discordant. Skipping read.')
            been_warned[_WARN_MISMATCH_VAL] = True
        elif warn_val == _WARN_OVRWRT_VAL:
            if not been_warned[_WARN_OVRWRT_VAL]:
                _warning_message(
                    'Basecalls exsit in specified slot for some reads. ' +
                    'Set --overwrite option to overwrite these basecalls.')
            been_warned[_WARN_OVRWRT_VAL] = True
        else:
            _warning_message('Invalid wanring code encountered.')

    return iter_added, been_warned

def _annotate_with_fastqs(
        fastq_fns, fast5s_read_ids, fastq_slot, fq_slot_prepped, num_processes,
        bc_grp_name, bc_subgrp_name, overwrite):
    if VERBOSE: _status_message('Annotating FAST5s with sequence from FASTQs.')
    fastq_rec_q = Queue(maxsize=_MAX_FASTQ_QUEUE_SIZE)
    # open a single process to read fastq files and feed the fastq record queue
    fq_feed_p = Process(target=_feed_seq_records_worker,
                        args=(fastq_fns, fastq_rec_q))
    fq_feed_p.start()

    # open fast5 annotation processes
    prog_q = Queue()
    warn_q = Queue()
    ann_args = (fastq_rec_q, fast5s_read_ids, fastq_slot, fq_slot_prepped,
                prog_q, warn_q, bc_grp_name, bc_subgrp_name, overwrite)
    ann_ps = []
    for p_id in range(num_processes):
        p = Process(target=_annotate_with_fastqs_worker, args=ann_args)
        p.start()
        ann_ps.append(p)

    if VERBOSE: bar = tqdm(total=len(fast5s_read_ids), smoothing=0)

    total_added_seqs = 0
    been_warned = dict((warn_code, False) for warn_code in _WARN_CODES)
    # process progress and warn queues until fastq filler process runs out of
    # files/records
    while fq_feed_p.is_alive():
        iter_added, been_warned = _get_ann_queues(prog_q, warn_q, been_warned)
        total_added_seqs += iter_added
        if VERBOSE: bar.update(iter_added)
        sleep(0.01)

    # put none records to trigger annotation processes to exit
    for _ in range(num_processes):
        fastq_rec_q.put(None)

    # process the rest of the records
    while any(p.is_alive() for p in ann_ps) or not prog_q.empty():
        iter_added, been_warned = _get_ann_queues(prog_q, warn_q, been_warned)
        total_added_seqs += iter_added
        if VERBOSE: bar.update(iter_added)
        sleep(0.01)
    if VERBOSE: bar.close()

    if VERBOSE: _status_message('Added sequences to a total of ' +
                                str(total_added_seqs) + ' reads.')

    return

def _prep_fastq_slot_worker(
        fast5_q, bc_grp, bc_subgrp, overwrite, read_ids_q, prog_q, warn_q):
    num_files_proc = 0
    been_warned_overwrite = False
    while not fast5_q.empty():
        try:
            fast5_fn = fast5_q.get(block=False)
        except queue.Empty:
            break

        num_files_proc += 1
        if num_files_proc % _PROC_UPDATE_INTERVAL == 0:
            prog_q.put(_PROC_UPDATE_INTERVAL)

        try:
            with h5py.File(fast5_fn) as fast5_data:
                try:
                    read_id = _prep_fast5_for_fastq(
                        fast5_data, bc_grp, bc_subgrp, overwrite)
                except NotImplementedError:
                    if not been_warned_overwrite:
                        been_warned_overwrite = True
                        warn_q.put(_WARN_OVRWRT_VAL)
                    continue
        except:
            continue
        if read_id is None:
            continue

        read_ids_q.put((read_id, fast5_fn))

    prog_q.put(num_files_proc % _PROC_UPDATE_INTERVAL)

    return

def _get_prep_queue(read_ids_q, prog_q, warn_q, fast5s_read_ids, been_warned):
    """
    Process all records from all fast5 prep queues
    """
    # only process up to _ITER_QUEUE_LIMIT items each iteration
    iter_processed = 0
    while True:
        try:
            read_id, fast5_fn = read_ids_q.get(block=False)
        except queue.Empty:
            break
        iter_processed += 1
        if iter_processed > _ITER_QUEUE_LIMIT: break

        if read_id in fast5s_read_ids:
            if not been_warned[_WARN_UNIQ_VAL]:
                _warning_message(
                    'Multiple FAST5 files contain the same read identifiers. ' +
                    'Ensure that FAST5 files are from a single run.')
            been_warned[_WARN_UNIQ_VAL] = True
            continue
        fast5s_read_ids[read_id] = fast5_fn

    while True:
        try:
            warn_val = warn_q.get(block=False)
        except queue.Empty:
            break
        if warn_val == _WARN_OVRWRT_VAL:
            if not been_warned[_WARN_OVRWRT_VAL]:
                _warning_message(
                    'Basecalls exsit in specified slot for some reads. ' +
                    'Set --overwrite option to overwrite these basecalls.')
            been_warned[_WARN_OVRWRT_VAL] = True
        else:
            _warning_message('Invalid wanring code encountered.')

    iter_prog = 0
    while True:
        try:
            iter_prog += prog_q.get(block=False)
        except queue.Empty:
            break

    return fast5s_read_ids, iter_prog, been_warned

def _get_read_ids_and_prep_fastq_slot(
        fast5s_dir, bc_grp, bc_subgrp, overwrite, num_processes):
    """
    Extract read id from /Raw group and prep fastq slots for annotation with
    associated FASTQ files.
    """
    if VERBOSE: _status_message('Getting read filenames.')
    fast5_fns = get_files_list(fast5s_dir)
    num_fast5s = len(fast5_fns)
    fast5_q = Queue()
    for fast5_fn in fast5_fns:
        fast5_q.put(fast5_fn)

    if VERBOSE: _status_message(
            'Preparing reads and extracting read identifiers.')
    read_ids_q = Queue()
    prog_q = Queue()
    warn_q = Queue()
    prep_args = (fast5_q, bc_grp, bc_subgrp, overwrite, read_ids_q,
                 prog_q, warn_q)
    prep_ps = []
    for p_id in range(num_processes):
        p = Process(target=_prep_fastq_slot_worker, args=prep_args)
        p.start()
        prep_ps.append(p)

    fast5s_read_ids = {}
    # Warn non-unique read_ids in directory
    been_warned = dict((warn_code, False) for warn_code in _WARN_CODES_PREP)
    if VERBOSE: bar = tqdm(total=num_fast5s, smoothing=0)
    while any(p.is_alive() for p in prep_ps):
        fast5s_read_ids, iter_prog, been_warned = _get_prep_queue(
            read_ids_q, prog_q, warn_q, fast5s_read_ids, been_warned)
        if VERBOSE: bar.update(iter_prog)
        sleep(0.01)

    fast5s_read_ids, iter_prog, been_warned = _get_prep_queue(
        read_ids_q, prog_q, warn_q, fast5s_read_ids, been_warned)
    if VERBOSE: bar.update(iter_prog)
    if VERBOSE: bar.close()

    return fast5s_read_ids

def _parse_sequencing_summary_files(fast5s_dir, seq_summary_fns):
    if VERBOSE: _status_message('Getting read filenames.')
    full_fast5_fns = {}
    # walk through directory structure searching for fast5 files
    for root, _, fns in os.walk(fast5s_dir):
        for fn in fns:
            if not fn.endswith('.fast5'): continue
            full_fast5_fns[fn] = os.path.join(root, fn)

    if VERBOSE: _status_message('Parsing sequencing summary files.')
    fast5s_read_ids = {}
    been_warned = False
    for seq_summary_fn in seq_summary_fns:
        with open(seq_summary_fn) as fp:
            try:
                header_fields = fp.readline().split()
                fn_field = next(i for i, h_field in enumerate(header_fields)
                                if re.match(_SEQ_SUMMARY_FN_FIELD, h_field))
                id_field = next(i for i, h_field in enumerate(header_fields)
                                if re.match(_SEQ_SUMMARY_ID_FIELD, h_field))
            except:
                _warning_message(
                    'Could not extract header information for sequencing ' +
                    'summary file: ' + seq_summary_fn)
                continue
            try:
                for line in fp:
                    rec_fields = line.split()
                    rec_short_fn = rec_fields[fn_field]
                    try:
                        rec_full_fn = full_fast5_fns[rec_short_fn]
                    except KeyError:
                        if not been_warned:
                            _warning_message(
                                'Some records from sequencing summaries ' +
                                'do not appear to have a matching file.')
                        been_warned = True
                        continue
                    # convert filename to full filename and link to read id
                    fast5s_read_ids[rec_fields[id_field]] = rec_full_fn
            except:
                _warning_message(
                    'Error parsing records for sequencing ' +
                    'summary file: ' + seq_summary_fn)

    return fast5s_read_ids


###################################
###### Filter Main Functions ######
###################################

def _clear_filters_main(args):
    for fast5s_dir in args.fast5_basedirs:
        clear_filters(fast5s_dir, args.corrected_group)

    return

def _filter_stuck_main(args):
    obs_filter = parse_obs_filter(args.obs_per_base_filter)
    for fast5s_dir in args.fast5_basedirs:
        filter_reads_for_stuck(fast5s_dir, args.corrected_group, obs_filter)

    return

def _filter_coverage_main(args):
    if not 0 < args.percent_to_filter < 100:
        _error_message_and_exit(
            '--percent-to-filter must be between 0 and 100.')

    for fast5s_dir in args.fast5_basedirs:
        filter_reads_for_coverage(
            fast5s_dir, args.corrected_group, args.percent_to_filter / 100.0)

    return

def _filter_q_score_main(args):
    if not 0 < args.q_score < 40:
        _error_message_and_exit('--q-score must be between 0 and 40.')

    for fast5s_dir in args.fast5_basedirs:
        filter_reads_for_qscore(
            fast5s_dir, args.basecall_group, args.corrected_group, args.q_score)

    return

def _filter_signal_matching_main(args):
    if not 0 < args.signal_matching_score < 10:
        _error_message_and_exit(
            '--signal-matching-score must be between 0 and 10.')

    for fast5s_dir in args.fast5_basedirs:
        filter_reads_for_signal_matching(
            fast5s_dir, args.corrected_group, args.signal_matching_score)

    return

def _filter_genome_pos_main(args):
    include_regs = parse_genome_regions(args.include_regions)

    for fast5s_dir in args.fast5_basedirs:
        filter_reads_for_genome_pos(
            fast5s_dir, args.corrected_group, include_regs)

    return

def _filter_main(args):
    global VERBOSE
    VERBOSE = not args.quiet

    if args.action_command == 'clear_filters':
        _clear_filters_main(args)
    elif args.action_command == 'genome_locations':
        _filter_genome_pos_main(args)
    elif args.action_command == 'stuck':
        _filter_stuck_main(args)
    elif args.action_command == 'level_coverage':
        _filter_coverage_main(args)
    elif args.action_command == 'q_score':
        _filter_q_score_main(args)
    elif args.action_command == 'raw_signal_matching':
        _filter_signal_matching_main(args)
    else:
        _error_message_and_exit('Invalid Tombo filter command.')

    return


##################################
###### Annotate FAST5s Main ######
##################################

def _annotate_reads_with_fastq_main(args):
    global VERBOSE
    VERBOSE = not args.quiet

    fast5s_basedir = (
        args.fast5_basedir if args.fast5_basedir.endswith('/') else
        args.fast5_basedir + '/')
    if args.sequencing_summary_filenames:
        fast5s_read_ids = _parse_sequencing_summary_files(
            fast5s_basedir, args.sequencing_summary_filenames)
        fq_slot_prepped = False
    else:
        fast5s_read_ids = _get_read_ids_and_prep_fastq_slot(
            fast5s_basedir, args.basecall_group, args.basecall_subgroup,
            args.overwrite, args.processes)
        fq_slot_prepped = True
    fastq_slot = '/'.join(('/Analyses', args.basecall_group,
                           args.basecall_subgroup))
    _annotate_with_fastqs(
        args.fastq_filenames, fast5s_read_ids, fastq_slot, fq_slot_prepped,
        args.processes, args.basecall_group, args.basecall_subgroup,
        args.overwrite)

    return


if __name__ == '__main__':
    raise NotImplementedError(
        'This is a module. See commands with `tombo -h`')
