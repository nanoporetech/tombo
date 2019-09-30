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
np.seterr(all='raise')

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
from ._default_parameters import PHRED_BASE, DNA_SAMP_TYPE, RNA_SAMP_TYPE
from ._c_helper import (
    c_new_mean_stds, c_new_means, c_valid_cpts_w_cap,
    c_valid_cpts_w_cap_t_test)
from ._c_dynamic_programming import (
    c_banded_traceback, c_adaptive_banded_forward_pass)


# list of classes/functions to include in API
__all__ = [
    'readData', 'intervalData', 'TomboReads',
    'TomboMotif', 'parse_motif_descs',
    'resquiggleParams', 'startClipParams',
    'resquiggleResults', 'alignInfo', 'genomeLocation',
    'sequenceData', 'dpResults', 'scaleValues', 'Fasta', 'seqSampleType',
    'get_seq_sample_type', 'get_raw_read_slot',
    'get_single_slot_read_centric', 'get_multiple_slots_read_centric']


VERBOSE = True

# single base conversion for motifs
SINGLE_LETTER_CODE = {
    'A':'A', 'C':'C', 'G':'G', 'T':'T', 'B':'[CGT]',
    'D':'[AGT]', 'H':'[ACT]', 'K':'[GT]', 'M':'[AC]',
    'N':'[ACGT]', 'R':'[AG]', 'S':'[CG]', 'V':'[ACG]',
    'W':'[AT]', 'Y':'[CT]'}
INVALID_BASES = re.compile('[^ACGT]')
INVALID_BASE_RUNS = re.compile('[^ACGT]+')


###############################
###### Custom TomboError ######
###############################

class TomboError(Exception):
    pass


######################################
###### Cython Wrapper Functions ######
######################################

# can't import TomboError in cython so wrap these functions in python
def valid_cpts_w_cap(*args, **kwargs):
    try:
        valid_cpts = c_valid_cpts_w_cap(*args, **kwargs)
    except NotImplementedError as e:
        raise TomboError(unicode(e))
    valid_cpts.sort()
    return valid_cpts

def valid_cpts_w_cap_t_test(*args, **kwargs):
    try:
        valid_cpts = c_valid_cpts_w_cap_t_test(*args, **kwargs)
    except NotImplementedError as e:
        raise TomboError(unicode(e))
    valid_cpts.sort()
    return valid_cpts

def banded_traceback(*args, **kwargs):
    try:
        return c_banded_traceback(*args, **kwargs)
    except NotImplementedError as e:
        raise TomboError(unicode(e))

def adaptive_banded_forward_pass(*args, **kwargs):
    try:
        return c_adaptive_banded_forward_pass(*args, **kwargs)
    except NotImplementedError as e:
        raise TomboError(unicode(e))


################################
###### Global Namedtuples ######
################################

class alignInfo(namedtuple(
        'alignInfo',
        ('ID', 'Subgroup', 'ClipStart', 'ClipEnd',
         'Insertions', 'Deletions', 'Matches', 'Mismatches'))):
    """Information from genomic read alignment

    Args:
        ID (str): read identifier
        Subgroup (str): sub-group location containing read information (e.g. 'BaseCalled_template')
        ClipStart (int): number of bases clipped from beginning of basecalls
        ClipEnd (int): number of bases clipped from end of basecalls
        Insertions (int): number of inserted bases in alignment
        Deletions (int): number of delected bases in alignment
        Matches (int): number of matched bases in alignment
        Mismatches (int): number of mis-matched bases in alignment
    """

# TODO convert rna to rev_sig
class readData(namedtuple('readData', (
        'start', 'end', 'filtered', 'read_start_rel_to_raw', 'strand', 'fn',
        'corr_group', 'rna', 'sig_match_score', 'mean_q_score', 'read_id'))):
    """Nanopore read meta-data

    Example::

        r_data = tombo_helper.readData(
            start=0, end=1000, filtered=False, read_start_rel_to_raw=100,
            strand='+', fn='path/to/read.fast5',
            corr_group='RawGenomeCorrected_000', rna=False,
            sig_match_score=1.0, mean_q_score=10.0)

        # output is a list of readData objects
        reads_index = tombo_helper.TomboReads(['test_data/native_reads',])
        cs_reads = reads_index.get_cs_reads('chr20', '+')

    Args:
        start (int): 0-based start mapped position
        end (int): 0-based open interval end mapped position
        filtered (bool): is this read filtered?
        read_start_rel_to_raw (int): start (in raw signal vector) of assigned bases
        strand (str): read mapped strand ('+', '-' or None)
        fn (str): raw read filename.
        corr_group (str): --corrected-group slot specified in re-squiggle command
        rna (bool): is this read RNA?
        sig_match_score (float): signal matching score. Default: None
        mean_q_score (float): mean basecalling q-score. Default: None
        read_id (str): read identifier Default: None
    """
# set default values for sig_match_score, q_score and read_id
readData.__new__.__defaults__ = (None, None, None)

class scaleValues(namedtuple(
        'scaleValues', (
            'shift', 'scale', 'lower_lim', 'upper_lim', 'outlier_thresh'))):
    """Signal normaliztion scaling parameters. For details see https://nanoporetech.github.io/tombo/resquiggle.html#signal-normalization

    Args:
        shift (float): shift scaling parameter
        scale (float): scale scaling parameter
        lower_lim (float): lower windsorizing threshold
        upper_lim (float): upper windsorizing threshold
        outlier_thresh (float): outlier threshold used to define `lower_lim` and `upper_lim`
    """

class resquiggleParams(namedtuple(
        'resquiggleParams',
        ('match_evalue', 'skip_pen', 'bandwidth', 'max_half_z_score',
         'running_stat_width', 'min_obs_per_base', 'raw_min_obs_per_base',
         'mean_obs_per_event', 'z_shift', 'stay_pen', 'use_t_test_seg',
         'band_bound_thresh', 'start_bw', 'start_save_bw', 'start_n_bases'))):
    """Re-squiggle parameters

    Args:
        match_evalue (float): expected value for matching event to sequence
        skip_pen (float): penalty for skipped sequence position
        bandwidth (int): adaptive bandwidth
        max_half_z_score (float): windsorize half z-scores above this value
        running_stat_width (int): running neighboring window width for segmentation scoring
        min_obs_per_base (int): minimum observations per genomic base
        mean_obs_per_event (int): mean number of raw obs. per event during segmentation
        z_shift (float): amount to shift z-scores for DP (derived from match_evalue)
        stay_pen (float): stay penalty for DP (derived from match_evalue)
        use_t_test_seg (bool): use t-test segmentation criterion (default: raw neighboring window difference)
        band_bound_thresh (int): bandwidth boundary threshold for determining if a read potentially left the adaptive bandwidth
        start_bw (int): bandwidth for read start identification
        start_save_bw (int): save bandwidth for read start identification (if start_bw_fails)
        start_n_bases (int): number of genomic bases to use for read start identification
    """
# set default values for start params
resquiggleParams.__new__.__defaults__ = (None, None, None)

class trimRnaParams(namedtuple(
        'trimRnaParams',
        ('moving_window_size', 'min_running_values', 'thresh_scale',
         'max_raw_obs'))):
    """Parameters to trim RNA adapters.
    """

class stallParams(namedtuple(
        'stallParams',
        ('window_size', 'threshold', 'min_consecutive_obs', 'edge_buffer',
         # percentile stall method params
         'lower_pctl', 'upper_pctl',
         # mean windows stall method params
         'mini_window_size', 'n_windows'))):
    """Parameters to identify RNA stalls
    """
# set method specific params to None
stallParams.__new__.__defaults__ = (None,) * 4

class startClipParams(namedtuple(
        'startClipParams',
        ('bandwidth', 'num_genome_bases'))):
    """Parameters to identify read start using bases clipped from mapping

    Args:
        bandwidth (int): bandwidth
        num_genome_bases (int): number of genome bases
    """

class resquiggleResults(namedtuple(
        'resquiggleResults',
        ('align_info', 'genome_loc', 'genome_seq', 'mean_q_score',
         'raw_signal', 'channel_info', 'read_start_rel_to_raw', 'segs',
         'scale_values', 'sig_match_score', 'norm_params_changed',
         'start_clip_bases', 'stall_ints'))):
    """Re-squiggle results

    Args:
        align_info (:class:`tombo.tombo_helper.alignInfo`): read alignment information
        genome_loc (:class:`tombo.tombo_helper.genomeLocation`): genome mapping location
        genome_seq (str): mapped genome sequence
        mean_q_score (float): mean basecalling q-score
        raw_signal (np.array::np.float64): raw signal (i.e. un-segmented; signal may be normalized) (optional)
        channel_info (:class:`tombo.tombo_helper.channelInfo`): channel information (optional)
        read_start_rel_to_raw (int): read start within raw signal (optional)
        segs (np.array::np.int64): relative raw signal segment positions (optional)
        scale_values (:class:`tombo.tombo_helper.scaleValues`): signal normalization scale values (optional)
        sig_match_score (float): expected to observed signal score (see :class:`tombo.tombo_stats.get_read_seg_score`; optional)
        norm_params_changed (bool): were scale parameters updated (optional)
        start_clip_bases (str): mapping clipped bases from start of read (optional)
        stall_ints (list): list of rna stall locations within raw signal (optional)
    """
# set default None values for when just mapping results are included
resquiggleResults.__new__.__defaults__ = (None,) * 9

class dpResults(namedtuple(
        'dpResults', ('read_start_rel_to_raw', 'segs', 'ref_means', 'ref_sds',
                      'genome_seq'))):
    """Dynamic programming results

    Args:
        read_start_rel_to_raw (int): read start within raw signal
        segs (np.array::np.int64): relative raw signal segment positions
        ref_means (np.array::np.float64): expected signal levels
        ref_sds (np.array::np.float64): expected SD of signal levels
        genome_seq (str): mapped genome sequence
    """

class genomeLocation(namedtuple('genomeLocation', ('Start', 'Strand', 'Chrom'))):
    """Genomic location

    Args:
        Start (int): 0-based genomic start position
        Strand (str): Strand (should be '+' or '-')
        Chrom (str): Chromosome
    """

class sequenceData(namedtuple('sequenceData', ('seq', 'id', 'mean_q_score'))):
    """Read sequence data from FASTQ record

    Args:
        seq (str): read sequence
        id (str): read id (extracted from FASTQ record line)
        mean_q_score (float): mean q-score
    """

class channelInfo(namedtuple(
        'channelInfo', ('offset', 'range', 'digitisation',
                        'number', 'sampling_rate'))):
    """Read channel information

    Args:
        offset (float): offset parameter
        range (float): range parameter
        digitisation (float): digitisation parameter
        number (int): channel number
        sampling_rate (int): number of raw samples per second
    """

class regionStats(namedtuple(
        'regionStats', ('reg_frac_standard_base', 'reg_poss', 'chrm', 'strand',
                        'start', 'reg_cov', 'ctrl_cov', 'valid_cov'))):
    """Region statistics

    Args:
        reg_frac_standard_base (np.array::np.float64): fraction of standard bases
        reg_poss (np.array::np.int64): positions for reported fractions
        chrm (str): chromosome name
        strand (str): strand (should be '+' or '-')
        start (int): 0-based region start
        reg_cov (np.array::np.int64): region read depth
        ctrl_cov (np.array::np.int64): region control sample read depth
        valid_cov (np.array::np.int64): region valid (tested) read depth
    """

class groupStats(namedtuple(
        'groupStats', ('reg_stats', 'reg_poss', 'chrm', 'strand',
                       'start', 'reg_cov', 'ctrl_cov'))):
    """Region statistics

    Args:
        reg_stats (np.array::np.float64): statistic for group comparison
        reg_poss (np.array::np.int64): positions for reported fractions
        chrm (str): chromosome name
        strand (str): strand (should be '+' or '-')
        start (int): 0-based region start
        reg_cov (np.array::np.int64): region read depth
        ctrl_cov (np.array::np.int64): region control sample read depth
    """

class seqSampleType(namedtuple(
        'seqSampleType', ('name', 'rev_sig'))):
    """Description of a sequencing sample type

    Args:
        name (str): name of the sequencing sample
        rev_sig (bool): Is the raw signal reversed (3' to 5')
    """


######################################
###### Various Helper Functions ######
######################################

def status_message(message, indent=False):
    pre_str = '\t' if indent else ''
    sys.stderr.write(pre_str + strftime('[%H:%M:%S] ') + message + '\n')
    sys.stderr.flush()
    return

def warning_message(message):
    sys.stderr.write(
        '*' * 20 + ' WARNING ' + '*' * 20 + '\n\t' +
        message + '\n')
    sys.stderr.flush()
    return

def error_message_and_exit(message):
    sys.stderr.write(
        '*' * 20 + ' ERROR ' + '*' * 20 + '\n\t' +
        message + '\n')
    sys.exit()
    return

def resolve_path(fn_path):
    """Helper function to resolve relative and linked paths that might
    give other packages problems.
    """
    return os.path.realpath(os.path.expanduser(fn_path))

COMP_BASES = dict(zip(map(ord, 'ACGT'), map(ord, 'TGCA')))
def comp_seq(seq):
    """Complement DNA sequence
    """
    return seq.translate(COMP_BASES)
def rev_comp(seq):
    """Reverse complement DNA sequence
    """
    return seq.translate(COMP_BASES)[::-1]

def invalid_seq(seq):
    return bool(INVALID_BASES.search(seq))

U_TO_T = {ord('U'):ord('T')}
def rev_transcribe(seq):
    """Convert U bases to T
    """
    return seq.translate(U_TO_T)

def get_mean_q_score(read_q):
    if sys.version_info[0] > 2:
        return np.mean([q_val - PHRED_BASE
                        for q_val in read_q.encode('ASCII')])
    return np.mean([ord(q_val) - PHRED_BASE
                    for q_val in read_q.encode('ASCII')])

def get_chrm_sizes(reads_index, ctrl_reads_index=None):
    """Get covered chromosome sizes from a set of reads
    """
    strand_chrm_sizes = defaultdict(list)
    for (chrm, strand), cs_read_cov in reads_index:
        try:
            strand_chrm_sizes[chrm].append(max(
                    r_data.end for r_data in cs_read_cov))
        except ValueError:
            continue
    if ctrl_reads_index is not None:
        for (chrm, strand), cs_read_cov in ctrl_reads_index:
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
    """Parse genome location strings and convert to 0-based coordinates
    """
    parsed_locs = []
    for chrm_pos_strand in genome_locs:
        # strip off any quotes and return up to the first 3 values
        split_vals = chrm_pos_strand.replace('"', '').replace(
            "'", "").split(':')[:3]
        # default to plus strand if not specified
        if len(split_vals) == 1:
            error_message_and_exit(
                'Invalid genome location provided: ' + chrm_pos_strand +
                '\n\t\tTry adding quotation marks around specified genome ' +
                'locations (especially for sequence identifiers with ' +
                'special characters).')
        elif len(split_vals) == 2:
            parsed_locs.append((
                split_vals[0], int(split_vals[1]) - 1, default_strand))
        else:
            parsed_locs.append((
                split_vals[0], int(split_vals[1]) - 1, split_vals[2]))

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
                raise TomboError('Invalid region text provided.')
        except:
            error_message_and_exit(
                'Invalid [--include-region] format.')

        parsed_regs[chrm].append(reg_pos)

    parsed_regs = dict(parsed_regs)
    for chrm in include_whole_chrms:
        parsed_regs[chrm] = None

    return parsed_regs

def parse_locs_file(locs_fn):
    """Parse BED files containing genomic locations (assumes single base
    locations, so end coordinate is ignored).
    """
    n_added, n_failed = 0, 0
    raw_locs = defaultdict(set)
    with open(locs_fn) as locs_fp:
        for line in locs_fp:
            try:
                chrm, pos, _, _, _, strand = line.split()[:6]
                # bed specs indicate 0-based start so don't shift here
                pos = int(pos)
                raw_locs[(chrm, strand)].add(pos)
                n_added += 1
            except:
                n_failed += 1
                continue

    if n_failed > 0:
        if n_added > 0:
            warning_message((
                'Some provided locations ({}) were incorrectly formatted. ' +
                'Ensure that provided file ({}) is in BED format.').format(
                    n_failed, locs_fn))
        else:
            warning_message((
                'All provided locations from {} were incorrectly formatted. ' +
                'Ensure that provided file is in BED format.').format(
                    locs_fn))

    return dict((cs, np.array(sorted(cs_poss)))
                for cs, cs_poss in raw_locs.items())

def parse_obs_filter(obs_filter):
    """Parse observations per base formatted filtering
    """
    if len(obs_filter) < 1:
        return None

    # parse obs_filter
    try:
        obs_filter = [list(map(int, pctl_nobs.split(':')))
                      for pctl_nobs in obs_filter]
    except:
        raise TomboError('Invalid format for observation filter')

    if any(pctl < 0 or pctl > 100 for pctl in map(itemgetter(0), obs_filter)):
       error_message_and_exit('Invalid percentile value.')

    return obs_filter

def get_seq_kmers(seq, kmer_width, rev_strand=False):
    """Compute expected signal levels for a sequence from a reference model

    Args:
        seq (str): genomic seqeunce to be converted to expected signal levels
        kmer_width (int): k-mer width
        rev_strand (bool): provided sequence is from reverse strand (so flip return order to genome forward direction)
    """
    seq_kmers = [seq[i:i + kmer_width]
                 for i in range(len(seq) - kmer_width + 1)]
    # get stat lookups from seq on native strand then flip if rev_strand
    if rev_strand:
        seq_kmers = seq_kmers[::-1]

    return seq_kmers

class TomboMotif(object):
    """Description of a sequence motif, including potentially modified position

    Attributes:
        raw_motif (str): input raw motif string
        motif_len (int): length of motif
        motif_pat (``re.compile``): python regular expression for this motif
        rev_comp_pat (``re.compile``): python regular expression for reverse complement of this motif
        is_palindrome (bool): is the motif palindromic (in the genomic not literary definition)
        mod_pos (int): modified base position within motif
        mod_base (str): modified base (should generally be A, C, G or T)

    .. automethod:: __init__
    """
    def _parse_motif(self, raw_motif, rev_comp_motif=False):
        conv_motif = ''.join(SINGLE_LETTER_CODE[letter]
                             for letter in raw_motif)
        if rev_comp_motif:
            # reverse complement and then flip any group brackets
            conv_motif = rev_comp(conv_motif).translate({
                ord('['):']', ord(']'):'['})
        return re.compile(conv_motif)

    def _compute_partial_patterns(self):
        """Compute patterns for partial matches that include the mod_pos
        at the start, end or within short sequences.

        Key into _partial_pats with:
            1) whether searching at start, end or within a short sequence
            2) length of the partial pattern

        Values are compiled partial pattern and mod_pos within pattern.
        Short patterns are lists.
        """
        self._partial_pats = {'start':{}, 'end':{}, 'short':{}}
        for offset in range(self.mod_pos - 1):
            self._partial_pats['start'][
                self.motif_len - offset - 1] = (self._parse_motif(
                    self.raw_motif[offset + 1:]), self.mod_pos - offset - 1)
        for offset in range(self.motif_len - self.mod_pos):
            self._partial_pats['end'][
                self.motif_len - offset - 1] = (self._parse_motif(
                    self.raw_motif[:-(offset + 1)]), self.mod_pos)
        for short_len in range(1, self.motif_len):
            self._partial_pats['short'][short_len] = [
                (self._parse_motif(self.raw_motif[offset:offset + short_len]),
                 self.mod_pos - offset)
                for offset in range(
                        max(0, self.mod_pos - short_len),
                        min(self.motif_len - short_len + 1, self.mod_pos))]
        return

    def __init__(self, raw_motif, mod_pos=None):
        """Parse string motif

        Args:
            raw_motif (str): sequence motif. supports IUPAC single letter codes (use T for RNA).
            mod_pos (int): 1-based position of modified base within the motif
        """
        # TODO convert mod_pos to 0-based coordinate
        # (1-based is much more error prone)
        invalid_chars = re.findall(
            '[^' + ''.join(SINGLE_LETTER_CODE) + ']', raw_motif)
        if len(invalid_chars) > 0:
            error_message_and_exit(
                'Invalid characters in motif: ' + ', '.join(invalid_chars))

        # basic motif parsing
        self.raw_motif = raw_motif
        self.motif_len = len(raw_motif)
        self.motif_pat = self._parse_motif(raw_motif)
        self.rev_comp_pat = self._parse_motif(raw_motif, True)

        self.is_palindrome = self.motif_pat == self.rev_comp_pat

        # parse modified position from motif if provided
        self.mod_pos = mod_pos
        if mod_pos is None:
            self.mod_base = None
        else:
            assert 0 < mod_pos <= self.motif_len
            self.mod_base = raw_motif[mod_pos - 1]
            if INVALID_BASES.match(self.mod_base):
                warning_message(
                    'Provided modified position is not a single base, which ' +
                    'is likely an error. Specified modified base is one of: ' +
                    ' '.join(SINGLE_LETTER_CODE[self.mod_base][1:-1]))
            self._compute_partial_patterns()

    def __repr__(self):
        return '\n'.join(('Raw Motif:\t' + self.raw_motif,
                          'Mod Position:\t' + str(self.mod_pos),
                          'Motif Pattern:\t' + str(self.motif_pat),
                          'Rev Comp Pattern:\t' + str(self.rev_comp_pat)))

    def matches_seq(self, seq):
        """Does the motif match provided sequence (including mod_pos within seq)?

        Including partial matches at beginning and end that include mod_pos.
        """
        # check matches to start of sequence
        for start_len in range(1, min(len(seq) + 1, self.motif_len)):
            try:
                start_pat, start_mod_pos = self._partial_pats[
                    'start'][start_len]
            except KeyError:
                continue
            if start_pat.match(seq[:start_len]):
                return True

        # check central sequence overlaps
        if len(seq) < self.motif_len:
            for short_pat, mod_pos in self._partial_pats['short'][len(seq)]:
                if short_pat.match(seq):
                    return True
        else:
            if self.motif_pat.search(seq):
                return True

        # check end of seq matches
        for end_len in range(1, min(len(seq) + 1, self.motif_len)):
            try:
                end_pat, end_mod_pos = self._partial_pats['end'][end_len]
            except KeyError:
                continue
            if end_pat.match(seq[-end_len:]):
                return True

        return False

    def find_mod_poss(self, seq):
        """Find all mod-base positions within the sequence.

        Including partial matches at beginning and end that include mod_pos.
        """
        seq_mod_poss = set()
        # check matches to start of sequence
        for start_len in range(1, min(len(seq) + 1, self.motif_len)):
            try:
                start_pat, start_mod_pos = self._partial_pats[
                    'start'][start_len]
            except KeyError:
                continue
            if start_pat.match(seq[:start_len]):
                seq_mod_poss.add(start_mod_pos)

        # check central sequence overlaps
        if len(seq) < self.motif_len:
            for short_pat, short_mod_pos in self._partial_pats[
                    'short'][len(seq)]:
                if short_pat.match(seq):
                    seq_mod_poss.add(short_mod_pos)
        else:
            for motif_match in self.motif_pat.finditer(seq):
                seq_mod_poss.add(motif_match.start() + self.mod_pos)

        # check end of seq matches
        for end_len in range(1, min(len(seq) + 1, self.motif_len)):
            try:
                end_pat, end_mod_pos = self._partial_pats['end'][end_len]
            except KeyError:
                continue
            if end_pat.match(seq[-end_len:]):
                seq_mod_poss.add(len(seq) - end_len + end_mod_pos)

        return sorted(seq_mod_poss)


def parse_motif_descs(stat_motif_descs):
    """Parse string motif descriptions as defined by ``tombo plot roc --motif-descriptions``

    Args:
        stat_motif_descs (str): string motifs description (see ``tombo plot roc --motif-descriptions``)

    Returns:
        list of tuples with :class:`tombo.tombo_helper.TomboMotif` and motif/modification names
    """
    parsed_motif_descs = []
    try:
        for motif_desc in stat_motif_descs.split('::'):
            raw_motif, mod_pos, mod_name = motif_desc.split(':')
            motif = TomboMotif(raw_motif, int(mod_pos))
            parsed_motif_descs.append((motif, mod_name))
    except:
        error_message_and_exit(
            'Invalid motif decriptions format. Format descriptions as: ' +
            '"motif:mod_pos:name[::motif2:mod_pos2:name2...]".')

    return parsed_motif_descs


###########################
###### FASTA Parsing ######
###########################

def get_rec_names(fasta_fn):
    with io.open(fasta_fn) as fasta_fp:
        all_rec_ids = [line.replace(">","").split()[0] for line in fasta_fp
                       if line.startswith('>')]

    return all_rec_ids

class Fasta(object):
    """Fasta file sequence format wrapper class. Will load faidx via ``pyfaidx`` package if installed, else the fasta will be loaded into memory for sequence extraction.

    .. automethod:: __init__
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
        """Load a fasta

        Args:
            fasta_fn (str): path to fasta file
            dry_run (bool): when pyfaidx is not installed, don't actually read sequence into memory.
            force_in_mem (bool): force genome to be loaded into memory even if pyfaidx is installed allowing on-disk access
            assume_dna_bases (bool): skip check for DNA or RNA bases (default: False)
        """
        self.fasta_fn = resolve_path(fasta_fn)
        self.has_rna_bases = False
        try:
            if force_in_mem: raise ImportError
            import pyfaidx
            self.has_pyfaidx = True
            try:
                self.index = pyfaidx.Faidx(self.fasta_fn)
            except UnicodeDecodeError:
                error_message_and_exit(
                    'FASTA file does not appear to be formatted correctly.')
        except ImportError:
            self.has_pyfaidx = False
            if not dry_run:
                self.index = self._load_in_mem()
        if not dry_run:
            self.has_rna_bases = (assume_dna_base or
                                  self._index_contains_uridines())

    def get_seq(self, chrm, start=None, end=None, error_end=True):
        """Extract sequence from a specific genomic region. Note if provided, start and end must both be provided or they will be ignored.

        Args:

            chrm (str): chromosome name
            start (int): 0-based start position
            end (int): 0-based open-interval end position
            error_end (bool): raise an error when the region requested extends beyond the chromosome (default: True)

        Returns:
            Genomic sequence requested. Sequence is converted to RNA sequence if applicable.
        """
        if self.has_pyfaidx:
            if not (start or end):
                r_seq = self.index.fetch(
                    chrm, 1, self.index.index[chrm].rlen).seq.upper()
            elif (start < 0 or start > self.index.index[chrm].rlen or (
                    error_end and (
                        end < 0 or end > self.index.index[chrm].rlen))):
                raise TomboError(
                    'Encountered invalid genome sequence request.')
            else:
                r_seq = self.index.fetch(chrm, start + 1, end).seq.upper()
        else:
            if (start is not None and (
                    start < 0 or start > len(self.index[chrm]))) or (
                        error_end and end is not None and
                        (end < 0 or
                         end > len(self.index[chrm]))):
                raise TomboError(
                    'Encountered invalid genome sequence request.')
            r_seq = self.index[chrm][start:end].upper()

        if self.has_rna_bases:
            r_seq = rev_transcribe(r_seq)

        return r_seq

    def iter_chrms(self):
        """Iterate over chromosome names
        """
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
    """Determine if a read is RNA or DNA

    Args:
        fast5_data (`h5py.File`)

    Returns:
        Boolean, incdicating whether the read appears to be RNA or DNA

    Note:
        This function uses the read meta-data and so non-standard processing pipelines may not prodcue expected values.
    """
    # check both experiment type and kit slots for "rna"
    exp_type = fast5_data['UniqueGlobalKey/context_tags'].attrs.get(
        'experiment_type')
    try:
        exp_type = exp_type.decode()
        # remove the word internal since it contains rna.
        exp_type = exp_type.replace('internal', '')
    except (AttributeError, TypeError):
        pass

    exp_kit = fast5_data['UniqueGlobalKey/context_tags'].attrs.get(
        'experiment_kit')
    try:
        exp_kit = exp_kit.decode()
        # remove the word internal since it contains rna.
        exp_kit = exp_kit.replace('internal', '')
    except (AttributeError, TypeError):
        pass

    if exp_type is None and exp_kit is None:
        return False

    return (
        (exp_type is not None and re.search('rna', exp_type) is not None) or
        (exp_kit is not None and re.search('rna', exp_kit) is not None))

def is_sample_rna(reads_index=None, fast5_fns=None, n_reads=50):
    """Determine if a set of reads are RNA or DNA from a small sample. Must provide either reads_index or fast5_fns.

    Args:
        reads_index (:class:`tombo.tombo_helper.TomboReads`)
        fast5_fns (list): list of fast5 read filename
        n_reads (int): number of reads to check (default: 50)

    Returns:
        False if any read does not appear to be an RNA read (see :class:`tombo.tombo_helper.is_read_rna`) else return True.
    """
    proc_reads = 0
    if reads_index is not None:
        for r_data in reads_index.iter_reads():
            if not r_data.rna:
                return False
            proc_reads += 1
            if proc_reads >= n_reads:
                break
    elif fast5_fns is not None:
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
    else:
        raise TomboError(
            'Must provide either reads_index or fast5_fns to ' +
            'determine is sample is RNA.')

    return True

def get_seq_sample_type(fast5_data=None, reads_index=None, fast5_fns=None,
                        num_reads=50):
    """Get the sequencing sample type from a single read or set of reads

    Args:
        fast5_data (`h5py.File`): open read h5py File object
        reads_index (:class:`tombo.tombo_helper.TomboReads`)
        fast5_fns (list): FAST5 read filenames
        num_reads (int): sample of reads to check for sample type
    """
    if fast5_data is None and reads_index is None and fast5_fns is None:
        raise TomboError('Must provide either fast5_data, reads_index, or ' +
                         'fast5_fns to determine sequencing sample type.')
    if fast5_data is None:
        return seqSampleType(RNA_SAMP_TYPE, True) if is_sample_rna(
            reads_index, fast5_fns, num_reads) else seqSampleType(
                DNA_SAMP_TYPE, False)
    return seqSampleType(RNA_SAMP_TYPE, True) if is_read_rna(
        fast5_data) else seqSampleType(DNA_SAMP_TYPE, False)


############################
###### Lock Functions ######
############################

def get_lock_fn(fast5s_dir):
    """Get filename for the lock file to indicate that this directory is currently being processed. This file should be saved to be deleted later.
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
    """Get all fast5 files recursively below this directory
    """
    all_fast5s = []
    # walk through directory structure searching for fast5 files
    for root, _, fns in os.walk(fast5s_dir):
        for fn in fns:
            if not fn.endswith('.fast5'): continue
            all_fast5s.append(os.path.join(root, fn))

    return all_fast5s

def clear_tombo_locks(lock_fns):
    """Clear all lock files
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
    """Get all fast5 files recursively below this directory and add a Tombo lock file to indicate that this directory is currently being re-squiggled
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
                error_message_and_exit(ignore_locks_mess)
            lock_fns.append(lock_fn)
            # create empty file indicating this directory is locked
            open(lock_fn, 'w').close()

            for fn in fns:
                if not fn.endswith('.fast5'): continue
                all_fast5s.append(os.path.join(root, fn))
    except:
        clear_tombo_locks(lock_fns)
        error_message_and_exit(
            'Unexpected error during file enumeration. Check that you have ' +
            'write permission within the specified [fast5_basedir].')

    return all_fast5s, lock_fns

def get_raw_read_slot(fast5_data):
    """Get the raw read slot from this FAST5 read file

    Args:
        fast5_data (`h5py.File`): open FAST5 read file object

    Example::

        all_aw_signal = get_raw_read_slot(fast5_data)['Signal'][:]

    Returns:
        The HDF5 group slot containing the raw signal data.
    """
    try:
        raw_read_slot = next(iter(fast5_data['/Raw/Reads'].values()))
    except KeyError:
        raise TomboError(
            'Raw data is not found in /Raw/Reads/Read_[read#]. Note that ' +
            'Tombo does not support multi-fast5 format.')

    return raw_read_slot

class TomboReads(object):
    """A set of reads with associated meta-data from re-squiggle processing

    .. automethod:: __init__
    """
    def _get_index_fn(self, fast5s_dir):
        """Get the filename for the requested directory and corrected group
        """
        # if directory comes with trailing slash, remove for processing
        if fast5s_dir.endswith('/'):
            fast5s_dir = fast5s_dir[:-1]
        split_dir = os.path.split(fast5s_dir)
        return os.path.join(split_dir[0], "." + split_dir[1] +
                            "." + self.corr_grp + '.tombo.index')


    # index building and writing class functions
    def _prep_for_writing(self, fast5s_dirs):
        assert len(fast5s_dirs) == 1, (
            'Must provide only a single FAST5 base directory when ' +
            'openning for writing.')
        fast5s_dir = fast5s_dirs[0]
        fast5s_dir = (fast5s_dir if fast5s_dir.endswith('/') else
                      fast5s_dir + '/')
        index_fn = self._get_index_fn(fast5s_dir)
        self.fast5s_dirs[fast5s_dir] = index_fn
        if os.path.exists(index_fn): os.remove(index_fn)

        # open default dict to fill with readData lists by (chrm, strand)
        self.reads_index = defaultdict(list)

        return

    def add_read_data(self, chrm, strand, read_data):
        """Add read data to the index

        Args:
            chrm (str): chromosome name
            strand (str): strand ('+' or '-')
            read_data (:class:`tombo.tombo_helper.readData`): read information
        """
        self.reads_index[(chrm, strand)].append(read_data)
        return

    def replace_index(self, new_reads_index):
        """Replace current reads index

        Args:
            new_reads_index (dict): dictionary with (chrm, strand) pointing to lists of :class:`tombo.tombo_helper.readData` objects
        """
        if sum(len(x) for x in new_reads_index.values()) == 0:
            raise TomboError('Cannot replace with an empty index.')

        self.reads_index = new_reads_index
        return

    def write_index_file(self):
        """Write index file

        Note:

            Must be an index from only a single ``fast5s_dir``

            Index filename will be: `.[fast5s_dir].[corr_grp].tombo.index`
            Note that this is a hidden file (starts with a ".")
        """
        status_message('Saving Tombo reads index to file.')
        assert len(self.fast5s_dirs) == 1, (
            'Cannot write index for TomboReads created from more than ' +
            'one base directory.')
        basedir, index_fn = next(iter(self.fast5s_dirs.items()))
        try:
            import cPickle as pickle
        except:
            import pickle

        index_data = defaultdict(list)
        for chrm_strand, cs_reads in self:
            for rd in cs_reads:
                # clip the basedir off the FAST5 filename in case later
                # functions are called from another relative path and
                # split corr_grp and bc_subgrp for easier filtering
                index_data[chrm_strand].append((
                    re.sub(basedir, '', rd.fn, 1), rd.start, rd.end,
                    rd.read_start_rel_to_raw, rd.corr_group.split('/')[0],
                    rd.corr_group.split('/')[-1], rd.filtered, rd.rna,
                    rd.sig_match_score, rd.mean_q_score, rd.read_id))

        with io.open(index_fn, 'wb') as index_fp:
            # note protocol 2 for py2/3 compatibility
            pickle.dump(dict(index_data), index_fp, protocol=2)

        return


    # Index parsing class functions
    def _parse_fast5s_wo_index(self, wo_index_fast5s_dirs):
        """Parse re-squiggled reads data from a list of fast5 directories
        """
        def get_read_data(read_fn, fast5_data, bc_subgrp):
            read_id = get_raw_read_slot(fast5_data).attrs.get('read_id')
            corr_data = fast5_data[
                '/'.join(('/Analyses', self.corr_grp, bc_subgrp))]
            rna = corr_data.attrs.get('rna')
            rna = False if rna is None else rna

            align_data = dict(corr_data['Alignment'].attrs.items())
            read_start_rel_to_raw = corr_data['Events'].attrs.get(
                'read_start_rel_to_raw')
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
                self.corr_grp + '/' + bc_subgrp, rna, read_id=read_id)


        files = [fn for fast5s_dir in wo_index_fast5s_dirs
                 for fn in get_files_list(fast5s_dir)]
        dir_reads_index = defaultdict(list)
        for read_fn in files:
            try:
                with h5py.File(read_fn, 'r') as fast5_data:
                    i_bc_subgrps = (
                        fast5_data['/Analyses/' + self.corr_grp].keys()
                        if self.bc_subgrps is None else self.bc_subgrps)
                    for bc_subgrp in i_bc_subgrps:
                        chrm, strand, r_data = get_read_data(
                            read_fn, fast5_data, bc_subgrp)
                        dir_reads_index[(chrm, strand)].append(r_data)
            except:
                # ignore errors and process all reads that don't error
                continue

        return dict(dir_reads_index)

    def _load_index_data(self, fast5s_dir):
        try:
            import cPickle as pickle
        except:
            import pickle
        with io.open(self.fast5s_dirs[fast5s_dir], 'rb') as index_fp:
            raw_index_data = pickle.load(index_fp)

        # determine the index type used. Index information was added around
        # version 1.3 making the index data 2 elements longer.
        # so check which one is used here.
        try:
            num_index_vals = len(next(iter(raw_index_data.values()))[0])
        except StopIteration:
            raise TomboError('Tombo index file appears to be empty')
        if num_index_vals == 8:
            def convert_r_data(from_base_fn, start, end, rsrtr,
                               c_grp, s_grp, filtered, rna):
                return readData(start, end, filtered, rsrtr, strand,
                                os.path.join(fast5s_dir, from_base_fn),
                                self.corr_grp + '/' + s_grp, rna)
        elif num_index_vals == 10:
            def convert_r_data(
                    from_base_fn, start, end, rsrtr, c_grp, s_grp, filtered, rna,
                    sig_match_score, mean_q_score):
                return readData(start, end, filtered, rsrtr, strand,
                                os.path.join(fast5s_dir, from_base_fn),
                                self.corr_grp + '/' + s_grp, rna,
                                sig_match_score, mean_q_score)
        elif num_index_vals == 11:
            def convert_r_data(
                    from_base_fn, start, end, rsrtr, c_grp, s_grp, filtered, rna,
                    sig_match_score, mean_q_score, read_id):
                return readData(start, end, filtered, rsrtr, strand,
                                os.path.join(fast5s_dir, from_base_fn),
                                self.corr_grp + '/' + s_grp, rna,
                                sig_match_score, mean_q_score, read_id)
        else:
            raise TomboError('Invalid Tombo index file.')

        dir_reads_index = {}
        for (chrm, strand), cs_raw_data in raw_index_data.items():
            cs_data = [convert_r_data(*r_data) for r_data in cs_raw_data]
            # don't add chrm/strand if all reads are filtered
            if len(cs_data) > 0:
                dir_reads_index[(chrm, strand)] = cs_data

        return dir_reads_index

    def _parse_fast5s_w_index(self, fast5s_dir):
        """Use index file to parse information about a set of reads
        """
        try:
            curr_dir_reads_index = self._load_index_data(fast5s_dir)
        except UnicodeDecodeError:
            warning_message(
                'Invalid Tombo index file.\n\t\tThis occurs most often ' +
                'when the re-squiggle command was completed using a Tombo ' +
                'build against a different python version (2 or 3).')
            raise TomboError

        if not self.remove_filtered and self.bc_subgrps is None:
            return curr_dir_reads_index

        filt_dir_reads_index = {}
        for (chrm, strand), cs_raw_data in curr_dir_reads_index.items():
            cs_data = [
                rd for rd in cs_raw_data
                if rd.corr_group.split('/')[0] == self.corr_grp and
                (self.bc_subgrps is None or
                 rd.corr_group.split('/')[-1] in self.bc_subgrps) and
                (not self.remove_filtered or not rd.filtered)]
            # don't add chrm/strand if all reads are filtered
            if len(cs_data) > 0:
                filt_dir_reads_index[(chrm, strand)] = cs_data

        return filt_dir_reads_index

    def _merge_dir_indices(self, w_index_ri, wo_index_ri):
        """Merge coverage from serveral parsed sets of data
        """
        all_raw_reads_index = w_index_ri + [wo_index_ri,]
        reads_index = defaultdict(list)
        for chrm_strand in set([cs for d_ri in all_raw_reads_index
                                for cs in d_ri]):
            for dir_reads_index in all_raw_reads_index:
                if chrm_strand not in dir_reads_index: continue
                reads_index[chrm_strand].extend(dir_reads_index[chrm_strand])

        return dict(reads_index)

    def _parse_fast5s(self, fast5s_dirs):
        wo_index_dirs = []
        w_index_covs = []
        warn_index = False
        # determine if index exists for each directory and load appropriately
        for fast5s_dir in fast5s_dirs:
            fast5s_dir = (fast5s_dir if fast5s_dir.endswith('/') else
                          fast5s_dir + '/')
            self.fast5s_dirs[fast5s_dir] = self._get_index_fn(fast5s_dir)
            if os.path.exists(self.fast5s_dirs[fast5s_dir]):
                try:
                    w_index_covs.append(self._parse_fast5s_w_index(fast5s_dir))
                except TomboError:
                    warning_message(
                        'Failed to parse tombo index file for ' + fast5s_dir +
                        ' directory. Creating temporary index from ' +
                        'FAST5 files.')
                    wo_index_dirs.append(fast5s_dir)
            else:
                if not warn_index:
                    warning_message(
                        'Tombo index file does not exist for one or more ' +
                        'directories.\n\t\tIf --skip-index was not set for ' +
                        're-squiggle command, ensure that the specified ' +
                        'directory is the same as for the re-squiggle command.')
                    warn_index = True
                wo_index_dirs.append(fast5s_dir)

        wo_index_cov = self._parse_fast5s_wo_index(wo_index_dirs)
        self.reads_index = self._merge_dir_indices(w_index_covs, wo_index_cov)

        return

    def __init__(self, fast5s_basedirs, corrected_group='RawGenomeCorrected_000',
                 basecall_subgroups=None, for_writing=False, remove_filtered=True,
                 sample_name=None):
        """Parse data from a list of re-squiggle fast5 directories

        Args:
            fast5s_basedirs (list): fast5 base directories, which have been processed by ``tombo resquiggle``
            corrected_group (str): Analysis slot containing the re-squiggle information (optional; default: 'RawGenomeCorrected_000')
            basecall_subgroups (list): basecall subgroups (optional; default: Process all basecall subgroups)
            for_writing (bool): open TomboReads to write index (optional; default: Open for reading and parse re-squiggled reads)
            remove_filtered (bool): remove filtered reads as indicated by index file (optional; default: True)
            sample_name (str): for verbose output
        """
        if VERBOSE and not for_writing:
            status_mess = ('Parsing Tombo index file(s).'
                           if sample_name is None else
                           'Parsing ' + sample_name + ' Tombo index file(s).')
            status_message(status_mess)
        assert isinstance(fast5s_basedirs, list), (
            'fast5s_basedirs must be a list.')

        self.fast5s_dirs = {}
        self.corr_grp = corrected_group
        self.bc_subgrps = basecall_subgroups
        self.sample_name = sample_name
        self.remove_filtered = remove_filtered
        self.for_writing = for_writing
        self.coverage = None
        if self.for_writing:
            self._prep_for_writing(fast5s_basedirs)
        else:
            self._parse_fast5s(fast5s_basedirs)

        return

    def _compute_coverage(self):
        if VERBOSE: status_message('Calculating read coverage.')
        self.coverage = {}
        for chrm_strand, cs_reads in self:
            if len(cs_reads) == 0: continue
            cs_coverage = np.zeros(max(r_data.end for r_data in cs_reads),
                                   dtype=np.int64)
            for r_data in cs_reads:
                cs_coverage[r_data.start:r_data.end] += 1
            self.coverage[chrm_strand] = cs_coverage
        return self

    def _add_coverages(self, ctrl_reads_index):
        merged_cov = {}
        # if self or control reads don't have coverage, compute it
        if self.coverage is None:
            self._compute_coverage()
        if ctrl_reads_index.coverage is None:
            ctrl_reads_index._compute_coverage()

        for chrm_strand, ctrl_cs_cov in ctrl_reads_index.coverage.items():
            if chrm_strand in self.coverage:
                self_cs_cov = self.coverage[chrm_strand]
                # copy longer array and add shorter array coverage
                if self_cs_cov.shape[0] > ctrl_cs_cov.shape[0]:
                    merged_cs_cov = self_cs_cov.copy()
                    merged_cs_cov[:ctrl_cs_cov.shape[0]] += ctrl_cs_cov
                else:
                    merged_cs_cov = ctrl_cs_cov.copy()
                    merged_cs_cov[:self_cs_cov.shape[0]] += self_cs_cov
            else:
                merged_cs_cov = ctrl_cs_cov.copy()
            merged_cov[chrm_strand] = merged_cs_cov

        return merged_cov

    def iter_coverage_regions(self, ctrl_reads_index=None):
        """Get genome coverage for a set of reads

        Args:

            ctrl_reads_index (:class:`tombo.tombo_helper.TomboReads`): a second set of tombo reads to add for coverage

        Returns:
            Yeilds (chrm, strand, cs_cov, cs_cov_starts) indicating read coverage levels

            cs_cov and cs_cov_starts are numpy arrays containing the coverage level and the 0-based genomic start positions of those intervals
        """
        if VERBOSE: status_message('Calculating read coverage regions.')
        if self.coverage is None:
            self._compute_coverage()

        coverage = (self.coverage if ctrl_reads_index is None else
                    self._add_coverages(ctrl_reads_index))
        for (chrm, strand), cs_cov in coverage.items():
            cs_cov_starts = np.concatenate([
                [0,], np.where(np.diff(cs_cov))[0] + 1,
                [cs_cov.shape[0],]])
            cs_cov = cs_cov[cs_cov_starts[:-1]]
            yield chrm, strand, cs_cov, cs_cov_starts

        return

    def iter_cov_regs(
            self, cov_thresh, region_size=None, ctrl_reads_index=None):
        """Iterate over regions with coverage greater than or equal to cov_thresh.

        If region_size is provided, regions are rounded to the nearest region_sized windows and only the region start is yielded (e.g. chrm, strand, start). If not provided (chrm, strand, start, end) is yielded.
        """
        def round_reg_start(x):
            return int(region_size * np.floor(x / float(region_size)))
        def round_reg_end(x):
            return int(region_size * np.ceil(x / float(region_size)))


        for chrm, strand, cov, starts in self.iter_coverage_regions(
                ctrl_reads_index):
            curr_reg_start = -1
            valid_cov = np.where(np.diff(np.concatenate([
                [False,], np.greater_equal(cov, cov_thresh), [False,]])))[0]
            for cov_start_i, cov_end_i in zip(valid_cov[:-1], valid_cov[1:]):
                cov_start, cov_end = starts[cov_start_i], starts[cov_end_i]
                if region_size is None:
                    yield chrm, strand, cov_start, cov_end
                    continue
                for reg_start in range(round_reg_start(cov_start),
                                       round_reg_end(cov_end), region_size):
                    if reg_start != curr_reg_start:
                        yield chrm, strand, reg_start
                        curr_reg_start = reg_start

        return

    def get_all_cs(self):
        """Get list of all (chromosome, strand) stored in index.
        """
        return list(self.reads_index.keys())

    def is_empty(self):
        """Are any reads stored in the index?
        """
        for cs_reads in self.reads_index.values():
            if len(cs_reads) > 0:
                return False
        return True

    def __contains__(self, chrm_strand):
        return chrm_strand in self.reads_index

    def __iter__(self):
        self._iter = iter(self.reads_index.items())
        return self

    def __next__(self):
        return next(self._iter)

    # for python2 compatibility
    def next(self):
        return self.__next__()

    def iter_reads(self):
        """Iterate over reads stored in the index
        """
        for _, cs_reads in self:
            for rd in cs_reads:
                yield rd
        return

    def get_cs_reads(self, chrm, strand, invalid_return=[]):
        """Extract the list of reads stored in a specified chromosome and strand

        Args:

            chrm (str): chromosome name
            strand (str): strand ('+' or '-')
            invalid_return: value to return for invalid (chrm, strand)
        """
        if (chrm, strand) not in self.reads_index:
            return invalid_return
        return self.reads_index[(chrm, strand)]

    def _get_strand_spec_cov(self, chrm, pos, strand, invalid_return):
        if (chrm, strand) not in self.coverage:
            return invalid_return
        if pos >= self.coverage[(chrm, strand)].shape[0]:
            return invalid_return
        return self.coverage[(chrm, strand)][pos]

    def get_coverage(self, chrm, pos, strand=None, invalid_return=0):
        """Get coverage at specified position

        Args:
            chrm (str): chromosome name
            pos (int): 0-based genomic position
            strand (str): interval strand ('+', '-' or None). Default: None (max over both strands)
            invalid__return: return value for invalid (bad chrm, strand or pos beyond coverage) position. Default: 0
        """
        if self.coverage is None:
            self._compute_coverage()
        if strand is None:
            try:
                return max(
                    self._get_strand_spec_cov(chrm, pos, '+', invalid_return),
                    self._get_strand_spec_cov(chrm, pos, '-', invalid_return))
            # with None invalid_return value could be max over 2 None's
            except TypeError:
                return invalid_return
        return self._get_strand_spec_cov(chrm, pos, strand, invalid_return)

    def get_cs_coverage(self, chrm, strand, invalid_return=None):
        """Extract coverage levels over a specified chromosome and strand

        Args:

            chrm (str): chromosome name
            strand (str): strand ('+' or '-')
            invalid_return: value to return for invalid (chrm, strand)

        Returns:
            numpy array (numpy.int64) containing read coverage levels
        """
        if self.coverage is None:
            self._compute_coverage()
        if (chrm, strand) not in self.coverage:
            return invald_return
        return self.coverage[(chrm, strand)]

    def iter_cs_coverage(self):
        """Iterate over coverage levels across all stored (chrm, strands)
        """
        if self.coverage is None:
            self._compute_coverage()
        return self.coverage.items()


###########################################
###### Events Table Access Functions ######
###########################################

def get_multiple_slots_read_centric(r_data, slot_names, corr_grp=None):
    """Extract multiple read-centric slot_names from this read's Events table

    Args:

        fast5_data (:class:`tombo.tombo_helper.readData` or an open ``h5py.File`` object)
        slot_name (str): slot from which to extract data (valid values: norm_mean, norm_stdev, start, length, base)
        corr_grp (str): corrected group slot from which to extract data (default: use value from ``fast5_data`` object; required for ``h5py.File``)

    Returns:
        A tuple of numpy arrays specified by the ``slot_names``
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
        r_event_data = r_data[event_slot_name][:]
        if do_close: r_data.close()
    except:
        # probably truncated file or events don't exist
        return [None,] * len(slot_names)

    return [r_event_data[slot_name] for slot_name in slot_names]

def get_single_slot_read_centric(r_data, slot_name, corr_grp=None):
    """Extract read-centric slot_name from this read's Events table

    Args:

        fast5_data (:class:`tombo.tombo_helper.readData` or an open ``h5py.File`` object)
        slot_name (str): slot from which to extract data (valid values: norm_mean, norm_stdev, start, length, base)
        corr_grp (str): corrected group slot from which to extract data (default: use value from ``fast5_data`` object; required for ``h5py.File``)
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
    """Extract genome-centric slot_name from this read's Events table
    """
    r_slot_values = get_single_slot_read_centric(r_data, slot_name)
    if r_slot_values is None:
        return None

    if r_data.strand == '-':
        r_slot_values = r_slot_values[::-1]

    return r_slot_values

def get_mean_slot_genome_centric(cs_reads, chrm_len, slot_name):
    """Get the mean over all reads at each covered genomic location for this slots value
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

def iter_mean_slot_values(
        reads_index, chrm_sizes, slot_name, ctrl_reads_index=None):
    """Iterate through chromosomes and strands yielding mean slots values over
    all reads at each covered genomic location.

    Generator returns chrmosome, strand, cs_mean_values,
    cs_ctrl_mean_values tuples (4 return values).
    """
    # ignore divide by zero errors that occur where there is no
    # coverage. Need to correct nan values after subtracting two sets of
    # coverage so leave as nan for now
    old_err_settings = np.seterr(all='ignore')
    for chrm, strand in [(c, s) for c in sorted(chrm_sizes) for s in ('+', '-')]:
        if ctrl_reads_index is None:
            if (chrm, strand) not in reads_index: continue
            cs_mean_values = get_mean_slot_genome_centric(
                reads_index.get_cs_reads(chrm, strand),
                chrm_sizes[chrm], slot_name)
            yield chrm, strand, cs_mean_values, None
        else:
            cs_mean_values, ctrl_cs_mean_values = None, None
            if (chrm, strand) in reads_index:
                cs_mean_values = get_mean_slot_genome_centric(
                    reads_index.get_cs_reads(chrm, strand),
                    chrm_sizes[chrm], slot_name)
            if (chrm, strand) in ctrl_reads_index:
                ctrl_cs_mean_values = get_mean_slot_genome_centric(
                    ctrl_reads_index.get_cs_reads(chrm, strand),
                    chrm_sizes[chrm], slot_name)
            if cs_mean_values is None and ctrl_cs_mean_values is None: continue
            yield chrm, strand, cs_mean_values, ctrl_cs_mean_values

    _ = np.seterr(**old_err_settings)

    return

def get_largest_signal_differences(
        reads_index, ctrl_reads_index, num_regions, num_bases):
    chrm_sizes = get_chrm_sizes(reads_index, ctrl_reads_index)

    all_largest_diff_poss = []
    for chrm, strand, cs_sig_means, ctrl_cs_sig_means in iter_mean_slot_values(
            reads_index, chrm_sizes, 'norm_mean', ctrl_reads_index):
        if cs_sig_means is None or ctrl_cs_sig_means is None: continue
        chrm_diffs = np.nan_to_num(np.abs(cs_sig_means - ctrl_cs_sig_means))
        chrm_max_diff_regs = np.argsort(chrm_diffs)[::-1][:num_regions]
        all_largest_diff_poss.extend((
            chrm_diffs[pos], max(pos - int(num_bases / 2.0), 0),
            chrm, strand) for pos in chrm_max_diff_regs)

    return sorted(all_largest_diff_poss, reverse=True)[:num_regions]

def get_signal_differences(reads_index, ctrl_reads_index):
    """Helper function to compute all signal differences
    """
    chrm_sizes = get_chrm_sizes(reads_index, ctrl_reads_index)

    all_diffs =  {}
    for chrm, strand, cs_sig_means, ctrl_cs_sig_means in iter_mean_slot_values(
            reads_index, chrm_sizes, 'norm_mean', ctrl_reads_index):
        if cs_sig_means is None or ctrl_cs_sig_means is None: continue
        all_diffs[(chrm, strand)] = np.nan_to_num(
            cs_sig_means - ctrl_cs_sig_means)

    return all_diffs


####################################
###### Genomic Interval Class ######
####################################

class intervalData(object):
    """Genome/transcriptome interval information

    Example::

        int_data = tombo_helper.intervalData(
            chrm='chr20', start=10000, end=10100, strand='+')

    Note:
        All intervalData functions return the object in order to allow chaining of interval commands. (e.g. ``int_data.add_reads(reads_index).add_seq()``)

    .. automethod:: __init__
    """
    def __setattr__(self, name, value):
        if name == 'chrm':
            if not isinstance(value, unicode):
                raise TypeError(name + ' must be a string')
        elif name in ('start', 'end'):
            if not (isinstance(value, int) or isinstance(value, np.integer)):
                raise TypeError(name + ' must be an int')
        elif name == 'strand':
            if value not in (None, '+', '-'):
                raise TypeError('strand must be either None, "+", or "-"')
        elif name in ('reg_id', 'reg_text', 'seq'):
            if value is not None and not isinstance(value, unicode):
                raise TypeError(name + ' must be a string')
        elif name == 'reads':
            if value is not None and not isinstance(value, list):
                raise TypeError('reads must be a list')
        else:
            raise ValueError(name + ' is an invalid attribute for intervalData')
        super(intervalData, self).__setattr__(name, value)

    def __init__(self, chrm, start, end, strand=None, reg_id=None,
                 reg_text='', reads=None, seq=None):
        """Initialize a genome/transcriptome interval object

        Args:
            chrm (str): chromosome name.
            start (int): 0-based start position.
            end (int): 1-based (or open interval) end position.
            strand (str): interval strand ('+', '-' or None). Default: None
            reg_id (str): used to keep track of multiple intervals. Default: None
            reg_text (str): region description text. Default: ''
            reads (list): list of readData values overlapping region. Default: None
            seq (str): region genomic sequence. Default: None
        """
        self.chrm = unicode(chrm)
        self.start = start
        self.end = end
        self.strand = unicode(strand) if strand is not None else None
        self.reg_id = unicode(reg_id) if reg_id is not None else None
        self.reg_text = unicode(reg_text)
        self.reads = reads
        self.seq = seq

    def update(self, **kwargs):
        """Update slots specified in keyword arguments (kwargs)
        """
        for k, v in kwargs.items():
            self.__setattr__(k, v)
        # return self to allow chaining and auto-return
        return self

    def __repr__(self):
        # convert reads and seq if they are too long
        self_vars = vars(self).copy()
        if self_vars['reads'] is not None:
            self_vars['reads'] = '<{:d} reads>'.format(len(self_vars['reads']))
        if self_vars['seq'] is not None and len(self_vars['seq']) > 50:
            self_vars['seq'] = '<{:d} bases of sequence>'.format(
                len(self_vars['seq']))
        self_vars['reg_text'] = '"' + self_vars['reg_text'] + '"'
        return '<tombo.tombo_helper.intervalData object> :\n' + '\n'.join(
            "{:>15} : {}".format(k, str(v)) for k, v in self_vars.items())

    def copy(self, include_reads=True):
        """Create a copy of this interval. Useful when adding sets of reads from multiple samples to compare.

        Args:
            include_reads (bool): include current reads slot in the new object (default: True)
        """
        cp_reads = self.reads if include_reads else None
        return type(self)(self.chrm, self.start, self.end, self.strand,
                          self.reg_id, self.reg_text, cp_reads, self.seq)

    def merge(self, other_reg):
        """Create a copy of this interval with the reads from this interval and `other_reg`

        Args:
            other_reg (:class:`tombo.tombo_helper.intervalData`): a second region to merge with this interval

        Note:
            Aside from reads, all other attributes will be taken from this interval
        """
        merged_reg_data = self.copy()
        self_reads = self.reads if self.reads is not None else []
        other_reads = other_reg.reads if other_reg.reads is not None else []
        return merged_reg_data.update(reads=self_reads + other_reads)

    def expand_interval(self, expand_width):
        """Expand this interval by the specified amount (effects only the start and stop attributes; NOT seq or reads)

        Args:

            expand_width (int): amount by which to expand the interval
        """
        self.update(start=self.start - expand_width,
                    end=self.end + expand_width)
        return self

    def add_reads(self, reads_index, require_full_span=False):
        """Add reads overlapping  this interval

        Args:
            reads_index (:class:`tombo.tombo_helper.TomboReads`): reads index
            require_full_span (bool): require that reads span then entire interval (default: include all reads overlapping the interval)
        """
        def get_c_s_data(strand):
            if require_full_span:
                def read_intersects_interval(r_start, r_end):
                    return r_start <= self.start and r_end >= self.end
            else:
                def read_intersects_interval(r_start, r_end):
                    return not (r_start >= self.end or r_end <= self.start)

            # get all reads intersecting the interval
            return [
                r_data for r_data in reads_index.get_cs_reads(self.chrm, strand)
                if read_intersects_interval(r_data.start, r_data.end)]


        # get all reads that overlap this interval
        # note that this includes partial overlaps as these contribute
        # to coverage and other statistics so can't really restrict to
        # full coverage as previous versions of code did
        if self.strand is not None:
            return self.update(reads=get_c_s_data(self.strand))

        # if strand is None, get data from both strands
        return self.update(reads=get_c_s_data('+') + get_c_s_data('-'))

    def _update_seq(self, r_data, reg_base_data):
        """Update the sequence for the region based on this read
        """
        read_bases = get_single_slot_read_centric(r_data, 'base')
        if read_bases is None:
            warning_message(
                'Unable to extract data from read. Potentially corrupted file ' +
                'or invalid Tombo index file for this directory.')
            return reg_base_data, max(0, r_data.start - self.start)
        r_seq = b''.join(read_bases).decode()

        if r_data.strand == '-':
            r_seq = rev_comp(r_seq)

        # if read starts before the interval
        if r_data.start <= self.start:
            r_end_overlap = r_data.end - self.start
            # if read covers the whole interval
            if r_data.end > self.end:
                r_end_clip = r_data.end - self.end
                reg_base_data = r_seq[-r_end_overlap:-r_end_clip]
                return reg_base_data, len(reg_base_data)
            # end of read overlaps beginning of interval
            reg_base_data[:r_end_overlap] = r_seq[-r_end_overlap:]
            return reg_base_data, r_end_overlap
        # read doesn't cover the beginning of region
        if r_data.end > self.end:
            # beginning of read covers to the end of the region
            r_begin_overlap = self.end - r_data.start
            reg_base_data[-r_begin_overlap:] = r_seq[:r_begin_overlap]
            return reg_base_data, len(reg_base_data)
        # first read is completely contained in the interval
        r_len = r_data.end - r_data.start
        r_int_start = r_data.start - self.start
        reg_base_data[r_int_start:r_int_start + r_len] = r_seq
        return reg_base_data, r_int_start + r_len

    def add_seq(self, genome_index=None, error_end=True):
        """Extract the forward strand genomic sequence for an interval from reads or genome_index if provided

        Args:
            genome_index (:class:`tombo.tombo_helper.Fasta`) Tombo FASTA sequence object
        """
        if genome_index is not None:
            return self.update(seq=genome_index.get_seq(
                self.chrm, self.start, self.end, error_end=error_end))

        # handle case where no read overlaps whole region
        # let each read contibute its sequence and fill the rest
        # with dashes
        reg_base_data = ['-'] * (self.end - self.start)
        if self.reads is None or len(self.reads) == 0:
            return self.update(seq=''.join(reg_base_data))
        # get region sequence by moving through reads that
        # cover the region, but only extract seqeunce from the
        # (close to) minimal reads
        s_reg_reads = sorted(self.reads, key=lambda r: (r.start, r.end))
        # begin filling sequence with first (by start pos) read
        reg_base_data, curr_cov_pos = self._update_seq(
            s_reg_reads.pop(0), reg_base_data)
        # if there was only one read return now
        if len(s_reg_reads) == 0 or curr_cov_pos >= self.end - self.start:
            return self.update(seq=''.join(reg_base_data))

        # get next read (by start pos)
        curr_read = s_reg_reads.pop(0)
        for next_read in s_reg_reads:
            # once the next read start passes the region covered thus far
            # add the sequence for the saved curr_read to the reg sequence
            if next_read.start >= curr_cov_pos:
                # add read with curr longest end position to the region seq
                reg_base_data, curr_cov_pos = self._update_seq(
                    curr_read, reg_base_data)
                curr_read = next_read
                # if the whole interval is covered return the sequence
                if curr_cov_pos >= self.end - self.start:
                    return self.update(seq=''.join(reg_base_data))
                continue
            if next_read.end > curr_read.end:
                curr_read = next_read

        reg_base_data, _ = self._update_seq(curr_read, reg_base_data)

        return self.update(seq=''.join(reg_base_data))

    def get_base_levels(self, read_rows=False, num_reads=None):
        """Extract base levels from this interval.

        Args:
            read_rows (bool): return array with reads as rows (default: row columns)
            num_reads (int): maximal number of reads to output. randomly downsample if more reads are present (default: output all reads)

        Return:
            `np.array` containing read mean levels with rows corresponding to interval position and columns to individual reads (or reverse if `read_rows`)
        """
        def get_read_reg_events(r_data):
            r_means = get_single_slot_genome_centric(r_data, 'norm_mean')
            if r_means is None: return None
            if r_data.start > self.start and r_data.end < self.end:
                # handle reads that are contained in a region
                # create region with nan values
                r_reg_means = np.full(self.end - self.start, np.NAN)
                r_reg_means[r_data.start - self.start:
                            r_data.end - self.start] = r_means
            elif r_data.start > self.start:
                # handle reads that start in middle of region
                start_overlap = self.end - r_data.start
                # create region with nan values
                r_reg_means = np.full(self.end - self.start, np.NAN)
                r_reg_means[-start_overlap:] = r_means[:start_overlap]
            elif r_data.end < self.end:
                # handle reads that end inside region
                end_overlap = r_data.end - self.start
                # create region with nan values
                r_reg_means = np.full(self.end - self.start, np.NAN)
                r_reg_means[:end_overlap] = r_means[-end_overlap:]
            else:
                r_reg_means = r_means[
                    self.start - r_data.start:self.end - r_data.start]

            return r_reg_means


        if self.reads is None or len(self.reads) == 0:
            raise TomboError(
                'Must annotate region with reads ' +
                '(see `TomboInterval.add_reads`) to extract base levels.')
        if num_reads is not None:
            np.random.shuffle(self.reads)
        reg_events = []
        for r_data in self.reads:
            if self.strand is not None and r_data.strand != self.strand:
                continue
            r_means = get_read_reg_events(r_data)
            if r_means is None: continue
            reg_events.append(r_means)
            if num_reads is not None and len(reg_events) >= num_reads:
                break

        if read_rows:
            return np.row_stack(reg_events)
        return np.column_stack(reg_events)


def filter_empty_regions(plot_intervals):
    num_filt = sum(len(p_int.reads) == 0 for p_int in plot_intervals)
    plot_intervals = [p_int for p_int in plot_intervals if len(p_int.reads) > 0]
    if len(plot_intervals) == 0:
        error_message_and_exit('No reads in any selected regions.')
    if VERBOSE and num_filt > 0:
        warning_message('Some selected regions contain no reads.')

    return plot_intervals

def get_unique_intervals(plot_intervals, covered_poss=None, num_regions=None):
    # unique genomic regions filter
    uniq_p_intervals = []
    used_intervals = defaultdict(set)
    for reg_data in plot_intervals:
        # could have significant region immediately next to
        # beginning/end of reads
        interval_poss = list(range(reg_data.start, reg_data.end))
        if reg_data.start not in used_intervals[(
                reg_data.chrm, reg_data.strand)] and (
                    covered_poss is None or all(
                        pos in covered_poss[(reg_data.chrm, reg_data.strand)]
                        for pos in interval_poss)):
            uniq_p_intervals.append(reg_data)
            used_intervals[(reg_data.chrm, reg_data.strand)].update(
                interval_poss)
        if num_regions is not None and len(uniq_p_intervals) >= num_regions:
            break

    return uniq_p_intervals


###########################################
###### Special Data Access Functions ######
###########################################

def get_channel_info(fast5_data):
    """Get channel information for a read
    """
    try:
        fast5_attrs = fast5_data['UniqueGlobalKey/channel_id'].attrs
    except KeyError:
        raise TomboError("No channel_id group in HDF5 file. " +
                         "Probably mux scan HDF5 file.")

    try:
        channel_info = channelInfo(
            fast5_attrs.get('offset'), fast5_attrs.get('range'),
            fast5_attrs.get('digitisation'), fast5_attrs.get('channel_number'),
            fast5_attrs.get('sampling_rate').astype(np.int64))
    except KeyError:
        raise TomboError("Channel info parameters not available.")

    return channel_info

def get_raw_signal(r_data, int_start, int_end):
    """Extract raw signal from where this read overlaps a particular genomic region
    """
    with h5py.File(r_data.fn, 'r') as fast5_data:
        # retrieve shift and scale computed in correction script
        corr_subgrp = fast5_data['/Analyses/' + r_data.corr_group]
        event_starts = corr_subgrp['Events']['start']
        events_end = event_starts[-1] + corr_subgrp['Events']['length'][-1]
        segs = np.concatenate([event_starts, [events_end,]])

        scale_values = scaleValues(
            corr_subgrp.attrs.get('shift'), corr_subgrp.attrs.get('scale'),
            corr_subgrp.attrs.get('lower_lim'),
            corr_subgrp.attrs.get('upper_lim'),
            corr_subgrp.attrs.get('outlier_threshold'))
        all_sig = get_raw_read_slot(fast5_data)['Signal'][:]

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
    """Parse correction data from an event resquiggled read
    """
    try:
        with h5py.File(r_data.fn, 'r') as fast5_data:
            corr_grp = fast5_data['/Analyses/' + r_data.corr_group]
            events_grp = corr_grp['Events']
            event_starts = events_grp['start']
            events_end = event_starts[-1] + events_grp['length'][-1]
            new_segs = np.concatenate([event_starts, [events_end,]])

            raw_grp = get_raw_read_slot(fast5_data)
            read_id = raw_grp.attrs.get('read_id')
            try:
                read_id = read_id.decode()
            except (AttributeError, TypeError):
                pass
            signal_data = raw_grp['Signal'][:]

            raw_offset = events_grp.attrs.get('read_start_rel_to_raw')
            scale_values = scaleValues(*[
                corr_grp.attrs.get(attr_name) for attr_name in (
                    'shift', 'scale', 'lower_lim', 'upper_lim',
                    'outlier_threshold')])

            old_segs = corr_grp['Alignment/read_segments'][:]
            old_align_vals = list(map(
                lambda x: x.decode(),
                corr_grp['Alignment/read_alignment'][:]))
            new_align_vals = list(map(
                lambda x: x.decode(),
                corr_grp['Alignment/genome_alignment'][:]))
    except:
        return None

    if r_data.rna:
        signal_data = signal_data[::-1]

    return (read_id, signal_data, raw_offset, scale_values, old_segs,
            old_align_vals, new_align_vals, events_end, new_segs)

def get_all_read_data(r_data):
    """Extract most relevant read data from this read
    """
    try:
        with h5py.File(r_data.fn, 'r') as fast5_data:
            # note that it's more efficient to try to access the slot
            # and except the error that check if the slot exists first
            corr_subgrp = fast5_data['/Analyses/' + r_data.corr_group]
            algn_subgrp = dict(corr_subgrp['Alignment'].attrs.items())
            event_data = corr_subgrp['Events'][:]
            r_attrs = dict(corr_subgrp.attrs.items())
            all_sig = get_raw_read_slot(fast5_data)['Signal'][:]
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
            r_attrs.get('norm_type'), r_attrs.get('outlier_threshold'),
            genomeLocation(
                algn_subgrp['mapped_start'], algn_subgrp['mapped_strand'],
                algn_subgrp['mapped_chrom']))

def get_reads_events(cs_reads):
    """Extract read base levels split by genomic position

    Args:

        cs_reads (list): a list of reads from a single (chromosome, strand)

    Returns:
        A dictionary with 0-based genomic positions pointing to read signal levels
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
        if read_means.shape[0] != r_data.end - r_data.start: continue
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


###################################
###### FAST5 Write Functions ######
###################################

def prep_fast5(fast5_fn, corr_grp, overwrite, in_place,
               bc_grp=None, return_fp=False):
    """Prepare a read for re-squiggle processing (This deletes old re-squiggle info for this read)
    """
    def try_close_prep_err(fast5_data, err_str):
        try:
            fast5_data.close()
        except:
            pass
        # is_tombo_error = True
        return err_str, fast5_fn, True

    # several checks to prepare the FAST5 file for correction before
    # processing to save compute
    if not in_place:
        return ('Not currently implementing new hdf5 file writing',
                fast5_fn, True)
    # check that the file is writeable before trying to correct
    if not os.access(fast5_fn, os.W_OK):
        return 'FAST5 file is not writable', fast5_fn, True

    try:
        # create group to store data
        fast5_data = h5py.File(fast5_fn, 'r+')
        try:
            analyses_grp = fast5_data['/Analyses']
        except:
            return try_close_prep_err(
                fast5_data,
                'Base calls not found in FAST5 (see `tombo preprocess`)')
        try:
            # check that the requested basecalls group exsists
            if bc_grp is not None:
                analyses_grp[bc_grp]
        except:
            return try_close_prep_err(
                fast5_data,
                'Base calls not found in FAST5 (see `tombo preprocess`)')

        try:
            corr_grp_ptr = analyses_grp[corr_grp]
            if not overwrite:
                return try_close_prep_err(
                    fast5_data, "Tombo data exists in [--corrected-group] " +
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
        return 'Error opening or writing to fast5 file', fast5_fn, True

    if return_fp:
        return fast5_data

    try:
        fast5_data.close()
    except:
        return 'Error closing fast5 file', fast5_fn, True

    return

def write_error_status(fn, corr_grp, bc_subgrp, error_text):
    """Write error message for a read into the FAST5 file
    """
    with h5py.File(fn, 'r+') as fast5_data:
        analysis_grp = fast5_data['/Analyses']
        corr_grp = analysis_grp[corr_grp]
        if bc_subgrp is not None:
            # add subgroup matching subgroup from original basecalls
            corr_subgrp = corr_grp.create_group(bc_subgrp)
            corr_subgrp.attrs['status'] = error_text
        else:
            corr_grp.attrs['status'] = error_text

    return

def write_new_fast5_group(
        fast5_data, corr_grp_slot, rsqgl_res, norm_type, compute_sd,
        alignVals=None, old_segs=None, rna=False):
    """Write new fast5 group with re-squiggle data
    """
    try:
        # compute event data before accessing fast5 file
        if compute_sd:
            norm_means, norm_stds = c_new_mean_stds(
                rsqgl_res.raw_signal, rsqgl_res.segs)
        else:
            norm_means = c_new_means(rsqgl_res.raw_signal, rsqgl_res.segs)
            norm_stds = repeat(np.NAN)

        # had to shift to names formats numpy array specification due to
        # python2 numpy unicode issues. See discussion here:
        # https://github.com/numpy/numpy/issues/2407
        event_data = np.array(
            list(zip(norm_means, norm_stds,
                     rsqgl_res.segs[:-1], np.diff(rsqgl_res.segs),
                     list(rsqgl_res.genome_seq))),
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
        raise TomboError('Error computing new events')

    do_close = False
    if not isinstance(fast5_data, h5py.File):
        try:
            fast5_data = h5py.File(fast5_data, 'r+')
            do_close = True
        except:
            raise TomboError(
                'Error opening file for new group writing. This should ' +
                'have been caught during the alignment phase. Check that ' +
                'there are no other tombo processes or processes ' +
                'accessing these HDF5 files running simultaneously.')

    try:
        analysis_grp = fast5_data['/Analyses']
        corr_grp = analysis_grp[corr_grp_slot]
        # add subgroup matching subgroup from original basecalls
        corr_subgrp = corr_grp.create_group(rsqgl_res.align_info.Subgroup)
        corr_subgrp.attrs['status'] = 'success'
        # TODO change to rev_sig
        corr_subgrp.attrs['rna'] = rna
        if rsqgl_res.sig_match_score is not None:
            corr_subgrp.attrs[
                'signal_match_score'] = rsqgl_res.sig_match_score
        corr_subgrp.attrs['shift'] = rsqgl_res.scale_values.shift
        corr_subgrp.attrs['scale'] = rsqgl_res.scale_values.scale
        corr_subgrp.attrs['norm_type'] = norm_type
        if rsqgl_res.scale_values.lower_lim is not None:
            corr_subgrp.attrs[
                'lower_lim'] = rsqgl_res.scale_values.lower_lim
        if rsqgl_res.scale_values.upper_lim is not None:
            corr_subgrp.attrs[
                'upper_lim'] = rsqgl_res.scale_values.upper_lim
        if rsqgl_res.scale_values.outlier_thresh is not None:
            corr_subgrp.attrs[
                'outlier_threshold'] = rsqgl_res.scale_values.outlier_thresh

        # store alignment statistics
        corr_alignment = corr_subgrp.create_group('Alignment')
        corr_alignment.attrs['mapped_start'] = rsqgl_res.genome_loc.Start
        corr_alignment.attrs[
            'mapped_end'] = rsqgl_res.genome_loc.Start + len(rsqgl_res.segs) - 1
        corr_alignment.attrs[
            'mapped_strand'] = rsqgl_res.genome_loc.Strand
        corr_alignment.attrs['mapped_chrom'] = rsqgl_res.genome_loc.Chrom

        if rsqgl_res.align_info is not None:
            corr_alignment.attrs[
                'clipped_bases_start'] = rsqgl_res.align_info.ClipStart
            corr_alignment.attrs[
                'clipped_bases_end'] = rsqgl_res.align_info.ClipEnd
            corr_alignment.attrs[
                'num_insertions'] = rsqgl_res.align_info.Insertions
            corr_alignment.attrs[
                'num_deletions'] = rsqgl_res.align_info.Deletions
            corr_alignment.attrs[
                'num_matches'] = rsqgl_res.align_info.Matches
            corr_alignment.attrs[
                'num_mismatches'] = rsqgl_res.align_info.Mismatches

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
            'read_start_rel_to_raw'] = rsqgl_res.read_start_rel_to_raw
    except:
        raise TomboError(
            'Error writing resquiggle information back into fast5 file.')

    if do_close:
        try:
            fast5_data.close()
        except:
            raise TomboError(
                'Error closing fast5 file after writing resquiggle information.')

    return


if __name__ == '__main__':
    sys.stderr.write('This is a module. See commands with `tombo -h`')
    sys.exit(1)
