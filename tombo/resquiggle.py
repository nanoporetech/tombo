import os, sys

import re
import h5py
import Queue
import pkg_resources

import numpy as np
np.seterr(all='raise')
import multiprocessing as mp

from time import sleep
from subprocess import call
from tempfile import NamedTemporaryFile
from collections import defaultdict, namedtuple

# import tombo modules/functions
import tombo_stats as ts
import tombo_helper as th

from c_helper import c_new_means, c_valid_cpts_w_cap_t_test
from model_resquiggle import forward_pass, traceback
from dynamic_programming import c_reg_z_scores, c_banded_forward_pass, \
    c_banded_traceback, c_base_z_scores, c_adaptive_banded_forward_pass

VERBOSE = False

_PROFILE_RSQGL = False
_PROFILE_ALIGN = False

_DEBUG_FIT = False
_DEBUG_FULL = False
_DEBUG_MIDDLE = False
_NUM_DEBUG_ENDS = 250

# allow this many times the alignment batch size into the queue of
# reads to be resquiggled
ALIGN_BATCH_MULTIPLIER = 5
PROGRESS_INTERVAL = 500

# table containing the number of observations per event and the running
# difference width for different conditions this will need to be updated
# particularly for RNA motor updates most likely
OBS_PER_EVENT_TABLE = {'rna':(60, 20), 'dna':(7, 3)}

readInfo = namedtuple(
    'readInfo',
    ('ID', 'Subgroup', 'ClipStart', 'ClipEnd',
     'Insertions', 'Deletions', 'Matches', 'Mismatches'))
mapperData = namedtuple('mapperData', ('exe', 'type', 'index'))
# set default index to None
mapperData.__new__.__defaults__ = (None,)

M5_FIELDS = (
    'qName', 'qLength', 'qStart', 'qEnd', 'qStrand',
    'tName', 'tLength', 'tStart', 'tEnd', 'tStrand',
    'score', 'numMatch', 'numMismatch', 'numIns', 'numDel',
    'mapQV', 'qAlignedSeq', 'matchPattern', 'tAlignedSeq')
SAM_FIELDS = (
    'qName', 'flag', 'rName', 'pos', 'mapq',
    'cigar', 'rNext', 'pNext', 'tLen', 'seq', 'qual')
CIGAR_PAT = re.compile('(\d+)([MIDNSHP=X])')


############################################
########## Debug output functions ##########
############################################

def _write_middle_debug(z_scores, fwd_pass, band_event_starts,
                       debug_fp, reg_id, debug_num_seq=_NUM_DEBUG_ENDS,
                       short=False):
    fwd_pass = fwd_pass[1:]
    if fwd_pass.shape[0] < debug_num_seq:
        debug_num_seq = fwd_pass.shape[0]
    debug_end_start = len(band_event_starts) - debug_num_seq
    debug_fp.write('\n'.join(
        '\t'.join(map(str, (band_pos + band_event_starts[seq_pos], seq_pos,
                            score, str(reg_id) + 'z_begin')))
        for seq_pos, s_data in enumerate(z_scores[:debug_num_seq])
        for band_pos, score in enumerate(s_data)) + '\n')
    debug_fp.write('\n'.join(
        '\t'.join(map(str, (band_pos + band_event_starts[seq_pos], seq_pos,
                            score, str(reg_id) + 'fwd_begin')))
        for seq_pos, s_data in enumerate(fwd_pass[:debug_num_seq])
        for band_pos, score in enumerate(s_data)) + '\n')
    if short: return

    debug_fp.write('\n'.join(
        '\t'.join(map(str, (
            band_pos + band_event_starts[debug_end_start + seq_pos], seq_pos,
            score, str(reg_id) + 'z_end')))
        for seq_pos, s_data in enumerate(z_scores[-debug_num_seq:])
        for band_pos, score in enumerate(s_data)) + '\n')
    debug_fp.write('\n'.join(
        '\t'.join(map(str, (
            band_pos + band_event_starts[debug_end_start + seq_pos], seq_pos,
            score, str(reg_id) + 'fwd_end')))
        for seq_pos, s_data in enumerate(fwd_pass[-debug_num_seq:])
        for band_pos, score in enumerate(s_data)) + '\n')

    return

def _write_full_debug(fwd_pass_move, band_event_starts, top_max_pos,
                     z_scores, debug_fp, failed_fp, reg_id, final_score):
    read_tb = c_banded_traceback(fwd_pass_move, band_event_starts, top_max_pos)
    prev_event_pos = read_tb[0]
    band_poss = []
    event_scores = []
    for seq_pos, event_pos in enumerate(read_tb[1:]):
        seq_e_poss = range(prev_event_pos, event_pos)
        seq_band_poss = [e_pos - band_event_starts[seq_pos]
                         for e_pos in seq_e_poss]
        band_poss.extend(seq_band_poss)
        event_scores.extend([z_scores[seq_pos][b_pos] for b_pos in seq_band_poss])
        prev_event_pos = event_pos

    debug_fp.write('\n'.join(
        '\t'.join(map(str, (e_pos, b_pos, e_score, str(reg_id))))
        for e_pos, (b_pos, e_score) in enumerate(zip(
                band_poss, event_scores))) + '\n')
    fail_str = (('Failed ' if final_score < 0 else 'Pass ') + str(final_score) +
                ' ' + str(float(final_score) / len(read_tb)))
    failed_fp.write(fail_str + '\t' + str(reg_id) + '\n')

    return

def _write_tb_debug(fwd_pass_move, band_event_starts, top_max_pos,
                   debug_fp, reg_id, debug_num_seq=_NUM_DEBUG_ENDS):
    read_tb = c_banded_traceback(fwd_pass_move, band_event_starts, top_max_pos)
    debug_fp.write('\n'.join(
        '\t'.join(map(str, (event_pos, seq_pos,
                            str(reg_id) + 'fwd_begin')))
        for seq_pos, event_pos in enumerate(read_tb[:debug_num_seq])) + '\n')
    debug_fp.write('\n'.join(
        '\t'.join(map(str, (event_pos, seq_pos,
                            str(reg_id) + 'fwd_end')))
        for seq_pos, event_pos in enumerate(read_tb[-debug_num_seq:])) + '\n')
    return


#################################################
########## Raw Signal Re-squiggle Code ##########
#################################################

def get_model_fit_segs(
        segs, norm_signal, r_ref_means, r_ref_sds, running_diff_width,
        del_fix_window, max_new_cpts, extra_sig_factor=1.1):
    """
    Find new segments at skipped bases during dynamic programming segmentation.

    :param all_del_ranges: List of start stop tuple ranges of deletion locations
    :param segs: current Read segment locations
    :param norm_signal: Normalized read siganl
    :param r_ref_means: Read refererence means from genomic sequence
    :param r_ref_sds: Read refererence standard deviations from genomic sequence
    :param running_diff_width: Width of moving neighboring windows over which to compute segment locations
    :param del_fix_window: amount to extend skipped base windows
    :param max_new_cpts: Maximum new changepoints to find
    :param extra_sig_factor: Amount of extra signal to require in order to perform signal space re-squiggle

    :returns: New segments with skipped bases resolved
    """
    def get_deletion_ranges():
        all_del_ranges = []
        for del_pos in np.where(np.diff(segs) == 0)[0]:
            if (all_del_ranges and
                del_pos <= all_del_ranges[-1][1] + del_fix_window + 2):
                all_del_ranges[-1] = (all_del_ranges[-1][0],
                                      del_pos + del_fix_window + 1)
            else:
                all_del_ranges.append((del_pos - del_fix_window,
                                       del_pos + del_fix_window + 1))
        if len(all_del_ranges) == 0:
            return

        # potentially trim first and last windows
        if all_del_ranges[0][0] < 0:
            all_del_ranges[0] = (0, all_del_ranges[0][1])
        if all_del_ranges[-1][1] > len(segs) - 1:
            all_del_ranges[-1] = (all_del_ranges[-1][0], len(segs) - 1)

        if max_new_cpts is not None and max([
                end - start for start, end in  all_del_ranges]) > max_new_cpts:
            raise NotImplementedError, (
                'Too many changepoints required for re-squiggle algorithm')

        return all_del_ranges

    all_del_ranges = get_deletion_ranges()
    if all_del_ranges is None:
        return segs

    for start, end in all_del_ranges:
        n_events = end - start
        sig_start, sig_end = segs[start], segs[end]
        sig_len = sig_end - sig_start
        if sig_len <= (n_events * running_diff_width) * extra_sig_factor:
            raise NotImplementedError, (
                'Too little signal around event-aligned genomic deletion')

        # since there are no read starts to start from
        pseudo_starts = np.linspace(0, sig_len, n_events + 1, dtype=np.int32)
        # find signal space z-scores mapping without real banding by allowing
        # entire window to be searched
        reg_z_scores = c_reg_z_scores(
            norm_signal[sig_start:sig_end], r_ref_means[start:end],
            r_ref_sds[start:end], pseudo_starts,
            0, n_events, n_events, running_diff_width)
        reg_fwd_scores = forward_pass(reg_z_scores, running_diff_width)
        # perform signal based scoring segmentation
        #  - it is ~60X faster than base space
        reg_segs = traceback(reg_fwd_scores, running_diff_width) + sig_start
        segs[start+1:end] = reg_segs

    if np.diff(segs).min() < 1:
        raise NotImplementedError, (
            'New segments include zero length events')
    if segs[0] < 0:
        raise NotImplementedError, (
            'New segments start with negative index')
    if segs[-1] > norm_signal.shape[0]:
        raise NotImplementedError, (
            'New segments end past raw signal values')

    return segs



##########################################################
########## Standard banding dynamic programming ##########
##########################################################

def _get_masked_event_mapping(
        event_means, r_ref_means, r_ref_sds,
        mapped_start_offset, mapped_end_offset, skip_pen, stay_pen, z_shift,
        bandwidth, band_boundary_thresh=5, score_thresh=0.0,
        mask_fill_z_score=-10, mask_bases=50, end_event_gap=0,
        reg_id=None, debug_fps=None):
    """
    Perform banded dynamic programming sequence to event alignment forcing
    the path to start and end at the previously discovered locations.
    This is performed by masking the z-scores outside a "cone" extended
    mask_bases from the beginning and end of the middle of the read.
    """
    half_bandwidth = bandwidth / 2
    seq_len = r_ref_means.shape[0]
    events_len = event_means.shape[0]

    # check if the mapped start and end positions are too close to the end of
    # the events array and extend the bandwidth window if so
    band_events_start_pos = (
        0 if half_bandwidth <= mapped_start_offset else
        mapped_start_offset - half_bandwidth)
    band_events_end_pos = (
        events_len - bandwidth if half_bandwidth <= mapped_end_offset
        else events_len - half_bandwidth - mapped_end_offset)
    band_event_starts = np.linspace(
        band_events_start_pos, band_events_end_pos, seq_len).astype(np.int32)

    # figure out how many bases on each end contain masking to only compute
    # masked z-scores here
    start_mask_seq_len = max(
        mask_bases, next(i for i, bes in enumerate(band_event_starts)
                         if bes >= mapped_start_offset))
    end_mask_seq_len = max(
        mask_bases, next(
            i for i, bes in enumerate(band_event_starts[::-1])
            if bes + bandwidth <= events_len - mapped_end_offset))
    assert start_mask_seq_len + end_mask_seq_len < seq_len, (
        'Invalid masking encountered in dynamic sequence to events mapping')

    # get masked z-scores at the beginning of the read
    mask_start_pos = np.linspace(
        mapped_start_offset + 1 + end_event_gap,
        band_event_starts[mask_bases - 1] + bandwidth,
        mask_bases).astype(np.int32)
    def get_start_mask_z_score(seq_pos, event_pos):
        start_mask_len = max(mapped_start_offset - event_pos, 0)
        end_mask_len = (0 if seq_pos >= mask_bases else
                        bandwidth - (mask_start_pos[seq_pos] - event_pos))
        event_vals = event_means[event_pos + start_mask_len:
                                 event_pos + bandwidth - end_mask_len]
        b_z_scores = c_base_z_scores(
            event_vals, r_ref_means[seq_pos], r_ref_sds[seq_pos])
        masked_z_scores = np.concatenate([
            [mask_fill_z_score] * start_mask_len, b_z_scores,
            [mask_fill_z_score] * end_mask_len])
        return masked_z_scores
    start_scores = [get_start_mask_z_score(seq_pos, event_pos)
                    for seq_pos, event_pos in enumerate(
                            band_event_starts[:start_mask_seq_len])]

    # now the same for the end masked positions
    mask_end_pos = np.linspace(
        events_len - mapped_end_offset - end_event_gap - 1,
        band_event_starts[-mask_bases],
        mask_bases).astype(np.int32)
    def get_end_mask_z_score(seq_pos, event_pos):
        start_mask_len = (
            0 if seq_len - seq_pos - 1 >= mask_bases else
            mask_end_pos[seq_len - seq_pos - 1] - event_pos)
        end_mask_len = max(
            event_pos + bandwidth - events_len + mapped_end_offset - 1, 0)
        event_vals = event_means[event_pos + start_mask_len:
                                 event_pos + bandwidth - end_mask_len]
        b_z_scores = c_base_z_scores(
            event_vals, r_ref_means[seq_pos], r_ref_sds[seq_pos])
        masked_z_scores = np.concatenate([
            [mask_fill_z_score] * start_mask_len, b_z_scores,
            [mask_fill_z_score] * end_mask_len])
        return masked_z_scores
    end_scores = [get_end_mask_z_score(
        seq_pos + seq_len - end_mask_seq_len, event_pos)
                  for seq_pos, event_pos in enumerate(
                          band_event_starts[-end_mask_seq_len:])]

    # compute middle z_scores, combine and shift
    unmasked_z_scores = [
        c_base_z_scores(event_means[event_pos:event_pos + bandwidth],
                        r_ref_means[seq_pos + start_mask_seq_len],
                        r_ref_sds[seq_pos + start_mask_seq_len])
        for seq_pos, event_pos in enumerate(band_event_starts[
                start_mask_seq_len:-end_mask_seq_len])]
    shifted_z_scores = z_shift + np.row_stack(
        start_scores + unmasked_z_scores + end_scores)
    fwd_pass, fwd_pass_move = c_banded_forward_pass(
        shifted_z_scores, band_event_starts, skip_pen, stay_pen)

    # perform traceback
    top_max_pos = np.argmax(fwd_pass[-1,:])

    if _DEBUG_FULL:
        _write_full_debug(fwd_pass_move, band_event_starts, top_max_pos,
                         shifted_z_scores, debug_fps[2], debug_fps[3], reg_id,
                         fwd_pass[-1,top_max_pos])
    if _DEBUG_MIDDLE:
        _write_middle_debug(shifted_z_scores, fwd_pass, band_event_starts,
                           debug_fps[0], reg_id)
        _write_tb_debug(fwd_pass_move, band_event_starts, top_max_pos,
                       debug_fps[1], reg_id)

    if fwd_pass[-1,top_max_pos] < score_thresh:
        raise NotImplementedError, (
            'No valid path found through raw signal of long read')

    read_tb = c_banded_traceback(fwd_pass_move, band_event_starts, top_max_pos,
                                 band_boundary_thresh)
    start_trim_i = 0
    while read_tb[start_trim_i] < 0:
        read_tb[start_trim_i] = 0
        start_trim_i += 1
    end_trim_i = 1
    while read_tb[-end_trim_i] > events_len:
        read_tb[-end_trim_i] = events_len
        end_trim_i += 1

    # TODO: barrow code from debug_full to add strict filter for reads
    # with significant portion (500 events) not matching model well

    return read_tb

def _get_mapping_ends(
        event_means, r_ref_means, r_ref_sds, skip_pen, stay_pen, z_shift,
        seq_window, event_window, score_thresh, reg_id=None, debug_fps=None):
    if event_means.shape[0] < event_window:
        raise NotImplementedError, (
            'Read too short for eventless start/end discovery')
    if r_ref_means.shape[0] < seq_window:
        raise NotImplementedError, (
            'Genomic mapping too short for eventless start/end discovery')

    # banded z-scores (moving up one event per base for start/end discovery
    start_z_scores = z_shift - np.row_stack([
        np.abs(event_means[seq_pos:seq_pos + event_window] -
               r_ref_means[seq_pos]) / r_ref_sds[seq_pos]
        for seq_pos in range(seq_window)])
    start_band_event_starts = np.arange(seq_window, dtype=np.int32)

    start_fwd_pass, start_fwd_pass_move = c_banded_forward_pass(
        start_z_scores, start_band_event_starts, skip_pen, stay_pen)

    #print '\n'.join(
    #    '\t'.join(map(str, (pos, np.max(fwd_pass[pos-1,:]), 'start', reg_id)))
    #    for pos in [100,200,300,400])

    # find max along the top and right edges to start traceback
    top_max_pos = np.argmax(start_fwd_pass[-1,:])
    if start_fwd_pass[-1,top_max_pos] < score_thresh:
        raise NotImplementedError, (
            'No valid path found through start of raw signal')

    # perform traceback
    start_tb = c_banded_traceback(
        start_fwd_pass_move, start_band_event_starts, top_max_pos)


    # Now identify the end of the read
    n_events = event_means.shape[0]
    n_bases = r_ref_means.shape[0]
    end_band_event_starts = np.arange(n_events - event_window - seq_window,
                                      n_events - event_window, dtype=np.int32)
    end_z_scores = z_shift - np.row_stack([
        np.abs(event_means[end_band_event_starts[seq_pos]:
                           end_band_event_starts[seq_pos] + event_window] -
               r_ref_means[n_bases - seq_window + seq_pos]) / r_ref_sds[
                   n_bases - seq_window + seq_pos]
        for seq_pos in range(seq_window)])

    end_fwd_pass, end_fwd_pass_move = c_banded_forward_pass(
        end_z_scores, end_band_event_starts, skip_pen, stay_pen)
    #print '\n'.join(
    #    '\t'.join(map(str, (pos, np.max(fwd_pass[pos-1,:]), 'end', reg_id)))
    #    for pos in [100,200,300,400])

    # find max along the top and right edges to start traceback
    top_max_pos = np.argmax(end_fwd_pass[-1,:])
    if end_fwd_pass[-1,top_max_pos] < score_thresh:
        raise NotImplementedError, 'No valid path found through end of raw signal'

    # perform traceback
    end_tb = c_banded_traceback(
        end_fwd_pass_move, end_band_event_starts, top_max_pos)

    return start_tb[0], end_tb[-1]

def get_short_read_event_mapping(
        event_means, r_ref_means, r_ref_sds, skip_pen, stay_pen, z_shift,
        score_thresh, reg_id=None, debug_fps=None):
    """
    Perform banded dynamic programming sequence to event alignment without masking

    :param event_means: Numpy array with read base means
    :param r_ref_means: Numpy array with read reference means
    :param r_ref_sds: Numpy array with read reference standard deviations
    :param skip_pen: Penalty applied to skipped genomic bases
    :param stay_pen: Penalty applied to stay states (should shift to 0 expected value)
    :param z_shift: Shift z-scores by this amount (includes matching positive expected value)
    :param score_thresh: Threshold for a read to pass

    :returns: Event to sequence mapping for full length of short read
    """
    seq_len = r_ref_means.shape[0]
    events_len = event_means.shape[0]

    # create read starts in order to clip just the corners of the full events to
    # seqeunce matrix
    mask_len = min(seq_len, events_len) / 4
    band_event_starts = np.concatenate([
        np.zeros(seq_len - mask_len * 2),
        np.linspace(0, mask_len, mask_len * 2)]).astype(np.int32)
    bandwidth = events_len - mask_len

    shifted_z_scores = z_shift - np.row_stack([
        np.abs(event_means[event_pos:event_pos + bandwidth]
               - r_ref_means[seq_pos]) / r_ref_sds[seq_pos]
        for seq_pos, event_pos in enumerate(band_event_starts)])
    fwd_pass, fwd_pass_move = c_banded_forward_pass(
        shifted_z_scores, band_event_starts, skip_pen, stay_pen)

    # perform traceback
    top_max_pos = np.argmax(fwd_pass[-1,:])
    if _DEBUG_FULL:
        _write_full_debug(fwd_pass_move, band_event_starts, top_max_pos,
                         shifted_z_scores, debug_fps[2], debug_fps[3], reg_id,
                         fwd_pass[-1,top_max_pos])
    if _DEBUG_MIDDLE:
        _write_middle_debug(shifted_z_scores, fwd_pass, band_event_starts,
                           debug_fps[0], reg_id, short=True)
        _write_tb_debug(fwd_pass_move, band_event_starts, top_max_pos,
                       debug_fps[1], reg_id)

    if fwd_pass[-1,top_max_pos] < score_thresh:
        raise NotImplementedError, (
            'No valid path found through raw signal of short read')
    read_tb = c_banded_traceback(fwd_pass_move, band_event_starts, top_max_pos)

    return read_tb

def _find_base_assignment(
        norm_signal, min_base_obs, num_events, tb_model,
        genome_seq, genome_loc, skip_pen, stay_pen, z_shift, bandwidth,
        ends_bandwidth=2000, ends_seq_window=300, ends_thresh=0.0,
        reg_id=None, debug_fps=None):
    # get events before clipping
    valid_cpts = th.get_valid_cpts(norm_signal, min_base_obs, num_events)
    #valid_cpts = c_valid_cpts_w_cap_t_test(
    #                 norm_signal, min_base_obs, num_events)
    valid_cpts.sort()
    event_means = c_new_means(norm_signal, valid_cpts)
    kmer_ref, kmer_width, upstrm_bases, dnstrm_bases = tb_model

    r_ref_means, r_ref_sds = map(np.array, zip(*[
        kmer_ref[kmer] for kmer in [''.join(bs) for bs in zip(*[
            genome_seq[i:] for i in range(kmer_width)])]]))
    # trim genome seq to match model-able positions
    genome_seq = genome_seq[upstrm_bases:-dnstrm_bases]
    if genome_loc.Strand == '+':
        genome_loc = th.genomeLoc(
            genome_loc.Start + upstrm_bases, '+', genome_loc.Chrom)
    else:
        genome_loc = th.genomeLoc(
            genome_loc.Start + dnstrm_bases, '-', genome_loc.Chrom)

    # for short reads, just search the whole read with a larger bandwidth
    if (event_means.shape[0] < (ends_bandwidth + ends_seq_window) * 2 or
        r_ref_means.shape[0] < ends_seq_window * 2):
        seq_events = get_short_read_event_mapping(
            event_means, r_ref_means, r_ref_sds,
            skip_pen, stay_pen, z_shift, sr_bandwidth, ends_thresh,
            reg_id=reg_id, debug_fps=debug_fps)
        seq_segs = valid_cpts[seq_events]
        read_start_rel_to_raw = seq_segs[0]
        seq_segs = seq_segs - read_start_rel_to_raw
        return (seq_segs, r_ref_means, r_ref_sds, read_start_rel_to_raw,
                genome_seq, genome_loc)

    # identify the start and end of the read within the signal using a larger
    # bandwidth
    mapped_start, mapped_end = _get_mapping_ends(
        event_means, r_ref_means, r_ref_sds, skip_pen, stay_pen, z_shift,
        ends_seq_window, ends_bandwidth, ends_thresh,
        reg_id=reg_id, debug_fps=debug_fps)

    # get number of events to clip and how far into the events the
    # discovered start is located
    hald_bandwidth = bandwidth / 2
    if mapped_start < hald_bandwidth:
        events_start_clip = 0
        mapped_start_offset = mapped_start
    else:
        events_start_clip = mapped_start - hald_bandwidth
        mapped_start_offset = hald_bandwidth

    if mapped_end + hald_bandwidth > event_means.shape[0]:
        events_end_clip = event_means.shape[0]
        mapped_end_offset = event_means.shape[0] - mapped_end
    else:
        events_end_clip = mapped_end + hald_bandwidth
        mapped_end_offset = hald_bandwidth

    # now find full sequence to events path using a smaller bandwidth
    event_means = event_means[events_start_clip:events_end_clip]
    valid_cpts = valid_cpts[events_start_clip:events_end_clip + 1]
    read_tb = _get_masked_event_mapping(
        event_means, r_ref_means, r_ref_sds,
        mapped_start_offset, mapped_end_offset,
        skip_pen, stay_pen, z_shift, bandwidth,
        reg_id=reg_id, debug_fps=debug_fps)
    seq_segs = valid_cpts[read_tb]
    read_start_rel_to_raw = seq_segs[0]
    seq_segs = seq_segs - read_start_rel_to_raw

    return (seq_segs, r_ref_means, r_ref_sds, read_start_rel_to_raw,
            genome_seq, genome_loc)



##########################################################
########## Adaptive banding dynamic programming ##########
##########################################################

def get_masked_start_fwd_pass(
        event_means, r_ref_means, r_ref_sds, mapped_start_offset,
        skip_pen, stay_pen, z_shift, bandwidth, events_per_base,
        mask_fill_z_score=-10, mask_bases=50, end_event_gap=0,
        reg_id=None, debug_fps=None):
    """
    Perform banded dynamic programming sequence to event alignment forcing
    the path to start and end at the previously discovered locations.
    This is performed by masking the z-scores outside a "cone" extended
    mask_bases from the beginning and end of the middle of the read.

    :param event_means: Numpy array with read base means
    :param r_ref_means: Numpy array with read reference means
    :param r_ref_sds: Numpy array with read reference standard deviations
    :param mapped_start_offset: Previously identified start of genomic sequence within events
    :param skip_pen: Penalty applied to skipped genomic bases
    :param stay_pen: Penalty applied to stay states (should shift to 0 expected value)
    :param z_shift: Shift z-scores by this amount (includes matching positive expected value)
    :param bandwidth: Bandwidth over which to search for sequence to event mapping
    :param events_per_base: Average events per base for the start mapping

    :returns: Event to sequence mapping for start of read including forward pass scores, forward pass move
        values, band starts within the events vector and z-scores
    """
    half_bandwidth = bandwidth / 2

    # check if the mapped start position is too close to the end of
    # the events array and extend the bandwidth window if so
    band_events_start_pos = (
        0 if half_bandwidth <= mapped_start_offset else
        mapped_start_offset - half_bandwidth)

    tmp_seq_len = max(half_bandwidth, mask_bases,
                      int((half_bandwidth + 1) / events_per_base)) + 1
    band_event_starts = np.linspace(
        band_events_start_pos,
        band_events_start_pos + (tmp_seq_len  * events_per_base),
        tmp_seq_len).astype(np.int32)
    mask_seq_len = max(
        mask_bases, next(i + 2 for i, bes in enumerate(band_event_starts)
                         if bes >= mapped_start_offset))
    band_event_starts = band_event_starts[:mask_seq_len]

    # get masked z-scores at the beginning of the read
    mask_start_pos = np.linspace(
        mapped_start_offset + 1 + end_event_gap,
        band_event_starts[mask_bases - 1] + bandwidth,
        mask_bases).astype(np.int32)
    def get_start_mask_z_score(seq_pos, event_pos):
        start_mask_len = max(mapped_start_offset - event_pos, 0)
        end_mask_len = (0 if seq_pos >= mask_bases else
                        bandwidth - (mask_start_pos[seq_pos] - event_pos))
        event_vals = event_means[event_pos + start_mask_len:
                                 event_pos + bandwidth - end_mask_len]
        b_z_scores = c_base_z_scores(
            event_vals, r_ref_means[seq_pos], r_ref_sds[seq_pos])
        masked_z_scores = np.concatenate([
            [mask_fill_z_score] * start_mask_len, b_z_scores,
            [mask_fill_z_score] * end_mask_len])
        return masked_z_scores
    shifted_z_scores = z_shift + np.row_stack([
        get_start_mask_z_score(seq_pos, event_pos)
        for seq_pos, event_pos in enumerate(band_event_starts)])
    fwd_pass, fwd_pass_move = c_banded_forward_pass(
        shifted_z_scores, band_event_starts, skip_pen, stay_pen)

    return fwd_pass, fwd_pass_move, band_event_starts, shifted_z_scores

def get_mapping_start(
        event_means, r_ref_means, r_ref_sds, skip_pen, stay_pen, z_shift,
        seq_window, bandwidth, score_thresh, reg_id=None, debug_fps=None):
    """
    Perform banded dynamic programming sequence to event alignment through
    The beginning of an read to identify the start of genome sequence to event matching

    :param event_means: Numpy array with read base means
    :param r_ref_means: Numpy array with read reference means
    :param r_ref_sds: Numpy array with read reference standard deviations
    :param skip_pen: Penalty applied to skipped genomic bases
    :param stay_pen: Penalty applied to stay states (should shift to 0 expected value)
    :param z_shift: Shift z-scores by this amount (includes matching positive expected value)
    :param seq_window: Number of genomic bases to search over for the start of the read
    :param bandwidth: Bandwidth over which to search for sequence to event mapping
    :param score_thresh: Threshold for a read to pass

    :returns: Start of seqeunce to event alignment and the mean events_per_base through the queried portion of a read
    """
    if event_means.shape[0] < bandwidth:
        raise NotImplementedError, (
            'Read too short for eventless start/end discovery')
    if r_ref_means.shape[0] < seq_window:
        raise NotImplementedError, (
            'Genomic mapping too short for eventless start/end discovery')

    # banded z-scores (moving up one event per base for start/end discovery
    start_z_scores = z_shift - np.row_stack([
        np.abs(event_means[seq_pos:seq_pos + bandwidth] -
               r_ref_means[seq_pos]) / r_ref_sds[seq_pos]
        for seq_pos in range(seq_window)])
    start_band_event_starts = np.linspace(
        0, seq_window, seq_window).astype(np.int32)

    np.arange(seq_window, dtype=np.int32)

    start_fwd_pass, start_fwd_pass_move = c_banded_forward_pass(
        start_z_scores, start_band_event_starts, skip_pen, stay_pen)

    # find max along the top and right edges to start traceback
    top_max_pos = np.argmax(start_fwd_pass[-1,:])
    if start_fwd_pass[-1,top_max_pos] < score_thresh:
        # TODO: Add iterative start search with set number of iterations
        # for reads that start further into the read.
        raise NotImplementedError, (
            'No valid path found through start of raw signal')

    # perform traceback
    start_tb = c_banded_traceback(
        start_fwd_pass_move, start_band_event_starts, top_max_pos)
    # compute the average events per base to use for the start forward pass
    events_per_base = float(start_tb[-1] - start_tb[0]) / len(start_tb)

    return start_tb[0], events_per_base

def find_adaptive_base_assignment(
        norm_signal, min_base_obs, num_events, tb_model,
        genome_seq, genome_loc, skip_pen, stay_pen, z_shift, bandwidth, rna,
        start_bandwidth=2000, start_seq_window=500, start_thresh=0.0,
        band_boundary_thresh=5, reg_id=None, debug_fps=None):
    """
    Perform banded dynamic programming sequence to event alignment by first identifying the start of the
    sequence to event matching and then performing banded matching through the whole read

    :param norm_signal: Numpy array with normalized read signal
    :param min_base_obs: Minimum number of raw observations per base
    :param num_events: Number of events to identify in this read
    :param tb_model: A Tombo model
    :param genome_seq: Genomic sequence for this read
    :param genome_loc: Mapped genomic location for this read
    :param skip_pen: Penalty applied to skipped genomic bases
    :param stay_pen: Penalty applied to stay states (should shift to 0 expected value)
    :param z_shift: Shift z-scores by this amount (includes matching positive expected value)
    :param bandwidth: Bandwidth over which to search for sequence to event mapping
    :param rna: Is this an RNA read

    :returns: Start of seqeunce to event alignment and the mean events_per_base through the queried portion of a read
    """
    # get events
    # for RNA evenly smaller events could be detrimental to the fit
    # so perform slower segmentation which does not allow small events
    if rna:
        valid_cpts = c_valid_cpts_w_cap_t_test(
            norm_signal, min_base_obs, num_events)
    else:
        valid_cpts = th.get_valid_cpts(norm_signal, min_base_obs, num_events)
    valid_cpts.sort()
    event_means = c_new_means(norm_signal, valid_cpts)
    kmer_ref, kmer_width, upstrm_bases, dnstrm_bases = tb_model

    r_ref_means, r_ref_sds = map(np.array, zip(*[
        kmer_ref[kmer] for kmer in [''.join(bs) for bs in zip(*[
            genome_seq[i:] for i in range(kmer_width)])]]))
    # trim genome seq to match model-able positions
    genome_seq = genome_seq[upstrm_bases:-dnstrm_bases]
    seq_len = len(genome_seq)
    if genome_loc.Strand == '+':
        genome_loc = th.genomeLoc(
            genome_loc.Start + upstrm_bases, '+', genome_loc.Chrom)
    else:
        genome_loc = th.genomeLoc(
            genome_loc.Start + dnstrm_bases, '-', genome_loc.Chrom)

    # for short reads, just search the whole read with a larger bandwidth
    if (event_means.shape[0] < start_bandwidth + start_seq_window or
        seq_len < start_seq_window):
        seq_events = get_short_read_event_mapping(
            event_means, r_ref_means, r_ref_sds, skip_pen, stay_pen,
            z_shift, start_thresh, reg_id=reg_id, debug_fps=debug_fps)
        seq_segs = valid_cpts[seq_events]
        read_start_rel_to_raw = seq_segs[0]
        seq_segs = seq_segs - read_start_rel_to_raw
        return (seq_segs, r_ref_means, r_ref_sds, read_start_rel_to_raw,
                genome_seq, genome_loc)

    # identify the start and end of the read within the signal using a larger
    # bandwidth
    mapped_start, events_per_base = get_mapping_start(
        event_means, r_ref_means, r_ref_sds, skip_pen, stay_pen, z_shift,
        start_seq_window, start_bandwidth, start_thresh,
        reg_id=reg_id, debug_fps=debug_fps)

    # get number of events to clip and how far into the events the
    # discovered start is located
    half_bandwidth = bandwidth / 2
    if mapped_start < half_bandwidth:
        events_start_clip = 0
        mapped_start_offset = mapped_start
    else:
        events_start_clip = mapped_start - half_bandwidth
        mapped_start_offset = half_bandwidth

    # now find full sequence to events path using a smaller bandwidth
    event_means = event_means[events_start_clip:]
    valid_cpts = valid_cpts[events_start_clip:]
    (start_fwd_pass, start_fwd_pass_move,
     start_event_starts, start_z_scores) = get_masked_start_fwd_pass(
         event_means, r_ref_means, r_ref_sds,
         mapped_start_offset, skip_pen, stay_pen, z_shift, bandwidth,
         events_per_base, reg_id=reg_id, debug_fps=debug_fps)
    start_seq_len = start_event_starts.shape[0]
    fwd_pass = np.row_stack([
        start_fwd_pass,
        np.empty((seq_len - start_seq_len, bandwidth))])
    fwd_pass_move = np.row_stack([
        start_fwd_pass_move,
        np.empty((seq_len - start_seq_len, bandwidth), dtype=np.int32)])
    band_event_starts = np.concatenate(
        [start_event_starts,
         np.empty(seq_len - start_seq_len, dtype=np.int32)])
    #fwd_pass[start_seq_len+1:,:] = np.NAN
    #fwd_pass_move[start_seq_len+1:,:] = np.NAN
    #band_event_starts[start_seq_len:] = np.NAN

    if _DEBUG_FULL or _DEBUG_MIDDLE:
        shifted_z_scores = c_adaptive_banded_forward_pass(
            fwd_pass, fwd_pass_move, band_event_starts, event_means,
            r_ref_means, r_ref_sds, z_shift, skip_pen, stay_pen,
            start_seq_len, True)
        shifted_z_scores = np.row_stack([start_z_scores, shifted_z_scores])
    else:
        c_adaptive_banded_forward_pass(
            fwd_pass, fwd_pass_move, band_event_starts, event_means,
            r_ref_means, r_ref_sds, z_shift, skip_pen, stay_pen,
            start_seq_len)

    top_max_pos = np.argmax(fwd_pass[-1,:])
    if _DEBUG_FULL:
        _write_full_debug(fwd_pass_move, band_event_starts, top_max_pos,
                         shifted_z_scores, debug_fps[2], debug_fps[3], reg_id,
                         fwd_pass[-1,top_max_pos])
    if _DEBUG_MIDDLE:
        _write_middle_debug(shifted_z_scores, fwd_pass, band_event_starts,
                           debug_fps[0], reg_id)
        _write_tb_debug(fwd_pass_move, band_event_starts, top_max_pos,
                       debug_fps[1], reg_id)

    if fwd_pass[-1,top_max_pos] < start_thresh:
        raise NotImplementedError, (
            'No valid path found through raw signal of long read')
    read_tb = c_banded_traceback(
        fwd_pass_move, band_event_starts, top_max_pos, band_boundary_thresh)

    start_trim_i = 0
    while read_tb[start_trim_i] < 0:
        read_tb[start_trim_i] = 0
        start_trim_i += 1
    end_trim_i = 1
    events_len = event_means.shape[0]
    while read_tb[-end_trim_i] > events_len:
        read_tb[-end_trim_i] = events_len
        end_trim_i += 1

    seq_segs = valid_cpts[read_tb]
    read_start_rel_to_raw = seq_segs[0]
    seq_segs = seq_segs - read_start_rel_to_raw

    return (seq_segs, r_ref_means, r_ref_sds, read_start_rel_to_raw,
            genome_seq, genome_loc)



###################################################
########## Resquiggle initial processing ##########
###################################################

def resquiggle_read(
        fast5_fn, genome_seq, tb_model, outlier_thresh, genome_loc, read_info,
        basecall_group, corrected_group, compute_sd, skip_pen, stay_pen, z_shift,
        bandwidth, obs_filter, del_fix_window=5, min_event_to_seq_ratio=1.1,
        max_new_cpts=None, in_place=True, skip_index=False, reg_id=None,
        debug_fps=None, const_scale=None):
    """
    Perform banded dynamic programming sequence to event alignment for this read

    :param fast5_fn: Filename for a read
    :param genome_seq: Genomic sequence for this read
    :param tb_model: A Tombo model
    :param outlier_thresh: Outlier threshold for raw signal normalization
    :param genome_loc: Mapped genomic location for this read
    :param read_info: A read info named tuple for this read
    :param basecall_group: The basecalled read group to analyze
    :param corrected_group: The tombo corrected group to write results
    :param compute_sd: Should SD computations be computed and saved to file
    :param skip_pen: Penalty applied to skipped genomic bases
    :param stay_pen: Penalty applied to stay states (should shift to 0 expected value)
    :param z_shift: Shift z-scores by this amount (includes matching positive expected value)
    :param bandwidth: Bandwidth over which to search for sequence to event mapping
    :param obs_filter: Obervations per base filter to apply for filtered slot in FAST5
    """
    # errors should not happen here since these slotes were checked
    # in alignment function, but old zombie processes might cause
    # problems here
    try:
        fast5_data = h5py.File(fast5_fn, 'r')
        channel_info = th.get_channel_info(fast5_data)

        # extract raw data for this read
        all_raw_signal = fast5_data['/Raw/Reads/'].values()[0]['Signal'].value

        rna = th.is_read_rna(fast5_data)

        fast5_data.close()
    except:
        #raise
        raise NotImplementedError, (
            'Error opening file for re-squiggle. This should have ' +
            'been caught during the alignment phase. Check that there ' +
            'are no other tombo processes or processes accessing ' +
            'these HDF5 files running simultaneously')

    if rna:
        all_raw_signal = all_raw_signal[::-1]

    mean_obs_per_base, running_diff_width = (
        OBS_PER_EVENT_TABLE['rna'] if rna else OBS_PER_EVENT_TABLE['dna'])
    num_events = max(all_raw_signal.shape[0] / mean_obs_per_base,
                     int(len(genome_seq) * min_event_to_seq_ratio))
    if num_events / bandwidth > len(genome_seq):
        raise NotImplementedError, 'Too much raw signal for short mapped sequence'
    # normalize signal
    if const_scale is not None:
        norm_signal, scale_values = th.normalize_raw_signal(
            all_raw_signal, 0, all_raw_signal.shape[0],
            'median_const_scale', channel_info, outlier_thresh,
            const_scale=const_scale)
    else:
        norm_signal, scale_values = th.normalize_raw_signal(
            all_raw_signal, 0, all_raw_signal.shape[0],
            'median', channel_info, outlier_thresh)

    (segs, r_ref_means, r_ref_sds, read_start_rel_to_raw,
     genome_seq, genome_loc) = find_adaptive_base_assignment(
         norm_signal, running_diff_width, num_events, tb_model,
         genome_seq, genome_loc, skip_pen, stay_pen, z_shift, bandwidth, rna,
         reg_id=reg_id, debug_fps=debug_fps)
    norm_signal = norm_signal[read_start_rel_to_raw:
                              read_start_rel_to_raw + segs[-1]]

    # identify all stretches of genomic deletions within del_fix_window
    # to be fixed.
    segs = get_model_fit_segs(
        segs, norm_signal, r_ref_means, r_ref_sds,
        running_diff_width, del_fix_window, max_new_cpts)

    if segs.shape[0] != len(genome_seq) + 1:
        raise ValueError, ('Aligned sequence does not match number ' +
                           'of segments produced')

    # Output for testing/visualization of event-less re-squiggle
    if _DEBUG_FIT:
        norm_means = c_new_means(norm_signal, segs)
        #r_mean_z = np.mean([np.abs((b_m - b_ref_m) / b_ref_s)
        #                    for b_m, b_ref_m, b_ref_s in
        #                    zip(norm_means, r_ref_means, r_ref_sds)])
        #if r_mean_z > 1:
        #    print fast5_fn
        with open('eventless_testing.model.txt', 'w') as fp:
            fp.write('Position\tMean\tSD\n' + '\n'.join(
                '\t'.join(map(str, (pos, p_mean, p_std)))
                for pos, (p_mean, p_std) in enumerate(zip(
                        r_ref_means, r_ref_sds))) + '\n')
        with open('eventless_testing.seq.txt', 'w') as fp:
            fp.write('Base\tPosition\tSignalMean\n' + '\n'.join(
                '\t'.join(map(str, (b, pos, p_mean))) for pos, (b, p_mean) in
                enumerate(zip(genome_seq, norm_means))) + '\n')
        Position, Signal = [], []
        for base_i, (b_start, b_end) in enumerate(zip(segs[:-1], segs[1:])):
            Position.extend(
                base_i + np.linspace(0, 1, b_end - b_start, endpoint=False))
            Signal.extend(norm_signal[b_start:b_end])
        with open('eventless_testing.signal.txt', 'w') as fp:
            fp.write('Position\tSignal\n' + '\n'.join(
                '\t'.join(map(str, (pos, sig)))
                for pos, sig in zip(Position, Signal)) + '\n')

    if in_place:
        # create new hdf5 file to hold new read signal
        th.write_new_fast5_group(
            fast5_fn, genome_loc, read_start_rel_to_raw, segs, genome_seq,
            norm_signal, scale_values, corrected_group, read_info.Subgroup,
            'median', outlier_thresh, compute_sd, None, read_info, None)
    else:
        # create new hdf5 file to hold corrected read events
        pass

    if not skip_index:
        return th.prep_index_data(
            fast5_fn, genome_loc, read_start_rel_to_raw, segs,
            corrected_group, read_info.Subgroup, rna, obs_filter)

    return

def _resquiggle_worker(
        basecalls_q, progress_q, failed_reads_q, index_q, basecall_group,
        corrected_group, tb_model_fn, outlier_thresh, compute_sd, skip_pen,
        match_evalue, bandwidth, obs_filter, const_scale):
    num_processed = 0
    skip_index = index_q is None
    if not skip_index: proc_index_data = []
    debug_fps = None
    if _DEBUG_MIDDLE or _DEBUG_FULL:
        score_fp = open('debug_event_align.txt', 'w')
        score_fp.write('EventPos\tSeqPos\tScore\tRegion\n')
        tb_fp = open('debug_event_align.traceback.txt', 'w')
        tb_fp.write('EventPos\tSeqPos\tRegion\n')
        full_fit_fp = open('debug_event_align.full_fit.txt', 'w')
        full_fit_fp.write('EventPos\tBandPos\tEventScore\tRegion\n')
        full_failed_fp = open('debug_event_align.full_failed.txt', 'w')
        full_failed_fp.write('DidFail\tRegion\n')
        debug_fps = [score_fp, tb_fp, full_fit_fp, full_failed_fp]

    # parse tombo model (ignore alt_base and model_name)
    kmer_ref, upstrm_bases, _, _ = ts.parse_tombo_model(tb_model_fn)
    kmer_width = len(next(kmer_ref.iterkeys()))
    dnstrm_bases = kmer_width - upstrm_bases - 1
    tb_model = (kmer_ref, kmer_width, upstrm_bases, dnstrm_bases)

    # get dynamic programming parameters
    z_shift, stay_pen = ts.get_dynamic_prog_params(match_evalue)

    while True:
        try:
            fast5_fn, sgs_align_data = basecalls_q.get(block=False)
            # None values placed in queue when all files have
            # been processed
            if fast5_fn is None: break
        except Queue.Empty:
            sleep(1)
            continue

        num_processed += 1
        if num_processed % int(PROGRESS_INTERVAL / 5) == 0:
            progress_q.put(int(PROGRESS_INTERVAL / 5))
        # process different read subgroups sequentially so that the same
        # file is never open simultaneously
        for genome_seq, genome_loc, read_info in sgs_align_data:
            try:
                index_data = resquiggle_read(
                    fast5_fn, genome_seq, tb_model, outlier_thresh,
                    genome_loc, read_info, basecall_group, corrected_group,
                    compute_sd, skip_pen, stay_pen, z_shift, bandwidth,
                    obs_filter, skip_index=skip_index, reg_id=num_processed,
                    debug_fps=debug_fps, const_scale=const_scale)
                if not skip_index:
                    proc_index_data.append(index_data)
                    if index_data[1][6]:
                        failed_reads_q.put((
                            'Read filtered by observation per base ' +
                            'thresholds (revert with `tombo clear_filters`)',
                            read_info.Subgroup + th.FASTA_NAME_JOINER + fast5_fn))
            except Exception as e:
                # uncomment to identify mysterious errors
                #raise
                try:
                    th.write_error_status(
                        fast5_fn, corrected_group, read_info.Subgroup, str(e))
                except:
                    pass
                failed_reads_q.put((
                    str(e), read_info.Subgroup + th.FASTA_NAME_JOINER + fast5_fn))

    if not skip_index: index_q.put(proc_index_data)

    return

if _PROFILE_RSQGL:
    _resquiggle_wrapper = _resquiggle_worker
    def _resquiggle_worker(*args):
        import cProfile
        cProfile.runctx('_resquiggle_wrapper(*args)', globals(), locals(),
                        filename='resquiggle_eventless_main.prof')
        return


############################################
########## Genomic Alignment Code ##########
############################################

def clip_m5_alignment(alignVals, start, strand, chrm):
    """
    Clip hard and soft clipped bases from an m5 format alignment
    """
    # clip read to first matching bases
    start_clipped_read_bases = 0
    start_clipped_genome_bases = 0
    start_clipped_align_bases = 0
    r_base, g_base = alignVals[0]
    while r_base == '-' or g_base == '-':
        start_clipped_read_bases += int(r_base != '-')
        start_clipped_genome_bases += int(g_base != '-')
        start_clipped_align_bases += 1
        r_base, g_base = alignVals[start_clipped_align_bases]

    end_clipped_read_bases = 0
    end_clipped_genome_bases = 0
    end_clipped_align_bases = 0
    r_base, g_base = alignVals[-1]
    while r_base == '-' or g_base == '-':
        end_clipped_read_bases += int(r_base != '-')
        end_clipped_genome_bases += int(g_base != '-')
        end_clipped_align_bases += 1
        r_base, g_base = alignVals[-1 * (end_clipped_align_bases + 1)]

    alignVals = alignVals[start_clipped_align_bases:]
    if end_clipped_align_bases > 0:
        alignVals = alignVals[:-1*end_clipped_align_bases]

    if strand == '+' and start_clipped_genome_bases > 0:
        genome_loc = th.genomeLoc(
            start + start_clipped_genome_bases, '+', chrm)
    elif strand == '-' and end_clipped_genome_bases > 0:
        genome_loc = th.genomeLoc(
            start + end_clipped_genome_bases, '-', chrm)
    else:
        genome_loc = th.genomeLoc(start, strand, chrm)

    return alignVals, start_clipped_read_bases, \
        end_clipped_read_bases, genome_loc

def parse_m5_record(r_m5_record, read_id, bc_subgroup):
    """
    Parse a single m5 formatted alignment
    """
    if r_m5_record['tStrand'] != '+':
        raise NotImplementedError, (
            'Mapping indicates negative strand reference mapping')

    if r_m5_record['qStrand'] == "+":
        alignVals = zip(r_m5_record['qAlignedSeq'],
                        r_m5_record['tAlignedSeq'])
    else:
        alignVals = zip(th.rev_comp(r_m5_record['qAlignedSeq']),
                        th.rev_comp(r_m5_record['tAlignedSeq']))

    alignVals, start_clipped_bases, end_clipped_bases, genome_loc \
        = clip_m5_alignment(
            alignVals, int(r_m5_record['tStart']),
            r_m5_record['qStrand'], r_m5_record['tName'])
    tSeq = ''.join(zip(*alignVals)[1])

    # TOOD compute indel/match/mismatch counts
    read_info = readInfo(
        read_id, bc_subgroup, start_clipped_bases, end_clipped_bases,
        0, 0, 0, 0)

    return tSeq, genome_loc, read_info

def parse_m5_output(align_output, batch_reads_data):
    """
    Parse a batch of m5 formatted alignments
    """
    alignments = dict(
        (read_fn_sg, None) for read_fn_sg in batch_reads_data.keys())
    for line in align_output:
        r_m5_record = dict(zip(M5_FIELDS, line.strip().split()))
        if len(r_m5_record) != len(M5_FIELDS):
            continue
        # store the alignment if none is stored for this read or
        # if this read has the highest map quality thus far
        qName = r_m5_record['qName'].replace(th.FN_SPACE_FILLER, ' ')
        if alignments[qName] is None or (
                int(alignments[qName]['score']) > int(r_m5_record['score'])
                and int(r_m5_record['score']) < 255):
            alignments[qName] = r_m5_record

    batch_align_failed_reads = []
    batch_align_data = []
    for read_fn_sg, r_m5_record in alignments.iteritems():
        bc_subgroup, read_fn = read_fn_sg.split(th.FASTA_NAME_JOINER)
        read_id = batch_reads_data[read_fn_sg][1]
        if r_m5_record is None:
            batch_align_failed_reads.append(
                ('Alignment not produced.', read_fn_sg))
        else:
            try:
                batch_align_data.append((read_fn, parse_m5_record(
                    r_m5_record, read_id, bc_subgroup)))
            except Exception as e:
                batch_align_failed_reads.append((str(e), read_fn_sg))

    return batch_align_failed_reads, batch_align_data

def parse_sam_record(
        r_sam_record, genome_index, read_id, bc_subgroup,
        skip_align_stats=False):
    """
    Parse a single of sam formatted alignment
    """
    def parse_cigar(strand):
        # parse cigar string
        cigar = [
            (int(reg_len), reg_type) for reg_len, reg_type in
            CIGAR_PAT.findall(r_sam_record['cigar'])]
        if len(cigar) < 1:
            raise RuntimeError, 'Invalid cigar string produced'

        if strand == '-':
            cigar = cigar[::-1]

        return cigar

    def get_just_tseq(cigar, strand):
        start_clipped_bases = 0
        end_clipped_bases = 0
        # handle clipping elements (H and S)
        if cigar[0][1] == 'H':
            start_clipped_bases += cigar[0][0]
            cigar = cigar[1:]
        if cigar[-1][1] == 'H':
            end_clipped_bases += cigar[-1][0]
            cigar = cigar[:-1]
        if cigar[0][1] == 'S':
            start_clipped_bases += cigar[0][0]
            cigar = cigar[1:]
        if cigar[-1][1] == 'S':
            end_clipped_bases += cigar[-1][0]
            cigar = cigar[:-1]

        tLen = sum([reg_len for reg_len, reg_type in cigar
                    if reg_type in 'MDN=X'])
        tSeq = genome_index[r_sam_record['rName']][
            int(r_sam_record['pos']) - 1:
            int(r_sam_record['pos']) + tLen - 1]
        if strand == '-': tSeq = th.rev_comp(tSeq)

        # check that cigar starts and ends with matched bases
        while cigar[0][1] not in 'M=X':
            if cigar[0][1] in 'IP':
                tSeq = tSeq[cigar[0][0]:]
            else:
                start_clipped_bases += cigar[0][0]
            cigar = cigar[1:]
        while cigar[-1][1] not in 'M=X':
            if cigar[-1][1] in 'IP':
                tSeq = tSeq[:-cigar[-1][0]]
            else:
                end_clipped_bases += cigar[0][0]
            cigar = cigar[:-1]

        return tSeq, start_clipped_bases, end_clipped_bases

    def get_qseq(cigar, strand):
        # record clipped bases and remove from query seq as well as cigar
        qSeq = r_sam_record['seq'] if strand == '+' else th.rev_comp(
            r_sam_record['seq'])
        start_clipped_bases = 0
        end_clipped_bases = 0
        # handle clipping elements (H and S)
        if cigar[0][1] == 'H':
            start_clipped_bases += cigar[0][0]
            cigar = cigar[1:]
        if cigar[-1][1] == 'H':
            end_clipped_bases += cigar[-1][0]
            cigar = cigar[:-1]
        if cigar[0][1] == 'S':
            start_clipped_bases += cigar[0][0]
            qSeq = qSeq[cigar[0][0]:]
            cigar = cigar[1:]
        if cigar[-1][1] == 'S':
            end_clipped_bases += cigar[-1][0]
            qSeq = qSeq[:-cigar[-1][0]]
            cigar = cigar[:-1]

        return qSeq, start_clipped_bases, end_clipped_bases, cigar

    def get_tseq(qSeq, start_clipped_bases, end_clipped_bases, cigar, strand):
        tLen = sum([reg_len for reg_len, reg_type in cigar
                    if reg_type in 'MDN=X'])
        tSeq = genome_index[r_sam_record['rName']][
            int(r_sam_record['pos']) - 1:
            int(r_sam_record['pos']) + tLen - 1]
        if strand == '-': tSeq = th.rev_comp(tSeq)

        # check that cigar starts and ends with matched bases
        while cigar[0][1] not in 'M=X':
            if cigar[0][1] in 'IP':
                tSeq = tSeq[cigar[0][0]:]
            else:
                qSeq = qSeq[cigar[0][0]:]
                start_clipped_bases += cigar[0][0]
            cigar = cigar[1:]
        while cigar[-1][1] not in 'M=X':
            if cigar[-1][1] in 'IP':
                tSeq = tSeq[:-cigar[-1][0]]
            else:
                qSeq = qSeq[:-cigar[-1][0]]
                end_clipped_bases += cigar[0][0]
            cigar = cigar[:-1]

        qLen = sum([reg_len for reg_len, reg_type in cigar
                    if reg_type in 'MIP=X'])
        assert len(qSeq) == qLen, 'Read sequence from SAM and ' + \
            'cooresponding cigar string do not agree.'

        return tSeq, qSeq, start_clipped_bases, end_clipped_bases, cigar

    def get_align_stats(tSeq, qSeq, cigar, strand):
        num_ins, num_del, num_match, num_mismatch = 0, 0, 0, 0
        tPos, qPos = 0, 0
        for reg_len, reg_type in cigar:
            if reg_type in 'M=X':
                num_reg_match = sum(
                    qBase == tBase for qBase, tBase in
                    zip(qSeq[qPos:qPos+reg_len],
                        tSeq[tPos:tPos+reg_len]))
                num_match += num_reg_match
                num_mismatch += reg_len - num_reg_match
                tPos += reg_len
                qPos += reg_len
            elif reg_type in 'IP':
                num_ins += reg_len
                qPos += reg_len
            else:
                num_del += reg_len
                tPos += reg_len

        return num_ins, num_del, num_match, num_mismatch

    strand = '-' if int(r_sam_record['flag']) & 0x10 else '+'
    cigar = parse_cigar(r_sam_record['cigar'])
    if skip_align_stats:
        # if alignment statistics are not requested, then only the template
        # (genome) sequence is required and can be parsed slightly more quickly
        # not command line option is available for this at the moment, but
        # could be easily added with this code present. The resquiggle command
        # has a lot of command line options and this seems a better default
        (tSeq, start_clipped_bases,
         end_clipped_bases) = get_just_tseq(cigar, strand)
        num_ins, num_del, num_match, num_mismatch = 0, 0, 0, 0
    else:
        qSeq, start_clipped_bases, end_clipped_bases, cigar = get_qseq(
            cigar, strand)
        tSeq, qSeq, start_clipped_bases, end_clipped_bases, cigar = get_tseq(
            qSeq, start_clipped_bases, end_clipped_bases, cigar, strand)
        num_ins, num_del, num_match, num_mismatch = get_align_stats(
            tSeq, qSeq, cigar, strand)

    read_info = readInfo(
        read_id, bc_subgroup, start_clipped_bases, end_clipped_bases,
        num_ins, num_del, num_match, num_mismatch)
    genome_loc = th.genomeLoc(
        int(r_sam_record['pos']) - 1, strand, r_sam_record['rName'])

    return tSeq, genome_loc, read_info

def parse_sam_output(align_output, batch_reads_data, genome_index):
    """
    Parse a batch of sam formatted alignment
    """
    # create dictionary with empty slot to each read
    alignments = dict(
        (read_fn_sg, None) for read_fn_sg in batch_reads_data.keys())
    for line in align_output:
        if line.startswith('@'): continue
        r_sam_record = dict(zip(SAM_FIELDS, line.strip().split()))
        if len(r_sam_record) < len(SAM_FIELDS): continue
        if r_sam_record['rName'] == '*': continue
        # store the alignment if none is stored for this read or
        # if this read has the highest map quality thus far
        qName = r_sam_record['qName'].replace(th.FN_SPACE_FILLER, ' ')
        if alignments[qName] is None or (
                int(alignments[qName]['mapq']) > int(r_sam_record['mapq'])
                and int(r_sam_record['mapq']) < 255):
            alignments[qName] = r_sam_record

    batch_align_failed_reads = []
    batch_align_data = []
    for read_fn_sg, r_sam_record in alignments.iteritems():
        bc_subgroup, read_fn = read_fn_sg.split(th.FASTA_NAME_JOINER)
        read_id = batch_reads_data[read_fn_sg][1]
        if r_sam_record is None:
            batch_align_failed_reads.append(
                ('Alignment not produced (if all reads failed ' +
                 'check for index files)', read_fn_sg))
        else:
            try:
                batch_align_data.append((read_fn, parse_sam_record(
                    r_sam_record, genome_index, read_id, bc_subgroup)))
            except Exception as e:
                #raise
                batch_align_failed_reads.append((str(e), read_fn_sg))

    return batch_align_failed_reads, batch_align_data

def _prep_graphmap_options(
        genome_fn, read_fn, out_fn, output_format, num_align_ps):
    return ['align', '-r', genome_fn, '-d', read_fn, '-o', out_fn,
            '-L', output_format, '-t', str(num_align_ps)]

def _prep_bwa_mem_options(genome_fn, read_fn, num_align_ps):
    return ['mem', '-x', 'ont2d', '-v', '1', '-t', str(num_align_ps),
            genome_fn, read_fn]

def _prep_minimap2_options(genome_fn, read_fn, num_align_ps, index_fn):
    mapper_genome = genome_fn if index_fn is None else index_fn
    return ['-ax', 'map-ont', '-t', str(num_align_ps), mapper_genome, read_fn]

def align_to_genome(batch_reads_data, genome_fn, mapper_data, genome_index,
                    num_align_ps, output_format='sam'):
    """
    Align a batch of reads to the provided genome
    """
    # prepare fasta text with batch reads
    batch_reads_fasta = ''
    for read_fn_sg, (read_seq, read_id) in batch_reads_data.iteritems():
        # note spaces aren't allowed in read names so replace with
        # vertical bars and undo to retain file names
        batch_reads_fasta += ">" + read_fn_sg.replace(' ', th.FN_SPACE_FILLER) + \
                             '\n' + read_seq + '\n'

    read_fp = NamedTemporaryFile(suffix='.fasta')
    read_fp.write(batch_reads_fasta)
    read_fp.flush()
    out_fp = NamedTemporaryFile()

    # optionally suppress output from mapper with devnull sink
    with open(os.devnull, 'w') as FNULL:
        if mapper_data.type == 'graphmap':
            mapper_options = _prep_graphmap_options(
                genome_fn, read_fp.name, out_fp.name,
                output_format, num_align_ps)
            stdout_sink = FNULL
        elif mapper_data.type == 'bwa_mem':
            mapper_options = _prep_bwa_mem_options(
                genome_fn, read_fp.name, num_align_ps)
            stdout_sink = out_fp
        elif mapper_data.type == 'minimap2':
            mapper_options = _prep_minimap2_options(
                genome_fn, read_fp.name, num_align_ps, mapper_data.index)
            stdout_sink = out_fp
        else:
            raise RuntimeError, 'Mapper not supported'

        try:
            exitStatus = call([mapper_data.exe,] + mapper_options,
                              stdout=stdout_sink, stderr=FNULL)
            out_fp.seek(0)
            align_output = out_fp.readlines()
            # close files here so that they persist until
            # after basecalling is finished
            read_fp.close()
            out_fp.close()
        except:
            # whole mapping call failed so all reads failed
            return ([(
                'Problem running/parsing genome mapper. ' +
                'Ensure you have a compatible version installed.' +
                'Potentially failed to locate BWA index files.',
                read_fn_sg) for read_fn_sg
                     in batch_reads_data.keys()], [])

    if output_format == 'sam':
        batch_parse_failed_reads, batch_align_data = parse_sam_output(
            align_output, batch_reads_data, genome_index)
    elif output_format == 'm5':
        batch_parse_failed_reads, batch_align_data = parse_m5_output(
            align_output, batch_reads_data)
    else:
        raise RuntimeError, 'Mapper output type not supported'

    return batch_parse_failed_reads, batch_align_data

def get_read_seq(fast5_fn, basecall_group, basecall_subgroup):
    """
    Extract the read sequence from the Fastq slot providing useful error messages
    """
    try:
        fast5_data = h5py.File(fast5_fn, 'r')
    except:
        raise NotImplementedError, (
            'Error opening file for alignment. This should have ' +
            'been caught during the FAST5 prep phase. Check that there ' +
            'are no other tombo processes or processes accessing ' +
            'these FAST5 files running simultaneously')

    try:
        read_seq = fast5_data[
            '/Analyses/' + basecall_group + '/' + basecall_subgroup +
            '/Fastq'].value.split('\n')[1]
    except:
        raise RuntimeError, ('Fastq slot not present in --basecall-group')

    try:
        read_data = fast5_data['/Raw/Reads/'].values()[0]
    except:
        raise RuntimeError, (
            'Raw data is not found in Raw/Reads/Read_[read#]')

    # looks like read_id attribute has been removed in some files and attribute
    # is not really necessary for tombo
    try:
        read_id = read_data.attrs['read_id']
    except:
        try:
            read_id = str(read_data.attrs['read_num'])
        except:
            read_id = str(np.random.randint(1000000000))

    try:
        if th.is_read_rna(fast5_data):
            read_seq = th.rev_transcribe(read_seq)
    except:
        raise RuntimeError, 'Error determining whether read is DNA or RNA'

    try:
        fast5_data.close()
    except:
        raise RuntimeError, (
            'Could not close FAST5 file. Possibly corrupted file')

    return read_seq, read_id

def align_and_parse(
        fast5s_to_process, genome_fn, mapper_data,
        genome_index, basecall_group, basecall_subgroups, num_align_ps):
    """
    Align and parse a batch of reads
    """
    batch_reads_data = {}
    batch_get_data_failed_reads = []
    for fast5_fn in fast5s_to_process:
        for bc_subgroup in basecall_subgroups:
            try:
                read_seq, read_id = get_read_seq(
                    fast5_fn, basecall_group, bc_subgroup)
                batch_reads_data[bc_subgroup + th.FASTA_NAME_JOINER +
                                 fast5_fn] = (read_seq, read_id)
            except Exception as e:
                # uncomment to identify mysterious errors
                #raise
                batch_get_data_failed_reads.append((
                    str(e), bc_subgroup + th.FASTA_NAME_JOINER + fast5_fn))

    batch_align_failed_reads, batch_align_data = align_to_genome(
        batch_reads_data, genome_fn, mapper_data,
        genome_index, num_align_ps)
    # regroup reads by filename (for 2D reads to be processed together
    # and avoid the same HDF5 file being opened simultaneuously)
    fn_batch_align_data = defaultdict(list)
    for fast5_fn, sg_align_data in batch_align_data:
        fn_batch_align_data[fast5_fn].append(sg_align_data)
    # uncomment to identify mysterious errors
    #print "Get data errors: " + str(batch_get_data_failed_reads)
    #print "Align read errors: " + str(batch_align_failed_reads)

    return (batch_get_data_failed_reads + batch_align_failed_reads,
            fn_batch_align_data)

def align_reads(
        fast5_batch, genome_fn, mapper_data, genome_index,
        basecall_group, basecall_subgroups, corrected_group,
        basecalls_q, overwrite, num_align_ps, in_place=True):
    """
    Prepare FAST5s and then align the extracted sequences
    """
    batch_prep_failed_reads = []
    fast5s_to_process = []
    for fast5_fn in fast5_batch:
        prep_result = th.prep_fast5(
            fast5_fn, corrected_group, overwrite, in_place, basecall_group)
        if prep_result is None:
            fast5s_to_process.append(fast5_fn)
        else:
            batch_prep_failed_reads.append(prep_result)

    batch_align_failed_reads, batch_align_data = align_and_parse(
        fast5s_to_process, genome_fn, mapper_data,
        genome_index, basecall_group, basecall_subgroups, num_align_ps)
    for fast5_fn, sgs_align_data in batch_align_data.iteritems():
        basecalls_q.put((fast5_fn, sgs_align_data))
    # uncomment to identify mysterious errors
    #print "Prep reads fail: " + str(batch_prep_failed_reads)
    #print "Align reads fail: " + str(batch_align_failed_reads)

    return batch_prep_failed_reads + batch_align_failed_reads

def _alignment_worker(
        fast5_q, basecalls_q, progress_q, failed_reads_q, genome_fn,
        mapper_data, basecall_group, basecall_subgroups,
        corrected_group, overwrite, num_align_ps):
    # this is only needed for sam output format (not m5)
    genome_index = th.parse_fasta(genome_fn)
    while not fast5_q.empty():
        try:
            fast5_batch = fast5_q.get(block=False)
        except Queue.Empty:
            break

        batch_failed_reads = align_reads(
            fast5_batch, genome_fn, mapper_data,
            genome_index, basecall_group, basecall_subgroups,
            corrected_group, basecalls_q, overwrite, num_align_ps)
        # if a read didn't fail here it will be counted in the resquiggle worker
        progress_q.put(len(batch_failed_reads))
        for failed_read in batch_failed_reads:
            try:
                sg_fn = failed_read[1].split(th.FASTA_NAME_JOINER)
                if len(sg_fn) == 2:
                    subgroup, fast5_fn = sg_fn
                else:
                    subgroup, fast5_fn = None, sg_fn
                th.write_error_status(
                    fast5_fn, corrected_group, subgroup, failed_read[0])
            except:
                pass
            failed_reads_q.put(failed_read)

    return

if _PROFILE_ALIGN:
    _alignment_wrapper = _alignment_worker
    def _alignment_worker(*args):
        import cProfile
        cProfile.runctx('_alignment_wrapper(*args)', globals(), locals(),
                        filename='resquiggle_align.prof')
        return

def resquiggle_all_reads(
        fast5_fns, genome_fn, mapper_data,
        basecall_group, basecall_subgroups, corrected_group, tb_model_fn,
        outlier_thresh, overwrite, align_batch_size, num_align_ps,
        align_threads_per_proc, num_resquiggle_ps, compute_sd, skip_index,
        skip_pen, match_evalue, bandwidth, obs_filter, const_scale):
    """
    Perform genomic alignment and event-less re-squiggle algorithm batched across reads
    """
    manager = mp.Manager()
    fast5_q = manager.Queue()
    # set maximum number of parsed basecalls to sit in the middle queue
    basecalls_q = manager.Queue(align_batch_size * ALIGN_BATCH_MULTIPLIER)
    failed_reads_q = manager.Queue()
    progress_q = manager.Queue()
    index_q = manager.Queue() if not skip_index else None
    num_reads = 0
    fast5_batch = []
    for fast5_fn in fast5_fns:
        num_reads += 1
        fast5_batch.append(fast5_fn)
        # put batches of reads in queue
        if num_reads % align_batch_size == 0:
            fast5_q.put(fast5_batch)
            fast5_batch = []
    if len(fast5_batch) > 0:
        fast5_q.put(fast5_batch)

    if tb_model_fn is None:
        tb_model_fn = ts.get_default_standard_ref_from_files(fast5_fns)

    align_args = (
        fast5_q, basecalls_q, progress_q, failed_reads_q, genome_fn,
        mapper_data, basecall_group, basecall_subgroups,
        corrected_group, overwrite, align_threads_per_proc)
    align_ps = []
    for p_id in xrange(num_align_ps):
        p = mp.Process(target=_alignment_worker, args=align_args)
        p.start()
        align_ps.append(p)

    rsqgl_args = (basecalls_q, progress_q, failed_reads_q, index_q,
                  basecall_group, corrected_group, tb_model_fn, outlier_thresh,
                  compute_sd, skip_pen, match_evalue, bandwidth, obs_filter,
                  const_scale)
    resquiggle_ps = []
    for p_id in xrange(num_resquiggle_ps):
        p = mp.Process(target=_resquiggle_worker, args=rsqgl_args)
        p.start()
        resquiggle_ps.append(p)

    if VERBOSE: sys.stderr.write(
            'Correcting ' + str(num_reads) + ' files with ' +
            str(len(basecall_subgroups)) + ' subgroup(s)/read(s) ' +
            'each (Will print a dot for each ' + str(PROGRESS_INTERVAL) +
            ' reads completed).\n')
    tot_num_rec_proc = 0
    failed_reads = defaultdict(list)
    all_index_data = []
    while any(p.is_alive() for p in align_ps):
        try:
            errorType, fn = failed_reads_q.get(block=False)
            failed_reads[errorType].append(fn)
        except Queue.Empty:
            try:
                num_rec_proc = progress_q.get(block=False)
                num_int_proc = (
                    ((tot_num_rec_proc + num_rec_proc) / PROGRESS_INTERVAL) -
                    (tot_num_rec_proc / PROGRESS_INTERVAL))
                if num_int_proc > 0:
                    sys.stderr.write('.' * num_int_proc)
                    sys.stderr.flush()
                tot_num_rec_proc += num_rec_proc
            except Queue.Empty:
                # don't need to check index queue since this is only
                # filled once the resquiggle proceses are completed which is
                # only after the alignment workers have finished
                sleep(1)
                continue

    # add None entried to basecalls_q to indicate that all reads have
    # been basecalled and processed
    for _ in xrange(num_resquiggle_ps):
        basecalls_q.put((None, None))

    while any(p.is_alive() for p in resquiggle_ps):
        try:
            errorType, fn = failed_reads_q.get(block=False)
            failed_reads[errorType].append(fn)
        except Queue.Empty:
            try:
                num_rec_proc = progress_q.get(block=False)
                num_int_proc = (
                    ((tot_num_rec_proc + num_rec_proc) / PROGRESS_INTERVAL) -
                    (tot_num_rec_proc / PROGRESS_INTERVAL))
                if num_int_proc > 0:
                    sys.stderr.write('.' * num_int_proc)
                    sys.stderr.flush()
                tot_num_rec_proc += num_rec_proc
            except Queue.Empty:
                if index_q is not None:
                    try:
                        proc_index_data = index_q.get(block=False)
                        all_index_data.extend(proc_index_data)
                    except Queue.Empty:
                        sleep(1)
                        continue

    # empty any entries left in queue after processes have finished
    while not failed_reads_q.empty():
        errorType, fn = failed_reads_q.get(block=False)
        failed_reads[errorType].append(fn)
    if index_q is not None:
        while not index_q.empty():
            proc_index_data = index_q.get(block=False)
            all_index_data.extend(proc_index_data)

    # print newline after read progress dots
    if VERBOSE: sys.stderr.write('\n')

    return dict(failed_reads), all_index_data

def parse_files(args):
    if VERBOSE: sys.stderr.write('Getting file list.\n')
    try:
        if not os.path.isdir(args.fast5_basedir):
            sys.stderr.write(
                '*' * 60 + '\nERROR: Provided [fast5-basedir] is ' +
                'not a directory.\n' + '*' * 60 + '\n')
            sys.exit()
        fast5_basedir = (
            args.fast5_basedir if args.fast5_basedir.endswith('/') else
            args.fast5_basedir + '/')
        files = th.get_files_list(fast5_basedir)
        if not args.skip_index:
            index_fn = th.get_index_fn(fast5_basedir, args.corrected_group)
            if os.path.exists(index_fn): os.remove(index_fn)
    except OSError:
        sys.stderr.write(
            '*' * 60 + '\nERROR: Reads base directory, a sub-directory ' +
            'or an old (hidden) index file does not appear to be ' +
            'accessible. Check directory permissions.\n' + '*' * 60 + '\n')
        sys.exit()
    if len(files) < 1:
        sys.stderr.write(
            '*' * 60 + '\nERROR: No files identified in the specified ' +
            'directory or within immediate subdirectories.\n' + '*' * 60 + '\n')
        sys.exit()

    return files, fast5_basedir, index_fn

def get_mapper_data(args):
    if all(map_exe is None for map_exe in (
            args.minimap2_executable, args.bwa_mem_executable,
            args.graphmap_executable)):
        sys.stderr.write(
            '*' * 60 + '\nERROR: Must provide either a ' + \
            'minimap2, graphmap or bwa-mem executable.\n' + '*' * 60 + '\n')
        sys.exit()
    if args.minimap2_executable is not None:
        mapper_data = mapperData(
            args.minimap2_executable, 'minimap2', args.minimap2_index)
    elif args.bwa_mem_executable is not None:
        mapper_data = mapperData(args.bwa_mem_executable, 'bwa_mem')
    else:
        mapper_data = mapperData(args.graphmap_executable, 'graphmap')

    return mapper_data

def eventless_resquiggle_main(args):
    """
    Main method for event-less resquiggle
    """
    global VERBOSE
    VERBOSE = not args.quiet

    if args.basecall_group == args.corrected_group:
        sys.stderr.write(
            '********** ERROR *********\n\t--basecall-group and ' +
            '--corrected-group must be different.\n')
        sys.exit()

    mapper_data = get_mapper_data(args)

    files, fast5_basedir, index_fn = parse_files(args)

    const_scale = None
    if args.fixed_scale is not None:
        const_scale = args.fixed_scale
    elif not args.fit_scale_per_read:
        const_scale = th.estimate_global_scale(files)

    outlier_thresh = args.outlier_threshold if (
        args.outlier_threshold > 0) else None

    # resolve processor and thread arguments
    num_proc = 2 if args.processes < 2 else args.processes
    align_threads_per_proc = int(num_proc / 2) \
                             if args.align_threads_per_process is None else \
                                args.align_threads_per_process
    num_resquiggle_ps = int(num_proc / 2) \
                        if args.resquiggle_processes is None \
                           else args.resquiggle_processes

    obs_filter = th.parse_obs_filter(args.obs_per_base_filter) \
                 if 'obs_per_base_filter' in args else None

    failed_reads, all_index_data = resquiggle_all_reads(
        files, args.genome_fasta, mapper_data,
        args.basecall_group, args.basecall_subgroups, args.corrected_group,
        args.tombo_model_filename, outlier_thresh, args.overwrite,
        args.alignment_batch_size, args.align_processes, align_threads_per_proc,
        num_resquiggle_ps, args.include_event_stdev, args.skip_index,
        args.skip_penalty, args.match_expected_value, args.bandwidth, obs_filter,
        const_scale)
    if not args.skip_index:
        th.write_index_file(all_index_data, index_fn, fast5_basedir)
    fail_summary = [(err, len(fns)) for err, fns in failed_reads.items()]
    if len(fail_summary) > 0:
        total_num_failed = sum(zip(*fail_summary)[1])
        sys.stderr.write('Failed reads summary (' + str(total_num_failed) +
                         ' total failed):\n' + '\n'.join(
                             "\t" + err + " :\t" + str(n_fns)
                             for err, n_fns in sorted(fail_summary)) + '\n')
    else:
        sys.stderr.write('All reads successfully re-squiggled!\n')
    if args.failed_reads_filename is not None:
        with open(args.failed_reads_filename, 'w') as fp:
            fp.write('\n'.join((
                err + '\t' + ', '.join(fns)
                for err, fns in failed_reads.items())) + '\n')

    return

def _args_and_main():
    import option_parsers
    eventless_resquiggle_main(
        option_parsers.get_eventless_resquiggle_parser().parse_args())
    return

if __name__ == '__main__':
    _args_and_main()
