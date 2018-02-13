from __future__ import division, unicode_literals, absolute_import

from builtins import int, range, dict, map, zip

import os
import io
import re
import sys
import mappy
import queue

# Future warning from cython in h5py
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
import h5py

import numpy as np
np.seterr(all='raise')
import multiprocessing as mp

from time import sleep
from operator import itemgetter
from collections import defaultdict

if sys.version_info[0] > 2:
    unicode = str

# import tombo modules/functions
from . import tombo_stats as ts
from . import tombo_helper as th

from ._default_parameters import (
    SEG_PARAMS_TABLE, ALGN_PARAMS_TABLE, EXTRA_SIG_FACTOR, MASK_FILL_Z_SCORE,
    MASK_BASES, START_BANDWIDTH, START_SEQ_WINDOW, BAND_BOUNDARY_THRESH,
    DEL_FIX_WINDOW, MAX_DEL_FIX_WINDOW, MIN_EVENT_TO_SEQ_RATIO, MAX_RAW_CPTS)

from .dynamic_programming import traceback, forward_pass
from .c_helper import (
    c_new_means, c_valid_cpts_w_cap, c_valid_cpts_w_cap_t_test)
from .c_dynamic_programming import (
    c_reg_z_scores, c_banded_forward_pass, c_banded_traceback,
    c_base_z_scores, c_adaptive_banded_forward_pass)

VERBOSE = False
PROGRESS_INTERVAL = 1000

_PROFILE_RSQGL = False

_DEBUG_FIT = False
_DEBUG_FULL = False
_DEBUG_MIDDLE = False
_DEBUG_PARAMS = False
_DRY_RUN = any((_DEBUG_PARAMS, _DEBUG_FIT, _DEBUG_FULL, _DEBUG_MIDDLE))
_NUM_DEBUG_ENDS = 250


###############################################
########## Read Segmentation Scoring ##########
###############################################

def get_read_seg_score(norm_signal, segs, r_ref_means, r_ref_sds):
    return np.mean([
        np.abs((b_m - b_ref_m) / b_ref_s)
        for b_m, b_ref_m, b_ref_s in
        zip(c_new_means(norm_signal, segs), r_ref_means, r_ref_sds)])


##################################
########## Debug Output ##########
##################################

def _write_middle_debug(z_scores, fwd_pass, band_event_starts,
                       debug_fp, reg_id, debug_num_seq=_NUM_DEBUG_ENDS,
                       short=False):
    fwd_pass = fwd_pass[1:]
    if fwd_pass.shape[0] < debug_num_seq:
        debug_num_seq = fwd_pass.shape[0]
    debug_end_start = len(band_event_starts) - debug_num_seq
    debug_fp.write('\n'.join(
        '\t'.join(map(str, (band_pos + band_event_starts[seq_pos], seq_pos,
                            score, unicode(reg_id) + 'z_begin')))
        for seq_pos, s_data in enumerate(z_scores[:debug_num_seq])
        for band_pos, score in enumerate(s_data)) + '\n')
    debug_fp.write('\n'.join(
        '\t'.join(map(str, (band_pos + band_event_starts[seq_pos], seq_pos,
                            score, unicode(reg_id) + 'fwd_begin')))
        for seq_pos, s_data in enumerate(fwd_pass[:debug_num_seq])
        for band_pos, score in enumerate(s_data)) + '\n')
    if short: return

    debug_fp.write('\n'.join(
        '\t'.join(map(str, (
            band_pos + band_event_starts[debug_end_start + seq_pos], seq_pos,
            score, unicode(reg_id) + 'z_end')))
        for seq_pos, s_data in enumerate(z_scores[-debug_num_seq:])
        for band_pos, score in enumerate(s_data)) + '\n')
    debug_fp.write('\n'.join(
        '\t'.join(map(str, (
            band_pos + band_event_starts[debug_end_start + seq_pos], seq_pos,
            score, unicode(reg_id) + 'fwd_end')))
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
        seq_band_poss = [e_pos - band_event_starts[seq_pos]
                         for e_pos in range(prev_event_pos, event_pos)]
        band_poss.extend(seq_band_poss)
        event_scores.extend([z_scores[seq_pos][b_pos]
                             for b_pos in seq_band_poss])
        prev_event_pos = event_pos

    debug_fp.write('\n'.join(
        '\t'.join(map(str, (e_pos, b_pos, e_score, unicode(reg_id))))
        for e_pos, (b_pos, e_score) in enumerate(zip(
                band_poss, event_scores))) + '\n')
    fail_str = (('Failed ' if final_score < 0 else 'Pass ') +
                unicode(final_score) + ' ' + unicode(final_score / len(read_tb)))
    failed_fp.write(fail_str + '\t' + unicode(reg_id) + '\n')

    return

def _write_tb_debug(fwd_pass_move, band_event_starts, top_max_pos,
                   debug_fp, reg_id, debug_num_seq=_NUM_DEBUG_ENDS):
    read_tb = c_banded_traceback(fwd_pass_move, band_event_starts, top_max_pos)
    debug_fp.write('\n'.join(
        '\t'.join(map(str, (event_pos, seq_pos,
                            unicode(reg_id) + 'fwd_begin')))
        for seq_pos, event_pos in enumerate(read_tb[:debug_num_seq])) + '\n')
    debug_fp.write('\n'.join(
        '\t'.join(map(str, (event_pos, seq_pos,
                            unicode(reg_id) + 'fwd_end')))
        for seq_pos, event_pos in enumerate(read_tb[-debug_num_seq:])) + '\n')
    return

def _write_fit_debug(
        norm_signal, segs, r_ref_means, r_ref_sds, genome_seq):
    norm_means = c_new_means(norm_signal, segs)
    with io.open('eventless_testing.model.txt', 'wt') as fp:
        fp.write('Position\tMean\tSD\n' + '\n'.join(
            '\t'.join(map(str, (pos, p_mean, p_std)))
            for pos, (p_mean, p_std) in enumerate(zip(
                    r_ref_means, r_ref_sds))) + '\n')
    with io.open('eventless_testing.seq.txt', 'wt') as fp:
        fp.write('Base\tPosition\tSignalMean\n' + '\n'.join(
            '\t'.join(map(str, (b, pos, p_mean))) for pos, (b, p_mean) in
            enumerate(zip(genome_seq, norm_means))) + '\n')
    Position, Signal = [], []
    for base_i, (b_start, b_end) in enumerate(zip(segs[:-1], segs[1:])):
        Position.extend(
            base_i + np.linspace(0, 1, b_end - b_start, endpoint=False))
        Signal.extend(norm_signal[b_start:b_end])
    with io.open('eventless_testing.signal.txt', 'wt') as fp:
        fp.write('Position\tSignal\n' + '\n'.join(
            '\t'.join(map(str, (pos, sig)))
            for pos, sig in zip(Position, Signal)) + '\n')

    return

def _write_params_debug(
        norm_signal, segs, r_ref_means, r_ref_sds,
        running_stat_width, min_obs_per_base, mean_obs_per_event,
        match_evalue, skip_pen, bandwidth, fast5_fn):
    mean_half_z_score = get_read_seg_score(
        norm_signal, segs, r_ref_means, r_ref_sds)
    sys.stdout.write(
        '\t'.join(map(str, (
            running_stat_width, min_obs_per_base, mean_obs_per_event,
            match_evalue, skip_pen, bandwidth, fast5_fn,
            mean_half_z_score))) + '\n')

    return

def _open_debug_fps():
    score_fp = io.open('debug_event_align.txt', 'wt')
    score_fp.write('EventPos\tSeqPos\tScore\tRegion\n')
    tb_fp = io.open('debug_event_align.traceback.txt', 'wt')
    tb_fp.write('EventPos\tSeqPos\tRegion\n')
    full_fit_fp = io.open('debug_event_align.full_fit.txt', 'wt')
    full_fit_fp.write('EventPos\tBandPos\tEventScore\tRegion\n')
    full_failed_fp = io.open('debug_event_align.full_failed.txt', 'wt')
    full_failed_fp.write('DidFail\tRegion\n')
    debug_fps = [score_fp, tb_fp, full_fit_fp, full_failed_fp]

    return debug_fps


############################################
########## Raw Signal Re-squiggle ##########
############################################

def get_model_fit_segs(
        segs, norm_signal, r_ref_means, r_ref_sds, min_obs_per_base,
        max_raw_cpts=None, del_fix_window=DEL_FIX_WINDOW,
        max_del_fix_window=MAX_DEL_FIX_WINDOW,
        extra_sig_factor=EXTRA_SIG_FACTOR):
    """
    Find new segments at skipped bases during dynamic programming segmentation.

    :param segs: current Read segment locations
    :param norm_signal: Normalized read siganl
    :param r_ref_means: Read refererence means from genomic sequence
    :param r_ref_sds: Read refererence standard deviations from genomic sequence
    :param min_obs_per_base: Minimum raw observations to assign to each base
    :param max_raw_cpts: Maximum new changepoints to find from raw signal
    :param del_fix_window: initial bases to extend skipped base windows
    :param max_del_fix_window: max bases to extend skipped base windows
    :param extra_sig_factor: Amount of extra signal to require in order to
        perform signal space re-squiggle

    :returns: New segments with skipped bases resolved
    """
    def merge_del_windows(all_del_windows):
        merged_del_windows = []
        for start, end in all_del_windows:
            if (len(merged_del_windows) > 0 and
                start < merged_del_windows[-1][1]):
                merged_del_windows[-1] = (merged_del_windows[-1][0], end)
            else:
                merged_del_windows.append((start, end))

        return merged_del_windows

    def window_too_small(start, end):
        n_events = end - start
        sig_start, sig_end = segs[start], segs[end]
        sig_len = sig_end - sig_start
        # windows are expanded by one base and the extra signal factor
        # to allow some room to search for best path
        return sig_len <= ((n_events + 1) * min_obs_per_base) * extra_sig_factor

    def expand_small_windows(all_del_windows):
        expanded_del_windows = []
        windows_expanded = False
        for start, end in all_del_windows:
            if window_too_small(start, end):
                windows_expanded = True
                start -= 1
                end += 1
            expanded_del_windows.append((start, end))

        return expanded_del_windows, windows_expanded

    def trim_del_window_ends(all_del_windows):
        # potentially trim first and last windows
        if all_del_windows[0][0] < 0:
            all_del_windows[0] = (0, all_del_windows[0][1])
        if all_del_windows[-1][1] > len(segs) - 1:
            all_del_windows[-1] = (all_del_windows[-1][0], len(segs) - 1)

        return all_del_windows

    def get_deletion_windows():
        # get initial windows around deletions/skipped bases
        all_del_windows = []
        for del_pos in np.where(np.diff(segs) == 0)[0]:
            if (len(all_del_windows) > 0 and
                del_pos < all_del_windows[-1][1] + del_fix_window):
                all_del_windows[-1] = (all_del_windows[-1][0],
                                       del_pos + del_fix_window + 1)
            else:
                all_del_windows.append((del_pos - del_fix_window,
                                        del_pos + del_fix_window + 1))
        if len(all_del_windows) == 0:
            return

        windows_expanded = False
        all_del_windows = merge_del_windows(all_del_windows)
        all_del_windows = trim_del_window_ends(all_del_windows)
        # expand small windows until there are no more or the max del window
        # expansions have been attempted.
        for _ in range(max_del_fix_window - del_fix_window):
            all_del_windows, windows_expanded = expand_small_windows(
                all_del_windows)
            if not windows_expanded: break
            all_del_windows = merge_del_windows(all_del_windows)
            all_del_windows = trim_del_window_ends(all_del_windows)

        if windows_expanded and any(
                window_too_small(start, end) for start, end in all_del_windows):
            raise NotImplementedError(
                'Not enough raw signal around potential genomic deletion(s)')

        if max_raw_cpts is not None and max([
                end - start for start, end in all_del_windows]) > max_raw_cpts:
            raise NotImplementedError(
                'Read contains too many potential genomic deletions')

        return all_del_windows


    all_del_windows = get_deletion_windows()
    if all_del_windows is None:
        return segs

    for start, end in all_del_windows:
        n_events = end - start
        sig_start, sig_end = segs[start], segs[end]
        sig_len = sig_end - sig_start

        # find signal space z-scores mapping without real banding by allowing
        # entire window to be searched (c_reg_z_scores will clip base search
        # windows to enforce min_obs_per_base)
        pseudo_starts = np.linspace(0, sig_len, n_events + 1, dtype=np.int64)
        reg_z_scores = c_reg_z_scores(
            norm_signal[sig_start:sig_end], r_ref_means[start:end],
            r_ref_sds[start:end], pseudo_starts,
            0, n_events, n_events, min_obs_per_base)
        reg_fwd_scores = forward_pass(reg_z_scores, min_obs_per_base)
        # perform signal based scoring segmentation
        #  - it is ~60X faster than base space
        reg_segs = traceback(reg_fwd_scores, min_obs_per_base) + sig_start
        assert reg_segs.shape[0] == end - start - 1
        segs[start+1:end] = reg_segs

    if np.diff(segs).min() < 1:
        raise NotImplementedError('New segments include zero length events')
    if segs[0] < 0:
        raise NotImplementedError('New segments start with negative index')
    if segs[-1] > norm_signal.shape[0]:
        raise NotImplementedError('New segments end past raw signal values')

    return segs


#####################################################
########## Static Band Dynamic Programming ##########
#####################################################

def get_short_read_event_mapping(
        event_means, r_ref_means, r_ref_sds, skip_pen, stay_pen, z_shift,
        reg_id=None, debug_fps=None):
    """
    Perform banded dynamic programming sequence to event alignment
    without masking

    :param event_means: Numpy array with read base means
    :param r_ref_means: Numpy array with read reference means
    :param r_ref_sds: Numpy array with read reference standard deviations
    :param skip_pen: Penalty applied to skipped genomic bases
    :param stay_pen: Penalty applied to stay states (should shift to 0
        expected value)
    :param z_shift: Shift z-scores by this amount (includes matching
        positive expected value)

    :returns: Event to sequence mapping for full length of short read
    """
    seq_len = r_ref_means.shape[0]
    events_len = event_means.shape[0]

    # create read starts in order to clip just the corners of the full events to
    # seqeunce matrix
    mask_len = min(seq_len, events_len) // 4
    band_event_starts = np.concatenate([
        np.zeros(seq_len - mask_len * 2),
        np.linspace(0, mask_len, mask_len * 2)]).astype(np.int64)
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

    read_tb = c_banded_traceback(fwd_pass_move, band_event_starts, top_max_pos)

    return read_tb


#######################################################
########## Adaptive Band Dynamic Programming ##########
#######################################################

def get_masked_start_fwd_pass(
        event_means, r_ref_means, r_ref_sds, mapped_start_offset,
        skip_pen, stay_pen, z_shift, bandwidth, events_per_base,
        mask_fill_z_score=MASK_FILL_Z_SCORE,
        mask_bases=MASK_BASES, reg_id=None, debug_fps=None):
    """
    Perform banded dynamic programming sequence to event alignment forcing
    the path to start and end at the previously discovered locations.
    This is performed by masking the z-scores outside a "cone" extended
    mask_bases from the beginning and end of the middle of the read.

    :param event_means: Numpy array with read base means
    :param r_ref_means: Numpy array with read reference means
    :param r_ref_sds: Numpy array with read reference standard deviations
    :param mapped_start_offset: Previously identified start of genomic
        sequence within events
    :param skip_pen: Penalty applied to skipped genomic bases
    :param stay_pen: Penalty applied to stay states (should shift to 0
        expected value)
    :param z_shift: Shift z-scores by this amount (includes matching positive
        expected value)
    :param bandwidth: Bandwidth over which to search for sequence to
        event mapping
    :param events_per_base: Average events per base for the start mapping

    :returns: Event to sequence mapping for start of read including forward
        pass scores, forward pass move
        values, band starts within the events vector and z-scores
    """
    half_bandwidth = bandwidth // 2

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
        tmp_seq_len).astype(np.int64)
    mask_seq_len = max(
        mask_bases, next(i + 2 for i, bes in enumerate(band_event_starts)
                         if bes >= mapped_start_offset))
    band_event_starts = band_event_starts[:mask_seq_len]

    # get masked z-scores at the beginning of the read
    mask_start_pos = np.linspace(
        mapped_start_offset + 1,
        band_event_starts[mask_bases - 1] + bandwidth,
        mask_bases).astype(np.int64)
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
        del b_z_scores
        return masked_z_scores
    shifted_z_scores = z_shift + np.row_stack([
        get_start_mask_z_score(seq_pos, event_pos)
        for seq_pos, event_pos in enumerate(band_event_starts)])
    fwd_pass, fwd_pass_move = c_banded_forward_pass(
        shifted_z_scores, band_event_starts, skip_pen, stay_pen)

    return fwd_pass, fwd_pass_move, band_event_starts, shifted_z_scores

def get_mapping_start(
        event_means, r_ref_means, r_ref_sds, skip_pen, stay_pen, z_shift,
        seq_window, bandwidth, norm_signal, valid_cpts, score_thresh,
        min_obs_per_base, reg_id=None, debug_fps=None):
    """
    Perform banded dynamic programming sequence to event alignment through
    The beginning of an read to identify the start of genome sequence to
    event matching

    :param event_means: Numpy array with read base means
    :param r_ref_means: Numpy array with read reference means
    :param r_ref_sds: Numpy array with read reference standard deviations
    :param skip_pen: Penalty applied to skipped genomic bases
    :param stay_pen: Penalty applied to stay states (should shift to 0
        expected value)
    :param z_shift: Shift z-scores by this amount (includes matching positive
        expected value)
    :param seq_window: Number of genomic bases to search over for the start of
        the read
    :param bandwidth: Bandwidth over which to search for sequence to
        event mapping
    :param norm_signal: Normalized raw signal vector
    :param valid_cpts: Segmentation positions within norm_signal
    :param score_thresh: Read mean half-normal signal segmentation score
        threshold

    :returns: Start position (0-based) of seqeunce to event alignment within
        events and the mean events_per_base through the queried portion of read
    """
    if event_means.shape[0] < bandwidth:
        raise NotImplementedError(
            'Read too short for eventless start/end discovery')
    if r_ref_means.shape[0] < seq_window:
        raise NotImplementedError(
            'Genomic mapping too short for eventless start/end discovery')

    # banded z-scores (moving up one event per base for start/end discovery
    start_z_scores = z_shift - np.row_stack([
        np.abs(event_means[seq_pos:seq_pos + bandwidth] -
               r_ref_means[seq_pos]) / r_ref_sds[seq_pos]
        for seq_pos in range(seq_window)])
    start_band_event_starts = np.linspace(
        0, seq_window, seq_window).astype(np.int64)

    np.arange(seq_window, dtype=np.int64)

    start_fwd_pass, start_fwd_pass_move = c_banded_forward_pass(
        start_z_scores, start_band_event_starts, skip_pen, stay_pen)

    # find max along the top and right edges to start traceback
    top_max_pos = np.argmax(start_fwd_pass[-1,:])

    # perform traceback
    start_tb = c_banded_traceback(
        start_fwd_pass_move, start_band_event_starts, top_max_pos)

    # check that read start mapping is valid to avoid wasted compute on
    # adaptive dp
    start_segs = valid_cpts[start_tb]
    start_sig = norm_signal[start_segs[0]:start_segs[-1]]
    start_segs = start_segs - start_segs[0]
    start_segs = get_model_fit_segs(
        start_segs, start_sig, r_ref_means[:seq_window],
        r_ref_sds[:seq_window], min_obs_per_base)
    if get_read_seg_score(
            start_sig, start_segs, r_ref_means[:seq_window],
            r_ref_sds[:seq_window]) > score_thresh:
        raise NotImplementedError(
            'Poor raw to expected signal matching at read start')

    # compute the average events per base to use for the start forward pass
    events_per_base = (start_tb[-1] - start_tb[0]) / len(start_tb)
    start_loc = start_tb[0]

    return start_loc, events_per_base

def find_adaptive_base_assignment(
        norm_signal, running_stat_width, min_obs_per_base, num_events, std_ref,
        genome_seq, genome_loc, skip_pen, stay_pen, z_shift, bandwidth, is_rna,
        score_thresh, start_bandwidth=START_BANDWIDTH,
        start_seq_window=START_SEQ_WINDOW,
        band_boundary_thresh=BAND_BOUNDARY_THRESH, reg_id=None, debug_fps=None):
    """
    Perform banded dynamic programming sequence to event alignment by first
    identifying the start of the sequence to event matching and then
    performing banded matching through the whole read

    :param norm_signal: Numpy array with normalized read signal
    :param running_stat_width: Width of neighboring windows over which to
        compute changepoint stats
    :param min_obs_per_base: Minimum number of raw observations per base
    :param num_events: Number of events to identify in this read
    :param std_ref: A TomboModel object
    :param genome_seq: Genomic sequence for this read
    :param genome_loc: Mapped genomic location for this read
    :param skip_pen: Penalty applied to skipped genomic bases
    :param stay_pen: Penalty applied to stay states (should shift to 0
        expected value)
    :param z_shift: Shift z-scores by this amount (includes matching positive
        expected value)
    :param bandwidth: Bandwidth over which to search for sequence to
        event mapping
    :param is_rna: Is this an RNA read
    :param score_thresh: Read mean half-normal segmentation score threshold

    :returns: Start of seqeunce to event alignment and the mean
        events_per_base through the queried portion of a read
    """
    # get events
    if is_rna:
        # RNA bases show consistent variable spread so use t-test segmentation
        valid_cpts = c_valid_cpts_w_cap_t_test(
            norm_signal, min_obs_per_base, running_stat_width, num_events)
    else:
        valid_cpts = c_valid_cpts_w_cap(
            norm_signal, min_obs_per_base, running_stat_width, num_events)
        #valid_cpts = th.get_valid_cpts(
        #    norm_signal, running_stat_width, num_events)
    valid_cpts.sort()
    event_means = c_new_means(norm_signal, valid_cpts)

    dnstrm_bases = std_ref.kmer_width - std_ref.central_pos - 1
    r_ref_means, r_ref_sds, _, _ = ts.get_ref_from_seq(genome_seq, std_ref)
    # trim genome seq to match model-able positions
    genome_seq = genome_seq[std_ref.central_pos:-dnstrm_bases]
    seq_len = len(genome_seq)
    if genome_loc.Strand == '+':
        genome_loc = genome_loc._replace(
            Start=genome_loc.Start + std_ref.central_pos)
    else:
        genome_loc = genome_loc._replace(Start=genome_loc.Start + dnstrm_bases)

    # for short reads, just search the whole read with a larger bandwidth
    if (event_means.shape[0] < start_bandwidth + start_seq_window or
        seq_len < start_seq_window):
        seq_events = get_short_read_event_mapping(
            event_means, r_ref_means, r_ref_sds, skip_pen, stay_pen,
            z_shift, reg_id=reg_id, debug_fps=debug_fps)
        seq_segs = valid_cpts[seq_events]
        read_start_rel_to_raw = seq_segs[0]
        seq_segs = seq_segs - read_start_rel_to_raw

        return (seq_segs, r_ref_means, r_ref_sds, read_start_rel_to_raw,
                genome_seq, genome_loc)

    # identify the start and end of the read within the signal using a larger
    # bandwidth
    mapped_start, events_per_base = get_mapping_start(
        event_means, r_ref_means, r_ref_sds, skip_pen, stay_pen, z_shift,
        start_seq_window, start_bandwidth, norm_signal, valid_cpts,
        score_thresh, min_obs_per_base, reg_id=reg_id,
        debug_fps=debug_fps)

    # get number of events to clip and how far into the events the
    # discovered start is located
    half_bandwidth = bandwidth // 2
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
    fwd_pass = np.empty((seq_len+1, bandwidth), dtype=np.float64)
    fwd_pass[:start_seq_len+1] = start_fwd_pass
    fwd_pass_move = np.empty((seq_len+1, bandwidth), dtype=np.int64)
    fwd_pass_move[:start_seq_len+1] = start_fwd_pass_move
    band_event_starts = np.empty((seq_len,), dtype=np.int64)
    band_event_starts[:start_seq_len] = start_event_starts
    #fwd_pass[start_seq_len+1:,:] = np.NAN
    #fwd_pass_move[start_seq_len+1:,:] = np.NAN
    #band_event_starts[start_seq_len:] = np.NAN

    if _DEBUG_FULL or _DEBUG_MIDDLE:
        rest_z_scores = c_adaptive_banded_forward_pass(
            fwd_pass, fwd_pass_move, band_event_starts, event_means,
            r_ref_means, r_ref_sds, z_shift, skip_pen, stay_pen,
            start_seq_len, MASK_FILL_Z_SCORE, True)
        shifted_z_scores = np.empty((seq_len, bandwidth), dtype=np.float64)
        shifted_z_scores[:start_seq_len] = start_z_scores
        shifted_z_scores[start_seq_len:] = rest_z_scores
    else:
        c_adaptive_banded_forward_pass(
            fwd_pass, fwd_pass_move, band_event_starts, event_means,
            r_ref_means, r_ref_sds, z_shift, skip_pen, stay_pen,
            start_seq_len, MASK_FILL_Z_SCORE)

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


######################################
########## Re-squiggle Read ##########
######################################

def resquiggle_read(
        fast5_data, fast5_fn, genome_seq, genome_loc, align_info, std_ref,
        outlier_thresh, bc_grp, corr_grp, bio_samp_type, compute_sd,
        seg_params, sig_aln_params, obs_filter, max_raw_cpts=MAX_RAW_CPTS,
        min_event_to_seq_ratio=MIN_EVENT_TO_SEQ_RATIO,
        in_place=True, skip_index=False, reg_id=None, debug_fps=None,
        const_scale=None):
    """
    Perform banded dynamic programming sequence to event alignment for this read

    :param fast5_data: Open h5py object containing read information
    :param fast5_fn: Relative path to filename for index creation
    :param genome_seq: Genomic sequence for this read
    :param genome_loc: Mapped genomic location named tuple for this read
    :param align_info: A alignInfo named tuple for this read
    :param std_ref: A TomboModel object
    :param outlier_thresh: Outlier threshold for raw signal normalization
    :param bc_grp: The basecalled read group to analyze
    :param corr_grp: The tombo corrected group to write results
    :param bio_samp_type: Biological sample type (either 'DNA' or 'RNA' or
        None to determine from read)
    :param compute_sd: Should SD computations be computed and saved to file
    :param seg_params: 3 segmenation parameters (mean_obs_per_event,
        running_stat_width and min_obs_per_base)
    :param sig_aln_params: Signal align parameters (match_evalue, skip_pen
        and bandwidth)
    :param obs_filter: Obervations per base filter to apply for filtered slot
        in FAST5
    """
    try:
        channel_info = th.get_channel_info(fast5_data)

        # extract raw data for this read
        all_raw_signal = th.get_raw_read_slot(fast5_data)['Signal'].value
    except:
        #raise
        raise NotImplementedError(
            'Channel or raw signal information not found in FAST5 file')

    # flip raw signal for re-squiggling
    is_rna = bio_samp_type == 'RNA'
    if is_rna:
        all_raw_signal = all_raw_signal[::-1]

    if sig_aln_params is None:
        match_evalue, skip_pen, bandwidth, score_thresh = ALGN_PARAMS_TABLE[
            bio_samp_type]
    else:
        # unpack signal alignment parameters
        match_evalue, skip_pen, bandwidth, score_thresh = sig_aln_params
        bandwidth = int(bandwidth)
    z_shift, stay_pen = ts.get_dynamic_prog_params(match_evalue)

    if seg_params is None:
        (running_stat_width, min_obs_per_base,
         mean_obs_per_event) = SEG_PARAMS_TABLE[bio_samp_type]
    else:
        (running_stat_width, min_obs_per_base,
         mean_obs_per_event) = seg_params

    # compute number of events to find
    # ensure at least a minimal number of events per mapped sequence are found
    num_events = max(all_raw_signal.shape[0] // mean_obs_per_event,
                     int(len(genome_seq) * min_event_to_seq_ratio))
    # ensure that there isn't *far* too much signal for the mapped sequence
    # i.e. one adaptive bandwidth per base is too much to find a good mapping
    if num_events / bandwidth > len(genome_seq):
        raise NotImplementedError('Too much raw signal for mapped sequence')
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
         norm_signal, running_stat_width, min_obs_per_base, num_events, std_ref,
         genome_seq, genome_loc, skip_pen, stay_pen, z_shift, bandwidth, is_rna,
         score_thresh, reg_id=reg_id, debug_fps=debug_fps)
    norm_signal = norm_signal[read_start_rel_to_raw:
                              read_start_rel_to_raw + segs[-1]]

    # identify all stretches of genomic deletions within del_fix_window
    # to be fixed.
    segs = get_model_fit_segs(
        segs, norm_signal, r_ref_means, r_ref_sds,
        min_obs_per_base, max_raw_cpts)

    if get_read_seg_score(
            norm_signal, segs, r_ref_means, r_ref_sds) > score_thresh:
        raise NotImplementedError('Poor raw to expected signal matching')
    if segs.shape[0] != len(genome_seq) + 1:
        raise ValueError('Aligned sequence does not match number ' +
                         'of segments produced')

    # Output for testing/visualization of event-less re-squiggle
    if _DEBUG_PARAMS:
        _write_params_debug(
            norm_signal, segs, r_ref_means, r_ref_sds,
            running_stat_width, min_obs_per_base, mean_obs_per_event,
            match_evalue, skip_pen, bandwidth, fast5_fn)
    if _DEBUG_FIT:
        _write_fit_debug(
            norm_signal, segs, r_ref_means, r_ref_sds, genome_seq)

    if in_place:
        if not _DRY_RUN:
            # write re-squiggle event assignment to the read FAST5 file
            th.write_new_fast5_group(
                fast5_data, genome_loc, read_start_rel_to_raw, segs,
                genome_seq, norm_signal, scale_values, corr_grp,
                align_info.Subgroup, 'median', outlier_thresh,
                compute_sd, align_info=align_info, rna=is_rna)
    else:
        # create new hdf5 file to hold corrected read events
        pass

    if not skip_index:
        return th.prep_index_data(
            fast5_fn, genome_loc, read_start_rel_to_raw,
            segs, corr_grp, align_info.Subgroup, is_rna, obs_filter)

    return


#######################################
########## Genomic Alignment ##########
#######################################

def get_read_seq(fast5_data, bc_grp, bc_subgrp, bio_samp_type):
    """
    Extract the read sequence from the Fastq slot providing useful error messages
    """
    try:
        fastq_raw_value = fast5_data[
            '/Analyses/' + bc_grp + '/' + bc_subgrp + '/Fastq'].value
    except:
        raise NotImplementedError('Fastq slot not present in --basecall-group')

    # depending on how fastq data was stored it may already be encoded
    # as unicode, so this would fail.
    try:
        fastq_raw_value = fastq_raw_value.decode()
    except (TypeError, AttributeError):
        pass

    read_seq = fastq_raw_value.split('\n')[1]

    read_data = th.get_raw_read_slot(fast5_data)

    # looks like read_id attribute has been removed in some files and attribute
    # is not really necessary for tombo
    try:
        read_id = read_data.attrs['read_id']
    except:
        try:
            read_id = unicode(read_data.attrs['read_num'])
        except:
            read_id = unicode(np.random.randint(1000000000))

    try:
        if bio_samp_type is None:
            bio_samp_type = 'RNA' if th.is_read_rna(fast5_data) else 'DNA'
    except:
        raise NotImplementedError('Cannot determine whether read is DNA or RNA')
    if bio_samp_type == 'RNA':
        read_seq = th.rev_transcribe(read_seq)

    return read_seq, read_id, bio_samp_type

def map_read(fast5_data, bc_grp, bc_subgrp, corr_grp,
             aligner, genome_index, bio_samp_type):
    read_seq, read_id, bio_samp_type = get_read_seq(
        fast5_data, bc_grp, bc_subgrp, bio_samp_type)
    try:
        alignment = next(aligner.map(str(read_seq)))
    except StopIteration:
        raise NotImplementedError('Alignment not produced')

    chrm = alignment.ctg
    # subtract one to put into 0-based index
    ref_start = alignment.r_st
    ref_end = alignment.r_en
    strand = '+' if alignment.strand == 1 else '-'
    num_match = alignment.mlen
    num_ins, num_del, num_aligned = 0, 0 ,0
    for op_len, op in alignment.cigar:
        if op == 1: num_ins += op_len
        elif op in (2,3): num_del += op_len
        elif op in (0,7,8): num_aligned += op_len
        elif op == 6: pass
        else:
            # soft and hard clipping are not reported in the
            # mappy cigar
            raise NotImplementedError('Invalid cigar operation')
    if strand == '+':
        start_clipped_bases = alignment.q_st
        end_clipped_bases = len(read_seq) - alignment.q_en
    else:
        start_clipped_bases = len(read_seq) - alignment.q_en
        end_clipped_bases = alignment.q_st

    genome_seq = genome_index.get_seq(chrm, ref_start, ref_end)
    if strand == '-':
        genome_seq = th.rev_comp(genome_seq)
    assert len(genome_seq) == ref_end - ref_start, (
        'Discordant mapped position and sequence')
    align_info = th.alignInfo(
        read_id, bc_subgrp, start_clipped_bases, end_clipped_bases,
        num_ins, num_del, num_match, num_aligned - num_match)
    genome_loc = th.genomeLoc(ref_start, strand, chrm)

    return genome_seq, genome_loc, align_info, bio_samp_type

def load_minimap_index(genome_fn, mm_index):
    if mm_index:
        aligner = mappy.Aligner(str(mm_index), preset=str('map-ont'))
    else:
        aligner = mappy.Aligner(str(genome_fn), preset=str('map-ont'))

    return aligner


########################################
########## Re-squiggle Worker ##########
########################################

def _resquiggle_worker(
        fast5_q, progress_q, failed_reads_q, index_q, bc_grp, bc_subgrps,
        corr_grp, genome_fn, mm_index, tb_model_fn,
        outlier_thresh, compute_sd, sig_aln_params, obs_filter,
        const_scale, bio_samp_type, seg_params, overwrite, in_place=True):
    num_processed = 0
    debug_fps = None
    if _DEBUG_MIDDLE or _DEBUG_FULL:
        debug_fps = _open_debug_fps()

    # create minimap2 aligner instance
    aligner = load_minimap_index(genome_fn, mm_index)
    genome_index = th.Fasta(genome_fn)

    # parse tombo model (ignore alt_base and model_name)
    std_ref = ts.TomboModel(tb_model_fn)

    while True:
        try:
            fast5_fn = fast5_q.get(block=False)
        except queue.Empty:
            break

        num_processed += 1
        if num_processed % int(PROGRESS_INTERVAL / 10) == 0:
            progress_q.put(int(PROGRESS_INTERVAL / 10))

        if _DRY_RUN:
            prep_result = h5py.File(fast5_fn, 'r+')
        else:
            # prep the fast5 file for writing
            prep_result = th.prep_fast5(
                fast5_fn, corr_grp, overwrite, in_place, bc_grp,
                return_fp=True)
        if isinstance(prep_result, h5py.File):
            fast5_data = prep_result
        else:
            failed_reads_q.put(prep_result)
            continue

        for bc_subgrp in bc_subgrps:
            try:
                (genome_seq, genome_loc, align_info,
                 bio_samp_type) = map_read(
                     fast5_data, bc_grp, bc_subgrp, corr_grp, aligner,
                     genome_index, bio_samp_type)
                index_data = resquiggle_read(
                    fast5_data, fast5_fn, genome_seq, genome_loc, align_info,
                    std_ref, outlier_thresh, bc_grp, corr_grp, bio_samp_type,
                    compute_sd, seg_params, sig_aln_params, obs_filter,
                    skip_index=index_q is None, reg_id=num_processed,
                    debug_fps=debug_fps, const_scale=const_scale)
                if index_q is not None:
                    index_q.put(index_data)
                    if index_data[1][6]:
                        failed_reads_q.put((
                            'Read filtered by observation per base ' +
                            'thresholds (revert with `tombo clear_filters`)',
                            bc_subgrp + ':::' + fast5_fn))
            except Exception as e:
                # uncomment to identify mysterious errors
                #raise
                try:
                    th.write_error_status(
                        fast5_fn, corr_grp, bc_subgrp, unicode(e))
                except:
                    pass
                failed_reads_q.put((
                    unicode(e), bc_subgrp + ':::' + fast5_fn))

            try:
                fast5_data.close()
            except:
                pass

    return


if _PROFILE_RSQGL:
    _resquiggle_wrapper = _resquiggle_worker
    def _resquiggle_worker(*args):
        import cProfile
        cProfile.runctx('_resquiggle_wrapper(*args)', globals(), locals(),
                        filename='resquiggle_eventless_main.prof')
        return


###########################################
########## Re-squiggle All Reads ##########
###########################################

def resquiggle_all_reads(
        fast5_fns, genome_fn, bc_grp, bc_subgrps, corr_grp, tb_model_fn,
        bio_samp_type, outlier_thresh, overwrite, num_ps, compute_sd, skip_index,
        sig_aln_params, obs_filter, const_scale, seg_params, mm_index):
    """
    Perform genomic alignment and event-less re-squiggle algorithm
    batched across reads
    """
    manager = mp.Manager()
    fast5_q = manager.Queue()
    failed_reads_q = manager.Queue()
    index_q = manager.Queue() if not skip_index else None
    progress_q = manager.Queue()
    for fast5_fn in fast5_fns:
        fast5_q.put(fast5_fn)

    rsqgl_args = (fast5_q, progress_q, failed_reads_q, index_q,
                  bc_grp, bc_subgrps, corr_grp, genome_fn,
                  mm_index, tb_model_fn, outlier_thresh,
                  compute_sd, sig_aln_params, obs_filter, const_scale,
                  bio_samp_type, seg_params, overwrite)
    resquiggle_ps = []
    for p_id in range(num_ps):
        p = mp.Process(target=_resquiggle_worker, args=rsqgl_args)
        p.start()
        resquiggle_ps.append(p)

    if VERBOSE: sys.stderr.write(
            'Correcting ' + unicode(len(fast5_fns)) + ' files with ' +
            unicode(len(bc_subgrps)) + ' subgroup(s)/read(s) ' +
            'each (Will print a dot for each ' + unicode(PROGRESS_INTERVAL) +
            ' reads completed).\n')
    tot_num_rec_proc = 0
    failed_reads = defaultdict(list)
    all_index_data = []
    while any(p.is_alive() for p in resquiggle_ps):
        try:
            errorType, fn = failed_reads_q.get(block=False)
            failed_reads[errorType].append(fn)
        except queue.Empty:
            try:
                num_rec_proc = progress_q.get(block=False)
                num_int_proc = (
                    ((tot_num_rec_proc + num_rec_proc) // PROGRESS_INTERVAL) -
                    (tot_num_rec_proc // PROGRESS_INTERVAL))
                if num_int_proc > 0:
                    sys.stderr.write('.' * num_int_proc)
                    sys.stderr.flush()
                tot_num_rec_proc += num_rec_proc
            except queue.Empty:
                if index_q is not None:
                    try:
                        r_index_data = index_q.get(block=False)
                        all_index_data.append(r_index_data)
                    except queue.Empty:
                        sleep(1)
                        continue

    # empty any entries left in queue after processes have finished
    while not failed_reads_q.empty():
        errorType, fn = failed_reads_q.get(block=False)
        failed_reads[errorType].append(fn)
    if index_q is not None:
        while not index_q.empty():
            r_index_data = index_q.get(block=False)
            all_index_data.append(r_index_data)

    # print newline after read progress dots
    if VERBOSE: sys.stderr.write('\n')

    return dict(failed_reads), all_index_data

def parse_files(args):
    if VERBOSE: sys.stderr.write('Getting file list.\n')
    try:
        if not os.path.isdir(args.fast5_basedir):
            th._error_message_and_exit(
                'Provided [fast5-basedir] is not a directory.')
        fast5_basedir = (
            args.fast5_basedir if args.fast5_basedir.endswith('/') else
            args.fast5_basedir + '/')
        files = th.get_files_list(fast5_basedir)
        if args.skip_index:
            index_fn = None
        else:
            index_fn = th.get_index_fn(fast5_basedir, args.corrected_group)
            if os.path.exists(index_fn): os.remove(index_fn)
    except OSError:
        th._error_message_and_exit(
            'Reads base directory, a sub-directory ' +
            'or an old (hidden) index file does not appear to be ' +
            'accessible. Check directory permissions.')
    if len(files) < 1:
        th._error_message_and_exit(
            'No files identified in the specified ' +
            'directory or within immediate subdirectories.')

    if not th.reads_contain_basecalls(
            files, args.basecall_group, num_reads=1000):
        th._error_message_and_exit(
            'Reads do not to contain basecalls. Check --basecall-group option ' +
            'if basecalls are stored in non-standard location or use ' +
            '`tombo annotate_raw_with_fastqs` to add basecalls from FASTQ ' +
            'files to raw FAST5 files.')

    return files, fast5_basedir, index_fn


###################################
########## Main Function ##########
###################################

def eventless_resquiggle_main(args):
    """
    Main method for event-less resquiggle
    """
    global VERBOSE
    VERBOSE = not args.quiet
    th.VERBOSE = VERBOSE
    ts.VERBOSE = VERBOSE

    if args.basecall_group == args.corrected_group:
        th._error_message_and_exit(
            '--basecall-group and --corrected-group must ' +
            'be different.')

    # check simple arguments for validity first
    outlier_thresh = args.outlier_threshold if (
        args.outlier_threshold > 0) else None

    obs_filter = th.parse_obs_filter(args.obs_per_base_filter) \
                 if 'obs_per_base_filter' in args else None
    aligner = load_minimap_index(args.genome_fasta, args.minimap2_index)
    if not aligner:
        th._error_message_and_exit(
            'Failed to load --genome-fasta or --minimap2-index for mapping.')
    del aligner

    # load genome once here to index it if using pyfaidx so it isn't built in
    # each process seperately.
    genome_index = th.Fasta(args.genome_fasta, dry_run=True)
    del genome_index

    files, fast5_basedir, index_fn = parse_files(args)

    tb_model_fn = args.tombo_model_filename
    bio_samp_type = args.bio_sample_type
    if tb_model_fn is None:
        tb_model_fn, bio_samp_type = ts.get_default_standard_ref_from_files(
            files, bio_samp_type)
    if not os.path.exists(tb_model_fn):
        th._error_message_and_exit('Invalid tombo model file provided.')

    const_scale = None
    if args.fixed_scale is not None:
        const_scale = args.fixed_scale
    elif not args.fit_scale_per_read:
        const_scale = th.estimate_global_scale(files)

    failed_reads, all_index_data = resquiggle_all_reads(
        files, args.genome_fasta,
        args.basecall_group, args.basecall_subgroups, args.corrected_group,
        tb_model_fn, bio_samp_type, outlier_thresh,
        args.overwrite, args.processes, args.include_event_stdev,
        args.skip_index, args.signal_align_parameters, obs_filter,
        const_scale, args.segmentation_parameters, args.minimap2_index)
    if not args.skip_index:
        th.write_index_file(all_index_data, index_fn, fast5_basedir)
    fail_summary = [(err, len(fns)) for err, fns in failed_reads.items()]
    if len(fail_summary) > 0:
        total_num_failed = sum(map(itemgetter(1), fail_summary))
        sys.stderr.write('Failed reads summary (' + unicode(total_num_failed) +
                         ' total failed):\n' + '\n'.join(
                             "\t" + err + " :\t" + unicode(n_fns)
                             for err, n_fns in sorted(fail_summary)) + '\n')
    else:
        sys.stderr.write('All reads successfully re-squiggled!\n')
    if args.failed_reads_filename is not None:
        with io.open(args.failed_reads_filename, 'wt') as fp:
            fp.write('\n'.join((
                err + '\t' + ', '.join(fns)
                for err, fns in failed_reads.items())) + '\n')

    return

def _args_and_main():
    import _option_parsers
    eventless_resquiggle_main(
        _option_parsers.get_eventless_resquiggle_parser().parse_args())
    return

if __name__ == '__main__':
    _args_and_main()
