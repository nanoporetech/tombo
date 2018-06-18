from __future__ import division, unicode_literals, absolute_import

from builtins import int, range, dict, map, zip

import os
import io
import re
import sys
import mappy
import queue
import threading

# Future warning from cython in h5py
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
import h5py

import numpy as np
np.seterr(all='raise')
import multiprocessing as mp

from tqdm import tqdm
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
    DEL_FIX_WINDOW, MAX_DEL_FIX_WINDOW, MIN_EVENT_TO_SEQ_RATIO, MAX_RAW_CPTS,
    PHRED_BASE, SHIFT_CHANGE_THRESH, SCALE_CHANGE_THRESH, SIG_MATCH_THRESH)

from .dynamic_programming import traceback, forward_pass
from .c_helper import (
    c_new_means, c_valid_cpts_w_cap, c_valid_cpts_w_cap_t_test)
from .c_dynamic_programming import (
    c_reg_z_scores, c_banded_forward_pass, c_banded_traceback,
    c_base_z_scores, c_adaptive_banded_forward_pass)

VERBOSE = False
PROC_UPDATE_INTERVAL = 100

_PROFILE_RSQGL = False

_DEBUG_FIT = False
_DEBUG_FULL = False
_DEBUG_MIDDLE = False
_DEBUG_PARAMS = False
_DRY_RUN = any((_DEBUG_PARAMS, _DEBUG_FIT, _DEBUG_FULL, _DEBUG_MIDDLE))
_NUM_DEBUG_ENDS = 250

MAX_QUEUE_SIZE = 1000

###############################################
########## Read Segmentation Scoring ##########
###############################################

def get_read_seg_score(r_means, r_ref_means, r_ref_sds):
    return np.mean([
        np.abs((b_m - b_ref_m) / b_ref_s)
        for b_m, b_ref_m, b_ref_s in zip(r_means, r_ref_means, r_ref_sds)])


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
    r_means = c_new_means(norm_signal, segs)
    mean_half_z_score = get_read_seg_score(r_means, r_ref_means, r_ref_sds)
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
        extra_sig_factor=EXTRA_SIG_FACTOR, max_half_z_score=None):
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
            0, n_events, n_events, min_obs_per_base,
            max_half_z_score=max_half_z_score)
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
        reg_id=None, debug_fps=None, max_half_z_score=None):
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

    shifted_z_scores = np.empty((band_event_starts.shape[0], bandwidth))
    for seq_pos, event_pos in enumerate(band_event_starts):
        if max_half_z_score is None:
            shifted_z_scores[seq_pos,:] = z_shift - np.abs(
                event_means[event_pos:event_pos + bandwidth]
                - r_ref_means[seq_pos]) / r_ref_sds[seq_pos]
        else:
            shifted_z_scores[seq_pos,:] = z_shift - np.minimum(
                max_half_z_score, np.abs(
                    event_means[event_pos:event_pos + bandwidth]
                    - r_ref_means[seq_pos]) / r_ref_sds[seq_pos])

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
        mask_bases=MASK_BASES, reg_id=None, debug_fps=None,
        max_half_z_score=None):
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
    assert event_means.shape[0] - mapped_start_offset >= bandwidth, (
        'Read sequence to signal matching starts too far into events for ' +
        'full adaptive assignment.')
    # if max_half_z_score is none set it to valid float for cython
    # z-score computation
    if max_half_z_score is None:
        do_winsorize_z = False
        max_half_z_score = 0.0
    else:
        do_winsorize_z = True

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
        band_events_start_pos + (tmp_seq_len * events_per_base),
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
            event_vals, r_ref_means[seq_pos], r_ref_sds[seq_pos],
            do_winsorize_z=do_winsorize_z, max_half_z_score=max_half_z_score)
        masked_z_scores = np.concatenate([
            [mask_fill_z_score] * start_mask_len, b_z_scores,
            [mask_fill_z_score] * end_mask_len])
        return masked_z_scores
    shifted_z_scores = np.empty((band_event_starts.shape[0], bandwidth))
    for seq_pos, event_pos in enumerate(band_event_starts):
        shifted_z_scores[seq_pos,:] = get_start_mask_z_score(seq_pos, event_pos)
    shifted_z_scores += z_shift
    fwd_pass, fwd_pass_move = c_banded_forward_pass(
        shifted_z_scores, band_event_starts, skip_pen, stay_pen)

    return fwd_pass, fwd_pass_move, band_event_starts, shifted_z_scores

def get_mapping_start(
        event_means, r_ref_means, r_ref_sds, skip_pen, stay_pen, z_shift,
        seq_window, bandwidth, norm_signal, valid_cpts,
        min_obs_per_base, reg_id=None, debug_fps=None, max_half_z_score=None):
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

    :returns: Start position (0-based) of seqeunce to event alignment within
        events and the mean events_per_base through the queried portion of read
    """
    if event_means.shape[0] < bandwidth:
        raise NotImplementedError('Read too short for start/end discovery')
    if r_ref_means.shape[0] < seq_window:
        raise NotImplementedError(
            'Genomic mapping too short for start/end discovery')

    # banded z-scores (moving up one event per base for start/end discovery
    start_z_scores = np.empty((seq_window, bandwidth))
    for seq_event_pos in range(seq_window):
        if max_half_z_score is None:
            start_z_scores[seq_event_pos,:] = z_shift - np.abs(
                event_means[seq_event_pos:seq_event_pos + bandwidth]
                - r_ref_means[seq_event_pos]) / r_ref_sds[seq_event_pos]
        else:
            start_z_scores[seq_event_pos,:] = z_shift - np.minimum(
                max_half_z_score, np.abs(
                    event_means[seq_event_pos:seq_event_pos + bandwidth]
                    - r_ref_means[seq_event_pos]) / r_ref_sds[seq_event_pos])
    start_band_event_starts = np.arange(seq_window, dtype=np.int64)

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
        r_ref_sds[:seq_window], min_obs_per_base,
        max_half_z_score=max_half_z_score)
    start_means = c_new_means(start_sig, start_segs)
    #if get_read_seg_score(start_means, r_ref_means[:seq_window],
    #                      r_ref_sds[:seq_window]) > sig_match_thresh:
    #    raise NotImplementedError(
    #        'Poor raw to expected signal matching at read start')

    # compute the average events per base to use for the start forward pass
    events_per_base = (start_tb[-1] - start_tb[0]) / len(start_tb)
    start_loc = start_tb[0]

    return start_loc, events_per_base

def find_adaptive_base_assignment(
        norm_signal, running_stat_width, min_obs_per_base, num_events, std_ref,
        genome_seq, genome_loc, skip_pen, stay_pen, z_shift, bandwidth, is_rna,
        start_bandwidth=START_BANDWIDTH, start_seq_window=START_SEQ_WINDOW,
        band_boundary_thresh=BAND_BOUNDARY_THRESH, reg_id=None, debug_fps=None,
        max_half_z_score=None):
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
        #valid_cpts = ts.get_valid_cpts(
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

    # for short reads, just search the whole read with an appropriate bandwidth
    if (event_means.shape[0] < start_bandwidth + start_seq_window or
        seq_len < max(start_seq_window, bandwidth / 2)):
        seq_events = get_short_read_event_mapping(
            event_means, r_ref_means, r_ref_sds, skip_pen, stay_pen,
            z_shift, reg_id=reg_id, debug_fps=debug_fps,
            max_half_z_score=max_half_z_score)
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
        min_obs_per_base, reg_id=reg_id,
        debug_fps=debug_fps, max_half_z_score=max_half_z_score)

    # get number of events to clip and how far into the events the
    # discovered start is located
    half_bandwidth = bandwidth // 2
    if mapped_start < half_bandwidth:
        events_start_clip = 0
        mapped_start_offset = mapped_start
    else:
        events_start_clip = mapped_start - half_bandwidth
        mapped_start_offset = half_bandwidth

    # process long enough reads that start too far into read for normal
    # adaptive processing just as with short reads
    if (event_means.shape[0] - mapped_start_offset -
        events_start_clip < bandwidth):
        seq_events = get_short_read_event_mapping(
            event_means, r_ref_means, r_ref_sds, skip_pen, stay_pen,
            z_shift, reg_id=reg_id, debug_fps=debug_fps,
            max_half_z_score=max_half_z_score)
        seq_segs = valid_cpts[seq_events]
        read_start_rel_to_raw = seq_segs[0]
        seq_segs = seq_segs - read_start_rel_to_raw

        return (seq_segs, r_ref_means, r_ref_sds, read_start_rel_to_raw,
                genome_seq, genome_loc)

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

    # if max_half_z_score is none set it to valid float for cython
    # z-score computation
    if max_half_z_score is None:
        do_winsorize_z = False
        max_half_z_score = 0.0
    else:
        do_winsorize_z = True

    if _DEBUG_FULL or _DEBUG_MIDDLE:
        rest_z_scores = c_adaptive_banded_forward_pass(
            fwd_pass, fwd_pass_move, band_event_starts, event_means,
            r_ref_means, r_ref_sds, z_shift, skip_pen, stay_pen,
            start_seq_len, MASK_FILL_Z_SCORE, do_winsorize_z, max_half_z_score,
            return_z_scores=True)
        shifted_z_scores = np.empty((seq_len, bandwidth), dtype=np.float64)
        shifted_z_scores[:start_seq_len] = start_z_scores
        shifted_z_scores[start_seq_len:] = rest_z_scores
    else:
        c_adaptive_banded_forward_pass(
            fwd_pass, fwd_pass_move, band_event_starts, event_means,
            r_ref_means, r_ref_sds, z_shift, skip_pen, stay_pen,
            start_seq_len, MASK_FILL_Z_SCORE, do_winsorize_z, max_half_z_score)

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
        all_raw_signal, channel_info, genome_seq, genome_loc,
        align_info, std_ref, outlier_thresh, corr_grp,
        bio_samp_type, seg_params, sig_aln_params,
        fast5_fn=None, max_raw_cpts=MAX_RAW_CPTS,
        min_event_to_seq_ratio=MIN_EVENT_TO_SEQ_RATIO, skip_index=False,
        reg_id=None, debug_fps=None, const_scale=None, skip_seq_scaling=False,
        scale_values=None, use_save_bandwith=False):
    """
    Perform banded dynamic programming sequence to event alignment for this read

    :param all_raw_signal: Vector containing raw (DAC) current signal values
    :param channel_info: Channel info containing info for signal normalization
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
    :param seg_params: 3 segmenation parameters (mean_obs_per_event,
        running_stat_width and min_obs_per_base)
    :param sig_aln_params: Signal align parameters (match_evalue, skip_pen
        bandwidth, save_bandwidth, signal_matching_threshold and windsorizor
        score)
    """
    # flip raw signal for re-squiggling
    is_rna = bio_samp_type == 'RNA'
    if is_rna:
        all_raw_signal = all_raw_signal[::-1]

    if sig_aln_params is None:
        (match_evalue, skip_pen, bandwidth, save_bandwidth,
         max_half_z_score) = ALGN_PARAMS_TABLE[bio_samp_type]
    else:
        # unpack signal alignment parameters
        (match_evalue, skip_pen, bandwidth, save_bandwidth,
         max_half_z_score) = sig_aln_params
        bandwidth = int(bandwidth)
        save_bandwidth = int(save_bandwidth)

    if use_save_bandwith:
        bandwidth = save_bandwidth
    z_shift, stay_pen = ts.get_dynamic_prog_params(match_evalue)

    if seg_params is None:
        (running_stat_width, min_obs_per_base,
         mean_obs_per_event) = SEG_PARAMS_TABLE[bio_samp_type]
    else:
        (running_stat_width, min_obs_per_base, mean_obs_per_event) = seg_params

    # compute number of events to find
    # ensure at least a minimal number of events per mapped sequence are found
    num_events = max(all_raw_signal.shape[0] // mean_obs_per_event,
                     int(len(genome_seq) * min_event_to_seq_ratio))
    # ensure that there isn't *far* too much signal for the mapped sequence
    # i.e. one adaptive bandwidth per base is too much to find a good mapping
    if num_events / bandwidth > len(genome_seq):
        raise NotImplementedError('Too much raw signal for mapped sequence')

    # normalize signal
    # note that channel_info is only used for pA normalization, which is not
    # available here. This option is retained here in case some channel
    # info should become useful in the future. The primary target for this is
    # the before median parameter.
    if scale_values is not None:
        norm_signal, scale_values = ts.normalize_raw_signal(
            all_raw_signal, 0, all_raw_signal.shape[0],
            scale_values=scale_values)
    elif const_scale is not None:
        norm_signal, scale_values = ts.normalize_raw_signal(
            all_raw_signal, 0, all_raw_signal.shape[0],
            'median_const_scale', channel_info, outlier_thresh,
            const_scale=const_scale)
    else:
        norm_signal, scale_values = ts.normalize_raw_signal(
            all_raw_signal, 0, all_raw_signal.shape[0],
            'median', channel_info, outlier_thresh)

    (segs, r_ref_means, r_ref_sds, read_start_rel_to_raw,
     genome_seq, genome_loc) = find_adaptive_base_assignment(
         norm_signal, running_stat_width, min_obs_per_base, num_events, std_ref,
         genome_seq, genome_loc, skip_pen, stay_pen, z_shift, bandwidth, is_rna,
         reg_id=reg_id, debug_fps=debug_fps, max_half_z_score=max_half_z_score)
    norm_signal = norm_signal[read_start_rel_to_raw:
                              read_start_rel_to_raw + segs[-1]]

    # identify all stretches of genomic deletions within del_fix_window
    # to be fixed.
    segs = get_model_fit_segs(
        segs, norm_signal, r_ref_means, r_ref_sds,
        min_obs_per_base, max_raw_cpts, max_half_z_score=max_half_z_score)

    if skip_seq_scaling:
        norm_params_changed = False
    else:
        (shift, scale, shift_corr_factor,
         scale_corr_factor) = ts.calc_kmer_fitted_shift_scale(
             scale_values.shift, scale_values.scale,
             c_new_means(norm_signal, segs), r_ref_means, method='theil_sen')
        scale_values = th.scaleValues(
            shift, scale, scale_values.lower_lim, scale_values.upper_lim)
        # re-normalize signal with new fitted parameters
        norm_signal = (norm_signal - shift_corr_factor) / scale_corr_factor
        # determine if normalization parameters changed enough to warrant
        # re-squiggling again
        norm_params_changed = (
            np.abs(shift_corr_factor) > SHIFT_CHANGE_THRESH or
            np.abs(scale_corr_factor - 1) > SCALE_CHANGE_THRESH)

    sig_match_score = get_read_seg_score(c_new_means(norm_signal, segs),
                                         r_ref_means, r_ref_sds)
    if segs.shape[0] != len(genome_seq) + 1:
        raise ValueError('Aligned sequence does not match number ' +
                         'of segments produced')

    # Output for testing/visualization of re-squiggle
    if _DEBUG_PARAMS:
        _write_params_debug(
            norm_signal, segs, r_ref_means, r_ref_sds,
            running_stat_width, min_obs_per_base, mean_obs_per_event,
            match_evalue, skip_pen, bandwidth, fast5_fn)
    if _DEBUG_FIT:
        _write_fit_debug(
            norm_signal, segs, r_ref_means, r_ref_sds, genome_seq)

    return (genome_loc, read_start_rel_to_raw, segs, genome_seq, norm_signal,
            scale_values, corr_grp, align_info, is_rna, sig_match_score,
            norm_params_changed)


#######################################
########## Genomic Alignment ##########
#######################################

def get_read_seq(fast5_data, bc_grp, bc_subgrp, bio_samp_type, q_score_thresh):
    """
    Extract the read sequence from the Fastq slot providing useful error
    messages
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

    s_fastq = fastq_raw_value.split('\n')
    read_seq, read_q = s_fastq[1], s_fastq[3]

    # compute read q-score
    if sys.version_info[0] > 2:
        mean_q_score = np.mean([q_val - PHRED_BASE
                               for q_val in read_q.encode('ASCII')])
    else:
        mean_q_score = np.mean([ord(q_val) - PHRED_BASE
                               for q_val in read_q.encode('ASCII')])
    if q_score_thresh is not None and mean_q_score < q_score_thresh:
        raise NotImplementedError('Read filtered by q-score.')

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

    return read_seq, read_id, bio_samp_type, mean_q_score

def map_read(fast5_data, bc_grp, bc_subgrp, corr_grp, aligner, bio_samp_type,
             map_thr_buf, q_score_thresh):
    read_seq, read_id, bio_samp_type, mean_q_score = get_read_seq(
        fast5_data, bc_grp, bc_subgrp, bio_samp_type, q_score_thresh)
    try:
        alignment = next(aligner.map(str(read_seq), buf=map_thr_buf))
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

    # extract genome sequence from mappy aligner
    genome_seq = aligner.seq(chrm, ref_start, ref_end)
    if sys.version_info[0] < 3:
        genome_seq = genome_seq.decode()
    if strand == '-':
        genome_seq = th.rev_comp(genome_seq)
    assert len(genome_seq) == ref_end - ref_start, (
        'Discordant mapped position and sequence')
    align_info = th.alignInfo(
        read_id, bc_subgrp, start_clipped_bases, end_clipped_bases,
        num_ins, num_del, num_match, num_aligned - num_match)
    genome_loc = th.genomeLoc(ref_start, strand, chrm)

    return genome_seq, genome_loc, align_info, bio_samp_type, mean_q_score

def _io_and_map_read(
        fast5_data, failed_reads_q, bc_subgrps, bc_grp, corr_grp, aligner,
        bio_samp_type, map_thr_buf, fast5_fn, num_processed, map_conn,
        outlier_thresh, compute_sd, obs_filter, index_q, q_score_thresh,
        sig_match_thresh):
    try:
        # extract channel and raw data for this read
        channel_info = th.get_channel_info(fast5_data)
        all_raw_signal = th.get_raw_read_slot(fast5_data)['Signal'].value
    except:
        failed_reads_q.put(
            ('Channel or raw signal information not found in FAST5 file',
             fast5_fn))
        return

    for bc_subgrp in bc_subgrps:
        try:
            # TODO warn if reads appear to switch bio sample type
            (genome_seq, genome_loc, align_info, bio_samp_type,
             mean_q_score) = map_read(
                 fast5_data, bc_grp, bc_subgrp, corr_grp, aligner,
                 bio_samp_type, map_thr_buf, q_score_thresh)
            if th.invalid_seq(genome_seq):
                raise NotImplementedError(
                    'Reference mapping contains non-canonical bases ' +
                    '(transcriptome reference cannot contain U bases)')
            map_conn.send([
                all_raw_signal, channel_info, fast5_fn, genome_seq,
                genome_loc, align_info, bio_samp_type, num_processed])
            # wait until re-squiggle returns
            read_failed, rsqgl_data = map_conn.recv()
            if read_failed:
                failed_reads_q.put((
                    rsqgl_data[0], bc_subgrp + ':::' + rsqgl_data[1]))
                continue

            # unpack data needed to write new event data
            # this is the return data from resquiggle_read
            (genome_loc, read_start_rel_to_raw, segs, genome_seq,
             norm_signal, scale_values, corr_grp, align_info,
             is_rna, sig_match_score) = rsqgl_data
            if not _DRY_RUN:
                # write re-squiggle event assignment to the read FAST5 file
                th.write_new_fast5_group(
                    fast5_data, genome_loc, read_start_rel_to_raw, segs,
                    genome_seq, norm_signal, scale_values, corr_grp,
                    align_info.Subgroup, 'median', outlier_thresh,
                    compute_sd, align_info=align_info, rna=is_rna,
                    sig_match_score=sig_match_score)

            if index_q is not None:
                # check that read passes reversible filters
                is_filtered = False
                if sig_match_score > sig_match_thresh:
                    failed_reads_q.put((
                        'Poor raw to expected signal matching ' +
                        '(revert with `tombo clear_filters`)',
                        bc_subgrp + ':::' + fast5_fn))
                    is_filtered = True
                elif obs_filter is not None:
                    base_lens = np.diff(segs)
                    is_filtered = any(np.percentile(base_lens, pctl) > thresh
                                      for pctl, thresh in obs_filter)
                    failed_reads_q.put((
                        'Read filtered by observation per base ' +
                        'thresholds (revert with `tombo clear_filters`)',
                        bc_subgrp + ':::' + fast5_fn))
                # prep and load data into index queue
                index_q.put(th.prep_index_data(
                    fast5_fn, genome_loc, read_start_rel_to_raw, segs,
                    corr_grp, align_info.Subgroup, is_rna, is_filtered,
                    sig_match_score, mean_q_score))
        except Exception as e:
            # uncomment to identify mysterious errors
            #map_conn.send(None)
            #raise
            try:
                th.write_error_status(
                    fast5_fn, corr_grp, bc_subgrp, unicode(e))
            except:
                pass
            failed_reads_q.put((
                unicode(e), bc_subgrp + ':::' + fast5_fn))

    return


#########################################
########## Re-squiggle Workers ##########
#########################################

def _resquiggle_worker(
        rsqgl_conns, std_ref, outlier_thresh, corr_grp, bio_samp_type,
        seg_params, sig_aln_params, skip_index, const_scale, skip_seq_scaling,
        max_scaling_iters):
    debug_fps = None
    if _DEBUG_MIDDLE or _DEBUG_FULL:
        debug_fps = _open_debug_fps()

    while len(rsqgl_conns) > 0:
        # get next active connection or wait for one to be ready
        try:
            conn_num, rsqgl_conn = next(
                (conn_num, rsqgl_conn)
                for conn_num, rsqgl_conn in enumerate(rsqgl_conns)
                if rsqgl_conn.poll())
        except StopIteration:
            sleep(0.1)
            continue

        try:
            map_info = rsqgl_conn.recv()
            if map_info is None:
                # this thread has finished the reads queue
                del rsqgl_conns[conn_num]
                continue

            (all_raw_signal, channel_info, fast5_fn, genome_seq, genome_loc,
             align_info, bio_samp_type, reg_id) = map_info

            rsqgl_data = resquiggle_read(
                all_raw_signal, channel_info, genome_seq, genome_loc,
                align_info, std_ref, outlier_thresh, corr_grp, bio_samp_type,
                seg_params, sig_aln_params, fast5_fn=fast5_fn,
                skip_index=skip_index, reg_id=reg_id, debug_fps=debug_fps,
                const_scale=const_scale, skip_seq_scaling=skip_seq_scaling)
            n_iters = 1
            while n_iters < max_scaling_iters and rsqgl_data[-1]:
                rsqgl_data = resquiggle_read(
                    all_raw_signal, channel_info, genome_seq, genome_loc,
                    align_info, std_ref, outlier_thresh, corr_grp, bio_samp_type,
                    seg_params, sig_aln_params, fast5_fn=fast5_fn,
                    skip_index=skip_index, reg_id=reg_id, debug_fps=debug_fps,
                    skip_seq_scaling=skip_seq_scaling, scale_values=rsqgl_data[5])
                n_iters += 1
        except Exception as e:
            try:
                rsqgl_data = resquiggle_read(
                    all_raw_signal, channel_info, genome_seq, genome_loc,
                    align_info, std_ref, outlier_thresh, corr_grp, bio_samp_type,
                    seg_params, sig_aln_params, fast5_fn=fast5_fn,
                    skip_index=skip_index, reg_id=reg_id, debug_fps=debug_fps,
                    const_scale=const_scale, skip_seq_scaling=skip_seq_scaling,
                    use_save_bandwith=True)
                n_iters = 1
                while n_iters < max_scaling_iters and rsqgl_data[-1]:
                    rsqgl_data = resquiggle_read(
                        all_raw_signal, channel_info, genome_seq, genome_loc,
                        align_info, std_ref, outlier_thresh, corr_grp,
                        bio_samp_type, seg_params, sig_aln_params,
                        fast5_fn=fast5_fn, skip_index=skip_index,
                        reg_id=reg_id, debug_fps=debug_fps,
                        skip_seq_scaling=skip_seq_scaling,
                        scale_values=rsqgl_data[5],
                        use_save_bandwith=True)
                    n_iters += 1
            except Exception as e:
                # uncomment to identify mysterious errors
                # added connection closing to avoid deadlocks here
                #for rsqgl_conn in rsqgl_conns:
                #    rsqgl_conn.send(None)
                #raise
                rsqgl_conn.send([True, [unicode(e), fast5_fn]])
                continue

        rsqgl_conn.send([False, rsqgl_data[:-1]])

    return

if _PROFILE_RSQGL:
    _resquiggle_wrapper = _resquiggle_worker
    def _resquiggle_worker(*args):
        import cProfile
        cProfile.runctx('_resquiggle_wrapper(*args)', globals(), locals(),
                        filename='resquiggle_main.prof')
        return

def _io_and_mappy_thread_worker(
        fast5_q, progress_q, failed_reads_q, index_q, bc_grp, bc_subgrps,
        corr_grp, aligner, outlier_thresh, compute_sd, sig_aln_params,
        sig_match_thresh, obs_filter, bio_samp_type, overwrite, map_conn,
        q_score_thresh):
    # get mappy aligner thread buffer
    map_thr_buf = mappy.ThreadBuffer()

    num_processed = 0
    while True:
        try:
            fast5_fn = fast5_q.get(block=False)
        except queue.Empty:
            sleep(0.1)
            continue

        if fast5_fn is None:
            # signal that all reads have been processed to child process
            map_conn.send(None)
            # update with all reads processed from this thread
            progress_q.put(num_processed % PROC_UPDATE_INTERVAL)
            break

        num_processed += 1
        if num_processed % PROC_UPDATE_INTERVAL == 0:
            progress_q.put(PROC_UPDATE_INTERVAL)

        if _DRY_RUN:
            prep_result = h5py.File(fast5_fn, 'r')
        else:
            # prep the fast5 file for writing
            prep_result = th.prep_fast5(
                fast5_fn, corr_grp, overwrite, True, bc_grp, return_fp=True)
        if isinstance(prep_result, h5py.File):
            fast5_data = prep_result
        else:
            failed_reads_q.put(prep_result)
            continue

        try:
            _io_and_map_read(
                fast5_data, failed_reads_q, bc_subgrps, bc_grp, corr_grp,
                aligner, bio_samp_type, map_thr_buf, fast5_fn,
                num_processed, map_conn, outlier_thresh, compute_sd,
                obs_filter, index_q, q_score_thresh, sig_match_thresh)
        finally:
            try:
                fast5_data.close()
            except:
                failed_reads_q.put(('Error closing fast5 file', fast5_fn))

    return


############################################
########## Multi-process Handling ##########
############################################

def _get_progress_queue(progress_q, prog_conn, max_value):
    if VERBOSE:
        th._status_message(
            'Re-squiggling reads (raw signal to genomic sequence alignment).')
        bar = tqdm(total=max_value, smoothing=0)

    tot_num_rec_proc = 0
    while True:
        try:
            iter_val = progress_q.get(block=False)
            tot_num_rec_proc += iter_val
            if VERBOSE: bar.update(iter_val)
        except queue.Empty:
            if prog_conn.poll():
                break
            sleep(0.1)
            continue

    if VERBOSE: bar.close()
    prog_conn.send(tot_num_rec_proc)

    return

def _get_failed_read_queue(failed_reads_q, failed_read_conn):
    failed_reads = defaultdict(list)
    # continue to process the failed reads queue until the end signal
    # is sent via the failed_read_conn
    while True:
        try:
            errorType, fn = failed_reads_q.get(block=False)
            failed_reads[errorType].append(fn)
        except queue.Empty:
            if failed_read_conn.poll():
                break
            sleep(0.1)
            continue

    # empty any entries left in queue after processes have finished
    while not failed_reads_q.empty():
        errorType, fn = failed_reads_q.get(block=False)
        failed_reads[errorType].append(fn)

    failed_read_conn.send(dict(failed_reads))

    return

def _get_index_queue(index_q, index_conn):
    all_index_data = []
    # continue to process the index queue until the end signal
    # is sent via the index_conn
    while True:
        try:
            r_index_data = index_q.get(block=False)
            all_index_data.append(r_index_data)
        except queue.Empty:
            if index_conn.poll():
                break
            sleep(0.1)
            continue

    # empty any entries left in queue after processes have finished
    while not index_q.empty():
        r_index_data = index_q.get(block=False)
        all_index_data.append(r_index_data)

    index_conn.send(all_index_data)

    return

def _fill_files_queue(fast5_q, fast5_fns, num_threads):
    for fast5_fn in fast5_fns:
        fast5_q.put(fast5_fn)
    for _ in range(num_threads):
        fast5_q.put(None)

    return

def resquiggle_all_reads(
        fast5_fns, aligner, bc_grp, bc_subgrps, corr_grp, std_ref,
        bio_samp_type, outlier_thresh, overwrite, num_ps, threads_per_proc,
        compute_sd, skip_index, sig_aln_params, sig_match_thresh, obs_filter,
        const_scale, seg_params, q_score_thresh, skip_seq_scaling,
        max_scaling_iters):
    """
    Perform genomic alignment and re-squiggle algorithm
    """
    fast5_q = mp.Queue(maxsize=MAX_QUEUE_SIZE)
    failed_reads_q = mp.Queue()
    index_q = mp.Queue(maxsize=MAX_QUEUE_SIZE) if not skip_index else None
    progress_q = mp.Queue()

    # open all multiprocessing pipes and queues before threading
    # as opening threads before all process are open seems to cause
    # a deadlock when some processes are started.
    # starting all multiprocess objects seems to fix this.
    files_p = mp.Process(target=_fill_files_queue,
                         args=(fast5_q, fast5_fns, num_ps * threads_per_proc))
    files_p.daemon = True
    files_p.start()

    map_conns = []
    rsqgl_ps = []
    for _ in range(num_ps):
        proc_rsqgl_conns = []
        for _ in range(threads_per_proc):
            # open mp pipe to communicate with re-squiggle process
            map_conn, rsqgl_conn = mp.Pipe()
            map_conns.append(map_conn)
            proc_rsqgl_conns.append(rsqgl_conn)
        # open re-squiggle process to void intensive processing hitting the GIL
        rsqgl_args = (
            proc_rsqgl_conns, std_ref, outlier_thresh, corr_grp, bio_samp_type,
            seg_params, sig_aln_params, index_q is None, const_scale,
            skip_seq_scaling, max_scaling_iters)
        rsqgl_process = mp.Process(target=_resquiggle_worker, args=rsqgl_args)
        rsqgl_process.daemon = True
        rsqgl_process.start()
        rsqgl_ps.append(rsqgl_process)

    # start queue getter processes
    main_prog_conn, prog_conn = mp.Pipe()
    prog_p = mp.Process(target=_get_progress_queue,
                        args=(progress_q, prog_conn, len(fast5_fns)))
    prog_p.daemon = True
    prog_p.start()
    # failed read queue getter
    main_failed_read_conn, failed_read_conn = mp.Pipe()
    failed_reads_p = mp.Process(target=_get_failed_read_queue,
                                args=(failed_reads_q, failed_read_conn))
    failed_reads_p.daemon = True
    failed_reads_p.start()
    # index queue getter
    if index_q is not None:
        main_index_conn, index_conn = mp.Pipe()
        index_p = mp.Process(target=_get_index_queue, args=(index_q, index_conn))
        index_p.daemon = True
        index_p.start()

    # now open mapping thread for each map connection created above
    resquiggle_ts = []
    for map_conn in map_conns:
        map_args = (fast5_q, progress_q, failed_reads_q, index_q, bc_grp,
                    bc_subgrps, corr_grp, aligner, outlier_thresh, compute_sd,
                    sig_aln_params, sig_match_thresh, obs_filter, bio_samp_type,
                    overwrite, map_conn, q_score_thresh)
        t = threading.Thread(target=_io_and_mappy_thread_worker,
                             args=map_args)
        t.daemon = True
        t.start()
        resquiggle_ts.append(t)

    # wait for all mapping and re-squiggling workers to finish
    files_p.join()
    for rsqgl_p in rsqgl_ps:
        rsqgl_p.join()
    for t in resquiggle_ts:
        t.join()

    # in a very unlikely case the progress queue could die while the
    # main process remains active and thus we would have a deadlock here
    if prog_p.is_alive():
        # send signal to getter queue to finish and return results
        main_prog_conn.send(True)
        # returns total number of processed reads if that is needed
        main_prog_conn.recv()
    main_failed_read_conn.send(True)
    failed_reads = main_failed_read_conn.recv()
    all_index_data = None
    if index_q is not None:
        main_index_conn.send(True)
        all_index_data = main_index_conn.recv()

    return failed_reads, all_index_data


###################################
########## Main Function ##########
###################################

def _parse_files_and_lock_dirs(args):
    if VERBOSE: th._status_message('Getting file list.')
    try:
        if not os.path.isdir(args.fast5_basedir):
            th._error_message_and_exit(
                'Provided [fast5-basedir] is not a directory.')
        fast5_basedir = (
            args.fast5_basedir if args.fast5_basedir.endswith('/') else
            args.fast5_basedir + '/')
        if args.skip_index:
            index_fn = None
        else:
            index_fn = th.get_index_fn(fast5_basedir, args.corrected_group)
            if os.path.exists(index_fn): os.remove(index_fn)

        files, lock_fns = th.get_files_list_and_lock_dirs(
            fast5_basedir, args.ignore_read_locks)
    except OSError:
        th._error_message_and_exit(
            'Reads base directory, a sub-directory ' +
            'or an old (hidden) index file does not appear to be ' +
            'accessible. Check directory permissions.')
    if len(files) < 1:
        th.clear_tombo_locks(lock_fns)
        th._error_message_and_exit(
            'No files identified in the specified ' +
            'directory or within immediate subdirectories.')

    if not th.reads_contain_basecalls(
            files, args.basecall_group, num_reads=1000):
        th.clear_tombo_locks(lock_fns)
        th._error_message_and_exit(
            'Reads do not to contain basecalls. Check --basecall-group ' +
            'option if basecalls are stored in non-standard location or ' +
            'use `tombo annotate_raw_with_fastqs` to add basecalls from ' +
            'FASTQ files to raw FAST5 files.')

    return files, fast5_basedir, index_fn, lock_fns

def _resquiggle_main(args):
    """
    Main method for resquiggle
    """
    global VERBOSE
    VERBOSE = not args.quiet
    th.VERBOSE = VERBOSE
    ts.VERBOSE = VERBOSE

    if args.print_advanced_arguments:
        from . import _option_parsers
        _option_parsers.print_advanced_resquiggle()
        sys.exit()

    if args.basecall_group == args.corrected_group:
        th._error_message_and_exit(
            '--basecall-group and --corrected-group must ' +
            'be different.')

    # check simple arguments for validity first
    outlier_thresh = args.outlier_threshold if (
        args.outlier_threshold > 0) else None
    obs_filter = th.parse_obs_filter(args.obs_per_base_filter) \
                 if 'obs_per_base_filter' in args else None

    if VERBOSE: th._status_message('Loading minimap2 reference.')
    # to be enabled when mappy genome sequence extraction bug is fixed
    aligner = mappy.Aligner(str(args.reference), preset=str('map-ont'))
    if not aligner:
        th._error_message_and_exit(
            'Failed to load reference genome FASTA for mapping.')

    # get files as late as possible in startup since it takes the longest
    # and so other errors can't happen after locks are written
    files, fast5_basedir, index_fn, lock_fns = _parse_files_and_lock_dirs(args)

    try:
        tb_model_fn = args.tombo_model_filename
        bio_samp_type = args.bio_sample_type
        if tb_model_fn is None:
            tb_model_fn, bio_samp_type = ts.get_default_standard_ref_from_files(
                files, bio_samp_type)
        else:
            bio_samp_type = 'RNA' if th.is_rna_from_files(files) else 'DNA'
        sig_match_thresh = args.signal_matching_score
        if sig_match_thresh is None:
            sig_match_thresh = SIG_MATCH_THRESH[bio_samp_type]
        if not os.path.exists(tb_model_fn):
            th._error_message_and_exit('Invalid tombo model file provided.')
        # parse tombo model
        std_ref = ts.TomboModel(tb_model_fn)

        const_scale = None
        if args.fixed_scale is not None:
            const_scale = args.fixed_scale
        elif args.fit_global_scale:
            const_scale = ts.estimate_global_scale(files)

        failed_reads, all_index_data = resquiggle_all_reads(
            files, aligner, args.basecall_group, args.basecall_subgroups,
            args.corrected_group, std_ref, bio_samp_type, outlier_thresh,
            args.overwrite, args.processes, args.threads_per_process,
            args.include_event_stdev, args.skip_index,
            args.signal_align_parameters, sig_match_thresh,
            obs_filter, const_scale, args.segmentation_parameters, args.q_score,
            args.skip_sequence_rescaling, args.max_scaling_iterations)
    finally:
        th.clear_tombo_locks(lock_fns)

    if not args.skip_index:
        th.write_index_file(all_index_data, index_fn, fast5_basedir)
    fail_summary = [(err, len(fns)) for err, fns in failed_reads.items()]
    if len(fail_summary) > 0:
        total_num_failed = sum(map(itemgetter(1), fail_summary))
        th._status_message(
            'Failed reads summary (' + unicode(total_num_failed) +
            ' total failed):\n' + '\n'.join(
                "\t" + err + " :\t" + unicode(n_fns)
                for err, n_fns in sorted(fail_summary)))
    else:
        if len(files) == len(all_index_data):
            th._status_message('All reads successfully re-squiggled!')
        else:
            th._status_message('Tombo appears to have failed unexpectedly.')
    if args.failed_reads_filename is not None:
        with io.open(args.failed_reads_filename, 'wt') as fp:
            fp.write('\n'.join((
                err + '\t' + ', '.join(fns)
                for err, fns in failed_reads.items())) + '\n')

    return

def _args_and_main():
    import _option_parsers
    resquiggle_main(
        _option_parsers.get_resquiggle_parser().parse_args())
    return

if __name__ == '__main__':
    _args_and_main()
