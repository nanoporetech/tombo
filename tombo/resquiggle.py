from __future__ import division, unicode_literals, absolute_import

from builtins import int, range, dict, map, zip

import os
import io
import re
import sys
import queue
import traceback
import threading

# pip allows tombo install without correct version of mappy, so check here
try:
    import mappy
    if sys.version_info[0] > 2:
        mappy.Aligner(os.path.devnull).seq('')
    else:
        mappy.Aligner(os.path.devnull).seq(b'')
except AttributeError:
    th.error_message_and_exit('Tombo requires mappy version >= 2.10.')

# Future warning from cython in h5py
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
import h5py

import numpy as np
np.seterr(all='raise')
import multiprocessing as mp

from tqdm import tqdm
from tqdm._utils import _term_move_up

from time import sleep
from operator import itemgetter
from collections import defaultdict
from pkg_resources import resource_string

if sys.version_info[0] > 2:
    unicode = str

# import tombo modules/functions
from . import tombo_stats as ts
from . import tombo_helper as th

from ._default_parameters import (
    EXTRA_SIG_FACTOR, MASK_FILL_Z_SCORE,
    MASK_BASES, DEL_FIX_WINDOW, MAX_DEL_FIX_WINDOW,
    MIN_EVENT_TO_SEQ_RATIO, MAX_RAW_CPTS, SHIFT_CHANGE_THRESH,
    SCALE_CHANGE_THRESH, SIG_MATCH_THRESH, DNA_SAMP_TYPE, RNA_SAMP_TYPE,
    USE_RNA_EVENT_SCALE, RNA_SCALE_NUM_EVENTS, RNA_SCALE_MAX_FRAC_EVENTS,
    START_CLIP_PARAMS, STALL_PARAMS, COLLAPSE_RNA_STALLS, COLLAPSE_DNA_STALLS)
START_CLIP_PARAMS = th.startClipParams(*START_CLIP_PARAMS)
DEFAULT_STALL_PARAMS = th.stallParams(**STALL_PARAMS)

from ._c_dynamic_programming import (
    c_reg_z_scores, c_banded_forward_pass, c_base_z_scores,
    c_base_forward_pass, c_base_traceback)


# list of classes/functions to include in API
__all__ = [
    'get_read_seq', 'map_read', 'resquiggle_read',
    'segment_signal', 'find_adaptive_base_assignment',
    'resolve_skipped_bases_with_raw', 'find_seq_start_in_events',
    'find_static_base_assignment']


VERBOSE = True

_PROFILE_RSQGL = False

# use (mapping) clipped bases at the start of read to identify start position
USE_START_CLIP_BASES = False

# experimental RNA adapter trimming
TRIM_RNA_ADAPTER = False


# text output debugging
_DEBUG_PARAMS = False
_DEBUG_BANDWIDTH = False
_DEBUG_START_BANDWIDTH = False

# plot output debugging
_DEBUG_DP_ENDS = False
_DEBUG_DP_START = False
_DEBUG_CLIP_START = False
# fit debug plot requires r cowplot package to be installed
_DEBUG_FIT = False
_DEBUG_START_CLIP_FIT = False
# raw signal re-squiggle DP
_DEBUG_RAW_DP = False

# don't plot more than one debug type at a time
assert sum((
    _DEBUG_DP_ENDS, _DEBUG_FIT, _DEBUG_START_CLIP_FIT,
    _DEBUG_DP_START, _DEBUG_CLIP_START, _DEBUG_RAW_DP)) <= 1
_DEBUG_PLOTTING = any((
    _DEBUG_FIT, _DEBUG_START_CLIP_FIT, _DEBUG_DP_ENDS, _DEBUG_DP_START,
    _DEBUG_CLIP_START, _DEBUG_RAW_DP))
_DRY_RUN = any((
    _DEBUG_PARAMS, _DEBUG_BANDWIDTH, _DEBUG_START_BANDWIDTH, _DEBUG_PLOTTING))

_UNEXPECTED_ERROR_FN = 'unexpected_tombo_errors.{}.err'
_MAX_NUM_UNEXP_ERRORS = 50
_MAX_QUEUE_SIZE = 1000


##################################
########## Debug Output ##########
##################################

def _write_params_debug(
        norm_signal, segs, r_ref_means, r_ref_sds, rsqgl_params, read_id):
    r_means = ts.compute_base_means(norm_signal, segs)
    mean_half_z_score = ts.get_read_seg_score(r_means, r_ref_means, r_ref_sds)
    sys.stdout.write(
        '\t'.join(map(str, (
            rsqgl_params.running_stat_width,
            rsqgl_params.min_obs_per_base,
            rsqgl_params.raw_min_obs_per_base,
            rsqgl_params.mean_obs_per_event,
            rsqgl_params.match_evalue,
            rsqgl_params.skip_pen,
            rsqgl_params.bandwidth, read_id,
            mean_half_z_score))) + '\n')

    return

def _debug_plot_dp(
        z_scores, fwd_pass, band_event_starts, fwd_pass_move, top_max_pos,
        reg_id='0', debug_num_seq=500, short=False):
    reg_id = unicode(reg_id)
    fwd_pass = fwd_pass[1:]
    if fwd_pass.shape[0] < debug_num_seq:
        debug_num_seq = fwd_pass.shape[0]
    debug_end_start = len(band_event_starts) - debug_num_seq

    event_poss, seq_poss, scores, regs = [], [], [], []
    for seq_pos, (s_z_data, s_f_data) in enumerate(zip(z_scores, fwd_pass)):
        b_e_start = band_event_starts[seq_pos]
        for band_pos, score in enumerate(s_z_data):
            event_poss.append(band_pos + b_e_start)
            seq_poss.append(seq_pos)
            scores.append(score)
            regs.append(reg_id + '_z_begin')
        for band_pos, score in enumerate(s_f_data):
            event_poss.append(band_pos + b_e_start)
            seq_poss.append(seq_pos)
            scores.append(score)
            regs.append(reg_id + '_fwd_begin')
        if seq_pos >= debug_num_seq:
            break
    if not short:
        for seq_pos, (s_z_data, s_f_data) in enumerate(zip(
                z_scores[::-1], fwd_pass[::-1])):
            end_seq_pos = debug_num_seq - seq_pos - 1
            b_e_start = band_event_starts[debug_end_start + end_seq_pos]
            for band_pos, score in enumerate(s_z_data):
                event_poss.append(band_pos + b_e_start)
                seq_poss.append(end_seq_pos)
                scores.append(score)
                regs.append(reg_id + '_z_end')
            for band_pos, score in enumerate(s_f_data):
                event_poss.append(band_pos + b_e_start)
                seq_poss.append(end_seq_pos)
                scores.append(score)
                regs.append(reg_id + '_fwd_end')
            if seq_pos >= debug_num_seq:
                break

    dpDat = r.DataFrame({
        'EventPos':r.IntVector(event_poss),
        'SeqPos':r.IntVector(seq_poss),
        'Score':r.FloatVector(scores),
        'Region':r.StrVector(regs)})

    event_poss, seq_poss, regs = [], [], []
    read_tb = th.banded_traceback(fwd_pass_move, band_event_starts, top_max_pos)
    for seq_pos, event_pos in enumerate(read_tb[:debug_num_seq]):
        event_poss.append(event_pos)
        seq_poss.append(seq_pos)
        regs.append(reg_id + '_fwd_begin')
    if not short:
        for seq_pos, event_pos in enumerate(read_tb[-debug_num_seq:]):
            event_poss.append(event_pos)
            seq_poss.append(seq_pos)
            regs.append(reg_id + '_fwd_end')

    tbDat = r.DataFrame({
        'EventPos':r.IntVector(event_poss),
        'SeqPos':r.IntVector(seq_poss),
        'Region':r.StrVector(regs)})

    r.r(resource_string(__name__, 'R_scripts/debugDP.R').decode())
    r.globalenv[str('plotDP')](dpDat, tbDat)

    return

def _debug_raw_dp(z_scores, fwd_pass, read_tb, sig_data, reg_id='0'):
    reg_id = unicode(reg_id)

    event_poss, seq_poss, r_z_scores, fwd_scores = [], [], [], []
    for seq_pos, ((s_z_data, (b_e_start, b_e_end)), s_f_data) in enumerate(zip(
            z_scores, map(itemgetter(0), fwd_pass))):
        for band_pos, score in enumerate(s_z_data):
            r_z_scores.append(score)
            event_poss.append(band_pos + b_e_start)
            seq_poss.append(seq_pos)
        for band_pos, score in enumerate(s_f_data):
            fwd_scores.append(score)

    zDat = r.DataFrame({
        'Score':r.FloatVector(r_z_scores),
        'EventPos':r.IntVector(event_poss),
        'SeqPos':r.IntVector(seq_poss),
        'Region':r.StrVector([reg_id,] * len(seq_poss))})
    fwdDat = r.DataFrame({
        'Score':r.FloatVector(fwd_scores),
        'EventPos':r.IntVector(event_poss),
        'SeqPos':r.IntVector(seq_poss),
        'Region':r.StrVector([reg_id,] * len(seq_poss))})

    event_poss, seq_poss = [0,], [0,]
    for seq_pos, event_pos in enumerate(read_tb):
        event_poss.append(event_pos - 1)
        seq_poss.append(seq_pos)
        event_poss.append(event_pos)
        seq_poss.append(seq_pos + 1)
    event_poss.append(z_scores[-1][1][1] - 1)
    seq_poss.append(read_tb.shape[0])

    tbDat = r.DataFrame({
        'EventPos':r.IntVector(event_poss),
        'SeqPos':r.IntVector(seq_poss),
        'Region':r.StrVector([reg_id,] * len(seq_poss))})

    sigDat = r.DataFrame({
        'Pos':r.IntVector(list(range(sig_data.shape[0]))),
        'Signal':r.FloatVector(sig_data)
    })

    r.r(resource_string(__name__, 'R_scripts/debugRawDP.R').decode())
    r.globalenv[str('plotRawDP')](zDat, fwdDat, tbDat, sigDat)

    return

def _debug_fit(
        fwd_pass_move, band_event_starts, top_max_pos, z_scores, reg_id,
        final_score, bandwidth, event_means, r_ref_means,
        running_window=501, static_bw=False):
    read_tb = th.banded_traceback(fwd_pass_move, band_event_starts, top_max_pos)
    prev_event_pos = read_tb[0]
    band_poss, event_scores, ref_means = [], [], []
    for seq_pos, event_pos in enumerate(read_tb[1:]):
        seq_band_poss = [e_pos - band_event_starts[seq_pos]
                         for e_pos in range(prev_event_pos, event_pos)]
        band_poss.extend(seq_band_poss)
        event_scores.extend([z_scores[seq_pos][b_pos]
                             for b_pos in seq_band_poss])
        ref_means.extend([r_ref_means[seq_pos] for _ in seq_band_poss])
        prev_event_pos = event_pos

    if _DEBUG_BANDWIDTH or _DEBUG_START_BANDWIDTH:
        half_bandwidth = bandwidth // 2
        min_bw_edge_buffer = (
            half_bandwidth - np.max(np.abs(
                np.array(band_poss) - half_bandwidth)) if _DEBUG_BANDWIDTH else
            max(band_poss))
        sys.stdout.write('{:d}\t{:d}\t{}\t{}\n'.format(
            bandwidth, min_bw_edge_buffer, final_score / len(read_tb), reg_id))
        sys.stdout.flush()
        if not (_DEBUG_FIT or _DEBUG_START_CLIP_FIT):
            return

    if len(event_scores) > running_window:
        # make sure window is odd
        if running_window % 2 != 1:
            running_window =+ 1
        half_window = running_window // 2
        score_cumsum = np.concatenate([[0,], np.cumsum(event_scores)])
        event_scores = np.concatenate([
            np.repeat(np.NAN, half_window),
            (score_cumsum[:-running_window] - score_cumsum[running_window:]) /
            running_window, np.repeat(np.NAN, half_window)])

    reg_name = (unicode(reg_id) + '__' + unicode(final_score) + '__' +
                unicode(final_score / len(read_tb)))

    fitDat = r.DataFrame({
        'EventPos':r.IntVector(range(len(band_poss))),
        'BandPos':r.IntVector(band_poss),
        'EventMean':r.FloatVector(event_means[read_tb[0]:read_tb[-1]]),
        'ModelMean':r.FloatVector(ref_means),
        'EventScore':r.FloatVector(event_scores),
        'Region':r.StrVector([reg_name,] * len(band_poss))})

    r.r(resource_string(__name__, 'R_scripts/debugFit.R').decode())
    r.globalenv[str('plotFit')](fitDat, bandwidth)

    return

def _open_debug_pdf():
    # import plotting modules
    from rpy2 import robjects as r
    global r
    from rpy2.robjects.packages import importr
    importr(str('ggplot2'))

    if _DEBUG_DP_ENDS:
        r.r('pdf("debug_event_align.pdf", height=4.5, width=6)')
    elif _DEBUG_CLIP_START:
        r.r('pdf("debug_event_align.clip_start.pdf", height=4.5, width=10)')
    elif _DEBUG_DP_START:
        r.r('pdf("debug_event_align.start.pdf", height=4.5, width=10)')
    elif _DEBUG_FIT:
        importr(str('cowplot'))
        r.r('pdf("debug_event_align.full_fit.pdf", width=15, height=5)')
    elif _DEBUG_START_CLIP_FIT:
        importr(str('cowplot'))
        r.r('pdf("debug_event_align.start_clip_fit.pdf", width=15, height=5)')
    elif _DEBUG_RAW_DP:
        importr(str('cowplot'))
        r.r('pdf("debug_event_align.raw_dp.pdf", width=11, height=7)')
    else:
        th.error_message_and_exit('Must specify which debug plot to open.')

    return

def _close_debug_pdf():
    r.r('dev.off()')
    return


############################################
########## Raw Signal Re-squiggle ##########
############################################

def raw_forward_pass(reg_z_scores, min_obs_per_base):
    # dynamic programming algorithm to find modeled signal to base assignment

    # fill banded path with cumulative probabilties from the previous signal
    # either in the current base or the previous base (left or diagonal left
    # from associated plotting)

    # get the first row data
    prev_b_data, (prev_b_start, prev_b_end) = reg_z_scores[0]
    prev_b_fwd_data = np.cumsum(prev_b_data)
    # store number of observations since last diagonal at each position
    # - forces forward pass to allow legal traceback paths while
    #   enforcing the minimum observations per base threshold
    # - should also help from optimization pushing poor fitting bases
    #   to assign only an observation or two
    # - will also use this data to traceback all reasonable paths
    prev_b_last_diag = np.ones(prev_b_end - prev_b_start,
                               dtype=np.int64) * min_obs_per_base
    # first row is just a cumsum since there is no previous row
    reg_fwd_scores = [(prev_b_fwd_data, prev_b_last_diag,
                       (prev_b_start, prev_b_end))]

    for b_data, (b_start, b_end) in reg_z_scores[1:]:
        b_fwd_data, prev_b_last_diag = c_base_forward_pass(
            b_data, b_start, b_end,
            prev_b_data, prev_b_start, prev_b_end,
            prev_b_fwd_data, prev_b_last_diag, min_obs_per_base)

        # consider storing data to form traceback in one go at the
        # end of this loop
        reg_fwd_scores.append((
            b_fwd_data, prev_b_last_diag, (b_start, b_end)))
        prev_b_data, prev_b_fwd_data, prev_b_start, prev_b_end = (
            b_data, b_fwd_data, b_start, b_end)

    return reg_fwd_scores

def raw_traceback(reg_fwd_scores, min_obs_per_base):
    # traceback along maximally likely path

    # initilize array to store new segments
    new_segs = np.empty(len(reg_fwd_scores) - 1, dtype=np.int64)
    # get first two bases of data for lookups
    curr_b_data, _, (curr_start, curr_end) = reg_fwd_scores[-1]
    next_b_data, _, (next_start, next_end) = reg_fwd_scores[-2]
    new_segs[-1] = c_base_traceback(
        curr_b_data, curr_start, next_b_data, next_start, next_end,
        curr_end - 1, min_obs_per_base)
    for base_pos in range(len(reg_fwd_scores) - 3, -1, -1):
        curr_b_data, curr_start = next_b_data, next_start
        next_b_data, _, (next_start, next_end) = reg_fwd_scores[base_pos]
        new_segs[base_pos] = c_base_traceback(
            curr_b_data, curr_start, next_b_data, next_start, next_end,
            new_segs[base_pos+1] - 1, min_obs_per_base)

    return new_segs

def resolve_skipped_bases_with_raw(
        dp_res, norm_signal, rsqgl_params,
        max_raw_cpts=MAX_RAW_CPTS, del_fix_window=DEL_FIX_WINDOW,
        max_del_fix_window=MAX_DEL_FIX_WINDOW,
        extra_sig_factor=EXTRA_SIG_FACTOR):
    """Perform dynamic-programming forward pass over raw signal z-scores. For details, see https://nanoporetech.github.io/tombo/resquiggle.html#resolve-skipped-bases

    Args:
        dp_res (:class:`tombo_helper.dpResults`): dynamic programming results
        norm_signal (`np.array::np.float64`): normalized raw siganl
        rsqgl_params (:class:`tombo.tombo_helper.resquiggleParams`): parameters for the re-squiggle algorithm
        max_raw_cpts (int): maximum new changepoints to find from raw signal (optional)
        del_fix_window (int): initial bases to extend skipped base windows (optional)
        max_del_fix_window (int): max bases to extend skipped base windows (optional)
        extra_sig_factor (float): amount of extra signal required to perform signal space re-squiggle (optional)

    Returns:
        ``np.array::np.int64`` containing new deletion resolved base raw signal start positions
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
        sig_start, sig_end = dp_res.segs[start], dp_res.segs[end]
        sig_len = sig_end - sig_start
        # windows are expanded by one base and the extra signal factor
        # to allow some room to search for best path
        return sig_len <= ((n_events + 1) *
                           rsqgl_params.raw_min_obs_per_base) * extra_sig_factor

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
        if all_del_windows[-1][1] > len(dp_res.segs) - 1:
            all_del_windows[-1] = (all_del_windows[-1][0], len(dp_res.segs) - 1)

        return all_del_windows

    def get_deletion_windows():
        # get initial windows around deletions/skipped bases
        all_del_windows = []
        for del_pos in np.where(np.diff(dp_res.segs) == 0)[0]:
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
            raise th.TomboError(
                'Not enough raw signal around potential genomic deletion(s)')

        if max_raw_cpts is not None and max([
                end - start for start, end in all_del_windows]) > max_raw_cpts:
            raise th.TomboError(
                'Read contains too many potential genomic deletions')

        return all_del_windows


    all_del_windows = get_deletion_windows()
    resolved_segs = dp_res.segs.copy()
    if all_del_windows is None:
        return resolved_segs

    for start, end in all_del_windows:
        n_events = end - start
        sig_start, sig_end = dp_res.segs[start], dp_res.segs[end]
        sig_len = sig_end - sig_start

        # find signal space z-scores mapping without real banding by allowing
        # entire window to be searched (c_reg_z_scores will clip base search
        # windows to enforce min_obs_per_base)
        pseudo_starts = np.linspace(0, sig_len, n_events + 1, dtype=np.int64)
        reg_z_scores = c_reg_z_scores(
            norm_signal[sig_start:sig_end], dp_res.ref_means[start:end],
            dp_res.ref_sds[start:end], pseudo_starts,
            0, n_events, n_events, rsqgl_params.raw_min_obs_per_base,
            max_half_z_score=rsqgl_params.max_half_z_score)
        reg_fwd_scores = raw_forward_pass(
            reg_z_scores, rsqgl_params.raw_min_obs_per_base)
        # perform signal based scoring segmentation
        #  - it is ~60X faster than base space
        reg_segs = raw_traceback(
            reg_fwd_scores, rsqgl_params.raw_min_obs_per_base) + sig_start
        if _DEBUG_RAW_DP:
            _debug_raw_dp(reg_z_scores, reg_fwd_scores, reg_segs - sig_start,
            norm_signal[sig_start:sig_end])
        if reg_segs.shape[0] != end - start - 1:
            raise th.TomboError('Invalid segmentation results.')
        resolved_segs[start+1:end] = reg_segs

    if np.diff(resolved_segs).min() < 1:
        raise th.TomboError('New segments include zero length events')
    if resolved_segs[0] < 0:
        raise th.TomboError('New segments start with negative index')
    if resolved_segs[-1] > norm_signal.shape[0]:
        raise th.TomboError('New segments end past raw signal values')

    return resolved_segs


#####################################################
########## Static Band Dynamic Programming ##########
#####################################################

def find_static_base_assignment(
        event_means, r_ref_means, r_ref_sds, rsqgl_params,
        reg_id=None):
    """Align expected (from genome sequence) signal levels to observed using a dynamic programming approach with static bandwidth start positions

    Args:
        event_means (`np.array::np.float64`): read base means
        r_ref_means (`np.array::np.float64`): expected base signal levels
        r_ref_sds (`np.array::np.float64`): expected base level SDs
        rsqgl_params (:class:`tombo.tombo_helper.resquiggleParams`): parameters for the re-squiggle algorithm

    Returns:
        `np.array::np.int64` containng event to sequence mapping for full length of short read
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
        if rsqgl_params.max_half_z_score is None:
            shifted_z_scores[seq_pos,:] = rsqgl_params.z_shift - np.abs(
                event_means[event_pos:event_pos + bandwidth]
                - r_ref_means[seq_pos]) / r_ref_sds[seq_pos]
        else:
            shifted_z_scores[seq_pos,:] = rsqgl_params.z_shift - np.minimum(
                rsqgl_params.max_half_z_score, np.abs(
                    event_means[event_pos:event_pos + bandwidth]
                    - r_ref_means[seq_pos]) / r_ref_sds[seq_pos])

    fwd_pass, fwd_pass_move = c_banded_forward_pass(
        shifted_z_scores, band_event_starts, rsqgl_params.skip_pen,
        rsqgl_params.stay_pen)

    # perform traceback
    top_max_pos = np.argmax(fwd_pass[-1,:])
    if _DEBUG_FIT:
        _debug_fit(fwd_pass_move, band_event_starts, top_max_pos,
                   shifted_z_scores, reg_id, fwd_pass[-1,top_max_pos],
                   bandwidth, event_means, r_ref_means, static_bw=True)
    if _DEBUG_DP_ENDS:
        _debug_plot_dp(shifted_z_scores, fwd_pass, band_event_starts,
                       fwd_pass_move, top_max_pos, reg_id, short=True)

    read_tb = th.banded_traceback(fwd_pass_move, band_event_starts, top_max_pos)

    return read_tb


#######################################################
########## Adaptive Band Dynamic Programming ##########
#######################################################

def _get_masked_start_fwd_pass(
        event_means, r_ref_means, r_ref_sds, mapped_start_offset,
        rsqgl_params, events_per_base, mask_fill_z_score=MASK_FILL_Z_SCORE,
        mask_bases=MASK_BASES):
    if event_means.shape[0] - mapped_start_offset < rsqgl_params.bandwidth:
        raise th.TomboError(
            'Read sequence to signal matching starts too far into events for ' +
            'full adaptive assignment')
    # if max_half_z_score is none set it to valid float for cython
    # z-score computation
    if rsqgl_params.max_half_z_score is None:
        do_winsorize_z = False
        rsqgl_params = rsqgl_params._replace(max_half_z_score = 0.0)
    else:
        do_winsorize_z = True

    half_bandwidth = rsqgl_params.bandwidth // 2

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
        band_event_starts[mask_bases - 1] + rsqgl_params.bandwidth,
        mask_bases).astype(np.int64)
    def get_start_mask_z_score(seq_pos, event_pos):
        start_mask_len = max(mapped_start_offset - event_pos, 0)
        end_mask_len = (
            0 if seq_pos >= mask_bases else
            rsqgl_params.bandwidth - (mask_start_pos[seq_pos] - event_pos))
        # if the end mask does not clip back to the end of the events table
        # then extend the end mask to do so
        if (event_pos + rsqgl_params.bandwidth - end_mask_len >
            event_means.shape[0]):
            end_mask_len = (event_pos + rsqgl_params.bandwidth -
                            event_means.shape[0])
        event_vals = event_means[
            event_pos + start_mask_len:
            event_pos + rsqgl_params.bandwidth - end_mask_len]
        b_z_scores = c_base_z_scores(
            event_vals, r_ref_means[seq_pos], r_ref_sds[seq_pos],
            do_winsorize_z=do_winsorize_z,
            max_half_z_score=rsqgl_params.max_half_z_score)
        masked_z_scores = np.concatenate([
            [mask_fill_z_score - rsqgl_params.z_shift] *
            start_mask_len, b_z_scores,
            [mask_fill_z_score - rsqgl_params.z_shift] * end_mask_len])
        # This should have been handled above by checking the end_clip_len,
        # but raise an error here in case
        if masked_z_scores.shape[0] != rsqgl_params.bandwidth:
            raise th.TomboError('Masked z-score contains too few events.')
        return masked_z_scores
    shifted_z_scores = np.empty((band_event_starts.shape[0],
                                 rsqgl_params.bandwidth))
    for seq_pos, event_pos in enumerate(band_event_starts):
        shifted_z_scores[seq_pos,:] = get_start_mask_z_score(seq_pos, event_pos)
    shifted_z_scores += rsqgl_params.z_shift
    fwd_pass, fwd_pass_move = c_banded_forward_pass(
        shifted_z_scores, band_event_starts, rsqgl_params.skip_pen,
        rsqgl_params.stay_pen)

    return fwd_pass, fwd_pass_move, band_event_starts, shifted_z_scores

def find_seq_start_in_events(
        event_means, r_ref_means, r_ref_sds, rsqgl_params,
        num_bases, num_events, seq_samp_type=None, reg_id=None):
    """Identify most probably start of expected levels within observed events

    Args:
        event_means (`np.array::np.float64`): event normalized raw signal means
        r_ref_means (`np.array::np.float64`): expected base signal levels
        r_ref_sds (`np.array::np.float64`): expected base level SDs
        rsqgl_params (:class:`tombo.tombo_helper.resquiggleParams`): parameters for the re-squiggle algorithm
        num_bases (int): number of bases to process
        num_events (int): number of events to process
        reg_id (str): debug

    Returns:
        1) event position (0-based) corresponding to the start of expected signal levels
        2) mean events per base identified from start of the read
    """
    if event_means.shape[0] < num_events + num_bases:
        raise th.TomboError('Read too short for start/end discovery')
    if r_ref_means.shape[0] < num_bases:
        raise th.TomboError(
            'Genomic mapping too short for start/end discovery')

    # banded z-scores (moving up one event per base for start/end discovery
    start_z_scores = np.empty((num_bases, num_events))
    for seq_event_pos in range(num_bases):
        if rsqgl_params.max_half_z_score is None:
            start_z_scores[seq_event_pos,:] = rsqgl_params.z_shift - np.abs(
                event_means[seq_event_pos:seq_event_pos + num_events]
                - r_ref_means[seq_event_pos]) / r_ref_sds[seq_event_pos]
        else:
            start_z_scores[seq_event_pos,:] = rsqgl_params.z_shift - np.minimum(
                rsqgl_params.max_half_z_score, np.abs(
                    event_means[seq_event_pos:seq_event_pos + num_events]
                    - r_ref_means[seq_event_pos]) / r_ref_sds[seq_event_pos])
    start_band_event_starts = np.arange(num_bases, dtype=np.int64)

    start_fwd_pass, start_fwd_pass_move = c_banded_forward_pass(
        start_z_scores, start_band_event_starts, rsqgl_params.skip_pen,
        rsqgl_params.stay_pen)

    # find max along the top and right edges to start traceback
    top_max_pos = np.argmax(start_fwd_pass[-1,:])
    if _DEBUG_DP_START:
        _debug_plot_dp(
            start_z_scores, start_fwd_pass, start_band_event_starts,
            start_fwd_pass_move, top_max_pos, reg_id=reg_id, short=True)
    if _DEBUG_START_BANDWIDTH:
        _debug_fit(
            start_fwd_pass_move, start_band_event_starts, top_max_pos,
            start_z_scores, reg_id, start_fwd_pass[-1, top_max_pos],
            num_events, event_means, r_ref_means)

    # perform traceback
    start_tb = th.banded_traceback(
        start_fwd_pass_move, start_band_event_starts, top_max_pos)
    if (seq_samp_type is not None and
        ts.score_valid_bases(start_tb, event_means, r_ref_means, r_ref_sds) >
        SIG_MATCH_THRESH[seq_samp_type.name]):
        raise th.TomboError(
            'Poor raw to expected signal matching in beginning of read.')

    # compute the average events per base to use for the start forward pass
    events_per_base = (start_tb[-1] - start_tb[0]) / len(start_tb)
    start_loc = start_tb[0]

    return start_loc, events_per_base

def _trim_traceback(read_tb, events_len):
    start_trim_i = 0
    while read_tb[start_trim_i] < 0:
        read_tb[start_trim_i] = 0
        start_trim_i += 1
    end_trim_i = 1
    while read_tb[-end_trim_i] > events_len:
        read_tb[-end_trim_i] = events_len
        end_trim_i += 1

    return read_tb

def find_seq_start_from_clip_basecalls(
        event_means, rsqgl_params, start_clip_bases, genome_seq, std_ref,
        num_genome_bases, reg_id=None):
    """Perform dynamic programming over clipped basecalls and genome sequence to identify the start of the genomic mapping

    Args:
        event_means (`np.array::np.float64`): normalized raw signal event means
        rsqgl_params (:class:`tombo.tombo_helper.resquiggleParams`): parameters for the re-squiggle algorithm
        start_clip_bases (str): read bases clipped from before ``genome_seq``
        genome_seq (str): genomic mapping sequence (inlcuding extra bases based on k-mer size)
        std_ref (:class:`tombo.tombo_stats.TomboModel`): canonical model
        num_genome_bases (int): genome sequence length used to identify start position (needed for traceback)
        reg_id (str): debug

    Returns:
        1) event position (0-based) corresponding to the start of expected signal levels from genome_seq
        2) mean events per base identified from start of the read
    """
    dnstrm_bases = std_ref.kmer_width - std_ref.central_pos - 1
    start_genome_seq = genome_seq[
        std_ref.central_pos:num_genome_bases + dnstrm_bases]
    start_seq = start_clip_bases + start_genome_seq
    r_ref_means, r_ref_sds = std_ref.get_exp_levels_from_seq(start_seq)
    seq_len = r_ref_means.shape[0]

    # now find full sequence to events path using a smaller bandwidth
    (start_fwd_pass, start_fwd_pass_move,
     start_event_starts, start_z_scores) = _get_masked_start_fwd_pass(
         event_means, r_ref_means, r_ref_sds, 0,
         rsqgl_params, (rsqgl_params.bandwidth // 2) / float(MASK_BASES))
    start_seq_len = start_event_starts.shape[0]
    fwd_pass = np.empty((seq_len+1, rsqgl_params.bandwidth), dtype=np.float64)
    fwd_pass[:start_seq_len+1] = start_fwd_pass
    fwd_pass_move = np.empty((seq_len+1, rsqgl_params.bandwidth), dtype=np.int64)
    fwd_pass_move[:start_seq_len+1] = start_fwd_pass_move
    band_event_starts = np.empty((seq_len,), dtype=np.int64)
    band_event_starts[:start_seq_len] = start_event_starts

    if rsqgl_params.max_half_z_score is None:
        do_winsorize_z = False
        rsqgl_params = rsqgl_params._replace(max_half_z_score = 0.0)
    else:
        do_winsorize_z = True

    if _DEBUG_CLIP_START or _DEBUG_START_CLIP_FIT:
        # save z-scores for debug plotting
        rest_z_scores = th.adaptive_banded_forward_pass(
            fwd_pass, fwd_pass_move, band_event_starts, event_means,
            r_ref_means, r_ref_sds,
            z_shift=rsqgl_params.z_shift,
            skip_pen=rsqgl_params.skip_pen, stay_pen=rsqgl_params.stay_pen,
            start_seq_pos=start_seq_len, mask_fill_z_score=MASK_FILL_Z_SCORE,
            do_winsorize_z=do_winsorize_z,
            max_half_z_score=rsqgl_params.max_half_z_score, return_z_scores=True)
        shifted_z_scores = np.empty((
            seq_len, rsqgl_params.bandwidth), dtype=np.float64)
        shifted_z_scores[:start_seq_len] = start_z_scores
        shifted_z_scores[start_seq_len:] = rest_z_scores
    else:
        th.adaptive_banded_forward_pass(
            fwd_pass, fwd_pass_move, band_event_starts, event_means,
            r_ref_means, r_ref_sds,
            z_shift=rsqgl_params.z_shift,
            skip_pen=rsqgl_params.skip_pen, stay_pen=rsqgl_params.stay_pen,
            start_seq_pos=start_seq_len, mask_fill_z_score=MASK_FILL_Z_SCORE,
            do_winsorize_z=do_winsorize_z,
            max_half_z_score=rsqgl_params.max_half_z_score)

    top_max_pos = np.argmax(fwd_pass[-1,:])
    if _DEBUG_CLIP_START:
        _debug_plot_dp(shifted_z_scores, fwd_pass, band_event_starts,
                       fwd_pass_move, top_max_pos, short=True, reg_id=reg_id)

    if _DEBUG_START_CLIP_FIT:
        _debug_fit(fwd_pass_move, band_event_starts, top_max_pos,
                   shifted_z_scores, reg_id, fwd_pass[-1,top_max_pos],
                   rsqgl_params.bandwidth, event_means, r_ref_means)

    read_tb = th.banded_traceback(
        fwd_pass_move, band_event_starts, top_max_pos,
        rsqgl_params.band_bound_thresh)
    # trim invalid traceback positions
    read_tb = _trim_traceback(read_tb, events_len=event_means.shape[0])

    assert len(start_clip_bases) >= std_ref.central_pos, (
        'Invalid start clip base processing.')
    start_loc = read_tb[len(start_clip_bases) - std_ref.central_pos]
    events_per_base = (read_tb[-1] - start_loc) / (
        len(start_genome_seq) - dnstrm_bases)

    return start_loc, events_per_base

def get_rel_raw_coords(valid_cpts, seq_events):
    """get raw coordinates relative to the start of the assigned signal
    """
    seq_segs = valid_cpts[seq_events]
    read_start_rel_to_raw = seq_segs[0]
    seq_segs = seq_segs - read_start_rel_to_raw
    return seq_segs, read_start_rel_to_raw

def find_adaptive_base_assignment(
        valid_cpts, event_means, rsqgl_params, std_ref, genome_seq,
        start_clip_bases=None, start_clip_params=START_CLIP_PARAMS,
        seq_samp_type=th.seqSampleType(DNA_SAMP_TYPE, False), reg_id=None):
    """Align expected (from genome sequence) signal levels to observed using a dynamic programming approach with adaptive bandwidth start positions

    Args:
        valid_cpts (`np.array::np.int64`): raw signal base start locations
        event_means (`np.array::np.float64`): event normalized raw signal means
        rsqgl_params (:class:`tombo.tombo_helper.resquiggleParams`): parameters for the re-squiggle algorithm
        std_ref (:class:`tombo.tombo_stats.TomboModel`): canonical model
        genome_seq (str): genome sequence (from mapping)
        start_clip_bases (str): mapping read clipped bases
        start_clip_params (:class:`tombo.tombo_helper.startClipParams`): start clip basecall params
        seq_samp_type (:class:`tombo.tombo_helper.seqSampleType`): sequencing sample type (default: DNA)
        reg_id (str): debug

    Returns:
        :class:`tombo.tombo_helper.dpResults`
    """
    def get_short_read_results(r_ref_means, r_ref_sds, genome_seq):
        seq_events = find_static_base_assignment(
            event_means, r_ref_means, r_ref_sds, rsqgl_params, reg_id=reg_id)
        seq_segs, read_start_rel_to_raw = get_rel_raw_coords(
            valid_cpts, seq_events)
        return th.dpResults(
            read_start_rel_to_raw=read_start_rel_to_raw, segs=seq_segs,
            ref_means=r_ref_means, ref_sds=r_ref_sds, genome_seq=genome_seq)

    def run_fwd_pass():
        # find full sequence to events path using a smaller bandwidth
        (start_fwd_pass, start_fwd_pass_move,
         start_event_starts, start_z_scores) = _get_masked_start_fwd_pass(
             event_means[events_start_clip:], r_ref_means, r_ref_sds,
             mapped_start_offset, rsqgl_params, events_per_base)
        start_seq_len = start_event_starts.shape[0]
        fwd_pass = np.empty((seq_len+1, rsqgl_params.bandwidth), dtype=np.float64)
        fwd_pass[:start_seq_len+1] = start_fwd_pass
        fwd_pass_move = np.empty((seq_len+1, rsqgl_params.bandwidth),
                                 dtype=np.int64)
        fwd_pass_move[:start_seq_len+1] = start_fwd_pass_move
        band_event_starts = np.empty((seq_len,), dtype=np.int64)
        band_event_starts[:start_seq_len] = start_event_starts

        # if max_half_z_score is none set it to valid float for cython
        # z-score computation
        if rsqgl_params.max_half_z_score is None:
            do_winsorize_z = False
            wa_rsqgl_params = rsqgl_params._replace(max_half_z_score = 0.0)
        else:
            do_winsorize_z = True
            wa_rsqgl_params = rsqgl_params

        shifted_z_scores = None
        if _DEBUG_FIT or _DEBUG_BANDWIDTH or _DEBUG_DP_ENDS:
            # save z-scores for debug plotting
            rest_z_scores = th.adaptive_banded_forward_pass(
                fwd_pass, fwd_pass_move, band_event_starts,
                event_means[events_start_clip:], r_ref_means, r_ref_sds,
                wa_rsqgl_params.z_shift, wa_rsqgl_params.skip_pen,
                wa_rsqgl_params.stay_pen, start_seq_len, MASK_FILL_Z_SCORE,
                do_winsorize_z, wa_rsqgl_params.max_half_z_score,
                return_z_scores=True)
            shifted_z_scores = np.empty((
                seq_len, wa_rsqgl_params.bandwidth), dtype=np.float64)
            shifted_z_scores[:start_seq_len] = start_z_scores
            shifted_z_scores[start_seq_len:] = rest_z_scores
        else:
            th.adaptive_banded_forward_pass(
                fwd_pass, fwd_pass_move, band_event_starts,
                event_means[events_start_clip:],
                r_ref_means, r_ref_sds, wa_rsqgl_params.z_shift,
                wa_rsqgl_params.skip_pen, wa_rsqgl_params.stay_pen, start_seq_len,
                MASK_FILL_Z_SCORE, do_winsorize_z,
                wa_rsqgl_params.max_half_z_score)

        return fwd_pass, fwd_pass_move, band_event_starts, shifted_z_scores

    def plot_debug(shifted_z_scores):
        if _DEBUG_FIT or _DEBUG_BANDWIDTH:
            _debug_fit(
                fwd_pass_move, band_event_starts, top_max_pos, shifted_z_scores,
                reg_id, fwd_pass[-1, top_max_pos], rsqgl_params.bandwidth,
                event_means[events_start_clip:], r_ref_means)
        if _DEBUG_DP_ENDS:
            _debug_plot_dp(shifted_z_scores, fwd_pass, band_event_starts,
                           fwd_pass_move, top_max_pos, reg_id)
        return


    # if start clip bases are provided, run "cliped bases" start
    # identification algorithm
    if (start_clip_bases is not None and
        len(genome_seq) > start_clip_params.num_genome_bases):
        if len(start_clip_bases) < std_ref.central_pos:
            mapped_start = len(start_clip_bases) * 2
            events_per_base = 2
        else:
            clip_params = rsqgl_params._replace(
                bandwidth=start_clip_params.bandwidth)
            mapped_start, events_per_base = find_seq_start_from_clip_basecalls(
                event_means, clip_params, start_clip_bases, genome_seq, std_ref,
                start_clip_params.num_genome_bases, reg_id=reg_id)

    dnstrm_bases = std_ref.kmer_width - std_ref.central_pos - 1
    r_ref_means, r_ref_sds = std_ref.get_exp_levels_from_seq(genome_seq)
    # trim genome seq to match model-able positions
    genome_seq = genome_seq[std_ref.central_pos:-dnstrm_bases]
    seq_len = len(genome_seq)
    if seq_len != r_ref_means.shape[0]:
        raise th.TomboError('Discordant reference and seqeunce lengths.')

    # if read was too short for start clip map start discovery, but need r_ref*
    if (start_clip_bases is not None and
        seq_len <= start_clip_params.num_genome_bases):
        return get_short_read_results(r_ref_means, r_ref_sds, genome_seq)

    if start_clip_bases is None:
        # for short reads, just search the whole read with an appropriate
        # bandwidth
        if (event_means.shape[0] < rsqgl_params.start_bw +
            rsqgl_params.start_n_bases or
            seq_len < rsqgl_params.start_n_bases):
            return get_short_read_results(r_ref_means, r_ref_sds, genome_seq)
        try:
            # identify the start of genomic sequence within raw signal
            mapped_start, events_per_base = find_seq_start_in_events(
                event_means, r_ref_means, r_ref_sds, rsqgl_params,
                rsqgl_params.start_n_bases, rsqgl_params.start_bw,
                seq_samp_type, reg_id=reg_id)
        except th.TomboError:
            if (event_means.shape[0] < rsqgl_params.start_save_bw +
                rsqgl_params.start_n_bases):
                return get_short_read_results(r_ref_means, r_ref_sds, genome_seq)
            # if smaller (and faster) bandwidth did not find a sufficiently
            # scoring raw signal mapping, try again with larger save_bandwidth
            # and don't check score (by not passing seq_samp_type)
            mapped_start, events_per_base = find_seq_start_in_events(
                event_means, r_ref_means, r_ref_sds, rsqgl_params,
                rsqgl_params.start_n_bases, rsqgl_params.start_save_bw,
                reg_id=reg_id)

    if events_per_base == 0:
        raise th.TomboError(
            'Very poor signal quality. Read likely includes open pore.')

    # get number of events to clip and how far into the events the
    # discovered start is located
    half_bandwidth = rsqgl_params.bandwidth // 2
    if mapped_start < half_bandwidth:
        events_start_clip = 0
        mapped_start_offset = mapped_start
    else:
        events_start_clip = mapped_start - half_bandwidth
        mapped_start_offset = half_bandwidth

    # process long enough reads that start too far into read for normal
    # adaptive processing just as with short reads
    if (int((half_bandwidth + 1) / events_per_base) >= r_ref_means.shape[0] or
        (event_means.shape[0] - mapped_start_offset -
         events_start_clip < rsqgl_params.bandwidth)):
        return get_short_read_results(r_ref_means, r_ref_sds, genome_seq)

    fwd_pass, fwd_pass_move, band_event_starts, shifted_z_scores = run_fwd_pass()

    # find position of last base at the maximal score
    top_max_pos = np.argmax(fwd_pass[-1,:])

    # plot debugging if requested
    plot_debug(shifted_z_scores)

    read_tb = th.banded_traceback(
        fwd_pass_move, band_event_starts, top_max_pos,
        rsqgl_params.band_bound_thresh)
    # trim invalid traceback positions
    read_tb = _trim_traceback(
        read_tb, events_len=event_means.shape[0] - events_start_clip)

    # get segment positions within raw signal vector
    seq_segs, read_start_rel_to_raw = get_rel_raw_coords(
        valid_cpts[events_start_clip:], read_tb)

    return th.dpResults(
        read_start_rel_to_raw=read_start_rel_to_raw, segs=seq_segs,
        ref_means=r_ref_means, ref_sds=r_ref_sds, genome_seq=genome_seq)


######################################
########## Re-squiggle Read ##########
######################################

def segment_signal(
        map_res, num_events, rsqgl_params, outlier_thresh=None, const_scale=None):
    """Normalize and segment raw signal as defined by `rsqgl_params` into `num_events`.

    Args:
        map_res (:class:`tombo.tombo_helper.resquiggleResults`): containing mapping results only (attributes after ``read_start_rel_to_raw`` will all be ``None``)
        num_events (int): number of events to process
        rsqgl_params (:class:`tombo.tombo_helper.resquiggleParams`): parameters for the re-squiggle algorithm
        outlier_thresh (float): windsorize signal greater than this value (optional)

    Returns:
        1) identified event positions (0-based)
        2) normalized raw signal
        3) scale values (:class:`tombo.tombo_helper.scaleValues`)
    """
    if rsqgl_params.use_t_test_seg:
        # RNA bases show consistent variable spread so use t-test segmentation
        valid_cpts = th.valid_cpts_w_cap_t_test(
            map_res.raw_signal.astype(np.float64), rsqgl_params.min_obs_per_base,
            rsqgl_params.running_stat_width, num_events)

        # remove cpts within stall locations
        if map_res.stall_ints is not None:
            valid_cpts = ts.remove_stall_cpts(map_res.stall_ints, valid_cpts)

        if map_res.scale_values is not None:
            norm_signal, new_scale_values = ts.normalize_raw_signal(
                map_res.raw_signal, scale_values=map_res.scale_values)
        elif const_scale is not None:
            norm_signal, new_scale_values = ts.normalize_raw_signal(
                map_res.raw_signal, norm_type='median_const_scale',
                outlier_thresh=outlier_thresh, const_scale=const_scale)
        else:
            if USE_RNA_EVENT_SCALE:
                scale_values = ts.get_scale_values_from_events(
                    map_res.raw_signal, valid_cpts, outlier_thresh,
                    num_events=RNA_SCALE_NUM_EVENTS,
                    max_frac_events=RNA_SCALE_MAX_FRAC_EVENTS)
            else:
                scale_values = None
            norm_signal, new_scale_values = ts.normalize_raw_signal(
                map_res.raw_signal, scale_values=scale_values)
    else:
        # normalize signal
        if map_res.scale_values is not None:
            norm_signal, new_scale_values = ts.normalize_raw_signal(
                map_res.raw_signal, scale_values=map_res.scale_values)
        elif const_scale is not None:
            norm_signal, new_scale_values = ts.normalize_raw_signal(
                map_res.raw_signal, norm_type='median_const_scale',
                outlier_thresh=outlier_thresh, const_scale=const_scale)
        else:
            norm_signal, new_scale_values = ts.normalize_raw_signal(
                map_res.raw_signal, norm_type='median',
                outlier_thresh=outlier_thresh)

        valid_cpts = th.valid_cpts_w_cap(
            norm_signal, rsqgl_params.min_obs_per_base,
            rsqgl_params.running_stat_width, num_events)
        # remove cpts within stall locations
        if map_res.stall_ints is not None:
            valid_cpts = ts.remove_stall_cpts(map_res.stall_ints, valid_cpts)

    return valid_cpts, norm_signal, new_scale_values

def resquiggle_read(
        map_res, std_ref, rsqgl_params, outlier_thresh=None,
        all_raw_signal=None, max_raw_cpts=MAX_RAW_CPTS,
        min_event_to_seq_ratio=MIN_EVENT_TO_SEQ_RATIO, const_scale=None,
        skip_seq_scaling=False,
        seq_samp_type=th.seqSampleType(DNA_SAMP_TYPE, False)):
    """Identify raw signal to genome sequence assignment, using adaptive banded dynamic programming

    Args:
        map_res (:class:`tombo.tombo_helper.resquiggleResults`): mapping results
        std_ref (:class:`tombo.tombo_stats.TomboModel`): canonical base model
        rsqgl_params (:class:`tombo.tombo_helper.resquiggleParams`): parameters for the re-squiggle algorithm
        outlier_thresh (float): windsorize signal greater than this value (optional)
        all_raw_signal (`np.array::np.int64`): raw data acquisition (DAC) current signal values (optional; default use value in map_res)
        max_raw_cpts (int): read will fail if more than `max_raw_cpts` must be found to produce a valid re-squiggle results (optional)
        min_event_to_seq_ratio (float): minimum event to sequence ratio (optional)
        const_scale (float): constant scale value (optional; may be deprecated)
        skip_seq_scaling (bool): skip sequence-based scaling step
        seq_samp_type (:class:`tombo.tombo_helper.seqSampleType`): sequencing sample type (default: DNA)

    Returns:
        :class:`tombo.tombo_helper.resquiggleResults` containing raw signal to genome sequence alignment (note that ``raw_signal`` now contains trimmed and normalized raw signal)
    """
    if all_raw_signal is not None:
        map_res = map_res._replace(raw_signal = all_raw_signal)
    if map_res.raw_signal is None:
        raise th.TomboError(
            'Must have raw signal in order to complete re-squiggle algorithm')

    # compute number of events to find
    # ensure at least a minimal number of events per mapped sequence are found
    num_mapped_bases = len(map_res.genome_seq) - std_ref.kmer_width + 1
    num_events = ts.compute_num_events(
        map_res.raw_signal.shape[0], num_mapped_bases,
        rsqgl_params.mean_obs_per_event, min_event_to_seq_ratio)
    # ensure that there isn't *far* too much signal for the mapped sequence
    # i.e. one adaptive bandwidth per base is too much to find a good mapping
    if num_events / rsqgl_params.bandwidth > num_mapped_bases:
        raise th.TomboError('Too much raw signal for mapped sequence')

    valid_cpts, norm_signal, new_scale_values = segment_signal(
        map_res, num_events, rsqgl_params, outlier_thresh, const_scale)
    event_means = ts.compute_base_means(norm_signal, valid_cpts)

    dp_res = find_adaptive_base_assignment(
        valid_cpts, event_means, rsqgl_params, std_ref, map_res.genome_seq,
        start_clip_bases=map_res.start_clip_bases,
        seq_samp_type=seq_samp_type, reg_id=map_res.align_info.ID)
    # clip raw signal to only part mapping to genome seq
    norm_signal = norm_signal[dp_res.read_start_rel_to_raw:
                              dp_res.read_start_rel_to_raw + dp_res.segs[-1]]

    # identify all stretches of genomic deletions within del_fix_window
    # to be fixed.
    segs = resolve_skipped_bases_with_raw(
        dp_res, norm_signal, rsqgl_params, max_raw_cpts)

    if skip_seq_scaling:
        norm_params_changed = False
    else:
        (shift, scale, shift_corr_factor,
         scale_corr_factor) = ts.calc_kmer_fitted_shift_scale(
             new_scale_values.shift, new_scale_values.scale,
             ts.compute_base_means(norm_signal, segs), dp_res.ref_means,
             method='theil_sen')
        new_scale_values = new_scale_values._replace(
            shift=shift, scale=scale, outlier_thresh=outlier_thresh)
        # re-normalize signal with new fitted parameters
        norm_signal = (norm_signal - shift_corr_factor) / scale_corr_factor
        # determine if normalization parameters changed enough to warrant
        # re-squiggling again
        norm_params_changed = (
            np.abs(shift_corr_factor) > SHIFT_CHANGE_THRESH or
            np.abs(scale_corr_factor - 1) > SCALE_CHANGE_THRESH)

    sig_match_score = ts.get_read_seg_score(
        ts.compute_base_means(norm_signal, segs),
        dp_res.ref_means, dp_res.ref_sds)
    if segs.shape[0] != len(dp_res.genome_seq) + 1:
        raise th.TomboError('Aligned sequence does not match number ' +
                            'of segments produced')

    # Output for testing/visualization of re-squiggle
    if _DEBUG_PARAMS:
        _write_params_debug(
            norm_signal, segs, dp_res.ref_means, dp_res.ref_sds,
            rsqgl_params, map_res.align_info.ID)

    return map_res._replace(
        read_start_rel_to_raw=dp_res.read_start_rel_to_raw, segs=segs,
        genome_seq=dp_res.genome_seq, raw_signal=norm_signal,
        scale_values=new_scale_values, sig_match_score=sig_match_score,
        norm_params_changed=norm_params_changed)


#######################################
########## Genomic Alignment ##########
#######################################

def get_read_seq(
        fast5_data, bc_grp='Basecall_1D_000', bc_subgrp='BaseCalled_template',
        seq_samp_type=th.seqSampleType(DNA_SAMP_TYPE, False), q_score_thresh=0):
    """Extract read sequence from the Fastq slot providing useful error messages

    Args:

        fast5_data (:class:`tombo.tombo_helper.readData`): read information
        bc_grp (str): group location containing read information (optional; default: 'Basecall_1D_000')
        bc_subgrp (str): sub-group location containing read information (optional; default: 'BaseCalled_template')
        seq_samp_type (:class:`tombo.tombo_helper.seqSampleType`): sequencing sample type (default: DNA)
        q_score_thresh (float): basecalling mean q-score threshold (optional; default: 0/no filtering)

    Returns:
        :class:`tombo.tombo_helper.sequenceData`
    """
    try:
        fastq_raw_value = fast5_data[
            '/Analyses/' + bc_grp + '/' + bc_subgrp + '/Fastq'][()]
    except KeyError:
        raise th.TomboError('Fastq slot not present in --basecall-group')

    # depending on how fastq data was stored it may already be encoded
    # as unicode, so this would fail.
    try:
        fastq_raw_value = fastq_raw_value.decode()
    except (TypeError, AttributeError):
        pass

    s_fastq = fastq_raw_value.split('\n')
    read_seq, read_q = s_fastq[1], s_fastq[3]

    # compute read q-score
    mean_q_score = th.get_mean_q_score(read_q)
    if q_score_thresh is not None and mean_q_score < q_score_thresh:
        raise th.TomboError('Read filtered by q-score.')

    read_data = th.get_raw_read_slot(fast5_data)

    # looks like read_id attribute has been removed in some files and attribute
    # is not really necessary for tombo
    try:
        read_id = read_data.attrs.get('read_id')
    except KeyError:
        try:
            read_id = unicode(read_data.attrs.get('read_num'))
        except KeyError:
            read_id = unicode(np.random.randint(1000000000))

    # only really here for the API
    if seq_samp_type is None:
        seq_samp_type = th.get_seq_sample_type(fast5_data)
    if seq_samp_type.name == RNA_SAMP_TYPE:
        read_seq = th.rev_transcribe(read_seq)

    return th.sequenceData(seq=read_seq, id=read_id, mean_q_score=mean_q_score)

def map_read(
        fast5_data, aligner, std_ref,
        seq_samp_type=th.seqSampleType(DNA_SAMP_TYPE, False),
        bc_grp='Basecall_1D_000', bc_subgrp='BaseCalled_template',
        map_thr_buf=None, q_score_thresh=0):
    """Extract read sequence from the Fastq slot providing useful error messages

    Args:
        fast5_data (:class:`tombo.tombo_helper.readData`): read information
        aligner (mappy.Aligner): aligner object
        std_ref (:class:`tombo.tombo_stats.TomboModel`): canonical model (in order to extract extended genomic sequence)
        seq_samp_type (:class:`tombo.tombo_helper.seqSampleType`): sequencing sample type (default: DNA)
        bc_grp (str): group location containing read information (optional; default: 'Basecall_1D_000')
        bc_subgrp (str): sub-group location containing read information (optional; default: 'BaseCalled_template')
        map_thr_buf (mappy.ThreadBuffer): mappy thread buffer object (optional; default: None)
        q_score_thresh (float): basecalling mean q-score threshold (optional; default: 0/no filtering)

    Returns:
        :class:`tombo.tombo_helper.resquiggleResults` containing valid mapping values (signal to sequence assignment attributes will be ``None``)
    """
    seq_data = get_read_seq(
        fast5_data, bc_grp, bc_subgrp, seq_samp_type, q_score_thresh)
    try:
        alignment = next(aligner.map(str(seq_data.seq), buf=map_thr_buf))
    except StopIteration:
        raise th.TomboError('Alignment not produced')

    chrm = alignment.ctg
    # subtract one to put into 0-based index
    ref_start = alignment.r_st
    ref_end = alignment.r_en
    strand = '+' if alignment.strand == 1 else '-'
    num_match = alignment.mlen
    num_ins, num_del, num_aligned = 0, 0, 0
    for op_len, op in alignment.cigar:
        if op == 1: num_ins += op_len
        elif op in (2,3): num_del += op_len
        elif op in (0,7,8): num_aligned += op_len
        elif op == 6: pass
        else:
            # soft and hard clipping are not reported in the
            # mappy cigar
            raise th.TomboError('Invalid cigar operation')

    # store number of clipped bases relative to read sequence
    if strand == '+':
        num_start_clipped_bases = alignment.q_st
        num_end_clipped_bases = len(seq_data.seq) - alignment.q_en
    else:
        num_start_clipped_bases = len(seq_data.seq) - alignment.q_en
        num_end_clipped_bases = alignment.q_st

    align_info = th.alignInfo(
        seq_data.id.decode(), bc_subgrp, num_start_clipped_bases,
        num_end_clipped_bases, num_ins, num_del, num_match,
        num_aligned - num_match)

    # extract genome sequence from mappy aligner
    # expand sequence to get model levels for all sites (need to handle new
    # sequence coordinates downstream)
    dnstrm_bases = std_ref.kmer_width - std_ref.central_pos - 1
    if ((seq_samp_type.name == RNA_SAMP_TYPE and strand == '+') or
        (seq_samp_type.name == DNA_SAMP_TYPE and strand == '-' and
         USE_START_CLIP_BASES) or
        (seq_samp_type.name == DNA_SAMP_TYPE and strand == '+' and
         not USE_START_CLIP_BASES)):
        if ref_start < std_ref.central_pos:
            ref_start = std_ref.central_pos
        ref_seq_start = ref_start - std_ref.central_pos
        ref_seq_end = ref_end + dnstrm_bases
    else:
        if ref_start < dnstrm_bases:
            ref_start = dnstrm_bases
        ref_seq_start = ref_start - dnstrm_bases
        ref_seq_end = ref_end + std_ref.central_pos
    genome_seq = aligner.seq(chrm, ref_seq_start, ref_seq_end)
    if genome_seq is None or genome_seq == '':
        raise th.TomboReads('Invalid mapping location')

    if sys.version_info[0] < 3:
        genome_seq = genome_seq.decode()
    if strand == '-':
        genome_seq = th.rev_comp(genome_seq)
    # discordant mapping to sequence extraction is due to reads mapping up to
    # the end of a seqeunce record (and don't need to carry around record lens),
    # so don't error on these discordant lengths here
    #if len(genome_seq) != ref_end - ref_start + std_ref.kmer_width - 1:
    #    raise th.TomboError('Discordant mapped position and sequence')
    genome_loc = th.genomeLocation(ref_start, strand, chrm)

    # store sequence at the end of the read without an adapter
    # for simpler read start identification (start of RNA genomic sequence
    # end of DNA genomic sequence)
    start_clip_bases = None
    if USE_START_CLIP_BASES:
        start_clip_bases = seq_data.seq[alignment.q_en:][::-1]

    return th.resquiggleResults(
        align_info=align_info, genome_loc=genome_loc, genome_seq=genome_seq,
        mean_q_score=seq_data.mean_q_score, start_clip_bases=start_clip_bases)

def _io_and_map_read(
        fast5_data, failed_reads_q, bc_subgrps, bc_grp, corr_grp, aligner,
        seq_samp_type, map_thr_buf, fast5_fn, num_processed, map_conn,
        outlier_thresh, compute_sd, obs_filter, index_q, q_score_thresh,
        sig_match_thresh, std_ref):
    try:
        # extract channel and raw data for this read
        channel_info = th.get_channel_info(fast5_data)
    except th.TomboError:
        # channel info is not needed currently, so just pass
        channel_info = None
    all_raw_signal = th.get_raw_read_slot(fast5_data)['Signal'][:]

    for bc_subgrp in bc_subgrps:
        try:
            map_res = map_read(
                fast5_data, aligner, std_ref, seq_samp_type,
                bc_grp, bc_subgrp, map_thr_buf, q_score_thresh)
            if th.invalid_seq(map_res.genome_seq):
                raise th.TomboError(
                    'Reference mapping contains non-canonical bases ' +
                    '(transcriptome reference cannot contain U bases)')
            map_res = map_res._replace(
                raw_signal=all_raw_signal, channel_info=channel_info)

            # send mapping data to _resquiggle_worker process
            map_conn.send([map_res, fast5_fn])
            # wait until re-squiggle returns
            read_failed, rsqgl_res = map_conn.recv()

            if read_failed:
                failed_reads_q.put((
                    rsqgl_res[0], bc_subgrp + ':::' + rsqgl_res[1],
                    rsqgl_res[2]))
                continue

            if not _DRY_RUN:
                # write re-squiggle event assignment to the read FAST5 file
                th.write_new_fast5_group(
                    fast5_data, corr_grp, rsqgl_res, 'median', compute_sd,
                    rna=seq_samp_type.rev_sig)

            if index_q is not None:
                # check that read passes reversible filters
                is_filtered = False
                if rsqgl_res.sig_match_score > sig_match_thresh:
                    failed_reads_q.put((
                        'Poor raw to expected signal matching ' +
                        '(revert with `tombo filter clear_filters`)',
                        bc_subgrp + ':::' + fast5_fn, True))
                    is_filtered = True
                elif obs_filter is not None:
                    base_lens = np.diff(rsqgl_res.segs)
                    is_filtered = any(np.percentile(base_lens, pctl) > thresh
                                      for pctl, thresh in obs_filter)
                    failed_reads_q.put((
                        'Read filtered by observation per base ' +
                        'thresholds (revert with `tombo filter clear_filters`)',
                        bc_subgrp + ':::' + fast5_fn, True))
                # prep and load data into index queue
                mapped_end = (
                    rsqgl_res.genome_loc.Start + len(rsqgl_res.segs) - 1)
                index_q.put((
                    rsqgl_res.genome_loc.Chrom,
                    rsqgl_res.genome_loc.Strand, th.readData(
                        rsqgl_res.genome_loc.Start, mapped_end, is_filtered,
                        rsqgl_res.read_start_rel_to_raw,
                        rsqgl_res.genome_loc.Strand, fast5_fn,
                        corr_grp + '/' + bc_subgrp, seq_samp_type.rev_sig,
                        rsqgl_res.sig_match_score, rsqgl_res.mean_q_score,
                        rsqgl_res.align_info.ID)))
        except th.TomboError as e:
            try:
                th.write_error_status(
                    fast5_fn, corr_grp, bc_subgrp, unicode(e))
            except:
                pass
            failed_reads_q.put((
                unicode(e), bc_subgrp + ':::' + fast5_fn, True))
        except Exception as e:
            # is_tombo_error = False
            failed_reads_q.put((
                traceback.format_exc(), bc_subgrp + ':::' + fast5_fn, False))

    return


#########################################
########## Re-squiggle Workers ##########
#########################################

def _resquiggle_worker(
        rsqgl_conns, std_ref, outlier_thresh, corr_grp, seq_samp_type,
        rsqgl_params, save_params, const_scale, skip_seq_scaling,
        max_scaling_iters):
    def run_rsqgl_iters(map_res, params, fast5_fn, all_raw_signal):
        rsqgl_res = resquiggle_read(
            map_res, std_ref, params, outlier_thresh,
            const_scale=const_scale, skip_seq_scaling=skip_seq_scaling,
            seq_samp_type=seq_samp_type)
        n_iters = 1
        while n_iters < max_scaling_iters and rsqgl_res.norm_params_changed:
            rsqgl_res = resquiggle_read(
                map_res._replace(scale_values=rsqgl_res.scale_values),
                std_ref, params, outlier_thresh, all_raw_signal=all_raw_signal,
                seq_samp_type=seq_samp_type)
            n_iters += 1
        return rsqgl_res

    def adjust_map_res(map_res):
        if seq_samp_type.name == RNA_SAMP_TYPE:
            if TRIM_RNA_ADAPTER:
                # trim DNA adapter off of RNA signal
                adapter_end = ts.trim_rna(map_res.raw_signal, rsqgl_params)
                # trim off adapter
                map_res = map_res._replace(
                    raw_signal=map_res.raw_signal[adapter_end:])

            # flip raw signal for re-squiggling
            map_res = map_res._replace(raw_signal=map_res.raw_signal[::-1])

        elif seq_samp_type.name == DNA_SAMP_TYPE and USE_START_CLIP_BASES:
            # flip raw signal, genome and start clip seqs for re-squiggling
            map_res = map_res._replace(
                raw_signal=map_res.raw_signal[::-1],
                genome_seq=map_res.genome_seq[::-1])

        if ((COLLAPSE_RNA_STALLS and seq_samp_type.name == RNA_SAMP_TYPE) or
            (COLLAPSE_DNA_STALLS and seq_samp_type.name == DNA_SAMP_TYPE)):
            map_res = map_res._replace(
                stall_ints=ts.identify_stalls(
                    map_res.raw_signal, DEFAULT_STALL_PARAMS))

        return map_res

    def adjust_rsqgl_res(rsqgl_res, all_raw_signal):
        if seq_samp_type.name == DNA_SAMP_TYPE and USE_START_CLIP_BASES:
            # flip raw signal and events back for storage in genome direction
            rev_rsrtr = (all_raw_signal.shape[0] -
                         rsqgl_res.read_start_rel_to_raw -
                         rsqgl_res.segs[-1])
            rev_segs = -1 * (rsqgl_res.segs[::-1] - rsqgl_res.segs[-1])
            rsqgl_res = rsqgl_res._replace(
                read_start_rel_to_raw=rev_rsrtr, segs=rev_segs,
                genome_seq=rsqgl_res.genome_seq[::-1],
                raw_signal=rsqgl_res.raw_signal[::-1])

        return rsqgl_res


    if _DEBUG_PLOTTING:
        _open_debug_pdf()

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

            map_res, fast5_fn = map_info

            map_res = adjust_map_res(map_res)
            # save un-normalized signal for later iterations
            all_raw_signal = map_res.raw_signal
            try:
                rsqgl_res = run_rsqgl_iters(
                    map_res, rsqgl_params, fast5_fn, all_raw_signal)
            # if the resquiggle read fails for any reason
            except:
                rsqgl_res = run_rsqgl_iters(
                    map_res, save_params, fast5_fn, all_raw_signal)
            rsqgl_res = adjust_rsqgl_res(rsqgl_res, all_raw_signal)
        except th.TomboError as e:
            rsqgl_conn.send([True, [unicode(e), fast5_fn, True]])
            continue
        except Exception as e:
            rsqgl_conn.send([True, [traceback.format_exc(), fast5_fn, False]])
            continue

        rsqgl_conn.send([False, rsqgl_res])

    if _DEBUG_PLOTTING:
        _close_debug_pdf()

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
        corr_grp, aligner, outlier_thresh, compute_sd, sig_match_thresh,
        obs_filter, seq_samp_type, overwrite, map_conn, q_score_thresh,
        std_ref):
    # increase update interval as more reads are provided
    proc_update_interval = 1
    def update_progress(num_processed, proc_update_interval):
        if num_processed % proc_update_interval == 0:
            progress_q.put(proc_update_interval)
            # increase update interval as more reads are processed
            if num_processed == 100:
                proc_update_interval = 10
            if num_processed == 1000:
                proc_update_interval = 100
        return proc_update_interval

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
            progress_q.put(num_processed % proc_update_interval)
            break

        num_processed += 1
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
            proc_update_interval = update_progress(
                num_processed, proc_update_interval)
            continue

        try:
            _io_and_map_read(
                fast5_data, failed_reads_q, bc_subgrps, bc_grp, corr_grp,
                aligner, seq_samp_type, map_thr_buf, fast5_fn,
                num_processed, map_conn, outlier_thresh, compute_sd,
                obs_filter, index_q, q_score_thresh, sig_match_thresh, std_ref)
        except th.TomboError as e:
            failed_reads_q.put((str(e), fast5_fn, True))
        finally:
            try:
                fast5_data.close()
            except:
                failed_reads_q.put(('Error closing fast5 file', fast5_fn, True))
        proc_update_interval = update_progress(
            num_processed, proc_update_interval)

    return


############################################
########## Multi-process Handling ##########
############################################

def _get_progress_fail_queues(
        progress_q, failed_reads_q, pf_conn, num_reads, failed_reads_fn,
        num_update_errors=0):
    def format_fail_summ(header, fail_summ=[], num_proc=0, num_errs=None):
        summ_errs = sorted(fail_summ)[::-1]
        if num_errs is not None:
            summ_errs = summ_errs[:num_errs]
            if len(summ_errs) < num_errs:
                summ_errs.extend([
                    (None, '') for _ in range(num_errs - len(summ_errs))])
        errs_str = '\n'.join(
            "{:8.1f}% ({:>7} reads)".format(100 * n_fns / float(num_proc),
                                            n_fns) + " : " + '{:<80}'.format(
                                                err)
            if (n_fns is not None and num_proc > 0) else
            '     -----' for n_fns, err in summ_errs)
        return '\n'.join((header, errs_str))


    if VERBOSE:
        th.status_message(
            'Re-squiggling reads (raw signal to genomic sequence alignment).')
        if num_update_errors > 0:
            # add lines for dynamic error messages
            sys.stderr.write(
                '\n'.join(['' for _ in range(num_update_errors + 2)]))
        bar = tqdm(total=num_reads, smoothing=0)
        if num_update_errors > 0:
            prog_prefix = ''.join(
                [_term_move_up(),] * (num_update_errors + 1)) + '\r'
            bar_update_header = (
                str(num_update_errors) + ' most common unsuccessful ' +
                'read types (approx. %):')
            # write failed read update header
            bar.write(prog_prefix + format_fail_summ(
                bar_update_header, num_errs=num_update_errors),
                      file=sys.stderr)

    tot_num_rec_proc, last_prog_update = 0, 0
    non_tombo_errors = []
    failed_reads = defaultdict(list)
    # continue to process the failed reads queue until the end signal
    # is sent via the failed_read_conn
    while True:
        try:
            tot_num_rec_proc += progress_q.get(block=False)
        except queue.Empty:
            # only update once the progress queue has emptied to make fewer
            # updates to the bar (can be annoying for cat'ing stderr)
            if VERBOSE and tot_num_rec_proc > last_prog_update:
                bar.update(tot_num_rec_proc - last_prog_update)
                last_prog_update = tot_num_rec_proc
            try:
                errorType, fn, is_tombo_error = failed_reads_q.get(block=False)
                if is_tombo_error:
                    failed_reads[errorType].append(fn)
                else:
                    failed_reads['Unexpected error'].append(fn)
                    if len(non_tombo_errors) < _MAX_NUM_UNEXP_ERRORS:
                        non_tombo_errors.append(fn + '\n:::\n' + errorType)
            except queue.Empty:
                # check if all reads are done signal was sent from main thread
                if pf_conn.poll():
                    break
                if VERBOSE and num_update_errors > 0:
                    bar.write(prog_prefix + format_fail_summ(
                        bar_update_header,
                        [(len(fns), err) for err, fns in failed_reads.items()],
                        tot_num_rec_proc, num_update_errors),
                              file=sys.stderr)
                sleep(0.5)
                continue

    # empty any entries left in queues after processes have finished
    while not progress_q.empty():
        tot_num_rec_proc += progress_q.get(block=False)
    if VERBOSE and tot_num_rec_proc > last_prog_update:
        bar.update(tot_num_rec_proc - last_prog_update)

    while not failed_reads_q.empty():
        errorType, fn, is_tombo_error = failed_reads_q.get(block=False)
        if is_tombo_error:
            failed_reads[errorType].append(fn)
        else:
            if len(non_tombo_errors) < _MAX_NUM_UNEXP_ERRORS:
                non_tombo_errors.append(fn + '\n:::\n' + errorType)

    if VERBOSE: bar.close()

    # close out failed read printout
    fail_summary = [(len(fns), err) for err, fns in failed_reads.items()]
    if VERBOSE:
        if len(non_tombo_errors) > 0:
            # add random value to filename in case multiple runs are made
            # from the same location
            unex_err_fn = _UNEXPECTED_ERROR_FN.format(
                np.random.randint(10000))
            th.warning_message(
                'Unexpected errors occured. See full error stack traces ' +
                'for first (up to) {0:d} errors in "{1}"'.format(
                    _MAX_NUM_UNEXP_ERRORS, unex_err_fn))
            with io.open(unex_err_fn, 'w') as fp:
                fp.write('\n\n'.join(non_tombo_errors) + '\n')

        if len(fail_summary) > 0:
            total_num_failed = sum(map(itemgetter(0), fail_summary))
            header = (
                'Final unsuccessful reads summary ' +
                '({:.1%} reads unsuccessfully processed; {} total reads):'.format(
                    float(total_num_failed) / num_reads, total_num_failed))
            th.status_message(format_fail_summ(header, fail_summary, num_reads))
        else:
            th.status_message('All reads successfully re-squiggled!')

    if failed_reads_fn is not None:
        with io.open(failed_reads_fn, 'wt') as fp:
            fp.write('\n'.join((
                err + '\t' + ', '.join(fns)
                for err, fns in failed_reads.items())) + '\n')

    pf_conn.send(tot_num_rec_proc)

    return

def _get_index_queue(index_q, index_conn, fast5s_dir, corr_grp):
    # open TomboReads object for storing and eventually writing the index data
    reads_index = th.TomboReads([fast5s_dir,], corr_grp, for_writing=True)
    # continue to process the index queue until the end signal
    # is sent via the index_conn
    while True:
        try:
            reads_index.add_read_data(*index_q.get(block=False))
        except queue.Empty:
            if index_conn.poll():
                break
            sleep(0.1)
            continue

    # empty any entries left in queue after processes have finished
    while not index_q.empty():
        reads_index.add_read_data(*index_q.get(block=False))

    # write index out to file
    reads_index.write_index_file()

    return

def _fill_files_queue(fast5_q, fast5_fns, num_threads):
    for fast5_fn in fast5_fns:
        fast5_q.put(fast5_fn)
    for _ in range(num_threads):
        fast5_q.put(None)

    return

def resquiggle_all_reads(
        fast5_fns, aligner, bc_grp, bc_subgrps, corr_grp, std_ref,
        seq_samp_type, outlier_thresh, overwrite, num_ps, threads_per_proc,
        compute_sd, skip_index, rsqgl_params, save_params, sig_match_thresh,
        obs_filter, const_scale, q_score_thresh, skip_seq_scaling,
        max_scaling_iters, failed_reads_fn, fast5s_basedir, num_update_errors):
    """Perform genomic alignment and re-squiggle algorithm
    """
    fast5_q = mp.Queue(maxsize=_MAX_QUEUE_SIZE)
    failed_reads_q = mp.Queue()
    index_q = mp.Queue(maxsize=_MAX_QUEUE_SIZE) if not skip_index else None
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
            proc_rsqgl_conns, std_ref, outlier_thresh, corr_grp, seq_samp_type,
            rsqgl_params, save_params, const_scale, skip_seq_scaling,
            max_scaling_iters)
        rsqgl_process = mp.Process(target=_resquiggle_worker, args=rsqgl_args)
        rsqgl_process.daemon = True
        rsqgl_process.start()
        rsqgl_ps.append(rsqgl_process)

    # failed read and progress queues getter
    main_pf_conn, pf_conn = mp.Pipe()
    pf_p = mp.Process(target=_get_progress_fail_queues,
                      args=(progress_q, failed_reads_q, pf_conn, len(fast5_fns),
                            failed_reads_fn, num_update_errors))
    pf_p.daemon = True
    pf_p.start()

    # index queue getter
    if index_q is not None:
        main_index_conn, index_conn = mp.Pipe()
        index_p = mp.Process(target=_get_index_queue, args=(
            index_q, index_conn, fast5s_basedir, corr_grp))
        index_p.daemon = True
        index_p.start()

    # now open mapping thread for each map connection created above
    resquiggle_ts = []
    for map_conn in map_conns:
        map_args = (fast5_q, progress_q, failed_reads_q, index_q, bc_grp,
                    bc_subgrps, corr_grp, aligner, outlier_thresh, compute_sd,
                    sig_match_thresh, obs_filter, seq_samp_type,
                    overwrite, map_conn, q_score_thresh, std_ref)
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

    # in a very unlikely case the progress/fail queue could die while the
    # main process remains active and thus we would have a deadlock here
    if pf_p.is_alive():
        # send signal to getter queue to finish and return results
        main_pf_conn.send(True)
        # returns total number of processed reads if that is needed
        main_pf_conn.recv()
        pf_p.join()

    if index_q is not None:
        main_index_conn.send(True)
        index_p.join()

    return


###################################
########## Main Function ##########
###################################

def _parse_files_and_lock_dirs(args):
    if VERBOSE: th.status_message('Getting file list.')
    try:
        if not os.path.isdir(args.fast5s_basedir):
            th.error_message_and_exit(
                'Provided [fast5-basedir] is not a directory.')
        fast5s_basedir = (
            args.fast5s_basedir if args.fast5s_basedir.endswith('/') else
            args.fast5s_basedir + '/')

        files, lock_fns = th.get_files_list_and_lock_dirs(
            fast5s_basedir, args.ignore_read_locks)
    except OSError:
        th.error_message_and_exit(
            'Reads base directory or a sub-directory does not appear to be ' +
            'accessible. Check directory permissions.')
    if len(files) < 1:
        th.clear_tombo_locks(lock_fns)
        th.error_message_and_exit(
            'No files identified in the specified ' +
            'directory or within immediate subdirectories.')

    if not th.reads_contain_basecalls(
            files, args.basecall_group, num_reads=1000):
        th.clear_tombo_locks(lock_fns)
        th.error_message_and_exit(
            'Reads do not to contain basecalls. Check --basecall-group ' +
            'option if basecalls are stored in non-standard location or ' +
            'use `tombo preprocess annotate_raw_with_fastqs` to add ' +
            'basecalls from FASTQ files to raw FAST5 files.')

    return files, fast5s_basedir, lock_fns

def _resquiggle_main(args):
    """Main method for resquiggle
    """
    if args.processes > 1 and _DEBUG_PLOTTING:
        th.error_message_and_exit(
            'Cannot run multiple processes and debugging.')
    if _DEBUG_PLOTTING:
        th.warning_message(
            'Producing de-bug plotting output. Can be very slow and should ' +
            'only be run on a small number of files.')
    if _DRY_RUN:
        th.warning_message(
            'Producing de-bug output. Not saving re-squiggle results.')
    if _DEBUG_BANDWIDTH or _DEBUG_START_BANDWIDTH:
        sys.stdout.write(
            'bandwidth\tmin_bw_edge_buffer\tmean_dp_score\tread_id\n')

    global VERBOSE
    VERBOSE = not args.quiet
    th.VERBOSE = VERBOSE
    ts.VERBOSE = VERBOSE

    if args.print_advanced_arguments:
        from . import _option_parsers
        _option_parsers.print_advanced_resquiggle()
        sys.exit()

    if args.basecall_group == args.corrected_group:
        th.error_message_and_exit(
            '--basecall-group and --corrected-group must ' +
            'be different.')

    # check simple arguments for validity first
    outlier_thresh = args.outlier_threshold
    if outlier_thresh is not None and outlier_thresh <= 0:
        outlier_thresh = None
    obs_filter = th.parse_obs_filter(args.obs_per_base_filter) \
                 if 'obs_per_base_filter' in args else None

    if VERBOSE: th.status_message('Loading minimap2 reference.')
    # to be enabled when mappy genome sequence extraction bug is fixed
    aligner = mappy.Aligner(
        str(args.reference), preset=str('map-ont'), best_n=1)
    if not aligner:
        th.error_message_and_exit(
            'Failed to load reference genome FASTA for mapping.')

    # get files as late as possible in startup since it takes the longest
    # and so other errors can't happen after locks are written
    files, fast5s_basedir, lock_fns = _parse_files_and_lock_dirs(args)

    try:
        seq_samp_type = None
        if args.seq_sample_type is not None:
            seq_samp_type = th.seqSampleType(RNA_SAMP_TYPE, True) \
                            if args.seq_sample_type == RNA_SAMP_TYPE else \
                               th.seqSampleType(DNA_SAMP_TYPE, False)
        if args.tombo_model_filename is not None and seq_samp_type is None:
            seq_samp_type = th.get_seq_sample_type(fast5_fns=files)
        # parse tombo model
        std_ref = ts.TomboModel(
            ref_fn=args.tombo_model_filename, seq_samp_type=seq_samp_type,
            fast5_fns=files)
        seq_samp_type = std_ref.seq_samp_type
        if seq_samp_type.name == DNA_SAMP_TYPE and USE_START_CLIP_BASES:
            std_ref = std_ref.reverse_sequence_copy()

        sig_match_thresh = args.signal_matching_score
        if sig_match_thresh is None:
            sig_match_thresh = SIG_MATCH_THRESH[seq_samp_type.name]

        const_scale = None
        if args.fixed_scale is not None:
            const_scale = args.fixed_scale
        elif args.fit_global_scale:
            const_scale = ts.estimate_global_scale(files)

        rsqgl_params = ts.load_resquiggle_parameters(
            seq_samp_type, args.signal_align_parameters,
            args.segmentation_parameters)
        save_params = ts.load_resquiggle_parameters(
            seq_samp_type, args.signal_align_parameters,
            args.segmentation_parameters, use_save_bandwidth=True)

        resquiggle_all_reads(
            files, aligner, args.basecall_group, args.basecall_subgroups,
            args.corrected_group, std_ref, seq_samp_type, outlier_thresh,
            args.overwrite, args.processes, args.threads_per_process,
            args.include_event_stdev, args.skip_index,
            rsqgl_params, save_params, sig_match_thresh,
            obs_filter, const_scale, args.q_score,
            args.skip_sequence_rescaling, args.max_scaling_iterations,
            args.failed_reads_filename, fast5s_basedir,
            args.num_most_common_errors)
    finally:
        th.clear_tombo_locks(lock_fns)

    return

if __name__ == '__main__':
    sys.stderr.write('This is a module. See commands with `tombo -h`')
    sys.exit(1)
