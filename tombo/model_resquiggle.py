import os, sys

import Queue

import numpy as np
np.seterr(all='raise')
import multiprocessing as mp

from time import sleep
from itertools import repeat
from collections import defaultdict

# import tombo functions
import tombo_stats as ts
import tombo_helper as th

from c_helper import c_new_means
from dynamic_programming import c_reg_z_scores, c_base_forward_pass, \
    c_base_traceback

VERBOSE = False

OPTIMIZE_RSQGL = False

# DEBUG should only be run with a single thread and on a single read
DEBUG_SIGNAL = False
DEBUG_BASE = False
assert not (DEBUG_SIGNAL and DEBUG_BASE)
PROGRESS_INTERVAL = 100

# debug print functions
def write_fwd_scores(score_fp, ld_fp, reg_fwd_scores, reg_id, iter_num=0):
    score_fp.write('\n'.join(
        '\t'.join(map(str, (base_i, pos, score, reg_id, iter_num)))
        for base_i, (b_data, b_last_diag, b_range) in enumerate(
                reg_fwd_scores) for pos, score in
        zip(range(b_range[0], b_range[1]), b_data)) + '\n')
    ld_fp.write('\n'.join(
        '\t'.join(map(str, (base_i, pos, last_count, reg_id, iter_num)))
        for base_i, (b_data, ld_data, b_range) in enumerate(reg_fwd_scores)
        for pos, last_count in zip(range(b_range[0], b_range[1]),
                                   ld_data)) + '\n')
    return
def write_z_scores(zs_fp, z_scores, reg_id, iter_num=0):
    zs_fp.write('\n'.join(
        '\t'.join(map(str, (base_i, pos, zscore, reg_id, iter_num)))
        for base_i, (b_data, b_range) in enumerate(z_scores)
        for pos, zscore in zip(range(b_range[0], b_range[1]), b_data)) + '\n')
    return
def write_segs(segs_fp, segs, last_seg, reg_id, iter_num=0):
    dbg_sig_new_segs = np.concatenate([[0], segs, [last_seg]])
    segs_fp.write('\n'.join(
        '\t'.join(map(str, (base_i, b_start, reg_id, iter_num))) + '\n' +
        '\t'.join(map(str, (base_i, b_end - 1, reg_id, iter_num)))
        for base_i, (b_start, b_end) in enumerate(
                zip(dbg_sig_new_segs[:-1], dbg_sig_new_segs[1:]))) + '\n')
    return
def write_sig(sig_fp, sig, reg_id, iter_num=0):
    sig_fp.write('\n'.join('\t'.join(map(str, (pos, sig_val, reg_id, iter_num)))
                           for pos, sig_val in enumerate(sig)) + '\n')
    return
def write_switch(s_fp, switch_points, reg_id, iter_num=0):
    s_fp.write('\n'.join('\t'.join(map(str, (base_i, sig_i, reg_id, iter_num)))
                         for base_i, sig_is in enumerate(switch_points)
                         for sig_i in sig_is) + '\n')
    return


def forward_pass(reg_z_scores, min_obs_per_base):
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
                               dtype=np.int32) * min_obs_per_base
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

def traceback(reg_fwd_scores, min_obs_per_base):
    # traceback along maximally likely path

    # initilize array to store new segments
    new_segs = np.empty(len(reg_fwd_scores) - 1, dtype=np.int32)
    # get first two bases of data for lookups
    curr_base_sig = 1
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

def get_best_event_path(reg_z_scores, b_switch_pnts, min_obs_per_base):
    # calc cummulative sums for more efficient region sum computations
    reg_cumm_z = [(np.cumsum(np.concatenate([[0], b_data])), b_start)
                  for b_data, (b_start, _) in reg_z_scores]

    def get_base_z_mean(base_cumsum, b_start, curr_pos, prev_pos, prev_sum):
        return ((base_cumsum[curr_pos - b_start] -
                 base_cumsum[prev_pos - b_start]) / (
                     curr_pos - prev_pos)) + prev_sum

    prev_b_poss = [([0],0)]
    for base_sps, (b_cumm_z, b_start) in zip(b_switch_pnts, reg_cumm_z):
        curr_b_poss = []
        for switch_point in base_sps:
            # loop over switch points from last base
            prev_path, prev_sum = prev_b_poss[0]
            curr_max_path = prev_path
            curr_max_sum = get_base_z_mean(
                b_cumm_z, b_start, switch_point, prev_path[-1], prev_sum)
            for prev_path, prev_sum in prev_b_poss[1:]:
                # if this path extends past this next potential switch point
                # move on to the next switch point
                if prev_path[-1] + min_obs_per_base > switch_point:
                    break
                sp_event_mean_z = get_base_z_mean(
                    b_cumm_z, b_start, switch_point, prev_path[-1], prev_sum)
                if sp_event_mean_z > curr_max_sum:
                    curr_max_path = prev_path
                    curr_max_sum = sp_event_mean_z
            curr_b_poss.append((curr_max_path + [switch_point], curr_max_sum))
        prev_b_poss = curr_b_poss

    # get max over the final base
    end_pos = reg_z_scores[-1][-1][-1]
    prev_path, prev_sum = prev_b_poss[0]
    curr_max_path = prev_path
    b_cumm_z, b_start = reg_cumm_z[-1]
    curr_max_sum = get_base_z_mean(
        b_cumm_z, b_start, end_pos, prev_path[-1], prev_sum)
    for prev_path, prev_sum in prev_b_poss[1:]:
        sp_event_mean_z = get_base_z_mean(
            b_cumm_z, b_start, end_pos, prev_path[-1], prev_sum)
        if sp_event_mean_z > curr_max_sum:
            curr_max_path = prev_path
            curr_max_sum = sp_event_mean_z

    return np.array(curr_max_path[1:], dtype=np.int32)

def traceback_until(
        reg_fwd_scores, start_base, seq_pos, b_switch_pnts,
        tb_b_ranges, min_obs_per_base):
    """ perform traceback from this poition to the orgin or when the path
        hits another previous path
    """
    # get base data to become curr_data in first iteration
    next_b_data, _, (next_start, next_end) = reg_fwd_scores[start_base]
    for base_pos in range(start_base - 1, -1, -1):
        curr_b_data, curr_start = next_b_data, next_start
        next_b_data, _, (next_start, next_end) = reg_fwd_scores[base_pos]
        seq_pos = c_base_traceback(
            curr_b_data, curr_start, next_b_data, next_start, next_end,
            seq_pos - 1, min_obs_per_base)
        # update switch points and covered positions
        b_switch_pnts[base_pos].add(seq_pos)
        if seq_pos < tb_b_ranges[base_pos+1][0]:
            tb_b_ranges[base_pos+1] = (seq_pos, tb_b_ranges[base_pos+1][1])
        elif seq_pos > tb_b_ranges[base_pos+1][1]:
            tb_b_ranges[base_pos+1] = (tb_b_ranges[base_pos+1][0], seq_pos)
        else:
            # hit an already covered path so rest of path to origin
            # has been seen already
            break

    return b_switch_pnts, tb_b_ranges

def find_all_tb_paths(reg_z_scores, reg_fwd_scores, global_tb, min_obs_per_base,
                      min_b_stay_run):
    # find all *reasonable* locations where a potential move from one
    # base to the next could occur. These are marked by a transition
    # from a "move down" state to "stay" state indicated by the last
    # maximum taken from the next base as opposed to staying in the
    # current base at each signal transition.
    req_locations = []
    for base_pos, (_, b_last_diag, (b_start, b_end)) in enumerate(
            reg_fwd_scores[1:]):
        move_states = b_last_diag == 1
        stay_states = np.logical_not(move_states)
        valid_shifts = [move_states[:-min_b_stay_run],
                        stay_states[min_b_stay_run:]]
        for stay_offset in range(1, min_obs_per_base):
            valid_shifts.append(stay_states[
                stay_offset:-(min_b_stay_run - stay_offset)])
        req_locations.extend(zip(
            repeat(base_pos + 1),
            b_start + np.where(np.logical_and.reduce(valid_shifts))[0]))

    # store identified switch points from one base to the next
    b_switch_pnts = [set([pos]) for pos in global_tb]
    # store ranges in each base currently covered by a traceback path
    # which indicates the termination of a traceback iteration as the rest
    # of the path is then determined
    # TODO: This may have to be a set of ranges instead of a single interval
    # for the whole base (if a gap exists between two paths, but seems
    # unlikely with a window of 3 original bases.
    tb_b_ranges = np.concatenate([[0], global_tb, [
        reg_fwd_scores[-1][-1][-1] + 1]])
    tb_b_ranges = zip(tb_b_ranges[:-1], tb_b_ranges[1:] - 1)
    for base_pos, seq_pos in req_locations:
        path_i = []
        # add this position as a switch point for this base
        b_switch_pnts[base_pos-1].add(seq_pos)
        # if this position is already covered by a path continue
        if (tb_b_ranges[base_pos][0] <= seq_pos <=
            tb_b_ranges[base_pos][1] - min_obs_per_base):
            continue

        b_switch_pnts, tb_b_ranges = traceback_until(
            reg_fwd_scores, base_pos - 1, seq_pos,
            b_switch_pnts, tb_b_ranges, min_obs_per_base)

        # if the position is after any traced-back region traceback to the
        # orgin and to the region end, which requires a new forward pass with
        # a trimmed reg_z_scores. This may need some heuristic to determine
        # if a region has any potential to save computation here.
        if seq_pos > tb_b_ranges[base_pos][1] - min_obs_per_base:
            # perform forward pass
            b_data, (b_start, b_end) = reg_z_scores[base_pos]
            prev_b_start = seq_pos
            clipped_bases = prev_b_start - b_start
            b_data = b_data[clipped_bases:]
            s_reg_z_scores = [(b_data, (prev_b_start, b_end))]
            # trim z-scores to disallow impossible paths
            for b_data, (b_start, b_end) in reg_z_scores[base_pos+1:]:
                if b_start < prev_b_start + min_obs_per_base:
                    b_data = b_data[prev_b_start + min_obs_per_base -
                                    b_start:]
                    b_start = prev_b_start + min_obs_per_base
                s_reg_z_scores.append((b_data, (b_start, b_end)))
                prev_b_start = b_start
            s_reg_fwd_scores = forward_pass(s_reg_z_scores, min_obs_per_base)
            s_new_segs = traceback(s_reg_fwd_scores, min_obs_per_base)
            # update covered region and identified switch points data
            for s_base_pos, seg_pos in enumerate(s_new_segs):
                b_switch_pnts[base_pos + s_base_pos].add(seg_pos)
                if seg_pos < tb_b_ranges[base_pos + s_base_pos][0]:
                    tb_b_ranges[base_pos + s_base_pos] = (
                        seg_pos, tb_b_ranges[base_pos + s_base_pos][1])
                elif seg_pos > tb_b_ranges[base_pos + s_base_pos][1]:
                    tb_b_ranges[base_pos + s_base_pos] = (
                        tb_b_ranges[base_pos + s_base_pos][0], seg_pos)

    # sort switch points for each base
    b_switch_pnts = [sorted(b_sp) for b_sp in b_switch_pnts]

    return b_switch_pnts

def get_region_model_segs(
        reg_start, reg_end, r_b_starts, r_sig, max_base_shift,
        r_ref_means, r_ref_sds, min_obs_per_base, debug_fps=None, reg_id=None,
        min_b_stay_run=3, base_space_scoring=False, iter_num=None):
    def signal_space_pass(reg_z_scores):
        reg_fwd_scores = forward_pass(reg_z_scores, min_obs_per_base)
        # perform signal based scoring segmentation
        #  - it is ~60X faster than base space
        new_segs = traceback(reg_fwd_scores, min_obs_per_base)
        if DEBUG_SIGNAL:
            write_sig(debug_fps[0], r_sig[r_b_starts[reg_start]:
                                          r_b_starts[reg_end]], reg_id)
            write_z_scores(debug_fps[1], reg_z_scores, reg_id)
            write_segs(debug_fps[2],
                       r_b_starts[reg_start+1:reg_end] - r_b_starts[reg_start],
                       reg_z_scores[-1][-1][-1], reg_id)
            write_fwd_scores(debug_fps[3], debug_fps[4], reg_fwd_scores, reg_id)
            write_segs(debug_fps[5], new_segs,
                       reg_z_scores[-1][-1][-1], reg_id)

        return new_segs

    def base_space_pass(reg_z_scores):
        reg_fwd_scores = forward_pass(reg_z_scores, min_obs_per_base)
        # perform global traceback for the region
        global_tb = traceback(reg_fwd_scores, min_obs_per_base)
        # perform base-space scoring to avoid regions being *skipped through*
        # b/c of the signal space scoring allowing lower penalty for these bases
        b_switch_pnts = find_all_tb_paths(
            reg_z_scores, reg_fwd_scores, global_tb, min_obs_per_base,
            min_b_stay_run)
        new_segs = get_best_event_path(
            reg_z_scores, b_switch_pnts, min_obs_per_base)
        if DEBUG_BASE:
            write_sig(debug_fps[0], r_sig[r_b_starts[reg_start]:
                                          r_b_starts[reg_end]], reg_id, iter_num)
            write_z_scores(debug_fps[1], reg_z_scores, reg_id, iter_num)
            write_segs(debug_fps[2],
                       r_b_starts[reg_start+1:reg_end] - r_b_starts[reg_start],
                       reg_z_scores[-1][-1][-1], reg_id, iter_num)
            write_fwd_scores(debug_fps[3], debug_fps[4], reg_fwd_scores,
                             reg_id, iter_num)
            write_segs(debug_fps[5], global_tb, reg_z_scores[-1][-1][-1],
                       reg_id, iter_num)
            write_segs(debug_fps[6], new_segs, reg_z_scores[-1][-1][-1],
                       reg_id, iter_num)
            write_switch(debug_fps[7], b_switch_pnts, reg_id, iter_num)

        return new_segs

    if (min_obs_per_base * (reg_end - reg_start - 1) >=
        r_b_starts[reg_end - 1] - r_b_starts[reg_start]):
        raise NotImplementedError, (
            'Not enough signal to correct poor fitting region.')

    reg_z_scores = c_reg_z_scores(
        r_sig, r_ref_means, r_ref_sds, r_b_starts,
        reg_start, reg_end, max_base_shift, min_obs_per_base)
    if not base_space_scoring:
        new_segs = signal_space_pass(reg_z_scores)
    else:
        new_segs = base_space_pass(reg_z_scores)

    return new_segs + r_b_starts[reg_start]

def filter_regions(signif_shift_regs, r_prev_new_segs, r_pp_segs):
    if r_pp_segs is None: return signif_shift_regs
    filtered_regs = []
    for r_start, r_end in signif_shift_regs:
        if not all(r_prev_new_segs[r_start:r_end] == r_pp_segs[r_start:r_end]):
            filtered_regs.append((r_start, r_end))

    return filtered_regs

def model_resquiggle_read(
        r_data, kmer_ref, kmer_width, upstrm_bases, dnstrm_bases, z_trans_lag,
        z_thresh, reg_context, base_reg_context, max_base_shift, b_max_base_shift,
        min_obs_per_base, base_space_iters, new_corr_grp, compute_sd,
        debug_fps=None):
    # should also get signal here
    all_read_data = th.get_all_read_data(r_data)
    if all_read_data is None:
        raise NotImplementedError, ('Error parsing data from FAST5 file.')
    (r_means, r_seq, r_sig, r_b_starts, scale_vals, norm_type, outlier_thresh,
     genome_loc) = all_read_data
    r_ref_means, r_ref_sds = zip(*[
        kmer_ref[kmer] for kmer in [''.join(bs) for bs in zip(*[
            r_seq[i:] for i in range(kmer_width)])]])
    # add upstream NANs so all data passed to model shifts is on the same
    # coordinate system. Note that the nan values will never be accessed
    # as the shift regions don't let a region extend beyond the non-nan
    # statistic values
    r_ref_means = np.concatenate((([np.NAN] * upstrm_bases), r_ref_means))
    r_ref_sds = np.concatenate((([np.NAN] * upstrm_bases), r_ref_sds))

    # add NAN values so that shifted regions will line up with original
    # base regions since kmer upstream and downstream positions can't be tested
    window_z = np.concatenate((
        [np.NAN] * upstrm_bases,
        ts.calc_window_z_transform(
            r_means[upstrm_bases:-dnstrm_bases], r_ref_means[upstrm_bases:],
            r_ref_sds[upstrm_bases:], z_trans_lag), [np.NAN] * dnstrm_bases))
    signif_shift_regs = ts.get_read_signif_shift_regions(
        window_z, z_thresh, reg_context)

    # first perform signal space scored model re-squiggle
    r_prev_new_segs = r_b_starts
    for reg_id, (reg_start, reg_end) in enumerate(signif_shift_regs):
        reg_new_segs = get_region_model_segs(
            reg_start, reg_end, r_b_starts, r_sig, max_base_shift,
            r_ref_means, r_ref_sds, min_obs_per_base, debug_fps, reg_id)
        r_prev_new_segs[reg_start+1:reg_end] = reg_new_segs
    r_pp_segs = None

    for iter_num in range(base_space_iters):
        # get new base mean signal values
        # note that reference means and sds don't change since they are based
        # on sequence (which is un-changed)
        r_means = c_new_means(r_sig, r_prev_new_segs)
        window_z = np.concatenate((
            [np.NAN] * upstrm_bases,
            ts.calc_window_z_transform(
                r_means[upstrm_bases:-dnstrm_bases], r_ref_means[upstrm_bases:],
                r_ref_sds[upstrm_bases:], z_trans_lag), [np.NAN] * dnstrm_bases))
        signif_shift_regs = ts.get_read_signif_shift_regions(
            window_z, z_thresh, base_reg_context)
        # filter regions that didn't change in the last round of
        # base-space reqsuiggle
        signif_shift_regs = filter_regions(
            signif_shift_regs, r_prev_new_segs, r_pp_segs)

        # then perform base space scored model re-squiggle on those regions still
        # not fitting the model well (potentially sub-optimal scoring regions)
        r_new_segs = r_prev_new_segs
        for reg_id, (reg_start, reg_end) in enumerate(signif_shift_regs):
            reg_new_segs = get_region_model_segs(
                reg_start, reg_end, r_prev_new_segs, r_sig, b_max_base_shift,
                r_ref_means, r_ref_sds, min_obs_per_base, debug_fps, reg_id,
                base_space_scoring=True, iter_num=iter_num)
            r_new_segs[reg_start+1:reg_end] = reg_new_segs

        r_pp_segs = r_prev_new_segs
        r_prev_new_segs = r_new_segs

    bc_subgrp = r_data.corr_group.split('/')[1]
    th.write_new_fast5_group(
        r_data.fn, genome_loc, r_data.read_start_rel_to_raw,
        r_new_segs, r_seq, r_sig, scale_vals, new_corr_grp,
        bc_subgrp, norm_type, outlier_thresh, compute_sd)

    return

def model_resquiggle_worker(
        reads_q, failed_reads_q, tb_model_fn, z_trans_lag, z_thresh,
        reg_context, base_reg_context, max_base_shift, b_max_base_shift,
        min_obs_per_base, base_space_iters, new_corr_grp, compute_sd,
        overwrite, in_place, corr_group):
    kmer_ref, upstrm_bases, _, _ = ts.parse_tombo_model(tb_model_fn)
    kmer_width = len(next(kmer_ref.iterkeys()))
    dnstrm_bases = kmer_width - upstrm_bases - 1

    if DEBUG_SIGNAL or DEBUG_BASE:
        sig_fp = open('debug_signal_space.signal.txt', 'w')
        sig_fp.write('SignalPos\tSignal\tRegion\tIteration\n')
        zscore_fp = open('debug_signal_space.window_z_scores.txt', 'w')
        zscore_fp.write('BasePos\tSignalPos\tZScore\tRegion\tIteration\n')
        origP_fp = open('debug_signal_space.window_orig_path.txt', 'w')
        origP_fp.write('BasePos\tSignalPos\tRegion\tIteration\n')
        tb_fp = open('debug_signal_space.window_traceback.txt', 'w')
        tb_fp.write('BasePos\tSignalPos\tpathVal\tRegion\tIteration\n')
        ld_fp = open('debug_signal_space.window_last_diag.txt', 'w')
        ld_fp.write('BasePos\tSignalPos\tLastDiagCount\tRegion\tIteration\n')
        sigMaxP_fp = open('debug_signal_space.window_signal_max_path.txt', 'w')
        sigMaxP_fp.write('BasePos\tSignalPos\tRegion\tIteration\n')
        maxP_fp = open('debug_signal_space.window_max_path.txt', 'w')
        maxP_fp.write('BasePos\tSignalPos\tRegion\tIteration\n')
        spP_fp = open('debug_signal_space.window_switch_points.txt', 'w')
        spP_fp.write('BasePos\tSignalPos\tRegion\tIteration\n')
        debug_fps = (sig_fp, zscore_fp, origP_fp, tb_fp, ld_fp, sigMaxP_fp,
                     maxP_fp, spP_fp)
    else:
        debug_fps = None

    num_processed = 0
    while True:
        try:
            fn_reads = reads_q.get(block=False)
        except Queue.Empty:
            break

        num_processed += 1
        if VERBOSE and num_processed % PROGRESS_INTERVAL == 0:
            sys.stderr.write('.')
            sys.stderr.flush()

        prep_result = th.prep_fast5(
            fn_reads[0].fn, new_corr_grp, overwrite, in_place, corr_group)
        if prep_result is not None:
            try:
                th.write_error_status(
                    prep_result[1], corr_group, None, prep_result[0])
            except:
                pass
            failed_reads_q.put(prep_result)
            continue

        for r_data in fn_reads:
            try:
                model_resquiggle_read(
                    r_data, kmer_ref, kmer_width, upstrm_bases, dnstrm_bases,
                    z_trans_lag, z_thresh, reg_context, base_reg_context,
                    max_base_shift, b_max_base_shift, min_obs_per_base,
                    base_space_iters, new_corr_grp, compute_sd, debug_fps)
            except Exception as e:
                # uncomment to identify mysterious errors
                #raise
                try:
                    subgrp = r_data.corr_group.split('/')[1]
                    th.write_error_status(r_data.fn, corr_group, subgrp, str(e))
                except:
                    pass
                failed_reads_q.put((
                    str(e), r_data.corr_group + th.FASTA_NAME_JOINER + r_data.fn))

    return

if OPTIMIZE_RSQGL:
    model_resquiggle_wrapper = model_resquiggle_worker
    def model_resquiggle_worker(*args):
        import cProfile
        cProfile.runctx('model_resquiggle_wrapper(*args)', globals(), locals(),
                    filename='model_requiggle.prof')
        return


def model_resquiggle(
        f5_dirs1, corr_group, bc_subgrps,
        tb_model_fn, z_trans_lag, p_value_thresh, reg_context, base_reg_context,
        max_base_shift, b_max_base_shift, min_obs_per_base, base_space_iters,
        compute_sd, new_corr_grp, num_processes, overwrite, in_place=True):
    z_thresh = ts.p_value_to_z_score(p_value_thresh)
    raw_read_coverage = th.parse_fast5s(
        f5_dirs1, corr_group, bc_subgrps, new_corr_grp)

    if tb_model_fn is None:
        tb_model_fn = ts.get_default_standard_ref_from_files(fast5_fns)

    # load reads into Queue
    manager = mp.Manager()
    reads_q = manager.Queue()
    failed_reads_q = manager.Queue()

    # group reads by filename so slot is not deleted in 2D reads
    fn_grouped_reads = defaultdict(list)
    for cs_reads in raw_read_coverage.itervalues():
        for r_data in cs_reads:
            fn_grouped_reads[r_data.fn].append(r_data)
    num_reads = 0
    for fn_reads in fn_grouped_reads.itervalues():
        reads_q.put(fn_reads)
        num_reads += 1

    mod_rsqgl_args = (
        reads_q, failed_reads_q, tb_model_fn, z_trans_lag, z_thresh,
        reg_context, base_reg_context, max_base_shift, b_max_base_shift,
        min_obs_per_base, base_space_iters, new_corr_grp, compute_sd,
        overwrite, in_place, corr_group)
    mod_rsqgl_ps = []
    for p_id in xrange(num_processes):
        p = mp.Process(target=model_resquiggle_worker, args=mod_rsqgl_args)
        p.start()
        mod_rsqgl_ps.append(p)

    if VERBOSE: sys.stderr.write(
            'Correcting ' + str(num_reads) + ' files with ' +
            str(len(bc_subgrps)) + ' subgroup(s)/read(s) ' +
            'each (Will print a dot for each ' + str(PROGRESS_INTERVAL) +
            ' reads completed).\n')
    failed_reads = defaultdict(list)
    while any(p.is_alive() for p in mod_rsqgl_ps):
        try:
            errorType, fn = failed_reads_q.get(block=False)
            failed_reads[errorType].append(fn)
        except Queue.Empty:
            sleep(1)
            continue
    while not failed_reads_q.empty():
        errorType, fn = failed_reads_q.get(block=False)
        failed_reads[errorType].append(fn)

    # print newline after read progress dots
    if VERBOSE: sys.stderr.write('\n')

    return dict(failed_reads)

def model_resquiggle_main(args):
    global VERBOSE
    VERBOSE = not args.quiet
    th.VERBOSE = VERBOSE
    ts.VERBOSE = VERBOSE

    # whether or not to skip SD calculation due to time
    compute_sd = args.include_event_stdev

    failed_reads = model_resquiggle(
        args.fast5_basedirs, args.corrected_group, args.basecall_subgroups,
        args.tombo_model_filename, args.stouffer_z_context,
        args.p_value_threshold, args.region_context,
        args.base_score_region_context, args.max_bases_shift,
        args.base_score_max_bases_shift, args.min_obs_per_base,
        args.base_score_iterations, compute_sd,
        args.new_corrected_group, args.processes, args.overwrite)

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


if __name__ == '__main__':
    raise NotImplementedError, (
        'This is a module. See commands with `tombo -h`')
