from __future__ import unicode_literals, absolute_import

from builtins import int, range, dict

import numpy as np
np.seterr(all='raise')

from .c_dynamic_programming import c_base_forward_pass, c_base_traceback

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

def traceback(reg_fwd_scores, min_obs_per_base):
    # traceback along maximally likely path

    # initilize array to store new segments
    new_segs = np.empty(len(reg_fwd_scores) - 1, dtype=np.int64)
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


if __name__ == '__main__':
    raise NotImplementedError(
        'This is a module. See commands with `tombo -h`')
