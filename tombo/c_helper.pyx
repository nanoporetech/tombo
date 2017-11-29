from numpy cimport ndarray

import numpy as np
cimport numpy as np
DTYPE = np.float64
ctypedef np.float64_t DTYPE_t
from libc.math cimport log

cdef extern from "math.h":
    double sqrt(double m)

def c_mean_std(ndarray[DTYPE_t] values):
    """
    More efficient method to get both mean and standard deviation
    """
    cdef DTYPE_t v_mean, v_var
    cdef int idx
    cdef int v_len = values.shape[0]
    v_mean = 0
    for idx in range(v_len):
        v_mean += values[idx]
    v_mean /= v_len
    v_var = 0
    for idx in range(v_len):
        v_var += (values[idx] - v_mean)**2
    return v_mean, sqrt(v_var / v_len)

def c_new_mean_stds(ndarray[DTYPE_t] norm_signal not None,
                    ndarray[int] new_segs not None):
    cdef int n_segs = new_segs.shape[0] - 1
    cdef ndarray[DTYPE_t] means_arr = np.empty(n_segs, dtype=DTYPE)
    cdef ndarray[DTYPE_t] stds_arr = np.empty(n_segs, dtype=DTYPE)
    cdef DTYPE_t curr_sum, curr_var, seg_mean
    cdef int idx, seg_idx, seg_len
    for idx in range(n_segs):
        seg_len = new_segs[idx + 1] - new_segs[idx]
        curr_sum = 0
        for seg_idx in range(new_segs[idx], new_segs[idx + 1]):
            curr_sum += norm_signal[seg_idx]
        seg_mean = curr_sum / seg_len
        means_arr[idx] = seg_mean
        curr_var = 0
        for seg_idx in range(new_segs[idx], new_segs[idx + 1]):
            curr_var += (norm_signal[seg_idx] - seg_mean)**2
        stds_arr[idx] = sqrt(curr_var / seg_len)
    return means_arr, stds_arr

def c_new_means(ndarray[DTYPE_t] norm_signal not None,
                ndarray[int] new_segs not None):
    cdef int n_segs = new_segs.shape[0] - 1
    cdef ndarray[DTYPE_t] means_arr = np.empty(n_segs, dtype=DTYPE)
    cdef DTYPE_t curr_sum
    cdef int idx, seg_idx
    for idx in range(n_segs):
        curr_sum = 0
        for seg_idx in range(new_segs[idx], new_segs[idx + 1]):
            curr_sum += norm_signal[seg_idx]
        means_arr[idx] = curr_sum / (new_segs[idx + 1] - new_segs[idx])
    return means_arr

def c_apply_outlier_thresh(
        ndarray[DTYPE_t] raw_signal, float lower_lim, float upper_lim):
    cdef int raw_size = raw_signal.shape[0]
    cdef ndarray[DTYPE_t] clipped_signal = np.empty(raw_size, dtype=np.float)
    cdef int pos
    cdef DTYPE_t pos_sig
    for pos in range(raw_size):
        pos_sig = raw_signal[pos]
        if pos_sig > upper_lim:
            clipped_signal[pos] = upper_lim
        elif pos_sig < lower_lim:
            clipped_signal[pos] = lower_lim
        else:
            clipped_signal[pos] = pos_sig
    return clipped_signal

def c_valid_cpts_w_cap(
        ndarray[DTYPE_t] raw_signal, int min_base_obs, int num_cpts):
    cdef ndarray[DTYPE_t] raw_cumsum = np.cumsum(
        np.concatenate([[0.0], raw_signal]))
    # get difference between all neighboring min_base_obs regions
    cdef ndarray[int] candidate_poss = np.argsort(np.abs(
        (2 * raw_cumsum[min_base_obs:-min_base_obs]) -
        raw_cumsum[:-2*min_base_obs] -
        raw_cumsum[2*min_base_obs:])).astype(np.int32)[::-1]

    cdef ndarray[int] cpts = np.empty(num_cpts, dtype=np.int32)
    cpts[0] = candidate_poss[0] + min_base_obs
    blacklist_pos = set(range(
        candidate_poss[0] - min_base_obs + 1, candidate_poss[0] + min_base_obs))
    cdef int cand_pos
    cdef int num_cands = candidate_poss.shape[0]
    cdef int cand_idx = 1
    cdef int added_cpts = 1
    while added_cpts < num_cpts:
        cand_pos = candidate_poss[cand_idx]
        if cand_pos not in blacklist_pos:
            cpts[added_cpts] = cand_pos + min_base_obs
            added_cpts += 1
            blacklist_pos.update(range(
                cand_pos - min_base_obs + 1, cand_pos + min_base_obs))
        cand_idx += 1
        if cand_idx >= num_cands:
            raise NotImplementedError, 'Fewer changepoints found than requested'

    return cpts

def c_valid_cpts(ndarray[DTYPE_t] raw_signal, int min_base_obs):
    cdef ndarray[DTYPE_t] raw_cumsum = np.cumsum(
        np.concatenate([[0.0], raw_signal]))
    # get difference between all neighboring min_base_obs regions
    cdef ndarray[int] candidate_poss = np.argsort(np.abs(
        (2 * raw_cumsum[min_base_obs:-min_base_obs]) -
        raw_cumsum[:-2*min_base_obs] -
        raw_cumsum[2*min_base_obs:])).astype(np.int32)[::-1]

    cpts = [candidate_poss[0]]
    blacklist_pos = set()
    cdef int pos
    for pos in candidate_poss[1:]:
        if pos not in blacklist_pos:
            cpts.append(pos)
            blacklist_pos.update(range(
                pos-min_base_obs+1, pos+min_base_obs+1))

    return np.array(cpts) + min_base_obs

def c_valid_cpts_w_cap_t_test(
        ndarray[DTYPE_t] raw_signal, int min_base_obs, int num_cpts):
    cdef int pos, idx
    cdef DTYPE_t pos_diff, m1, m2, var1, var2
    cdef int num_cands = raw_signal.shape[0] - (min_base_obs * 2)
    # note these will not actually be t-scores, but will be a monotonic transform
    # so the rank order will be the same
    cdef ndarray[DTYPE_t] t_scores = np.empty(num_cands, dtype=DTYPE)
    for pos in range(num_cands):
        # compute means
        m1 = 0
        for idx in range(min_base_obs):
            m1 += raw_signal[pos + idx]
        m1 /= min_base_obs
        m2 = 0
        for idx in range(min_base_obs):
            m2 += raw_signal[pos + min_base_obs + idx]
        m2 /= min_base_obs

        # compute sum of variances
        var1 = 0
        for idx in range(min_base_obs):
            pos_diff = raw_signal[pos + idx] - m1
            var1 += pos_diff * pos_diff
        var2 = 0
        for idx in range(min_base_obs):
            pos_diff = raw_signal[pos + min_base_obs + idx] - m2
            var2 += pos_diff * pos_diff

        if var1 + var2 == 0:
            t_scores[pos] = 0.0
        elif m1 > m2:
            t_scores[pos] = (m1 - m2) / sqrt(var1 + var2)
        else:
            t_scores[pos] = (m2 - m1) / sqrt(var1 + var2)

    cdef ndarray[int] candidate_poss = np.argsort(
        t_scores).astype(np.int32)[::-1]

    cdef ndarray[int] cpts = np.empty(num_cpts, dtype=np.int32)
    cpts[0] = candidate_poss[0] + min_base_obs
    blacklist_pos = set(range(
        candidate_poss[0] - min_base_obs + 1, candidate_poss[0] + min_base_obs))
    cdef int cand_pos
    cdef int cand_idx = 1
    cdef int added_cpts = 1
    while added_cpts < num_cpts:
        cand_pos = candidate_poss[cand_idx]
        if cand_pos not in blacklist_pos:
            cpts[added_cpts] = cand_pos + min_base_obs
            added_cpts += 1
            blacklist_pos.update(range(
                cand_pos - min_base_obs + 1, cand_pos + min_base_obs))
        cand_idx += 1
        if cand_idx >= num_cands:
            raise NotImplementedError, 'Fewer changepoints found than requested.'

    return cpts

def c_calc_llh_ratio(
        ndarray[DTYPE_t] reg_means,
        ndarray[DTYPE_t] reg_ref_means, ndarray[DTYPE_t] reg_ref_vars,
        ndarray[DTYPE_t] reg_alt_means, ndarray[DTYPE_t] reg_alt_vars):
    cdef float ref_z_sum, ref_log_var_sum, alt_z_sum, alt_log_var_sum
    ref_z_sum, ref_log_var_sum, alt_z_sum, alt_log_var_sum = (0.0, 0.0 ,0.0 ,0.0)
    cdef float ref_diff, alt_diff, log_lh_ratio
    cdef int idx
    cdef int reg_len = reg_means.shape[0]
    for idx in range(reg_len):
        ref_diff = reg_means[idx] - reg_ref_means[idx]
        ref_z_sum += (ref_diff * ref_diff) / reg_ref_vars[idx]
        ref_log_var_sum += log(reg_ref_vars[idx])
        alt_diff = reg_means[idx] - reg_alt_means[idx]
        alt_z_sum += (alt_diff * alt_diff) / reg_alt_vars[idx]
        alt_log_var_sum += log(reg_alt_vars[idx])

    log_lh_ratio = alt_z_sum + alt_log_var_sum - ref_z_sum - ref_log_var_sum

    return log_lh_ratio
