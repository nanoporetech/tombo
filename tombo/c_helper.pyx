import numpy as np
cimport numpy as np

DTYPE = np.float64
ctypedef np.float64_t DTYPE_t

DTYPE_INT = np.int64
ctypedef np.int64_t DTYPE_INT_t

from libc.math cimport log, exp

cdef extern from "math.h":
    double sqrt(double m)

def c_mean_std(np.ndarray[DTYPE_t] values):
    """
    More efficient method to get both mean and standard deviation
    """
    cdef DTYPE_t v_mean, v_var
    cdef DTYPE_INT_t idx
    cdef DTYPE_INT_t v_len = values.shape[0]
    v_mean = 0
    for idx in range(v_len):
        v_mean += values[idx]
    v_mean /= v_len
    v_var = 0
    for idx in range(v_len):
        v_var += (values[idx] - v_mean)**2
    return v_mean, sqrt(v_var / v_len)

def c_new_mean_stds(np.ndarray[DTYPE_t] norm_signal not None,
                    np.ndarray[DTYPE_INT_t] new_segs not None):
    cdef DTYPE_INT_t n_segs = new_segs.shape[0] - 1
    cdef np.ndarray[DTYPE_t] means_arr = np.empty(n_segs, dtype=DTYPE)
    cdef np.ndarray[DTYPE_t] stds_arr = np.empty(n_segs, dtype=DTYPE)
    cdef DTYPE_t curr_sum, curr_var, seg_mean
    cdef DTYPE_INT_t idx, seg_idx, seg_len
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

def c_new_means(np.ndarray[DTYPE_t] norm_signal not None,
                np.ndarray[DTYPE_INT_t] new_segs not None):
    cdef DTYPE_INT_t n_segs = new_segs.shape[0] - 1
    cdef np.ndarray[DTYPE_t] means_arr = np.empty(n_segs, dtype=DTYPE)
    cdef DTYPE_t curr_sum
    cdef DTYPE_INT_t idx, seg_idx
    for idx in range(n_segs):
        curr_sum = 0
        for seg_idx in range(new_segs[idx], new_segs[idx + 1]):
            curr_sum += norm_signal[seg_idx]
        means_arr[idx] = curr_sum / (new_segs[idx + 1] - new_segs[idx])
    return means_arr

def c_apply_outlier_thresh(
        np.ndarray[DTYPE_t] raw_signal, DTYPE_t lower_lim, DTYPE_t upper_lim):
    cdef DTYPE_INT_t raw_size = raw_signal.shape[0]
    cdef np.ndarray[DTYPE_t] clipped_signal = np.empty(raw_size, dtype=DTYPE)
    cdef DTYPE_INT_t pos
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
        np.ndarray[DTYPE_t] raw_signal, DTYPE_INT_t min_base_obs,
        DTYPE_INT_t running_stat_width, DTYPE_INT_t num_cpts):
    cdef np.ndarray[DTYPE_t] raw_cumsum = np.cumsum(
        np.concatenate([[0.0], raw_signal]))
    # get difference between all neighboring running_stat_width regions
    cdef np.ndarray[DTYPE_INT_t] candidate_poss = np.argsort(np.abs(
        (2 * raw_cumsum[running_stat_width:-running_stat_width]) -
        raw_cumsum[:-2*running_stat_width] -
        raw_cumsum[2*running_stat_width:])).astype(DTYPE_INT)[::-1]

    cdef np.ndarray[DTYPE_INT_t] cpts = np.empty(num_cpts, dtype=DTYPE_INT)
    cpts[0] = candidate_poss[0] + running_stat_width
    blacklist_pos = set(range(
        candidate_poss[0] - min_base_obs + 1, candidate_poss[0] + min_base_obs))
    cdef DTYPE_INT_t cand_pos
    cdef DTYPE_INT_t num_cands = candidate_poss.shape[0] - (
        2 * running_stat_width)
    cdef DTYPE_INT_t cand_idx = 1
    cdef DTYPE_INT_t added_cpts = 1
    while added_cpts < num_cpts:
        cand_pos = candidate_poss[cand_idx]
        if cand_pos not in blacklist_pos:
            cpts[added_cpts] = cand_pos + running_stat_width
            added_cpts += 1
            blacklist_pos.update(range(
                cand_pos - min_base_obs + 1, cand_pos + min_base_obs))
        cand_idx += 1
        if cand_idx >= num_cands:
            raise NotImplementedError, 'Fewer changepoints found than requested'

    return cpts

def c_valid_cpts(np.ndarray[DTYPE_t] raw_signal, DTYPE_INT_t min_base_obs,
                 DTYPE_INT_t running_stat_width):
    cdef np.ndarray[DTYPE_t] raw_cumsum = np.cumsum(
        np.concatenate([[0.0], raw_signal]))
    # get difference between all neighboring running_stat_width regions
    cdef np.ndarray[DTYPE_INT_t] candidate_poss = np.argsort(np.abs(
        (2 * raw_cumsum[running_stat_width:-running_stat_width]) -
        raw_cumsum[:-2*running_stat_width] -
        raw_cumsum[2*running_stat_width:])).astype(DTYPE_INT)[::-1]

    cpts = [candidate_poss[0]]
    blacklist_pos = set()
    cdef DTYPE_INT_t pos
    for pos in candidate_poss[1:]:
        if pos not in blacklist_pos:
            cpts.append(pos)
            blacklist_pos.update(range(
                pos-min_base_obs+1, pos+min_base_obs+1))

    return np.array(cpts) + running_stat_width

def c_valid_cpts_w_cap_t_test(
        np.ndarray[DTYPE_t] raw_signal, DTYPE_INT_t min_base_obs,
        DTYPE_INT_t running_stat_width, DTYPE_INT_t num_cpts):
    cdef DTYPE_INT_t pos, idx
    cdef DTYPE_t pos_diff, m1, m2, var1, var2
    cdef DTYPE_INT_t num_cands = raw_signal.shape[0] - (running_stat_width * 2)
    # note these will not actually be t-scores, but will be a monotonic transform
    # so the rank order will be the same
    cdef np.ndarray[DTYPE_t] t_scores = np.empty(num_cands, dtype=DTYPE)
    for pos in range(num_cands):
        # compute means
        m1 = 0
        for idx in range(running_stat_width):
            m1 += raw_signal[pos + idx]
        m1 /= running_stat_width
        m2 = 0
        for idx in range(running_stat_width):
            m2 += raw_signal[pos + running_stat_width + idx]
        m2 /= running_stat_width

        # compute sum of variances
        var1 = 0
        for idx in range(running_stat_width):
            pos_diff = raw_signal[pos + idx] - m1
            var1 += pos_diff * pos_diff
        var2 = 0
        for idx in range(running_stat_width):
            pos_diff = raw_signal[pos + running_stat_width + idx] - m2
            var2 += pos_diff * pos_diff

        if var1 + var2 == 0:
            t_scores[pos] = 0.0
        elif m1 > m2:
            t_scores[pos] = (m1 - m2) / sqrt(var1 + var2)
        else:
            t_scores[pos] = (m2 - m1) / sqrt(var1 + var2)

    cdef np.ndarray[DTYPE_INT_t] candidate_poss = np.argsort(
        t_scores).astype(DTYPE_INT)[::-1]

    cdef np.ndarray[DTYPE_INT_t] cpts = np.empty(num_cpts, dtype=DTYPE_INT)
    cpts[0] = candidate_poss[0] + running_stat_width
    blacklist_pos = set(range(
        candidate_poss[0] - min_base_obs + 1, candidate_poss[0] + min_base_obs))
    cdef DTYPE_INT_t cand_pos
    cdef DTYPE_INT_t cand_idx = 1
    cdef DTYPE_INT_t added_cpts = 1
    while added_cpts < num_cpts:
        cand_pos = candidate_poss[cand_idx]
        if cand_pos not in blacklist_pos:
            cpts[added_cpts] = cand_pos + running_stat_width
            added_cpts += 1
            blacklist_pos.update(range(
                cand_pos - min_base_obs + 1, cand_pos + min_base_obs))
        cand_idx += 1
        if cand_idx >= num_cands:
            raise NotImplementedError, 'Fewer changepoints found than requested.'

    return cpts

def c_calc_llh_ratio(
        np.ndarray[DTYPE_t] reg_means,
        np.ndarray[DTYPE_t] reg_ref_means, np.ndarray[DTYPE_t] reg_ref_vars,
        np.ndarray[DTYPE_t] reg_alt_means, np.ndarray[DTYPE_t] reg_alt_vars):
    cdef DTYPE_t ref_z_sum, ref_log_var_sum, alt_z_sum, alt_log_var_sum
    ref_z_sum, ref_log_var_sum, alt_z_sum, alt_log_var_sum = 0.0, 0.0 ,0.0 ,0.0
    cdef DTYPE_t ref_diff, alt_diff, log_lh_ratio
    cdef DTYPE_INT_t idx
    for idx in range(reg_means.shape[0]):
        ref_diff = reg_means[idx] - reg_ref_means[idx]
        ref_z_sum += (ref_diff * ref_diff) / reg_ref_vars[idx]
        ref_log_var_sum += log(reg_ref_vars[idx])

        alt_diff = reg_means[idx] - reg_alt_means[idx]
        alt_z_sum += (alt_diff * alt_diff) / reg_alt_vars[idx]
        alt_log_var_sum += log(reg_alt_vars[idx])

    log_lh_ratio = alt_z_sum + alt_log_var_sum - ref_z_sum - ref_log_var_sum

    return log_lh_ratio

def c_calc_llh_ratio_const_var(
        np.ndarray[DTYPE_t] reg_means, np.ndarray[DTYPE_t] reg_ref_means,
        np.ndarray[DTYPE_t] reg_alt_means, DTYPE_t const_var):
    cdef DTYPE_t ref_diff, alt_diff, running_llhr, obs_mean
    cdef DTYPE_INT_t idx
    for idx in range(reg_means.shape[0]):
        obs_mean = reg_means[idx]
        ref_diff = obs_mean - reg_ref_means[idx]
        alt_diff = obs_mean - reg_alt_means[idx]
        running_llhr += ((alt_diff * alt_diff) -
                         (ref_diff * ref_diff)) / const_var

    return running_llhr

def c_calc_scaled_llh_ratio_const_var(
        np.ndarray[DTYPE_t] reg_means, np.ndarray[DTYPE_t] reg_ref_means,
        np.ndarray[DTYPE_t] reg_alt_means, DTYPE_t const_var,
        DTYPE_t scale_factor, DTYPE_t density_height_factor,
        DTYPE_t density_height_power):
    """
    Scale log likelihood ratio with the normal distribution halfway
    between the 2 distributions.

    scale_factor - sets the spread of the value (2 makes peaks equal the normal
        density centers, but this is very sharp near the boundary between the
        reference and alternative densities
    density_height_factor - globally scales the height of the scores. Set to
        approximately match log likelihood scale.
    density_height_power - scales the density height proportional to the
        difference between the reference and alternate means. 0.5 makes all
        densities peak at the same value. Recommend values between 0 and 0.5
        so that more divergent reference and alternate densities contrbute more
        to the score.
    """
    cdef DTYPE_t running_scaled_lhr = 0.0
    cdef DTYPE_t ref_diff, alt_diff, ref_mean, alt_mean, scale_diff, \
        obs_mean, means_diff
    cdef DTYPE_INT_t idx
    for idx in range(reg_means.shape[0]):
        ref_mean = reg_ref_means[idx]
        alt_mean = reg_alt_means[idx]
        if ref_mean == alt_mean:
            continue
        obs_mean = reg_means[idx]
        scale_mean = (alt_mean + ref_mean) / 2

        ref_diff = obs_mean - ref_mean
        alt_diff = obs_mean - alt_mean
        scale_diff = obs_mean - scale_mean
        means_diff = alt_mean - ref_mean
        if means_diff < 0:
            means_diff = means_diff * -1

        running_scaled_lhr += exp(
            -(scale_diff * scale_diff) / (scale_factor * const_var)) * (
                (alt_diff * alt_diff) - (ref_diff * ref_diff)) / (
                    const_var * (means_diff ** density_height_power) *
                    density_height_factor)

    return running_scaled_lhr
