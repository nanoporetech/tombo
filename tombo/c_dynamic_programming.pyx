cimport cython

import numpy as np
cimport numpy as np

DTYPE = np.float64
ctypedef np.float64_t DTYPE_t

DTYPE_INT = np.int64
ctypedef np.int64_t DTYPE_INT_t

from libcpp cimport bool

def c_base_z_scores(np.ndarray[DTYPE_t] b_sig not None,
                    DTYPE_t ref_mean, DTYPE_t ref_sd):
    cdef DTYPE_INT_t n_sig = b_sig.shape[0]
    b_z_scores = np.empty(n_sig, dtype=DTYPE)
    cdef DTYPE_t b_pos_z_score
    cdef DTYPE_INT_t idx
    for idx in range(n_sig):
        b_pos_z_score = (b_sig[idx] - ref_mean) / ref_sd
        if b_pos_z_score > 0:
            # convert all z-scores to lower tail
            b_pos_z_score *= -1
        b_z_scores[idx] = b_pos_z_score
    return b_z_scores

def c_reg_z_scores(
        np.ndarray[DTYPE_t] r_sig not None,
        np.ndarray[DTYPE_t] r_ref_means not None,
        np.ndarray[DTYPE_t] r_ref_sds not None,
        np.ndarray[DTYPE_INT_t] r_b_starts not None,
        DTYPE_INT_t reg_start, DTYPE_INT_t reg_end,
        DTYPE_INT_t max_base_shift, DTYPE_INT_t min_obs_per_base):
    cdef DTYPE_INT_t base_i, b_sig_start, b_sig_end, prev_sig_start, \
        prev_sig_end, idx
    cdef DTYPE_INT_t reg_len = reg_end - reg_start
    cdef np.ndarray[DTYPE_INT_t] sig_starts = np.empty(reg_len, dtype=DTYPE_INT)
    prev_start_set = False
    cdef np.ndarray[DTYPE_INT_t] base_range = np.arange(
        reg_start, reg_end, dtype=DTYPE_INT)
    for idx in range(reg_len):
        base_i = base_range[idx]
        b_sig_start = r_b_starts[max(reg_start, base_i - max_base_shift)]
        # clip observations from the beginning of a base if there is no
        # possible traceback path through that location
        if (prev_start_set and
            b_sig_start < prev_sig_start + min_obs_per_base):
            b_sig_start = prev_sig_start + min_obs_per_base
        prev_start_set = True
        sig_starts[idx] = b_sig_start
        prev_sig_start = b_sig_start

    cdef np.ndarray[DTYPE_INT_t] sig_ends = np.empty(reg_len, dtype=DTYPE_INT)
    prev_end_set = False
    # clip positions from the end of each base
    for idx in range(reg_len):
        base_i = base_range[reg_len - idx - 1]
        b_sig_end = r_b_starts[min(reg_end - 1, base_i + max_base_shift + 1)]
        # clip observations from the end of a base if there is no
        # possible traceback path through that location
        if (prev_end_set and
            b_sig_end > prev_sig_end - min_obs_per_base):
            b_sig_end = prev_sig_end - min_obs_per_base
        prev_end_set = True
        sig_ends[reg_len - idx - 1] = b_sig_end
        prev_sig_end = b_sig_end

    reg_scores = []
    for idx in range(reg_len):
        base_i = base_range[idx]
        b_sig_start = sig_starts[idx]
        b_sig_end = sig_ends[idx]
        # z-score computation is far more efficient than p-values and
        # produces *very* similar results
        reg_scores.append((
            c_base_z_scores(r_sig[b_sig_start:b_sig_end],
                            r_ref_means[base_i], r_ref_sds[base_i]), (
                                b_sig_start-r_b_starts[reg_start],
                                b_sig_end-r_b_starts[reg_start])))

    return reg_scores

def c_base_forward_pass(
        np.ndarray[DTYPE_t] b_data not None,
        DTYPE_INT_t b_start, DTYPE_INT_t b_end,
        np.ndarray[DTYPE_t] prev_b_data not None,
        DTYPE_INT_t prev_b_start, DTYPE_INT_t prev_b_end,
        np.ndarray[DTYPE_t] prev_b_fwd_data not None,
        np.ndarray[DTYPE_INT_t] prev_b_last_diag not None,
        DTYPE_INT_t min_obs_per_base):
    cdef DTYPE_INT_t b_len = b_end - b_start
    # forward pass cumulative z-scores for this base
    cdef np.ndarray[DTYPE_t] b_fwd_data = np.empty(b_len, dtype=DTYPE)
    # store last diagonal move to pass on to next base
    cdef np.ndarray[DTYPE_INT_t] b_last_diag = np.empty(b_len, dtype=DTYPE_INT)
    # use cumsum as it is much more efficient than sums
    cdef np.ndarray[DTYPE_t] prev_b_data_cumsum = np.cumsum(prev_b_data)
    cdef DTYPE_INT_t pos, last_valid_diag_lag, pos_diag_val
    cdef DTYPE_t diag_score, stay_base_score, pos_score, fwd_value

    # add the diagonally below position value for the first possible
    # position in each base
    fwd_value = b_data[0] + prev_b_fwd_data[b_start - prev_b_start - 1]
    b_fwd_data[0] = fwd_value
    b_last_diag[0] = 1

    # some bases end at the same position (could change this by trimming earlier)
    for pos in range(b_start + 1, prev_b_end + 1):
        last_valid_diag_lag = 1
        while (prev_b_last_diag[pos - prev_b_start - last_valid_diag_lag] +
               last_valid_diag_lag <= min_obs_per_base):
            last_valid_diag_lag += 1

        # about 50% of the model_resqiuiggle algorithm is spent on this
        # part of this loop
        # sum of last valid fwd position and original z-scores through to the
        # position before current pos
        diag_score = prev_b_fwd_data[pos - prev_b_start - last_valid_diag_lag]
        if last_valid_diag_lag > 1:
            diag_score += prev_b_data_cumsum[pos - prev_b_start - 1] - \
                          prev_b_data_cumsum[pos - prev_b_start -
                                             last_valid_diag_lag]
        stay_base_score = b_fwd_data[pos - b_start - 1]
        if diag_score > stay_base_score:
            pos_score, pos_diag_val = diag_score, 1
        else:
            # stayed in this base, so add one to the last stayed in base count
            pos_score, pos_diag_val = (
                stay_base_score, b_last_diag[pos - b_start - 1] + 1)
        fwd_value = b_data[pos - b_start] + pos_score
        b_fwd_data[pos - b_start] = fwd_value
        b_last_diag[pos - b_start] = pos_diag_val

    cdef DTYPE_INT_t idx, curr_last_diag, reg_left_len
    if b_end > prev_b_end + 1:
        # perform C cumsum until the end of the base
        # note no possible allowed diagonal moves here
        fwd_value = b_fwd_data[prev_b_end - b_start]
        curr_last_diag = b_last_diag[prev_b_end - b_start]
        reg_left_len = b_end - prev_b_end - 1
        for idx in range(reg_left_len):
            fwd_value += b_data[idx + prev_b_end - b_start + 1]
            curr_last_diag += 1
            b_fwd_data[idx + prev_b_end - b_start + 1] = fwd_value
            b_last_diag[idx + prev_b_end - b_start + 1] = curr_last_diag

    return b_fwd_data, b_last_diag

def c_base_traceback(
        np.ndarray[DTYPE_t] curr_b_data not None, DTYPE_INT_t curr_start,
        np.ndarray[DTYPE_t] next_b_data not None,
        DTYPE_INT_t next_start, DTYPE_INT_t next_end,
        DTYPE_INT_t sig_start, DTYPE_INT_t min_obs_per_base):
    cdef DTYPE_INT_t curr_base_sig = 1
    cdef DTYPE_INT_t sig_pos
    for sig_pos in range(sig_start, -1, -1):
        curr_base_sig += 1
        # if there is not enough signal in the current base or the next base
        # hasn't started yet continue on to the next signal position
        if curr_base_sig <= min_obs_per_base or sig_pos - 1 >= next_end:
            continue
        # if the current base has ended or the next base has a better score
        if (sig_pos <= curr_start or
            next_b_data[sig_pos-next_start-1] >
            curr_b_data[sig_pos-curr_start-1]):
            return sig_pos


# Eventless re-squiggle dynamic programming algorithm
@cython.wraparound(False)
@cython.boundscheck(False)
def c_banded_forward_pass(
        np.ndarray[DTYPE_t, ndim=2] shifted_z_scores not None,
        np.ndarray[DTYPE_INT_t, ndim=1] event_starts not None,
        DTYPE_t skip_pen, DTYPE_t stay_pen):
    cdef DTYPE_INT_t n_bases = shifted_z_scores.shape[0]
    cdef DTYPE_INT_t bandwidth = shifted_z_scores.shape[1]
    cdef np.ndarray[DTYPE_t, ndim=2] fwd_pass = np.empty((
        n_bases + 1, bandwidth))
    cdef np.ndarray[DTYPE_INT_t, ndim=2] fwd_pass_tb = np.empty(
        (n_bases + 1, bandwidth), dtype=DTYPE_INT)
    # zero starts let the read start anywhere along the beginning
    # (for finding the read start)
    cdef DTYPE_INT_t idx
    for idx in range(bandwidth):
        fwd_pass[0, idx] = 0.0

    cdef DTYPE_INT_t max_from, band_pos, seq_pos, prev_b_pos
    cdef DTYPE_t max_score, pos_z_score, skip_score, diag_score

    for seq_pos in range(n_bases):
        # set first band position to skip score if the bands have the same start
        if seq_pos == 0 or event_starts[seq_pos] == event_starts[seq_pos-1]:
            fwd_pass[seq_pos + 1, 0] = fwd_pass[seq_pos, 0] - skip_pen
            fwd_pass_tb[seq_pos + 1, 0] = 1
        # else use the match score
        else:
            fwd_pass[seq_pos + 1, 0] = (
                fwd_pass[seq_pos, event_starts[seq_pos] -
                         event_starts[seq_pos-1] - 1] +
                shifted_z_scores[seq_pos, 0])
            fwd_pass_tb[seq_pos + 1, 0] = 2

        for band_pos in range(1, bandwidth):
            pos_z_score = shifted_z_scores[seq_pos, band_pos]
            prev_b_pos = (band_pos + event_starts[seq_pos] -
                          event_starts[seq_pos-1]
                          if seq_pos > 0 else band_pos)

            # first set to stay state
            max_score = fwd_pass[seq_pos+1, band_pos-1] - stay_pen + pos_z_score
            max_from = 0
            # then check diagonal score
            if prev_b_pos - 1 < bandwidth:
                diag_score = fwd_pass[seq_pos, prev_b_pos-1] + pos_z_score
                if diag_score > max_score:
                    max_score = diag_score
                    max_from = 2
                # finally check skip score (note nested check to save some ops)
                if prev_b_pos < bandwidth:
                    skip_score = fwd_pass[seq_pos, prev_b_pos] - skip_pen
                    if skip_score > max_score:
                        max_score = skip_score
                        max_from = 1

            fwd_pass[seq_pos + 1, band_pos] = max_score
            fwd_pass_tb[seq_pos + 1, band_pos] = max_from

    return fwd_pass, fwd_pass_tb

def c_banded_traceback(
        np.ndarray[DTYPE_INT_t, ndim=2] fwd_pass_tb not None,
        np.ndarray[DTYPE_INT_t] event_starts not None, DTYPE_INT_t band_pos,
        DTYPE_INT_t band_boundary_thresh=-1):
    # first row in fwd pass is a pseudo-row and does not represent a base
    cdef DTYPE_INT_t n_bases = fwd_pass_tb.shape[0] - 1
    cdef DTYPE_INT_t bandwidth = fwd_pass_tb.shape[1]
    cdef np.ndarray[DTYPE_INT_t] seq_poss = np.empty(
        n_bases + 1, dtype=DTYPE_INT)
    cdef DTYPE_INT_t curr_event_pos = band_pos + event_starts[n_bases - 1]
    # last position is the end of the current looking window which is the
    # passed value
    seq_poss[n_bases] = curr_event_pos + 1
    cdef DTYPE_INT_t curr_seq_pos
    for curr_seq_pos in range(n_bases, 0, -1):
        band_pos = curr_event_pos - event_starts[curr_seq_pos-1]
        # 0 indicates stay in the current base
        while fwd_pass_tb[curr_seq_pos, band_pos] == 0:
            band_pos -= 1
        # if diagonal position
        if fwd_pass_tb[curr_seq_pos, band_pos] == 2:
            band_pos -= 1
        if (band_boundary_thresh >= 0 and
            min(band_pos, bandwidth - band_pos - 1) < band_boundary_thresh):
            raise NotImplementedError, (
                'Read event to sequence alignment extends beyond --bandwidth')
        curr_event_pos = event_starts[curr_seq_pos-1] + band_pos
        seq_poss[curr_seq_pos-1] = curr_event_pos + 1

    return seq_poss

@cython.wraparound(False)
@cython.boundscheck(False)
def c_argmax(np.ndarray[DTYPE_t] vals):
    cdef DTYPE_t val
    cdef DTYPE_t max_val = vals[0]
    cdef DTYPE_INT_t pos
    cdef DTYPE_INT_t max_pos = 0

    for pos in range(1, vals.shape[0]):
        val = vals[pos]
        if val > max_val:
            max_val = val
            max_pos = pos
    return max_pos

@cython.wraparound(False)
@cython.boundscheck(False)
def c_adaptive_banded_forward_pass(
        np.ndarray[DTYPE_t, ndim=2] fwd_pass not None,
        np.ndarray[DTYPE_INT_t, ndim=2] fwd_pass_tb not None,
        np.ndarray[DTYPE_INT_t] event_starts not None,
        np.ndarray[DTYPE_t] event_means not None,
        np.ndarray[DTYPE_t] r_ref_means not None,
        np.ndarray[DTYPE_t] r_ref_sds not None,
        DTYPE_t z_shift, DTYPE_t skip_pen, DTYPE_t stay_pen,
        DTYPE_INT_t start_seq_pos, DTYPE_t mask_fill_z_score,
        bool return_z_scores=False):
    cdef DTYPE_INT_t n_bases = fwd_pass.shape[0] - 1
    cdef DTYPE_INT_t bandwidth = fwd_pass.shape[1]
    cdef DTYPE_INT_t half_bandwidth = bandwidth / 2
    cdef DTYPE_INT_t n_events = event_means.shape[0]

    cdef DTYPE_INT_t event_pos, seq_pos, prev_band_start, curr_band_start, \
        band_pos, prev_b_pos, max_from
    cdef DTYPE_t pos_z_score, ref_mean, ref_sd, max_score, skip_score, diag_score

    cdef np.ndarray[DTYPE_t] shifted_z_scores = np.empty(bandwidth, dtype=DTYPE)
    cdef np.ndarray[DTYPE_t, ndim=2] all_shifted_z_scores
    if return_z_scores:
        all_shifted_z_scores = np.empty((n_bases - start_seq_pos, bandwidth),
                                        dtype=DTYPE)
    for seq_pos in range(start_seq_pos, n_bases):
        # determine adaptive location for this sequence position
        prev_band_start = event_starts[seq_pos - 1]
        curr_band_start = prev_band_start + c_argmax(fwd_pass[seq_pos]) \
                          - half_bandwidth + 1
        if curr_band_start < prev_band_start:
            curr_band_start = prev_band_start
        if curr_band_start >= n_events:
            # if this isn't within one of the last sequence position
            # the read is forced to skip to the end and will likely
            # not end in a favorable alignment
            if seq_pos < n_bases - 2:
                raise NotImplementedError, (
                    'Adaptive signal to seqeunce alignment extended ' +
                    'beyond raw signal')
            curr_band_start = n_events - 1
        event_starts[seq_pos] = curr_band_start

        # get shifted z-scores for adaptive band
        ref_mean = r_ref_means[seq_pos]
        ref_sd = r_ref_sds[seq_pos]
        if curr_band_start + bandwidth <= n_events:
            for event_pos in range(curr_band_start, curr_band_start + bandwidth):
                pos_z_score = (event_means[event_pos] - ref_mean) / ref_sd
                if pos_z_score < 0:
                    pos_z_score = -pos_z_score
                shifted_z_scores[
                    event_pos - curr_band_start] = z_shift - pos_z_score
        else:
            # if this band extends beyond events array pad with low score
            # note in the unlikely event that the traceback path goes through
            # this region it will be clipped later
            for event_pos in range(curr_band_start, n_events):
                pos_z_score = (event_means[event_pos] - ref_mean) / ref_sd
                if pos_z_score < 0:
                    pos_z_score = -pos_z_score
                shifted_z_scores[
                    event_pos - curr_band_start] = z_shift - pos_z_score
            for event_pos in range(n_events - curr_band_start, bandwidth):
                shifted_z_scores[event_pos] = mask_fill_z_score
        if return_z_scores:
            all_shifted_z_scores[seq_pos - start_seq_pos,:] = shifted_z_scores

        # now perform dynamic programming fill for this seq position

        # set first band position to skip score if the bands have the same start
        if curr_band_start == prev_band_start:
            fwd_pass[seq_pos + 1, 0] = fwd_pass[seq_pos, 0] - skip_pen
            fwd_pass_tb[seq_pos + 1, 0] = 1
        # else use the match score
        else:
            fwd_pass[seq_pos + 1, 0] = fwd_pass[
                seq_pos, curr_band_start - prev_band_start - 1] + \
                shifted_z_scores[0]
            fwd_pass_tb[seq_pos + 1, 0] = 2

        # profiling shows that >60% of the time is spent here. Not
        # functionalized now due to function call overheads
        for band_pos in range(1, bandwidth):
            pos_z_score = shifted_z_scores[band_pos]
            prev_b_pos = band_pos + curr_band_start - prev_band_start

            # first set to stay state
            max_score = fwd_pass[seq_pos+1, band_pos-1] - stay_pen + pos_z_score
            max_from = 0
            # then check diagonal score
            if prev_b_pos - 1 < bandwidth:
                diag_score = fwd_pass[seq_pos, prev_b_pos-1] + pos_z_score
                if diag_score > max_score:
                    max_score = diag_score
                    max_from = 2
                # finally check skip score (note nested check to save some ops)
                if prev_b_pos < bandwidth:
                    skip_score = fwd_pass[seq_pos, prev_b_pos] - skip_pen
                    if skip_score > max_score:
                        max_score = skip_score
                        max_from = 1

            fwd_pass[seq_pos + 1, band_pos] = max_score
            fwd_pass_tb[seq_pos + 1, band_pos] = max_from

    if return_z_scores:
        return all_shifted_z_scores

    return
