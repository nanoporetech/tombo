from __future__ import division, unicode_literals, absolute_import

from builtins import int, range, dict, map, zip, object

import os
import io
import re
import sys
import queue
import random
import pkg_resources

# Future warning from cython in h5py
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
import h5py

import numpy as np
np.seterr(all='raise')

from tqdm import tqdm
from time import sleep
from copy import deepcopy
from operator import itemgetter
from scipy import stats, optimize
from collections import defaultdict
from scipy.spatial.distance import pdist
from multiprocessing import Process, Queue, Pipe
from numpy.lib.recfunctions import append_fields
from scipy.cluster.hierarchy import single, leaves_list
from itertools import repeat, product, count, combinations

if sys.version_info[0] > 2:
    unicode = str

# import tombo functions
from . import tombo_helper as th

from ._c_helper import (
    c_mean_std, c_apply_outlier_thresh, c_new_means, c_calc_llh_ratio,
    c_calc_llh_ratio_const_var, c_calc_scaled_llh_ratio_const_var,
    c_new_mean_stds, c_compute_running_pctl_diffs, c_compute_slopes)

from ._default_parameters import (
    SMALLEST_PVAL, STANDARD_MODELS, ALTERNATE_MODELS,
    MIN_KMER_OBS_TO_EST, ALT_EST_BATCH, MAX_KMER_OBS, NUM_DENS_POINTS,
    LLR_THRESH, SAMP_COMP_THRESH, DE_NOVO_THRESH, KERNEL_DENSITY_RANGE,
    ROC_PLOT_POINTS, NANOPOLISH_CENTRAL_POS, NUM_READS_FOR_SCALE,
    ROBUST_QUANTS, MAX_POINTS_FOR_THEIL_SEN, NUM_READS_TO_ADJUST_MODEL,
    OCLLHR_SCALE, OCLLHR_HEIGHT, OCLLHR_POWER, FM_OFFSET_DEFAULT,
    MOST_SIGNIF_NUM_BATCHES_DEFAULT, DNA_SAMP_TYPE, RNA_SAMP_TYPE,
    MEAN_PRIOR_CONST, SD_PRIOR_CONST, ALGN_PARAMS_TABLE, SEG_PARAMS_TABLE,
    MIN_EVENT_TO_SEQ_RATIO, COLLAPSE_RNA_STALLS, STALL_PARAMS,
    OUTLIER_THRESH,
    RNA_SCALE_NUM_EVENTS, RNA_SCALE_MAX_FRAC_EVENTS, USE_RNA_EVENT_SCALE)
DEFAULT_STALL_PARAMS=th.stallParams(**STALL_PARAMS)

# list of classes/functions to include in API
__all__ = [
    'TomboStats', 'ModelStats', 'LevelStats', 'PerReadStats',
    'TomboModel', 'AltModel', 'normalize_raw_signal', 'compute_base_means',
    'get_read_seg_score', 'calc_kmer_fitted_shift_scale',
    'load_resquiggle_parameters', 'compute_num_events']


VERBOSE = True

_PROFILE_SIGNIF = False
_PROFILE_SIGNIF_STATS_OUT = False
_PROFILE_SIGNIF_PER_READ = False
_PROFILE_EST_REF = False
_PROFILE_CENTER_REF = False
_PROFILE_ALT_EST = False
_PROFILE_MOTIF_ALT_EST = False

_DEBUG_EST_STD = False
_DEBUG_EST_BW = 0.05
_DEBUG_EST_NUM_KMER_SAVE = 500

STAT_BLOCKS_QUEUE_LIMIT = 5

DNA_BASES = ['A','C','G','T']

HALF_NORM_EXPECTED_VAL = stats.halfnorm.expect()
LOG10_E = np.log10(np.exp(1))

STANDARD_MODEL_NAME = 'standard'

SAMP_COMP_TXT = 'sample_compare'
DE_NOVO_TXT = 'de_novo'
ALT_MODEL_TXT = 'model_compare'
PER_READ_STATS = (SAMP_COMP_TXT, DE_NOVO_TXT, ALT_MODEL_TXT)

# TODO only ks-test currently implemented
KS_TEST_TXT = 'ks_test'
U_TEST_TXT = 'u_test'
T_TEST_TXT = 't_test'
KS_STAT_TEST_TXT = 'ks_stat_test'
U_STAT_TEST_TXT = 'u_stat_test'
T_STAT_TEST_TXT = 't_stat_test'
LEVEL_STATS_TXTS = (
    KS_TEST_TXT, U_TEST_TXT, T_TEST_TXT,
    KS_STAT_TEST_TXT, U_STAT_TEST_TXT, T_STAT_TEST_TXT)

ALT_MODEL_SEP_CHAR = '_'

NORM_TYPES = ('none', 'pA', 'pA_raw', 'median', 'robust_median',
              'median_const_scale')

# options specifying testing methods
# assume constant SD in model to save on computation
CONST_SD_MODEL = True

STAT_BLOCKS_H5_NAME = 'Statistic_Blocks'
MOST_SIGNIF_H5_NAME = 'Most_Significant_Stats'
COV_DAMP_COUNTS_H5_NAME = 'Cov_Damp_Counts'
COV_THRESH_H5_NAME = 'Cov_Threshold'

# turned off by default (and not accessible via command line so
# hardcoded for now)
DEFAULT_TRIM_RNA_PARAMS = th.trimRnaParams(
    moving_window_size=50, min_running_values=100,
    thresh_scale=0.7, max_raw_obs=40000)


#############################################
##### Pair-wise Distance and Clustering #####
#############################################

def transform_and_trim_stats(reg_stats, are_pvals, trim_value):
    if are_pvals:
        reg_stats = -np.log10(reg_stats)
        nan_r_stats = np.nan_to_num(reg_stats)
        reg_stats[nan_r_stats > trim_value] = trim_value
    else:
        nan_r_stats = np.nan_to_num(reg_stats)
        reg_stats[nan_r_stats > trim_value] = trim_value
        reg_stats[nan_r_stats < -trim_value] = -trim_value
    return reg_stats

def order_reads(reg_stats):
    """Compute order of reads based on log p-values or -log likelihood ratios
    """
    if reg_stats.shape[0] == 1:
        return [0]
    # get pairwise distances between reads
    # will get some empty slice means warnings, so suppress those
    #   (no seterr for this specific warning)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        r_dists = pdist(reg_stats, lambda u, v:
                        np.nanmean(np.sqrt(((u-v)**2))))
    r_dists[np.isnan(r_dists)] = np.nanmax(r_dists) + 1
    # then perform single/min linkage clustering and return the leaf order
    return leaves_list(single(r_dists))

def sliding_window_dist(sig_diffs1, sig_diffs2, slide_span, num_bases):
    """Compute distance over the minimum over a sliding window
    """
    return np.sqrt(min(np.sum(np.square(
        sig_diffs1[i1:i1+num_bases] - sig_diffs2[i2:i2+num_bases]))
                       for i1 in range((slide_span * 2) + 1)
                       for i2 in range((slide_span * 2) + 1)))

def euclidian_dist(sig_diffs1, sig_diffs2):
    """Compute Euclidean distance
    """
    return np.sqrt(np.sum(np.square(sig_diffs1 - sig_diffs2)))

def get_pairwise_dists(reg_sig_diffs, index_q, dists_q, slide_span=None):
    """Compute pairwise distances between a set of signal shifts
    """
    if slide_span > 0:
        num_bases=reg_sig_diffs[0].shape[0] - (slide_span * 2)
    while True:
        try:
            index = index_q.get(block=False)
        except queue.Empty:
            break

        if slide_span > 0:
            row_dists = np.array(
                [sliding_window_dist(
                    reg_sig_diffs[index], reg_sig_diffs[j],
                    slide_span, num_bases)
                 for j in range(0,index+1)] +
                [0 for _ in range(index+1, len(reg_sig_diffs))])
        else:
            row_dists = np.array(
                [euclidian_dist(reg_sig_diffs[index], reg_sig_diffs[j])
                 for j in range(0,index+1)] +
                [0 for _ in range(index+1, len(reg_sig_diffs))])
        dists_q.put((index, row_dists))

    return


##################################
###### Signal Normalization ######
##################################

def compute_base_means(all_raw_signal, base_starts):
    """Efficiently compute new base mean values from raw signal and base start
    positions

    Args:
        all_raw_signal (`np.array`): raw nanopore signal obervation values
        base_starts (`np.array::np.int32`): 0-based base start positions within
            raw signal

    Returns:
        `np.array::np.float64` containing base mean levels
    """
    return c_new_means(all_raw_signal.astype(np.float64), base_starts)

def get_scale_values_from_events(
        all_raw_signal, valid_cpts, outlier_thresh,
        num_events=None, max_frac_events=None):
    if num_events is not None or max_frac_events is not None:
        if (num_events is None or
            valid_cpts.shape[0] * max_frac_events < num_events):
            num_events = int(valid_cpts.shape[0] * max_frac_events)
        valid_cpts = valid_cpts[:num_events]
    event_means = compute_base_means(all_raw_signal, valid_cpts)
    read_med = np.median(event_means)
    read_mad = np.median(np.abs(event_means - read_med))
    lower_lim = -outlier_thresh
    upper_lim = outlier_thresh

    return th.scaleValues(
        shift=read_med, scale=read_mad,
        lower_lim=lower_lim, upper_lim=upper_lim, outlier_thresh=None)

def trim_rna(
        all_raw_signal, rsqgl_params, trim_rna_params=DEFAULT_TRIM_RNA_PARAMS):
    all_raw_signal = all_raw_signal[:trim_rna_params.max_raw_obs]
    num_events = np.int64(
        all_raw_signal.shape[0] // rsqgl_params.mean_obs_per_event)
    # get stdev over delta-mean events
    valid_cpts = th.valid_cpts_w_cap(
        all_raw_signal.astype(np.float64), rsqgl_params.min_obs_per_base,
        rsqgl_params.running_stat_width, num_events)
    _, window_sds = c_new_mean_stds(
        all_raw_signal.astype(np.float64), valid_cpts)

    # now moving window through this array
    n_windows = (window_sds.size - trim_rna_params.moving_window_size) + 1
    s_bytes = window_sds.strides[0]
    moving_window_sds = np.lib.stride_tricks.as_strided(
        window_sds,
        shape=(n_windows, trim_rna_params.moving_window_size),
        strides=(s_bytes, s_bytes)).mean(-1)
    thresh = moving_window_sds.mean() * trim_rna_params.thresh_scale

    n_windows = (moving_window_sds.size - trim_rna_params.min_running_values) + 1
    s_bytes = moving_window_sds.strides[0]
    running_mins = np.lib.stride_tricks.as_strided(
        moving_window_sds,
        shape=(n_windows, trim_rna_params.min_running_values),
        strides=(s_bytes, s_bytes)).min(-1)
    try:
        pos_index = next(i for i, v in enumerate(running_mins) if v > thresh)
    except StopIteration:
        return 0

    return valid_cpts[pos_index]

def identify_stalls(all_raw_signal, stall_params, return_metric=False):
    """Identify locations where bases have stalled in the pore. Two methods
    availble depending on parameters specified in stall_params.
    """
    def compute_running_mean_diffs():
        """Compute average difference between n_window neighboring window means
        each of size window_size.
        """
        moving_average = np.cumsum(all_raw_signal)
        moving_average[stall_params.mini_window_size:] = (
            moving_average[stall_params.mini_window_size:] -
            moving_average[:-stall_params.mini_window_size])
        moving_average = moving_average[
            stall_params.mini_window_size - 1:] / stall_params.mini_window_size

        # extract moving window averages at n_window offsets
        offsets = [moving_average[
            int(stall_params.mini_window_size * offset):
            int(-stall_params.mini_window_size * (
                stall_params.n_windows - offset - 1))]
                   for offset in range(stall_params.n_windows - 1)] + [
                           moving_average[int(stall_params.mini_window_size * (
                               stall_params.n_windows - 1)):],]
        # compute difference between all pairwise offset
        diffs = [np.abs(offsets[i] - offsets[j])
                 for i in range(stall_params.n_windows)
                 for j in range(i + 1, stall_params.n_windows)]

        # compute average over offset differences at each valid position
        diff_sums = diffs[0].copy()
        for diff_i in diffs:
            diff_sums += diff_i
        return diff_sums / len(diffs)


    # if the raw signal is too short to compute stall metrics
    if all_raw_signal.shape[0] < stall_params.window_size:
        if return_metric:
            return [], np.repeat(np.NAN, all_raw_signal.shape[0])
        return []

    # identify potentially stalled signal from either running window means
    # or running percentile difference methods
    stall_metric = np.empty(all_raw_signal.shape, all_raw_signal.dtype)
    stall_metric[:] = np.NAN
    start_offset = int(stall_params.window_size * 0.5)
    end_offset = (all_raw_signal.shape[0] - stall_params.window_size
                  + start_offset + 1)
    if (stall_params.lower_pctl is not None and
        stall_params.upper_pctl is not None):
        stall_metric[start_offset:end_offset] = c_compute_running_pctl_diffs(
            all_raw_signal, np.int64(stall_params.window_size),
            np.float64(stall_params.lower_pctl),
            np.float64(stall_params.upper_pctl))
    elif (stall_params.n_windows is not None and
          stall_params.mini_window_size is not None):
        assert (stall_params.window_size ==
                stall_params.mini_window_size * stall_params.n_windows)
        stall_metric[start_offset:end_offset] = compute_running_mean_diffs()
    else:
        th.TomboError(
            'Must provide method specific parameters for stall detection')

    # identify contiguous windows over threshold for minimal stretches
    with np.errstate(invalid='ignore'):
        stall_locs = np.where(np.diff(np.concatenate(
            [[False], stall_metric <= stall_params.threshold])))[0]
    if stall_metric[-1] <= stall_params.threshold:
        stall_locs = np.concatenate([stall_locs, [stall_metric.shape[0]]])
    stall_locs = stall_locs.reshape(-1,2)
    stall_locs = stall_locs[(np.diff(stall_locs) >
                             stall_params.min_consecutive_obs).flatten()]
    if stall_locs.shape[0] == 0:
        if return_metric:
            return [], stall_metric
        return []

    # expand windows out to region that gave result below threshold
    # since windows are centered (minus edge buffer)
    expand_width = (stall_params.window_size // 2) - stall_params.edge_buffer
    if expand_width > 0:
        stall_locs[:,0] -= expand_width
        stall_locs[:,1] += expand_width
        # collapse intervals that now overlap
        merged_stall_locs = []
        prev_int = stall_locs[0]
        for curr_int in stall_locs:
            if curr_int[0] > prev_int[1]:
                # add previous interval to all intervals
                merged_stall_locs.append(prev_int)
                prev_int = curr_int
            else:
                # extend previous interval since these overlap
                prev_int[1] = curr_int[1]
        merged_stall_locs.append(prev_int)
        stall_locs = merged_stall_locs

    if return_metric:
        return stall_locs, stall_metric
    return stall_locs

def calc_kmer_fitted_shift_scale(
        prev_shift, prev_scale, r_event_means, r_model_means,
        r_model_inv_vars=None, method='theil_sen'):
    """Use robust Theil-Sen estimator to compute fitted shift and scale
    parameters based on read sequence

    Args:
        prev_shift (float): previous shift parameter
        prev_scale (float): previous scale parameter
        r_ref_means (`np.array::np.float64`): expected base signal levels
        r_ref_sds (`np.array::np.float64`): expected base signal level sds
        r_model_inv_vars (`np.array::np.float64`): expected base signal level
            inverse variances for method of moments (`mom`) computation
        method (str): one of `theil_sen`, `robust`, or `mom`

    Returns:
        Sequence-fitted scaling parameters

        1) shift parameter (float)
        2) scale parameter (float)
        3) shift correction factor; multiply by ``prev_scale`` and add to ``prev_shift`` to get ``shift`` (float)
        4) scale correction factor; multiply by ``prev_scale`` to get ``scale`` (float)
    """
    if method == 'robust':
        def read_lad_objective(x):
            return np.sum(np.abs(((r_event_means - x[0]) / x[1]) -
                                 r_model_means))

        shift_corr_factor, scale_corr_factor = optimize.minimize(
            read_lad_objective, np.array([0,1]), method='nelder-mead',
            options={'xtol': 1e-8}).x
    elif method == 'theil_sen':
        def compute_slopes(r_event_means, r_model_means):
            # despite computing each diff twice this vectorized solution is
            # about 10X faster than a list comprehension approach
            delta_event = r_event_means[:, np.newaxis] - r_event_means
            delta_model = r_model_means[:, np.newaxis] - r_model_means
            return delta_model[delta_event > 0] / delta_event[delta_event > 0]

        n_points = r_model_means.shape[0]
        # potentially sample points for long reads (>1kb)
        if r_model_means.shape[0] > MAX_POINTS_FOR_THEIL_SEN:
            n_points = MAX_POINTS_FOR_THEIL_SEN
            samp_ind = np.random.choice(
                r_model_means.shape[0], n_points, replace=False)
            r_model_means = r_model_means[samp_ind]
            r_event_means = r_event_means[samp_ind]
        # compute Theil-Sen slope estimator
        slope = np.median(c_compute_slopes(r_event_means, r_model_means))
        inter = np.median(r_model_means - (slope * r_event_means))
        if slope == 0:
            raise th.TomboError('Read failed sequence-based signal ' +
                                're-scaling parameter estimation.')
        # convert to shift and scale parameters (e.g. (obs - shift) / scale)
        scale_corr_factor = 1 / slope
        shift_corr_factor = -inter / slope
    elif method == 'mom':
        model_mean_var = r_model_means * r_model_inv_vars
        # prep kmer model coefficient matrix for the k-mers from this read
        model_mean_var_sum = model_mean_var.sum()
        coef_mat = np.array((
            (r_model_inv_vars.sum(), model_mean_var_sum),
            (model_mean_var_sum, (model_mean_var * r_model_means).sum())))

        # prep dependent values from this reads true events
        r_event_var = r_event_means * r_model_inv_vars
        r_event_var_mean = r_event_var * r_model_means
        dep_vect = np.array((r_event_var.sum(), r_event_var_mean.sum()))

        shift_corr_factor, scale_corr_factor = np.linalg.solve(
            coef_mat, dep_vect)
    else:
        th.error_message_and_exit(
            'Invalid k-mer fitted normalization parameter method: ' + method +
            '\n\t\tValid methods are "robust" and "mom".')

    # apply shift and scale values fitted from kmer conditional model
    shift = prev_shift + (shift_corr_factor * prev_scale)
    scale = prev_scale * scale_corr_factor

    return shift, scale, shift_corr_factor, scale_corr_factor

def estimate_global_scale(fast5_fns, num_reads=NUM_READS_FOR_SCALE):
    if VERBOSE: th.status_message('Estimating global scale parameter.')
    np.random.shuffle(fast5_fns)
    read_mads = []
    if VERBOSE:
        bar = tqdm(total=num_reads, desc='Total reads processed', smoothing=0)
    for fast5_fn in fast5_fns:
        try:
            with h5py.File(fast5_fn, 'r') as fast5_data:
                all_sig = th.get_raw_read_slot(fast5_data)['Signal'][:]
            shift = np.median(all_sig)
            read_mads.append(np.median(np.abs(all_sig - shift)))
            if VERBOSE: bar.update(1)
        except:
            continue
        if len(read_mads) >= num_reads:
            break

    if VERBOSE: bar.close()
    if len(read_mads) == 0:
        th.error_message_and_exit(
            'No reads contain raw signal for ' +
            'global scale parameter estimation.')
    if len(read_mads) < num_reads:
        th.warning_message(
            'Few reads contain raw signal for global scale parameter ' +
            'estimation. Results may not be optimal.')

    return np.mean(read_mads)

def normalize_raw_signal(
        all_raw_signal, read_start_rel_to_raw=0, read_obs_len=None,
        norm_type='median', outlier_thresh=None, channel_info=None,
        scale_values=None, event_means=None, model_means=None,
        model_inv_vars=None, const_scale=None):
    """Apply scaling and windsorizing parameters to normalize raw signal.

    Args:
        all_raw_signal (`np.array`): raw nanopore signal obervation values
        read_start_rel_to_raw (int): amount of signal to trim from beginning of
            the signal (default: 0)
        read_obs_len (int): length of signal to process from
            `read_start_rel_to_raw` (default: full length)
        norm_type (str): normalization type (`median` (default), `none`,
            `pA_raw`, `pA`, `median_const_scale`, `robust_median`; ignored is
            ``scale_values`` provided)
        outlier_thresh (float): windsorizing threshold (MAD units; default: None)
        channel_info (:class:`tombo.tombo_helper.channelInfo`): channel
            information (optional; only for `pA` and `pA_raw`)
        scale_values (:class:`tombo.tombo_helper.scaleValues`): scaling values
            (optional)
        event_means (`np.array`): for `pA` fitted scaling parameters (optional)
        model_means (`np.array`): for `pA` fitted scaling parameters (optional)
        model_inv_vars (`np.array`): for `pA` fitted scaling parameters
            (optional)
        const_scale (float): global scale parameter (optional)

    Returns:
        Normalized signal and scaling parameters

        1) normalized signal observations (`np.array::np.float64`)
        2) :class:`tombo.tombo_helper.scaleValues`
    """
    if read_obs_len is None:
        read_obs_len = all_raw_signal.shape[0] - read_start_rel_to_raw
    if norm_type not in NORM_TYPES and scale_values is None:
        raise th.TomboError(
            'Normalization type ' + norm_type + ' is not a valid ' +
            'option and shift or scale parameters were not provided.')

    raw_signal = all_raw_signal[read_start_rel_to_raw:
                                read_start_rel_to_raw + read_obs_len]
    if scale_values is None:
        if norm_type == 'none':
            shift, scale = 0, 1
        elif norm_type in ('pA_raw', 'pA'):
            # correct raw signal as described here:
            #     https://community.nanoporetech.com
            #           /posts/squiggle-plot-for-raw-data
            shift, scale = (
                -1 * channel_info.offset,
                channel_info.digitisation / channel_info.range)
            if norm_type == 'pA':
                # perform k-mer model fitted correction as in
                # nanocorr/nanopolish/albacore(pre-RNN)
                shift, scale, _, _ = calc_kmer_fitted_shift_scale(
                    shift, scale, event_means, model_means, model_inv_vars,
                    method='mom')
        elif norm_type == 'median':
            shift = np.median(raw_signal)
            scale = np.median(np.abs(raw_signal - shift))
        elif norm_type == 'median_const_scale':
            assert const_scale is not None
            shift = np.median(raw_signal)
            scale = const_scale
        elif norm_type == 'robust_median':
            shift = np.mean(np.percentile(raw_signal, ROBUST_QUANTS))
            scale = np.median(np.abs(raw_signal - shift))
    else:
        shift = scale_values.shift
        scale = scale_values.scale

    norm_signal = (raw_signal - shift) / scale

    # windsorize the raw signal
    lower_lim, upper_lim = None, None
    if outlier_thresh is not None or scale_values is not None:
        if outlier_thresh is not None:
            read_med = np.median(norm_signal)
            read_mad = np.median(np.abs(norm_signal - read_med))
            lower_lim = read_med - (read_mad * outlier_thresh)
            upper_lim = read_med + (read_mad * outlier_thresh)
        else:
            lower_lim = scale_values.lower_lim
            upper_lim = scale_values.upper_lim
        # provided scale values could contain None winsorizing limits still
        if lower_lim is not None and upper_lim is not None:
            norm_signal = c_apply_outlier_thresh(
                norm_signal, lower_lim, upper_lim)

    return norm_signal, th.scaleValues(
        shift, scale, lower_lim, upper_lim, outlier_thresh)


#############################
##### Tombo Model Class #####
#############################

class TomboModel(object):
    """Load, store and access Tombo model attributes and sequence-based
            expected mean and standard deviation levels (median normalization
            only).

    .. automethod:: __init__
    """
    def _center_model(self, shift_corr_factor, scale_corr_factor):
        centered_means = {}
        for kmer, k_mean in self.means.items():
            centered_means[kmer] = (
                k_mean * scale_corr_factor) + shift_corr_factor

        self.means = centered_means

        return

    def _make_constant_sd(self):
        med_sd = np.median(list(self.sds.values()))
        self.sds = dict((kmer, med_sd) for kmer in self.sds)
        return

    def write_model(self, ref_fn):
        """Write TomboModel to specified file

        Args:
            ref_fn (str): filename to write TomboModel
        """
        # Explicity use btype string names for py3 compatiblity as well as
        # pickle-ability of numpy arrays for consistency. See discussion here:
        # https://github.com/numpy/numpy/issues/2407
        ref_for_file = np.array(
            [(kmer, self.means[kmer], self.sds[kmer]) for kmer in self.means],
            dtype=[(str('kmer'), 'S' + unicode(self.kmer_width)),
                   (str('mean'), 'f8'), (str('sd'), 'f8')])

        with h5py.File(ref_fn, 'w') as ref_fp:
            ref_fp.create_dataset('model', data=ref_for_file, compression="gzip")
            ref_fp.attrs['central_pos'] = self.central_pos
            ref_fp.attrs['model_name'] = STANDARD_MODEL_NAME

        return

    def _parse_tombo_model(self):
        """Parse a tombo model file
        """
        try:
            with h5py.File(self.ref_fn, 'r') as ref_fp:
                ref_raw = ref_fp['model'][:]
                central_pos = ref_fp.attrs.get('central_pos')
                model_name = ref_fp.attrs.get('model_name')

                try:
                    model_name = model_name.decode()
                except (AttributeError, TypeError):
                    pass

        except:
            th.error_message_and_exit(
                'Invalid Tombo model file provided: ' + unicode(self.ref_fn))

        if model_name != STANDARD_MODEL_NAME:
            th.error_message_and_exit('Non-standard Tombo model provided.')

        mean_ref = {}
        sd_ref = {}
        for kmer, kmer_mean, kmer_std in ref_raw:
            kmer = kmer.decode()
            mean_ref[kmer] = kmer_mean
            sd_ref[kmer] = kmer_std

        self.means = mean_ref
        self.sds = sd_ref
        self.central_pos = central_pos
        self.name = model_name

        return

    def _parse_text_model(self):
        """
        Parse a text model file (such as those from nanopolish)
        """
        try:
            mean_ref, sd_ref = {}, {}
            with io.open(self.ref_fn) as fp:
                for line in fp:
                    if line.startswith('#'): continue
                    try:
                        kmer, kmer_mean, kmer_sd = line.split()[:3]
                        kmer_mean, kmer_sd = map(float, (kmer_mean, kmer_sd))
                    except ValueError:
                        # header or other non-kmer field
                        continue
                    mean_ref[kmer] = kmer_mean
                    sd_ref[kmer] = kmer_sd
        except:
            th.error_message_and_exit('Invalid text pA model file provided: '
                                       + unicode(self.ref_fn))

        self.means = mean_ref
        self.sds = sd_ref
        self.central_pos = NANOPOLISH_CENTRAL_POS
        self.name = STANDARD_MODEL_NAME

        return

    def _load_model(self, kmer_ref, central_pos):
        mean_ref = {}
        sd_ref = {}
        for kmer, kmer_mean, kmer_std in kmer_ref:
            # reference may or may not be stored as a numpy array
            try:
                kmer = kmer.decode()
            except AttributeError:
                pass
            mean_ref[kmer] = kmer_mean
            sd_ref[kmer] = kmer_std

        self.means = mean_ref
        self.sds = sd_ref
        self.central_pos = central_pos
        self.name = STANDARD_MODEL_NAME

        return

    def _add_invvar(self):
        self.inv_var = {}
        for kmer, stdev in self.sds.items():
            self.inv_var[kmer] = 1 / (stdev * stdev)

        return

    def _get_default_standard_ref(self, reads_index):
        if th.is_sample_rna(reads_index=reads_index):
            if VERBOSE: th.status_message(
                    'Loading default canonical ***** RNA ***** model.')
            std_ref_fn = STANDARD_MODELS[RNA_SAMP_TYPE]
            self.seq_samp_type = th.seqSampleType(RNA_SAMP_TYPE, True)
        else:
            if VERBOSE: th.status_message(
                    'Loading default canonical ***** DNA ***** model.')
            self.seq_samp_type = th.seqSampleType(DNA_SAMP_TYPE, False)
            std_ref_fn = STANDARD_MODELS[DNA_SAMP_TYPE]
        # get full filename path with setuptools
        self.ref_fn = th.resolve_path(pkg_resources.resource_filename(
            'tombo', 'tombo_models/' + std_ref_fn))

        return

    def _get_default_standard_ref_from_files(self, fast5_fns):
        if th.is_sample_rna(fast5_fns=fast5_fns):
            if VERBOSE: th.status_message(
                    'Loading default canonical ***** RNA ***** model.')
            std_ref_fn = STANDARD_MODELS[RNA_SAMP_TYPE]
            self.seq_samp_type = th.seqSampleType(RNA_SAMP_TYPE, True)
        else:
            if VERBOSE: th.status_message(
                    'Loading default canonical ***** DNA ***** model.')
            std_ref_fn = STANDARD_MODELS[DNA_SAMP_TYPE]
            self.seq_samp_type = th.seqSampleType(DNA_SAMP_TYPE, False)
        # get full filename path with setuptools
        self.ref_fn = th.resolve_path(pkg_resources.resource_filename(
            'tombo', 'tombo_models/' + std_ref_fn))

        return

    def __init__(
            self, ref_fn=None, is_text_model=False, kmer_ref=None,
            central_pos=None, seq_samp_type=None, reads_index=None,
            fast5_fns=None, minimal_startup=True):
        """Initialize a Tombo k-mer model object

        Args:
            ref_fn (str): tombo model filename
            is_text_model (bool): `ref_fn` is text (e.g.
                https://github.com/nanoporetech/kmer_models/blob/master/r9.4_180mv_450bps_6mer/template_median68pA.model)
            kmer_ref (list): containing 3-tuples 1) k-mer 2) expected level 3)
                level SD
            central_pos (int): base within k-mer to assign signal (only
                applicable when `kmer_ref` is provided)
            seq_samp_type (:class:`tombo.tombo_helper.seqSampleType`):
                sequencing sample type (default: None)
            reads_index (:class:`tombo.tombo_helper.TomboReads`): For
                determining `seq_samp_type`
            fast5_fns (list): fast5 read filenames from which to extract read
                metadata. For determining `seq_samp_type`
            minimal_startup (bool): don't compute inverse variances (default
                True)

        Note:

            Order of priority for initialization when multiple model
                specifications are provided:
                1) `ref_fn`
                2) `kmer_ref` (requires `central_pos`)
                3) `seq_samp_type`
                4) `reads_index`
                5) `fast5_fns`

            Last 3 options load a default model file included with Tombo. Last
                2 determine the sample type from read metadata.
        """
        if ref_fn is not None:
            self.ref_fn = th.resolve_path(ref_fn)
            if is_text_model:
                self._parse_text_model()
            else:
                self._parse_tombo_model()
            self.seq_samp_type = seq_samp_type
        elif kmer_ref is not None:
            assert central_pos is not None, (
                'central_pos must be provided is TomboModel is loaded ' +
                'with a kmer_ref')
            self._load_model(kmer_ref, central_pos)
            self.seq_samp_type = seq_samp_type
        else:
            if seq_samp_type is not None:
                self.seq_samp_type = seq_samp_type
                self.ref_fn = th.resolve_path(pkg_resources.resource_filename(
                    'tombo', 'tombo_models/' + STANDARD_MODELS[
                        seq_samp_type.name]))
            elif reads_index is not None:
                self._get_default_standard_ref(reads_index)
            elif fast5_fns is not None:
                self._get_default_standard_ref_from_files(fast5_fns)
            else:
                th.error_message_and_exit(
                    'Must provide initialization method for TomboModel.')
            self._parse_tombo_model()

        self.kmer_width = len(next(k for k in self.means))

        self.inv_var = None
        if not minimal_startup:
            self._add_invvar()

    def reverse_sequence_copy(self):
        """Return a copy of model for processing sequence/signal in reverse
        (default models are all saved in genome sequence forward (5p to 3p)
        direction)
        """
        rev_model = deepcopy(self)
        rev_model.central_pos = self.kmer_width - self.central_pos - 1
        rev_model.means = dict((kmer[::-1], kmer_mean)
                               for kmer, kmer_mean in self.means.items())
        rev_model.sds = dict((kmer[::-1], kmer_sds)
                               for kmer, kmer_sds in self.sds.items())
        if self.inv_var is not None:
            rev_model.inv_var = dict(
                (kmer[::-1], kmer_inv_var)
                for kmer, kmer_inv_var in self.inv_var.items())

        return rev_model

    def get_exp_levels_from_seq(self, seq, rev_strand=False):
        """Compute expected signal levels for a sequence

        Args:
            seq (str): genomic seqeunce to be converted to expected signal
                levels
            rev_strand (bool): provided sequence is from reverse strand (so
                flip return order to genome forward direction)

        Note:
            Returned expected signal levels will be trimmed compared to the
            passed sequence based on the `std_ref.kmer_width` and
            `std_ref.central_pos`.

        Returns:
            Expected signal level references
            1) ref_means (`np.array::np.float64`) expected signal levels
            2) ref_sds (`np.array::np.float64`) expected signal level sds
        """
        seq_kmers = th.get_seq_kmers(seq, self.kmer_width, rev_strand)

        try:
            ref_means = np.array([self.means[kmer] for kmer in seq_kmers])
            ref_sds = np.array([self.sds[kmer] for kmer in seq_kmers])
        except KeyError:
            th.error_message_and_exit(
                'Invalid sequence encountered from genome sequence.')

        return ref_means, ref_sds

    def get_exp_levels_from_kmers(self, seq_kmers):
        """Look up expected signal levels for a list of k-mers

        Args:
            seq_kmers (list): list of k-mers (as returned from
                :class:`tombo.tombo_helper.get_seq_kmers`)

        Returns:
            Expected signal level references

            1) ref_means (`np.array::np.float64`) expected signal levels
            2) ref_sds (`np.array::np.float64`) expected signal level sds
        """
        try:
            ref_means = np.array([self.means[kmer] for kmer in seq_kmers])
            ref_sds = np.array([self.sds[kmer] for kmer in seq_kmers])
        except KeyError:
            th.error_message_and_exit(
                'Invalid sequence encountered from genome sequence.')

        return ref_means, ref_sds

    def get_exp_levels_from_seq_with_gaps(self, reg_seq, rev_strand):
        # loop over regions without valid sequence (non-ACGT)
        reg_ref_means, reg_ref_sds = (
            np.empty(len(reg_seq) - self.kmer_width + 1),
            np.empty(len(reg_seq) - self.kmer_width + 1))
        reg_ref_means[:] = np.NAN
        reg_ref_sds[:] = np.NAN
        prev_ibr_end = 0
        for inv_base_run_m in th.INVALID_BASE_RUNS.finditer(reg_seq):
            ibr_start, ibr_end = inv_base_run_m.start(), inv_base_run_m.end()
            # if valid region is too short continue
            if ibr_start - prev_ibr_end < self.kmer_width:
                prev_ibr_end = ibr_end
                continue
            subreg_ref_means, subreg_ref_sds = self.get_exp_levels_from_seq(
                reg_seq[prev_ibr_end:ibr_start])
            reg_ref_means[prev_ibr_end:
                          ibr_start - self.kmer_width + 1] = subreg_ref_means
            reg_ref_sds[prev_ibr_end:
                        ibr_start - self.kmer_width + 1] = subreg_ref_sds
            prev_ibr_end = ibr_end

        # if there is valid sequence at the end of a region include it here
        if prev_ibr_end <= len(reg_seq) - self.kmer_width:
            subreg_ref_means, subreg_ref_sds = self.get_exp_levels_from_seq(
                reg_seq[prev_ibr_end:])
            reg_ref_means[prev_ibr_end:] = subreg_ref_means
            reg_ref_sds[prev_ibr_end:] = subreg_ref_sds

        if rev_strand:
            reg_ref_means = reg_ref_means[::-1]
            reg_ref_sds = reg_ref_sds[::-1]

        return reg_ref_means, reg_ref_sds


class AltModel(object):
    """Load, store and access alternate-base Tombo model attributes and
    sequence-based expected mean and standard deviation levels (median
    normalization only).

    .. automethod:: __init__
    """
    def write_model(self, ref_fn):
        """Write AltModel to specified file

        Args:
            ref_fn (str): filename to write AltModel
        """
        # Explicity use btype string names for py3 compatiblity as well as
        # pickle-ability of numpy arrays for consistency. See discussion here:
        # https://github.com/numpy/numpy/issues/2407
        ref_for_file = np.array(
            [(kmer, pos, self.means[(kmer, pos)], self.sds[(kmer, pos)])
             for kmer, pos in self.means],
            dtype=[(str('kmer'), 'S' + unicode(self.kmer_width)),
                   (str('pos'), 'u4'), (str('mean'), 'f8'), (str('sd'), 'f8')])

        with h5py.File(ref_fn, 'w') as ref_fp:
            ref_fp.create_dataset(
                'model', data=ref_for_file, compression="gzip")
            ref_fp.attrs['central_pos'] = self.central_pos
            ref_fp.attrs['model_name'] = self.name
            ref_fp.attrs['alt_base'] = self.alt_base
            ref_fp.attrs['motif'] = self.motif.raw_motif
            ref_fp.attrs['mod_pos'] = self.motif.mod_pos

        return

    def _make_constant_sd(self):
        med_sd = np.median(list(self.sds.values()))
        self.sds = dict((kmer, med_sd) for kmer in self.sds)
        return

    def _parse_alt_model(self):
        """Parse an alternate-base Tombo model file
        """
        try:
            with h5py.File(self.ref_fn, 'r') as ref_fp:
                ref_raw = ref_fp['model'][:]
                central_pos = ref_fp.attrs.get('central_pos')
                model_name = ref_fp.attrs.get('model_name')

                try:
                    model_name = model_name.decode()
                except (AttributeError, TypeError):
                    pass

                alt_base = ref_fp.attrs.get('alt_base')
                try:
                    alt_base = alt_base.decode()
                except (AttributeError, TypeError):
                    pass

                raw_motif = ref_fp.attrs.get('motif')
                try:
                    raw_motif = raw_motif.decode()
                except (AttributeError, TypeError):
                    pass
                mod_pos = ref_fp.attrs.get('mod_pos')

        except:
            th.error_message_and_exit(
                'Invalid alternate-base tombo model file provided: '
                + unicode(self.ref_fn))

        mean_ref = {}
        sd_ref = {}
        for kmer, pos, kmer_mean, kmer_std in ref_raw:
            kmer = kmer.decode()
            mean_ref[(kmer, pos)] = kmer_mean
            sd_ref[(kmer, pos)] = kmer_std

        self.means = mean_ref
        self.sds = sd_ref
        self.central_pos = central_pos
        self.alt_base = alt_base
        self.name = model_name
        self.motif = th.TomboMotif(raw_motif, mod_pos)

        return

    def _load_alt_model(self, alt_ref, central_pos, alt_base, name, motif):
        mean_ref = {}
        sd_ref = {}
        for kmer, pos, kmer_mean, kmer_std in alt_ref:
            # reference may or may not be stored as a numpy array
            try:
                kmer = kmer.decode()
            except AttributeError:
                pass
            mean_ref[(kmer, pos)] = kmer_mean
            sd_ref[(kmer, pos)] = kmer_std

        self.means = mean_ref
        self.sds = sd_ref
        self.central_pos = central_pos
        self.alt_base = alt_base
        self.name = name
        if motif is None:
            self.motif = th.TomboMotif(self.alt_base, 1)
        else:
            assert (isinstance(motif, th.TomboMotif) and
                    motif.mod_pos is not None)
            self.motif = motif

        return

    def _add_invvar(self):
        self.inv_var = {}
        for kmer_pos, stdev in self.sds.items():
            self.inv_var[kmer_pos] = 1 / (stdev * stdev)

        return

    def __init__(
            self, ref_fn=None, kmer_ref=None, central_pos=None, alt_base=None,
            name=None, motif=None, minimal_startup=True):
        """Initialize a Tombo k-mer model object

        Args:
            ref_fn (str): tombo model filename
            kmer_ref (list): containing 3-tuples 1) k-mer 2) expected level 3)
                level SD
            central_pos (int): base within k-mer to assign signal (only
                applicable when `kmer_ref` is provided)
            alt_base (int): "swap" base within k-mer (only applicable when
                `kmer_ref` is provided)
            name (str): model name (only applicable when `kmer_ref` is provided)
            motif (:class:`tombo.tombo_helper.TomboMotif`): model motif
                (defaults to alt_base motif; only applicable when `kmer_ref` is
                provided)
            minimal_startup (bool): don't compute inverse variances (default
                True)

        Note:

            ``kmer_ref``, ``central_pos``, ``alt_base``, ``name`` and ``motif``
                are ignored if ``ref_fn`` is provided.
        """
        if ref_fn is not None:
            self.ref_fn = th.resolve_path(ref_fn)
            self._parse_alt_model()
        elif kmer_ref is not None:
            assert central_pos is not None and alt_base is not None, (
                'central_pos and alt_base must be provided if AltModel is ' +
                'loaded with a kmer_ref')
            self._load_alt_model(kmer_ref, central_pos, alt_base, name, motif)
        else:
            th.error_message_and_exit(
                'Must provide initialization method for AltModel.')

        self.kmer_width = len(next(kmer for kmer, pos in self.means))

        self.inv_var = None
        if not minimal_startup:
            self._add_invvar()

    def get_exp_level(self, kmer, pos):
        try:
            return self.means[(kmer, pos)]
        except KeyError:
            return np.NAN

    def get_exp_sd(self, kmer, pos):
        try:
            return self.sds[(kmer, pos)]
        except KeyError:
            return np.NAN

    def get_exp_levels_from_kmers(self, seq_kmers, rev_strand=False):
        """Look up expected alternate-base signal levels across a central base.
        Alternative base to test should be last base in the first k-mer and
        first base in the last k-mer

        Args:
            seq_kmers (list): list of k-mers the same length as the k-mer

        Returns:
            Expected signal level references

            1) ref_means (`np.array::np.float64`) expected signal levels
            2) ref_sds (`np.array::np.float64`) expected signal level sds
        """
        pos_range = (range(self.kmer_width) if rev_strand else
                     range(self.kmer_width - 1, -1, -1))
        try:
            ref_means = np.array([
                self.get_exp_level(kmer, pos) for kmer, pos in zip(
                    seq_kmers, pos_range)])
            ref_sds = np.array([
                self.get_exp_sd(kmer, pos) for kmer, pos in zip(
                    seq_kmers, pos_range)])
        except KeyError:
            th.error_message_and_exit(
                'Invalid sequence encountered from genome sequence.')

        return ref_means, ref_sds


############################
##### Model Estimation #####
############################

def check_valid_alt_models(alt_refs, std_ref):
    """Parse several alternative tombo model files
    """
    for alt_ref in alt_refs.values():
        if (std_ref.central_pos != alt_ref.central_pos or
            std_ref.kmer_width != alt_ref.kmer_width):
            th.warning_message(
                'Standard and ' + alt_ref.ref_fn + ' alternative base ' +
                'models must be estimated using the same k-mer positions.')
            continue
        if not isinstance(alt_ref, AltModel):
            th.warning_message(
                'Alternative model ' + alt_ref.ref_fn + ' appears to be a ' +
                'standard model and will not be processed.')
            continue

    return alt_refs

def _print_alt_models():
    alt_model_types = [tuple(mod_name.split(ALT_MODEL_SEP_CHAR))
                       for mod_name in ALTERNATE_MODELS.keys()]
    alt_seq_samps = ['',] + sorted(set(list(zip(*alt_model_types))[0]))
    alt_mods = list(set(list(zip(*alt_model_types))[1]))
    row_format ="{:<10}" * (len(alt_seq_samps)) + '\n'
    sys.stderr.write(row_format.format(*alt_seq_samps))
    for alt_mod in alt_mods:
        has_mod = [alt_mod,]
        for seq_samp in alt_seq_samps[1:]:
            has_mod.append(' X' if (seq_samp, alt_mod) in alt_model_types
                           else '')
        sys.stderr.write(row_format.format(*has_mod))

    return

def load_default_alt_ref(alt_name, seq_samp_type):
    try:
        alt_model_fn = ALTERNATE_MODELS[
            seq_samp_type.name + ALT_MODEL_SEP_CHAR + alt_name]
    except KeyError:
        alt_model_fn = None
    if alt_model_fn is not None:
        # get full filename path with setuptools
        alt_model_fn = pkg_resources.resource_filename(
            'tombo', 'tombo_models/' + alt_model_fn)
    if alt_model_fn is None or not os.path.isfile(alt_model_fn):
        th.warning_message(
            'Tombo default model for ' + alt_name + ' in ' +
            seq_samp_type.name + ' does not exists.')
        return None

    return AltModel(ref_fn=alt_model_fn)

def load_alt_refs(alt_model_fns, alt_names, reads_index, std_ref,
                  seq_samp_type=None):
    alt_refs = {}
    if alt_model_fns is not None:
        # load alternative models from filenames
        for alt_model_fn in alt_model_fns:
            alt_ref = AltModel(alt_model_fn)
            if alt_ref.name in alt_refs:
                th.warning_message(
                    alt_ref.name + ' alternative model found in more than ' +
                    'one model file. Ignoring: ' + alt_model_fn)
                continue
            alt_refs[alt_ref.name] = alt_ref
    else:
        # load alternative models from internal defaults
        if seq_samp_type is None:
            seq_samp_type = th.get_seq_sample_type(reads_index=reads_index)
        for alt_name in alt_names:
            alt_ref = load_default_alt_ref(alt_name, seq_samp_type)
            if alt_ref is None:
                continue
            alt_refs[alt_name] = alt_ref

    check_valid_alt_models(alt_refs, std_ref)

    return alt_refs

def load_valid_models(
        tb_model_fn, plot_default_stnd, alt_model_fn,
        plot_default_alt, reads_index, ctrl_fast5s_dirs=None):
    # if no model was requested
    if (tb_model_fn is None and not plot_default_stnd and
        alt_model_fn is None and not plot_default_alt):
        return None, None

    std_ref = TomboModel(ref_fn=tb_model_fn, reads_index=reads_index)
    if alt_model_fn is not None:
        alt_ref = AltModel(ref_fn=alt_model_fn)
    elif plot_default_alt is not None:
        seq_samp_type = std_ref.seq_samp_type
        if seq_samp_type is None:
            seq_samp_type = th.get_seq_sample_type(reads_index=reads_index)
        alt_ref = load_default_alt_ref(plot_default_alt, seq_samp_type)
    else:
        alt_ref = None

    if ctrl_fast5s_dirs is not None and tb_model_fn is not None:
        th.warning_message(
            'Both a second set of FAST5s and a tombo model were ' +
            'provided. Two samples with model plotting is not ' +
            'currently available. Models requested will be ignored.')

    return std_ref, alt_ref

def calc_med_sd(vals):
    """Helper function to compute median and standard deviation from a
    numpy array
    """
    return np.median(vals), np.std(vals)

def get_region_kmer_levels(
        reg_data, cov_thresh, upstrm_bases, dnstrm_bases, cs_cov_thresh,
        est_mean, region_size, motif=None, valid_poss=None):
    """Compute mean or median and standard deviation for each k-mer
    """
    if cs_cov_thresh is not None:
        # sample reads until requested mean depth of coverage is achieved
        cs_num_bases_thresh = region_size * cs_cov_thresh
        np.random.shuffle(reg_data.reads)
        cumm_num_bases = np.cumsum([
            max(r_data.end, reg_data.end) - min(r_data.start, reg_data.start)
            for r_data in reg_data.reads])
        try:
            cs_num_reads = next((i for i, v in enumerate(cumm_num_bases)
                                 if v >= cs_num_bases_thresh))
            reg_data.update(reads=reg_data.reads[:cs_num_reads])
        except StopIteration:
            # if threshold is not met use all reads from region
            pass
    # TODO convert this to the intervalData method get_base_levels function
    # involves a bit of logical refactoring below (which should be much simpler)
    base_events = th.get_reads_events(reg_data.reads)
    if len(base_events) == 0:
        return

    # get intervals within the region where coverage is high enough
    # for model estimation
    reg_cov = np.array([
        len(base_events[pos]) if pos in base_events else 0
        for pos in range(reg_data.start, reg_data.end)])
    cov_intervals = np.where(np.diff(np.concatenate(
        [[False], reg_cov > cov_thresh])))[0]
    if reg_cov[-1] > cov_thresh:
        # if valid coverage goes until end of coverage, add the last
        # interval endpoint
        cov_intervals = np.concatenate([cov_intervals, [region_size]])
    if cov_intervals.shape[0] <= 1:
        return
    cov_intervals = cov_intervals.reshape(-1,2)

    kmer_width = upstrm_bases + dnstrm_bases + 1
    if motif is None:
        reg_kmer_levels = dict(
            (''.join(kmer), [])
            for kmer in product(DNA_BASES, repeat=kmer_width))
    else:
        reg_kmer_levels = dict(
            ((''.join(kmer), i_offset - 1), [])
            for kmer in product(DNA_BASES, repeat=kmer_width)
            for i_offset in motif.find_mod_poss(''.join(kmer)))

    # upstream and downstream changes the sequence selection
    # depending on the strand
    # TODO this will miss motif hits at region boundaries when the motif is
    # longer than either end of the k-me relative to the central position
    bb, ab = (upstrm_bases, dnstrm_bases) if reg_data.strand == '+' else \
             (dnstrm_bases, upstrm_bases)
    for cov_start, cov_end in cov_intervals:
        # get region sequence from expanded region to include k-mer lookups
        int_seq = reg_data.copy().update(
            start=reg_data.start + cov_start - bb,
            end=reg_data.start + cov_end + ab).add_seq().seq
        int_len = cov_end - cov_start

        # get valid positions to include in extracted levels
        # for standard extraction fill with None relative position
        if valid_poss is None and motif is None:
            int_poss = zip(range(int_len), repeat(None))
        # for motif/position extraction record relative modified position
        else:
            # get modified base position relative to coverage interval
            if valid_poss is not None:
                if (reg_data.chrm, reg_data.strand) not in valid_poss: continue
                reg_mod_poss = (valid_poss[(reg_data.chrm, reg_data.strand)] -
                                reg_data.start - cov_start)
                reg_mod_poss = reg_mod_poss[
                    np.logical_and(np.greater_equal(reg_mod_poss, 0),
                                   np.less(reg_mod_poss, int_len))]
            else:
                if reg_data.strand == '+':
                    reg_mod_poss = [
                        m.start() + motif.mod_pos - 1 - bb
                        for m in motif.motif_pat.finditer(int_seq)
                        if 0 <= m.start() + motif.mod_pos - 1 - bb < int_len]
                else:
                    reg_mod_poss = [
                        m.start() + motif.motif_len - motif.mod_pos - bb
                        for m in motif.rev_comp_pat.finditer(int_seq)
                        if 0 <= m.start() + motif.motif_len - motif.mod_pos - bb
                        < int_len]

            # record relative mod positions within k-mers
            int_poss = [
                (mod_pos - i_offset + bb, i_offset if reg_data.strand == '+'
                 else kmer_width - i_offset - 1)
                for mod_pos in reg_mod_poss for i_offset in range(kmer_width)
                if 0 <= mod_pos - i_offset + bb < int_len]

        for pos, offset in int_poss:
            pos_kmer = int_seq[pos:pos + kmer_width]
            if reg_data.strand == '-':
                pos_kmer = th.rev_comp(pos_kmer)
            #print(pos_kmer, offset,
            #      pos + reg_data.start + cov_start, reg_data.strand,
            #      reg_data.start + cov_start - bb,
            #      reg_data.start + cov_end + ab)
            kmer_key = pos_kmer if offset is None else (pos_kmer, offset)
            try:
                if est_mean:
                    reg_kmer_levels[kmer_key].append(c_mean_std(
                        base_events[pos + reg_data.start + cov_start]))
                else:
                    reg_kmer_levels[kmer_key].append(calc_med_sd(
                        base_events[pos + reg_data.start + cov_start]))
            except KeyError:
                continue

    return reg_kmer_levels

def _est_kmer_model_worker(
        region_q, kmer_level_q, progress_q, reads_index, cov_thresh,
        upstrm_bases, dnstrm_bases, cs_cov_thresh, est_mean, region_size,
        motif, valid_poss):
    while True:
        try:
            chrm, strand, reg_start = region_q.get(block=False)
        except queue.Empty:
            # sometimes throws false empty error with get(block=False)
            sleep(0.1)
            if not region_q.empty(): continue
            break

        reg_data = th.intervalData(
            chrm=chrm, start=reg_start, end=reg_start + region_size,
            strand=strand).add_reads(reads_index)
        if len(reg_data.reads) == 0:
            progress_q.put(1)
            continue

        reg_kmer_levels = get_region_kmer_levels(
            reg_data, cov_thresh, upstrm_bases, dnstrm_bases,
            cs_cov_thresh, est_mean, region_size, motif, valid_poss)
        if reg_kmer_levels is not None:
            kmer_level_q.put(reg_kmer_levels)
        progress_q.put(1)

    return

if _PROFILE_EST_REF:
    _est_kmer_model_wrapper = _est_kmer_model_worker
    def _est_kmer_model_worker(*args):
        import cProfile
        cProfile.runctx('_est_kmer_model_wrapper(*args)', globals(), locals(),
                        filename='est_kmer_model.prof')
        return

def extract_kmer_levels(
        reads_index, region_size, cov_thresh, upstrm_bases, dnstrm_bases,
        cs_cov_thresh, est_mean=False, num_processes=1,
        motif=None, valid_poss=None):
    region_q = Queue()
    kmer_level_q = Queue()
    progress_q = Queue()
    num_regions = 0
    for chrm, strand, reg_start in reads_index.iter_cov_regs(
            cov_thresh, region_size):
        region_q.put((chrm, strand, reg_start))
        num_regions += 1

    est_args = (
        region_q, kmer_level_q, progress_q, reads_index, cov_thresh,
        upstrm_bases, dnstrm_bases, cs_cov_thresh, est_mean, region_size,
        motif, valid_poss)
    est_ps = []
    for p_id in range(num_processes):
        p = Process(target=_est_kmer_model_worker, args=est_args)
        p.start()
        est_ps.append(p)

    if VERBOSE:
        th.status_message('Extracting average k-mer levels.')
        bar = tqdm(total=num_regions, smoothing=0)
    all_reg_kmer_levels = []
    while any(p.is_alive() for p in est_ps):
        try:
            reg_kmer_levels = kmer_level_q.get(block=False)
            all_reg_kmer_levels.append(reg_kmer_levels)
        except queue.Empty:
            try:
                iter_proc = progress_q.get(block=False)
                if VERBOSE: bar.update(iter_proc)
            except queue.Empty:
                sleep(1)
                continue

    # clear rest of queue
    while not kmer_level_q.empty():
        reg_kmer_levels = kmer_level_q.get(block=False)
        all_reg_kmer_levels.append(reg_kmer_levels)
    if VERBOSE:
        while not progress_q.empty():
            iter_proc = progress_q.get(block=False)
            if VERBOSE: bar.update(iter_proc)
        bar.close()

    if len(all_reg_kmer_levels) == 0:
        th.error_message_and_exit(
            'No genomic positions contain --minimum-test-reads. Consider ' +
            'setting this option to a lower value.')

    return all_reg_kmer_levels

def tabulate_kmer_levels(all_reg_kmer_levels, min_kmer_obs):
    if VERBOSE: th.status_message('Tabulating k-mer model statistics.')
    all_kmer_mean_sds = []
    if _DEBUG_EST_STD:
        kmer_dens = []
        save_x = np.linspace(KERNEL_DENSITY_RANGE[0], KERNEL_DENSITY_RANGE[1],
                             _DEBUG_EST_NUM_KMER_SAVE)
    kmer_width = len(next(iter(all_reg_kmer_levels[0].keys())))
    for kmer in product(DNA_BASES, repeat=kmer_width):
        kmer = ''.join(kmer)
        try:
            kmer_levels = np.concatenate([
                reg_kmer_levels[kmer] for reg_kmer_levels in all_reg_kmer_levels
                if len(reg_kmer_levels[kmer]) > 0])
        except ValueError:
            th.error_message_and_exit(
                'At least one k-mer is not covered at any poitions by ' +
                '--minimum-test-reads.\n\t\tConsider fitting to a smaller ' +
                'k-mer via the --upstream-bases and --downstream-bases, ' +
                'or lowering --minimum-test-reads.\n\t\tNote that this may ' +
                'result in a lower quality model.')
        if kmer_levels.shape[0] < min_kmer_obs:
            min_obs = min(
                sum(len(reg_levs[''.join(kmer)])
                    for reg_levs in all_reg_kmer_levels)
                for kmer in product(DNA_BASES, repeat=kmer_width)
                if motif is None or motif.matches_seq(''.join(kmer)))
            th.error_message_and_exit(
                'K-mers represeneted in fewer observations than ' +
                'requested in the provided reads. Consider a shorter ' +
                'k-mer or providing more reads.\n\t' + unicode(min_obs) +
                ' observations found in least common kmer.')
        all_kmer_mean_sds.append((kmer, np.median(kmer_levels[:,0]),
                                  np.median(kmer_levels[:,1])))
        if _DEBUG_EST_STD:
            kmer_kde = stats.gaussian_kde(
                kmer_levels[:,0],
                bw_method=_DEBUG_EST_BW / kmer_levels[:,0].std(ddof=1))
            with np.errstate(under='ignore'):
                kmer_dens.append((kmer, kmer_kde.evaluate(save_x)))

    if _DEBUG_EST_STD:
        with io.open('debug_est_standard_ref.density.txt', 'wt') as fp:
            fp.write('Kmer\tSignal\tDensity\n')
            fp.write('\n'.join('\t'.join(map(str, (kmer, x, y)))
                               for kmer, dens_i in kmer_dens
                               for x, y in zip(save_x, dens_i)) + '\n')
    return all_kmer_mean_sds

# methods needed for (re-squiggle) segmentation are also needed here for
# RNA event-based scaling (for model scale centering)
def load_resquiggle_parameters(
        seq_samp_type, sig_aln_params=None, seg_params=None,
        use_save_bandwidth=False):
    """Load parameters for re-squiggle algorithm

    Args:
        seq_samp_type (:class:`tombo.tombo_helper.seqSampleType`): sequencing
            sample type
        sig_aln_params (tuple): signal alignment parameters (optional; default:
            load seq_samp_type defaults)
        seg_params (tuple): segmentation parameters (optional; default: load
            seq_samp_type defaults)
        use_save_bandwidth (bool): load larger "save" bandwidth

    Returns:
        :class:`tombo.tombo_helper.resquiggleParams`
    """
    if sig_aln_params is None:
        (match_evalue, skip_pen, bandwidth, save_bandwidth, max_half_z_score,
         band_bound_thresh, start_bw, start_save_bw,
         start_n_bases) = ALGN_PARAMS_TABLE[seq_samp_type.name]
    else:
        # unpack signal alignment parameters
        (match_evalue, skip_pen, bandwidth, save_bandwidth,
         max_half_z_score, band_bound_thresh, start_bw, start_save_bw,
         start_n_bases) = sig_aln_params
        bandwidth = int(bandwidth)
        save_bandwidth = int(save_bandwidth)
        band_bound_thresh = int(band_bound_thresh)
        start_bw = int(start_bw)
        start_save_bw = int(start_save_bw)
        start_n_bases = int(start_n_bases)
    if use_save_bandwidth:
        bandwidth = save_bandwidth

    if seg_params is None:
        (running_stat_width, min_obs_per_base, raw_min_obs_per_base,
         mean_obs_per_event) = SEG_PARAMS_TABLE[seq_samp_type.name]
    else:
        (running_stat_width, min_obs_per_base, raw_min_obs_per_base,
         mean_obs_per_event) = seg_params

    z_shift, stay_pen = get_dynamic_prog_params(match_evalue)

    rsqgl_params = th.resquiggleParams(
        match_evalue, skip_pen, bandwidth, max_half_z_score,
        running_stat_width, min_obs_per_base, raw_min_obs_per_base,
        mean_obs_per_event, z_shift, stay_pen,
        seq_samp_type.name == RNA_SAMP_TYPE, band_bound_thresh,
        start_bw, start_save_bw, start_n_bases)

    return rsqgl_params

def compute_num_events(
        signal_len, seq_len, mean_obs_per_event,
        min_event_to_seq_ratio=MIN_EVENT_TO_SEQ_RATIO):
    """Compute number of events to find for this read

    Args:
        signal_len (int): length of raw signal
        seq_len (int): length of sequence
        mean_obs_per_base (int): mean raw observations per genome base
        min_event_to_seq_ratio (float): minimum event to sequence ratio
            (optional)

    Returns:
        Number of events to find for this read
    """
    return max(signal_len // mean_obs_per_event,
               int(seq_len * min_event_to_seq_ratio))

def remove_stall_cpts(stall_ints, valid_cpts):
    if len(stall_ints) == 0:
        return valid_cpts

    # RNA data contains stall regions that can cause problems for
    # banded dynamic programming so they are removed here
    stall_int_iter = iter(stall_ints)
    curr_stall_int = next(stall_int_iter)
    non_stall_cpts = []
    # loop over valid cpts
    for i, cpt in enumerate(valid_cpts):
        # iterate through stall intervals until the current interval end
        # is greater than the cpt to check against
        while cpt > curr_stall_int[1]:
            try:
                curr_stall_int = next(stall_int_iter)
            except StopIteration:
                break
        if not (curr_stall_int[0] < cpt < curr_stall_int[1]):
            non_stall_cpts.append(i)

    return valid_cpts[non_stall_cpts]

def center_model_to_median_norm(
        reads_index, init_ref, max_reads=NUM_READS_TO_ADJUST_MODEL):
    upstrm_bases = init_ref.central_pos
    dnstrm_bases = init_ref.kmer_width - init_ref.central_pos - 1
    def get_event_scale_values(all_raw_signal, r_len):
        rsqgl_params = load_resquiggle_parameters(
            th.seqSampleType(RNA_SAMP_TYPE, True))

        num_events = compute_num_events(
            all_raw_signal.shape[0], r_len,
            rsqgl_params.mean_obs_per_event, MIN_EVENT_TO_SEQ_RATIO)
        valid_cpts = th.valid_cpts_w_cap_t_test(
            all_raw_signal.astype(np.float64), rsqgl_params.min_obs_per_base,
            rsqgl_params.running_stat_width, num_events)
        if COLLAPSE_RNA_STALLS:
            valid_cpts = remove_stall_cpts(
                identify_stalls(all_raw_signal, DEFAULT_STALL_PARAMS),
                valid_cpts)
        scale_values = get_scale_values_from_events(
            all_raw_signal, valid_cpts, OUTLIER_THRESH,
            num_events=RNA_SCALE_NUM_EVENTS,
            max_frac_events=RNA_SCALE_MAX_FRAC_EVENTS)
        return normalize_raw_signal(
            all_raw_signal, scale_values=scale_values)

    def get_read_corr_factors(r_data):
        with h5py.File(r_data.fn, 'r+') as fast5_data:
            all_raw_signal = th.get_raw_read_slot(fast5_data)['Signal'][:]
            event_starts, r_seq = th.get_multiple_slots_read_centric(
                fast5_data, ('start', 'base'), r_data.corr_group)

        if r_data.rna:
            all_raw_signal = all_raw_signal[::-1]
            if USE_RNA_EVENT_SCALE:
                norm_signal, scale_values = get_event_scale_values(
                    all_raw_signal, r_data.end - r_data.start)
            else:
                # use raw signal median normalization
                norm_signal, scale_values = normalize_raw_signal(
                    all_raw_signal)
        else:
            norm_signal, scale_values = normalize_raw_signal(all_raw_signal)

        event_starts = event_starts.astype(np.int64)
        rsrtr = r_data.read_start_rel_to_raw + event_starts[upstrm_bases]
        # since the last segment end wasn't extracted the events are already
        # clipped by one, so deal with event boundary clipping here
        if dnstrm_bases > 1:
            event_starts = event_starts[upstrm_bases:-(dnstrm_bases - 1)]
        else:
            assert dnstrm_bases == 1, (
                'Must have at least one upstream and downstream base for ' +
                'a Tombo model.')
            event_starts = event_starts[upstrm_bases:]
        # reset event starts to 0
        event_starts -= event_starts[0]

        norm_signal = norm_signal[rsrtr:rsrtr + event_starts[-1]]

        r_seq = b''.join(r_seq).decode()
        r_ref_means, _ = init_ref.get_exp_levels_from_seq(r_seq)

        (_, _, shift_corr_factor,
         scale_corr_factor) = calc_kmer_fitted_shift_scale(
             scale_values.shift, scale_values.scale,
             compute_base_means(norm_signal, event_starts), r_ref_means,
             method='theil_sen')

        return shift_corr_factor, scale_corr_factor


    th.status_message('Centering model to normalized signal')
    all_shift_corr_factors, all_scale_corr_factors = [], []
    all_reads = list(reads_index.iter_reads())
    np.random.shuffle(all_reads)
    for r_data in all_reads:
        try:
            r_shift_corr_factor, r_scale_corr_factor = get_read_corr_factors(
                r_data)
            all_shift_corr_factors.append(r_shift_corr_factor)
            all_scale_corr_factors.append(r_scale_corr_factor)
            if len(all_shift_corr_factors) >= max_reads:
                break
        except:
            continue

    if len(all_shift_corr_factors) < max_reads:
        if len(all_shift_corr_factors) == 0:
            th.error_message_and_exit(
                'No reads succcessfully processed for sequence-based ' +
                'normalization parameter re-fitting.')
        th.warning_message(
            'Fewer reads succcessfully processed for sequence-based ' +
            'normalization parameter re-fitting than requested.')

    # compute median shift and scale correction factors
    # scale parameter should be taken in log space, but median performs
    # the same computation
    med_shift_corr_factor = np.median(all_shift_corr_factors)
    med_scale_corr_factor = np.median(all_scale_corr_factors)

    th.status_message('Shift and scale adjustments to match model to ' +
                       'median normalization: ' + str(med_shift_corr_factor) +
                       " " + str(med_scale_corr_factor))
    init_ref._center_model(med_shift_corr_factor, med_scale_corr_factor)

    return init_ref

if _PROFILE_CENTER_REF:
    center_model_to_median_norm_wrapper = center_model_to_median_norm
    def center_model_to_median_norm(*args):
        import cProfile
        cProfile.runctx(
            'center_model_to_median_norm_wrapper(*args)', globals(), locals(),
            filename='center_kmer_model.prof')
        return

def estimate_kmer_model(
        fast5s_dirs, corr_grp, bc_subgrps,
        cov_thresh, upstrm_bases, dnstrm_bases, min_kmer_obs,
        kmer_specific_sd, cs_cov_thresh, est_mean, region_size, num_processes):
    """Estimate a standard tombo k-mer model
    """
    reads_index = th.TomboReads(fast5s_dirs, corr_grp, bc_subgrps)
    all_reg_kmer_levels = extract_kmer_levels(
        reads_index, region_size, cov_thresh, upstrm_bases, dnstrm_bases,
        cs_cov_thresh, est_mean, num_processes)

    all_kmer_mean_sds = tabulate_kmer_levels(all_reg_kmer_levels, min_kmer_obs)

    # adjust model to match median normalization best via Theil-Sen optimizer
    # fit this will increase the accuracy of median normalized re-squiggle
    # results and should reduce the need for (or number of) iterative
    # re-squiggle runs
    init_ref = TomboModel(kmer_ref=all_kmer_mean_sds, central_pos=upstrm_bases)

    centered_ref = center_model_to_median_norm(reads_index, init_ref)

    if not kmer_specific_sd:
        centered_ref._make_constant_sd()

    return centered_ref


########################################
##### Alternative Model Estimation #####
########################################

def _parse_base_levels_worker(
        reads_q, kmer_level_q, kmer_width, central_pos, completed_kmers):
    dnstrm_bases = kmer_width - central_pos - 1
    proc_kmer_levels = dict(
        (''.join(kmer), [])
        for kmer in product(DNA_BASES, repeat=kmer_width))
    while True:
        try:
            r_fn, corr_slot = reads_q.get(block=False)
        except queue.Empty:
            # sometimes throws false empty exception with get(block=False)
            sleep(0.01)
            if not reads_q.empty(): continue
            break

        with h5py.File(r_fn, 'r') as fast5_data:
            r_means, r_seq = th.get_multiple_slots_read_centric(
                fast5_data, ['norm_mean', 'base'], corr_slot)
        if r_means is None: continue
        r_seq = b''.join(r_seq).decode()
        for kmer, level in zip(
                (r_seq[i:i + kmer_width]
                 for i in range(len(r_seq) - kmer_width + 1)),
                r_means[central_pos:-dnstrm_bases]):
            if kmer in completed_kmers: continue
            proc_kmer_levels[kmer].append(level)

    kmer_level_q.put(proc_kmer_levels)

    return

def get_batch_kmer_levels(
        reads_q, kmer_level_q, all_reads, parse_levels_batch_size,
        std_ref, completed_kmers, num_processes):
    no_more_reads = False
    try:
        for _ in range(parse_levels_batch_size):
            r_data = next(all_reads)
            reads_q.put((r_data.fn, r_data.corr_group))
    except StopIteration:
        no_more_reads = True

    base_lev_args = (reads_q, kmer_level_q, std_ref.kmer_width,
                     std_ref.central_pos, completed_kmers)
    base_lev_ps = []
    for p_id in range(num_processes):
        p = Process(target=_parse_base_levels_worker, args=base_lev_args)
        p.start()
        base_lev_ps.append(p)

    batch_kmer_levels = []
    while any(p.is_alive() for p in base_lev_ps):
        try:
            proc_kmer_levels = kmer_level_q.get(block=False)
            batch_kmer_levels.append(proc_kmer_levels)
        except queue.Empty:
            sleep(1)
            continue
    while not kmer_level_q.empty():
        proc_kmer_levels = kmer_level_q.get(block=False)
        batch_kmer_levels.append(proc_kmer_levels)

    return batch_kmer_levels, no_more_reads

def parse_base_levels(
        all_reads, std_ref, parse_levels_batch_size, kmer_obs_thresh,
        max_kmer_obs, min_kmer_obs_to_est, num_processes):
    """
    Parse base levels and store grouped by k-mer
    """
    reads_q = Queue()
    kmer_level_q = Queue()

    all_kmer_levels = dict(
        (''.join(kmer), [])
        for kmer in product(DNA_BASES, repeat=std_ref.kmer_width))
    # store set of k-mers with enough observations to save on memory footprint
    # while filling more rare k-mers
    completed_kmers = set()
    if VERBOSE:
        all_reads_bar = tqdm(total=len(all_reads), smoothing=0,
                             desc='Number of total reads used', leave=True)
        min_kmer_bar = tqdm(total=kmer_obs_thresh, smoothing=0,
                            desc='K-mer with fewest observations', leave=True)
        curr_min_kmer_count = 0
    all_reads = iter(all_reads)
    n_batches = 0
    while True:
        batch_kmer_levels, no_more_reads = get_batch_kmer_levels(
            reads_q, kmer_level_q, all_reads, parse_levels_batch_size,
            std_ref, completed_kmers, num_processes)

        # only add observations for k-mers that have not been seen enough times
        # save memory for abundent k-mers
        batch_total_kmers = []
        for kmer in set(all_kmer_levels).difference(completed_kmers):
            all_kmer_levels[kmer].extend((
                kmer_level for proc_kmer_levels in batch_kmer_levels
                for kmer_level in proc_kmer_levels[kmer]))
            batch_total_kmers.append(len(all_kmer_levels[kmer]))
            if batch_total_kmers[-1] > max_kmer_obs:
                completed_kmers.add(kmer)
        if VERBOSE:
            if no_more_reads:
                all_reads_bar.update(all_reads_bar.total - all_reads_bar.n)
            else:
                all_reads_bar.update(parse_levels_batch_size)
            new_min_kmer_count = min(batch_total_kmers)
            min_kmer_bar.update(new_min_kmer_count - curr_min_kmer_count)
            curr_min_kmer_count = new_min_kmer_count
            sleep(0.1)
        else:
            curr_min_kmer_count = min(batch_total_kmers)
        if curr_min_kmer_count > kmer_obs_thresh or no_more_reads:
            break

        n_batches += 1

    if VERBOSE:
        all_reads_bar.close()
        min_kmer_bar.close()

    fewest_kmer_obs = min(len(kmer_levels) for kmer_levels in
                          all_kmer_levels.values())
    if fewest_kmer_obs < kmer_obs_thresh:
        if fewest_kmer_obs < min_kmer_obs_to_est:
            th.error_message_and_exit(
                'Too few minimal k-mer observations to continue to ' +
                'alternative estimation. Minimal k-mer has ' +
                unicode(fewest_kmer_obs) + ' total observations and ' +
                unicode(min_kmer_obs_to_est) +
                ' observations per k-mer are required.')
        th.warning_message(
            'Requested minimal k-mer observations not found in all reads. ' +
            'Continuing to estimation using a k-mer with ' +
            unicode(fewest_kmer_obs) + ' total observations')

    return all_kmer_levels

def write_kmer_densities_file(dens_fn, kmer_dens, save_x):
    with io.open(dens_fn, 'wt') as fp:
        fp.write('Kmer\tSignal\tDensity\n')
        fp.write('\n'.join('\t'.join(map(str, (kmer, x, y)))
                           for kmer, dens_i in kmer_dens.items()
                           for x, y in zip(save_x, dens_i)) + '\n')

    return

def parse_kmer_densities_file(dens_fn):
    kmer_dens_raw = defaultdict(list)
    with io.open(dens_fn) as dens_fp:
        # read in header
        dens_fp.readline()
        for line in dens_fp:
            kmer, _, dens_i = line.split()
            kmer_dens_raw[kmer].append(float(dens_i))

    kmer_dens = {}
    first_len = None
    for kmer, dens_i in kmer_dens_raw.items():
        if first_len is None: first_len = len(dens_i)
        if len(dens_i) != first_len:
            th.error_message_and_exit('Density file is valid.')
        kmer_dens[kmer] = np.array(dens_i)

    return kmer_dens

def est_kernel_density(
        reads_index, std_ref, kmer_obs_thresh, density_basename, save_x,
        kernel_dens_bw, num_processes, alt_or_stnd_name='alt',
        parse_levels_batch_size=ALT_EST_BATCH, max_kmer_obs=MAX_KMER_OBS,
        min_kmer_obs_to_est=MIN_KMER_OBS_TO_EST):
    all_reads = list(reads_index.iter_reads())
    np.random.shuffle(all_reads)
    base_levels = parse_base_levels(
        all_reads, std_ref, parse_levels_batch_size, kmer_obs_thresh,
        max_kmer_obs, min_kmer_obs_to_est, num_processes)

    if VERBOSE: th.status_message('Fitting kernel densities for k-mer levels.')
    kmer_dens = {}
    for kmer, norm_levels in base_levels.items():
        norm_levels = np.array(norm_levels)
        kmer_kde = stats.gaussian_kde(
            norm_levels, bw_method=kernel_dens_bw / norm_levels.std(ddof=1))
        with np.errstate(under='ignore'):
            kmer_dens[kmer] = kmer_kde.evaluate(save_x)

    if density_basename is not None:
        write_kmer_densities_file(
            density_basename + '.' + alt_or_stnd_name + '_density.txt',
            kmer_dens, save_x)

    return kmer_dens

def estimate_kmer_densities(
        fast5s_dirs, ctrl_fast5s_dirs, corr_grp, bc_subgrps,
        standard_ref_fn, seq_samp_type, kmer_obs_thresh, density_basename,
        kernel_dens_bw, save_x, num_processes):
    reads_index = th.TomboReads(fast5s_dirs, corr_grp, bc_subgrps)
    ctrl_reads_index = th.TomboReads(ctrl_fast5s_dirs, corr_grp, bc_subgrps)

    if VERBOSE: th.status_message('Parsing standard model file.')
    std_ref = TomboModel(ref_fn=standard_ref_fn, seq_samp_type=seq_samp_type,
                         reads_index=reads_index)

    if VERBOSE: th.status_message('Parsing base levels from alternative reads.')
    alt_dens = est_kernel_density(
        reads_index, std_ref, kmer_obs_thresh, density_basename,
        save_x, kernel_dens_bw, num_processes, 'alternate')
    if VERBOSE: th.status_message('Parsing base levels from standard reads.')
    std_dens = est_kernel_density(
        ctrl_reads_index, std_ref, kmer_obs_thresh, density_basename,
        save_x, kernel_dens_bw, num_processes, 'control')

    return alt_dens, std_dens, std_ref

def load_kmer_densities(
        alt_dens_fn, std_dens_fn, fast5s_dirs, corr_grp, bc_subgrps,
        std_ref_fn, seq_samp_type):
    if VERBOSE: th.status_message('Parsing standard model file.')
    reads_index = None
    if std_ref_fn is None and seq_samp_type is None:
        if fast5s_dirs is None:
            th.error_message_and_exit(
                'Must provide a FAST5s directory, a canonical model ' +
                'file or spcify the biological sample type.')
        reads_index = th.TomboReads(fast5s_dirs, corr_grp, bc_subgrps)
    std_ref = TomboModel(ref_fn=std_ref_fn, seq_samp_type=seq_samp_type,
                         reads_index=reads_index)

    if VERBOSE: th.status_message('Parsing density files.')
    alt_dens = parse_kmer_densities_file(alt_dens_fn)
    std_dens = parse_kmer_densities_file(std_dens_fn)
    num_dens_points = next(v for v in alt_dens.values()).shape[0]
    if num_dens_points != next(v for v in std_dens.values()).shape[0]:
        th.error_message_and_exit(
            'Alternative and standard density ' +
            'estimates do not correspond.')

    save_x = np.linspace(KERNEL_DENSITY_RANGE[0], KERNEL_DENSITY_RANGE[1],
                         num_dens_points)

    return alt_dens, std_dens, std_ref, save_x

def isolate_alt_density(
        alt_dens, std_dens, alt_base, alt_frac_pctl, std_ref, save_x):
    def calc_mean(dens):
        return np.average(save_x[dens>1e-10], weights=dens[dens>1e-10])


    # estimate density shift from k-mers without the alternate base
    no_alt_std_means, no_alt_mean_diffs = [], []
    for kmer in std_dens:
        if alt_base in kmer: continue
        no_alt_std_means.append(calc_mean(std_dens[kmer]))
        kmer_alt_mean = calc_mean(alt_dens[kmer])
        no_alt_mean_diffs.append(kmer_alt_mean - no_alt_std_means[-1])
    calc_offset = np.poly1d(np.polyfit(no_alt_std_means, no_alt_mean_diffs, 2))
    save_x_unit = save_x[1] - save_x[0]

    shifted_alt_dens = {}
    for kmer, kmer_alt_dens in alt_dens.items():
        est_offset = int(calc_offset(calc_mean(std_dens[kmer])) / save_x_unit)
        # if std density mean is estimated to be greater
        if est_offset < 0:
            # shift alt dens to the right
            shifted_alt_dens[kmer] = np.concatenate([
                [0.0,] * -est_offset, kmer_alt_dens[:est_offset]])
        else:
            # else shift alt dens to the left
            shifted_alt_dens[kmer] = np.concatenate([
                kmer_alt_dens[est_offset:], [0.0,] * est_offset])

    def get_peak_frac(kmer_std_dens, kmer_alt_dens):
        # find highest peak in standard density
        std_peak = np.argmax(kmer_std_dens)

        # find closest peak in alternative density
        alt_local_peaks = np.where(np.concatenate([
            [False,], np.logical_and(
                kmer_alt_dens[1:-1] > kmer_alt_dens[:-2],
                kmer_alt_dens[1:-1] > kmer_alt_dens[2:]) + [False,]]))[0]
        matched_alt_peak = alt_local_peaks[
            np.argmin(abs(alt_local_peaks - std_peak))]
        return kmer_alt_dens[matched_alt_peak] / kmer_std_dens[std_peak]


    # estimate the alternative base incorporation rate
    std_frac = np.percentile([
        get_peak_frac(std_dens[kmer], shifted_alt_dens[kmer])
        for kmer in std_dens if kmer.count(alt_base) == 1], alt_frac_pctl)
    if VERBOSE: th.status_message(
            'Alternative base incorporation rate estimate: ' +
            unicode(1 - std_frac))
    if std_frac >= 1:
        th.warning_message(
            'Alternative base incorporation rate ' +
            'estimate is approximately 0. Consider lowering ' +
            '--alt-fraction-percentile.')

    # get mean model SD. most models will be constant, but use mean in case
    model_sd = np.mean(list(std_ref.sds.values()))
    # subtract off fraction of standard density from alt density
    # to estimate mean of isolated alternative distribution
    alt_ref = []
    for kmer, std_level in std_ref.means.items():
        if kmer.count(alt_base) == 0:
            continue
        # assuming random incorporation the prortion of standard base
        # observations at this k-mer is the standard fraction raised
        # to the number of alt_base occurences in the k-mer
        kmer_std_frac = std_frac**kmer.count(alt_base)
        with np.errstate(under='ignore'):
            diff_dens = shifted_alt_dens[kmer] - (
                std_dens[kmer] * kmer_std_frac)
            diff_dens[diff_dens < 0] = 0
            alt_level = np.average(save_x, weights=diff_dens)
        # add alt mean for each alt base position within the kmer
        for m in re.finditer(alt_base, kmer):
            alt_ref.append((kmer, m.start(), alt_level, model_sd))

    alt_ref = AltModel(
        kmer_ref=alt_ref, central_pos=std_ref.central_pos, alt_base=alt_base)

    return alt_ref

def estimate_alt_model(
        fast5s_dirs, ctrl_fast5s_dirs, corr_grp, bc_subgrps,
        std_ref_fn, seq_samp_type, alt_base, alt_frac_pctl,
        kmer_obs_thresh, density_basename, kernel_dens_bw, alt_dens_fn,
        std_dens_fn, num_processes, num_dens_points=NUM_DENS_POINTS):
    """Estimate an alternative model from a sample with a single,
    known, randomly-incorporated alternative base
    """
    if alt_dens_fn is None or std_dens_fn is None:
        save_x = np.linspace(KERNEL_DENSITY_RANGE[0], KERNEL_DENSITY_RANGE[1],
                             num_dens_points)
        alt_dens, std_dens, std_ref = estimate_kmer_densities(
            fast5s_dirs, ctrl_fast5s_dirs, corr_grp, bc_subgrps,
            std_ref_fn, seq_samp_type, kmer_obs_thresh, density_basename,
            kernel_dens_bw, save_x, num_processes)
    else:
        alt_dens, std_dens, std_ref, save_x = load_kmer_densities(
            alt_dens_fn, std_dens_fn, fast5s_dirs, corr_grp, bc_subgrps,
            std_ref_fn, seq_samp_type)

    if VERBOSE: th.status_message('Isolating alternative base distribtuions.')
    # perform alternative density isolation algorithm
    alt_ref = isolate_alt_density(
        alt_dens, std_dens, alt_base, alt_frac_pctl, std_ref, save_x)

    return alt_ref

if _PROFILE_ALT_EST:
    _est_alt_wrapper = estimate_alt_model
    def estimate_alt_model(*args):
        import cProfile
        cProfile.runctx('_est_alt_wrapper(*args)', globals(), locals(),
                        filename='est_alt_model.prof')
        return None

def tabulate_mod_kmer_levels(all_reg_kmer_levels, min_kmer_obs, motif):
    if VERBOSE: th.status_message('Tabulating k-mer model statistics.')
    all_kmer_mean_sds = []
    if _DEBUG_EST_STD:
        kmer_dens = []
        save_x = np.linspace(KERNEL_DENSITY_RANGE[0], KERNEL_DENSITY_RANGE[1],
                             _DEBUG_EST_NUM_KMER_SAVE)
    kmer_width = len(next(iter(all_reg_kmer_levels[0].keys()))[0])
    for kmer, offset in [(kmer, offset - 1)
                         for kmer in product(DNA_BASES, repeat=kmer_width)
                         for offset in motif.find_mod_poss(''.join(kmer))]:
        kmer = ''.join(kmer)
        try:
            kmer_levels = np.concatenate([
                reg_kmer_levels[(kmer, offset)]
                for reg_kmer_levels in all_reg_kmer_levels
                if len(reg_kmer_levels[(kmer, offset)]) > 0])
        except ValueError:
            th.error_message_and_exit(
                'At least one modified k-mer is not covered at any poitions ' +
                'by --minimum-test-reads.\n\t\tConsider fitting to a smaller ' +
                'k-mer via the --upstream-bases and --downstream-bases, ' +
                'or lowering --minimum-test-reads.\n\t\tNote that this may ' +
                'result in a lower quality model.')
        if kmer_levels.shape[0] < min_kmer_obs:
            min_obs = min(
                sum(len(reg_levs[(''.join(kmer), i_offset - 1)])
                    for reg_levs in all_reg_kmer_levels)
                for kmer in product(DNA_BASES, repeat=kmer_width)
                for i_offset in motif.find_mod_poss(''.join(kmer)))
            th.error_message_and_exit(
                'K-mers represeneted in fewer observations than ' +
                'requested in the provided reads. Consider a shorter ' +
                'k-mer or providing more reads.\n\t' + unicode(min_obs) +
                ' observations found in least common kmer.')
        all_kmer_mean_sds.append((kmer, offset, np.median(kmer_levels[:,0]),
                                  np.median(kmer_levels[:,1])))
        if _DEBUG_EST_STD:
            kmer_kde = stats.gaussian_kde(
                kmer_levels[:,0],
                bw_method=_DEBUG_EST_BW / kmer_levels[:,0].std(ddof=1))
            with np.errstate(under='ignore'):
                kmer_dens.append((kmer, offset, kmer_kde.evaluate(save_x)))

    if _DEBUG_EST_STD:
        with io.open('debug_est_motif_alt_ref.density.txt', 'wt') as fp:
            fp.write('Kmer\tOffset\tSignal\tDensity\n')
            fp.write('\n'.join('\t'.join(map(str, (kmer, offset, x, y)))
                               for kmer, offset, dens_i in kmer_dens
                               for x, y in zip(save_x, dens_i)) + '\n')
    return all_kmer_mean_sds

def estimate_motif_alt_model(
        fast5s_dirs, corr_grp, bc_subgrps, motif_desc,
        upstrm_bases, dnstrm_bases, valid_locs_fn, min_kmer_obs,
        cov_thresh, cs_cov_thresh, region_size, num_processes):
    """Estimate a motif-centered alternate-base tombo k-mer model
    """
    reads_index = th.TomboReads(fast5s_dirs, corr_grp, bc_subgrps)
    try:
        raw_motif, mod_pos = motif_desc.split(":")
    except:
        th.error_message_and_exit('Invalid motif decription format.')
    motif = th.TomboMotif(raw_motif, int(mod_pos))
    valid_poss = None if valid_locs_fn is None else th.parse_locs_file(
        valid_locs_fn)

    all_reg_kmer_levels = extract_kmer_levels(
        reads_index, region_size, cov_thresh, upstrm_bases, dnstrm_bases,
        cs_cov_thresh, False, num_processes, motif, valid_poss)

    # process motif kmers with relative position stored
    all_mod_kmer_mean_sds = tabulate_mod_kmer_levels(
        all_reg_kmer_levels, min_kmer_obs, motif)

    alt_ref = AltModel(
        kmer_ref=all_mod_kmer_mean_sds, central_pos=upstrm_bases,
        alt_base=motif.mod_base, motif=motif)

    alt_ref._make_constant_sd()

    return alt_ref

if _PROFILE_MOTIF_ALT_EST:
    _est_motif_alt_wrapper = estimate_motif_alt_model
    def estimate_motif_alt_model(*args):
        import cProfile
        cProfile.runctx('_est_motif_alt_wrapper(*args)', globals(), locals(),
                        filename='est_motif_alt_model.prof')
        return None


####################################
##### Core Statistical Testing #####
####################################

def p_value_to_z_score(pvalue):
    """Helper function to convert p-value to z-score
    """
    return -stats.norm.ppf(pvalue)

def z_score_to_p_value(zscore):
    """Helper function to convert z-score to p-value
    """
    return stats.norm.cdf(zscore)

def correct_multiple_testing(pvals):
    """Use FDR Benjamini-Hochberg multiple testing correction
    """
    pvals = np.asarray(pvals)

    pvals_sortind = np.argsort(pvals)
    pvals_sorted = pvals[pvals_sortind]
    sortrevind = pvals_sortind.argsort()

    nobs = len(pvals)
    ecdffcator = np.arange(1, nobs + 1) / nobs
    # ignore underflow values
    with np.errstate(under='ignore'):
        pvals_corrected_raw = pvals_sorted / ecdffcator
    pvals_corrected = np.minimum.accumulate(
        pvals_corrected_raw[::-1])[::-1]
    del pvals_corrected_raw
    pvals_corrected[pvals_corrected>1] = 1

    return pvals_corrected[sortrevind]

def calc_vectorized_fm_pvals(split_pvals, filter_nan=True):
    """Compute Fisher's Method p-values in a vectorized fashion
    """
    if filter_nan:
        chi_stats = [np.sum(np.log(base_pvals[~np.isnan(base_pvals)])) * -2
                     for base_pvals in split_pvals]
        chi_shapes = [
            np.sum(~np.isnan(base_pvals)) * 2 for base_pvals in split_pvals]
    else:
        chi_stats = [np.sum(np.log(base_pvals)) * -2
                     for base_pvals in split_pvals]
        chi_shapes = [base_pvals.shape[0] * 2
                     for base_pvals in split_pvals]

    f_pvals = stats.chi2.sf(chi_stats, chi_shapes)
    return f_pvals

def calc_window_fishers_method(pvals, lag):
    """Compute Fisher's Method over a moving window across a set of p-values
    """
    assert lag > 0, 'Invalid p-value window provided.'
    width = (lag * 2) + 1
    if pvals.shape[-1] < width:
        raise th.TomboError("P-values vector too short for Fisher's Method " +
                            "window compuation.")
    with np.errstate(invalid='ignore'):
        pvals = np.maximum(pvals, SMALLEST_PVAL)
    log_sums = np.lib.stride_tricks.as_strided(
        np.log(pvals),
        shape=pvals.shape[:-1] + (pvals.shape[-1] - width + 1, width),
        strides=pvals.strides + (pvals.strides[-1],)).sum(-1)
    f_pvals = np.empty(pvals.shape)
    f_pvals[:] = np.NAN
    with np.errstate(invalid='ignore'):
        f_pvals[...,lag:-lag] = stats.chi2.sf(log_sums * -2, width * 2)

    return f_pvals

def calc_window_means(stats, lag):
    """Compute mean over a moving window across a set of statistics
    """
    assert lag > 0, 'Invalid window provided.'
    width = (lag * 2) + 1
    if stats.shape[-1] < width:
        raise th.TomboError(
            "Statistics vector too short for window mean compuation.")
    m_stats = np.empty(stats.shape)
    m_stats[:] = np.NAN
    m_stats[...,lag:-lag] = np.mean(np.lib.stride_tricks.as_strided(
        stats, shape=stats.shape[:-1] + (stats.shape[-1] - width + 1, width),
        strides=stats.strides + (stats.strides[-1],)), -1)

    return m_stats

def calc_window_z_transform(r_means, ref_means, ref_sds, lag):
    """Compute Stouffer's Z-transformation across a read
    """
    z_scores = np.abs(r_means - ref_means) / ref_sds
    width = (lag * 2) + 1
    window_z_trans = np.lib.stride_tricks.as_strided(
        z_scores,
        shape=(z_scores.shape[0] - (lag * 2), width),
        strides=z_scores.strides +
        (z_scores.strides[-1],)).sum(1) / np.sqrt(width)
    window_z_trans = np.concatenate([np.array([np.NAN] * lag),
                                     window_z_trans, np.array([np.NAN] * lag)])

    return window_z_trans

def calc_mann_whitney_z_score(samp1, samp2):
    """Compute Mann-Whitney z-scores comparing two samples
    """
    s1_len = samp1.shape[0]
    s2_len = samp2.shape[0]
    tot_len = s1_len + s2_len

    all_vals = np.concatenate([samp1, samp2])
    ranks = np.empty(tot_len, int)
    ranks[all_vals.argsort()] = np.arange(1, tot_len + 1)
    s1_ranks_sum = ranks[:s1_len].sum()
    #s2_ranks_sum = ranks[s1_len:].sum()

    u1 = s1_ranks_sum - (s1_len * (s1_len + 1)) / 2
    #u2 = s2_ranks_sum - (s2_len * (s2_len + 1)) / 2

    mu = s1_len * s2_len / 2
    rhou = np.sqrt(s1_len * s2_len * (s1_len + s2_len + 1) / 12)

    z = np.abs(u1 - mu) / rhou

    return z

def get_read_seg_score(r_means, r_ref_means, r_ref_sds):
    """Compute expected to observed signal matching score

    Args:
        r_means (`np.array::np.float64`): observed base signal levels
        r_ref_means (`np.array::np.float64`): expected base signal levels
        r_ref_sds (`np.array::np.float64`): expected base signal level sds

    Returns:
        Mean half z-score for observed versus expected signal levels
    """
    return np.mean(np.abs((r_means - r_ref_means) / r_ref_sds))

def score_valid_bases(read_tb, event_means, r_ref_means, r_ref_sds):
    """Compute expected to observed signal matching score for bases not deleted
    in dynamic programming

    Args:
        read_tb (`np.array::np.int32`): event changepoints
        r_means (`np.array::np.float64`): observed base signal levels
        r_ref_means (`np.array::np.float64`): expected base signal levels
        r_ref_sds (`np.array::np.float64`): expected base signal level sds

    Returns:
        Mean half z-score for observed versus expected signal levels (for valid
            bases)
    """
    valid_bases = np.where(np.diff(read_tb) != 0)[0]
    if valid_bases.shape[0] == 0:
        raise th.TomboError('Invalid path through read start')
    valid_ref_means, valid_ref_sds = (
        r_ref_means[valid_bases], r_ref_sds[valid_bases])
    base_means = np.array([event_means[b_start:b_end].mean()
                           for b_start, b_end in zip(read_tb[:-1], read_tb[1:])
                           if b_start != b_end])
    return get_read_seg_score(base_means, valid_ref_means, valid_ref_sds)

def get_dynamic_prog_params(match_evalue):
    """Compute dynamic programming shift parameters from an expected match
    expected value
    """
    z_shift = HALF_NORM_EXPECTED_VAL + match_evalue
    stay_pen = match_evalue
    return z_shift, stay_pen


##########################
##### Statistics I/O #####
##########################

def compute_auc(tp_rate, fp_rate):
    return np.sum(tp_rate[:-1] * (fp_rate[1:] - fp_rate[:-1]))

def compute_mean_avg_precison(tp_rate, precision):
    return np.sum(np.diff(np.concatenate([[0,], tp_rate, [1,]]),) *
                  np.concatenate([[0,], precision, [1,]])[:-1])

def compute_accuracy_rates(stat_has_mod, num_plot_points=ROC_PLOT_POINTS):
    """Given a list or numpy array of true/false values, function returns
    num_plot_point evenly spaced values along the true positive, false positive
    and precision arrays
    """
    tp_cumsum = np.cumsum(stat_has_mod)
    tp_rate = tp_cumsum / tp_cumsum[-1]
    fp_cumsum = np.cumsum(np.logical_not(stat_has_mod))
    fp_rate = fp_cumsum / fp_cumsum[-1]

    precision = tp_cumsum / np.arange(1, len(stat_has_mod) + 1, dtype=float)

    # trim to number of requested points
    tp_rate = tp_rate[np.linspace(0, tp_rate.shape[0] - 1,
                                  num_plot_points).astype(np.int64)]
    fp_rate = fp_rate[np.linspace(0, fp_rate.shape[0] - 1,
                                  num_plot_points).astype(np.int64)]
    precision = precision[np.linspace(0, precision.shape[0] - 1,
                                      num_plot_points + 1).astype(np.int64)][1:]

    return tp_rate, fp_rate, precision

def _compute_motif_stats(
        stats, motif_descs, genome_index,
        stats_per_block=None, total_stats_limit=None):
    all_motif_stats = dict(
        (mod_name, []) for mod_name in list(zip(*motif_descs))[1])
    before_bases = max((
        motif.mod_pos for motif in list(zip(*motif_descs))[0])) - 1
    after_bases = max((motif.motif_len - motif.mod_pos
                       for motif in list(zip(*motif_descs))[0]))
    total_num_stats = 0
    for chrm, strand, start, end, block_stats in stats:
        if strand == '+':
            seq_start = max(start - before_bases, 0)
            seq_end = end + after_bases
        else:
            seq_start = max(start - after_bases, 0)
            seq_end = end + before_bases

        reg_seq = genome_index.get_seq(
            chrm, seq_start, seq_end, error_end=False)
        # TODO potentially keep all mod sites when they are extremely rare
        # randomly sub-sample per-read stats here
        if (stats_per_block is not None and
            block_stats.shape[0] > stats_per_block):
            block_stats = block_stats[np.random.choice(
                block_stats.shape[0], stats_per_block, replace=False)]
        total_num_stats += block_stats.shape[0]
        for r_pos_stat in block_stats:
            # extract position sequence
            if strand == '+':
                r_pos_seq = reg_seq[
                    r_pos_stat['pos'] - seq_start - before_bases:
                    r_pos_stat['pos'] - seq_start + after_bases + 1]
            else:
                r_pos_seq = th.rev_comp(reg_seq[
                    r_pos_stat['pos'] - seq_start - after_bases:
                    r_pos_stat['pos'] - seq_start + before_bases + 1])

            # add statistic and whether the sequence matches each motif
            for motif, mod_name in motif_descs:
                if r_pos_seq[before_bases] != motif.mod_base: continue
                all_motif_stats[mod_name].append((
                    r_pos_stat[stats._stat_slot],
                    bool(motif.motif_pat.match(
                        r_pos_seq[before_bases - motif.mod_pos + 1:]))))

        if (total_stats_limit is not None and
            total_num_stats >= total_stats_limit):
            break

    return all_motif_stats

def _compute_ground_truth_stats(stats, ground_truth_locs):
    all_stats = []
    mod_locs, unmod_locs, mod_name = ground_truth_locs
    for chrm, strand, start, end, block_stats in stats:
        try:
            cs_mod_locs = mod_locs[(chrm, strand)]
        except KeyError:
            cs_mod_locs = np.array([])
        try:
            cs_unmod_locs = unmod_locs[(chrm, strand)]
        except KeyError:
            cs_unmod_locs = np.array([])
        block_mod_locs = cs_mod_locs[
            np.logical_and(np.greater_equal(cs_mod_locs, start),
                           np.less(cs_mod_locs, end))]
        block_unmod_locs = cs_unmod_locs[
            np.logical_and(np.greater_equal(cs_unmod_locs, start),
                           np.less(cs_unmod_locs, end))]
        block_valid_locs = block_stats[
            np.isin(block_stats['pos'], np.concatenate([
                block_mod_locs, block_unmod_locs]))]
        all_stats.extend(zip(
            block_valid_locs[stats._stat_slot],
            np.isin(block_valid_locs['pos'], block_mod_locs)))

    return {mod_name:all_stats}

def _compute_ctrl_motif_stats(
        stats, ctrl_stats, motif_descs, genome_index,
        stats_per_block=None, total_stats_limit=None):
    all_motif_stats = dict(
        (mod_name, []) for mod_name in list(zip(*motif_descs))[1])
    before_bases = max((
        motif.mod_pos for motif in list(zip(*motif_descs))[0])) - 1
    after_bases = max((motif.motif_len - motif.mod_pos
                       for motif in list(zip(*motif_descs))[0]))
    total_num_stats = 0
    for chrm, strand, start, end, block_stats in stats:
        if strand == '+':
            seq_start = max(start - before_bases, 0)
            seq_end = end + after_bases
        else:
            seq_start = max(start - after_bases, 0)
            seq_end = end + before_bases

        reg_seq = genome_index.get_seq(
            chrm, seq_start, seq_end, error_end=False)

        ctrl_block_stats = ctrl_stats.get_reg_stats(chrm, strand, start, end)
        for motif, mod_name in motif_descs:
            if strand == '+':
                mod_poss = np.array([
                    m.start() + motif.mod_pos - 1
                    for m in motif.motif_pat.finditer(reg_seq)]) + seq_start
            else:
                mod_poss = np.array([
                    m.start() + motif.motif_len - motif.mod_pos
                    for m in motif.rev_comp_pat.finditer(reg_seq)]) + seq_start

            # TODO possibly use stats_per_block here, but motif stats are
            # less common, so probably not as useful here
            for r_pos_stat in block_stats[
                    np.isin(block_stats['pos'], mod_poss)]:
                all_motif_stats[mod_name].append((
                    r_pos_stat[stats._stat_slot], True))
                total_num_stats += 1
            if ctrl_block_stats is not None:
                for r_pos_stat in ctrl_block_stats[
                        np.isin(ctrl_block_stats['pos'], mod_poss)]:
                    all_motif_stats[mod_name].append((
                        r_pos_stat[stats._stat_slot], False))
                    total_num_stats += 1

        if (total_stats_limit is not None and
            total_num_stats >= total_stats_limit):
            break

    return all_motif_stats

def calc_damp_fraction(cov_damp_counts, fracs, valid_cov):
    """Compute dampened fraction of un-modified reads using provided modified
    and un-modified pseudo-counts from cov_damp_counts

    See https://nanoporetech.github.io/tombo/text_output.html?highlight=dampened#text-output-browser-files for more details
    """
    damp_fracs = np.empty(fracs.shape[0])
    damp_fracs[:] = np.nan
    non_mod_counts = np.round(fracs * valid_cov)
    # compute dampened fraction of modified reads by adding psuedo-counts
    # to the modified and un-modified counts (equivalent to a beta prior
    # on the fraction estimation as a binomial variable)
    damp_fracs = (non_mod_counts + cov_damp_counts['unmod']) / (
        valid_cov + sum(list(cov_damp_counts.values())))

    return damp_fracs

class ModelStats(object):
    """Parse and retrieve relevant information from a standard (per-genomic
    base) Tombo statistics file.

    .. automethod:: __init__
    """
    # TODO add attributes
    def _parse_stats(self):
        if self.stats_fn is None or not os.path.isfile(self.stats_fn):
            th.error_message_and_exit(
                'Statistics file not provided or provided file does not exist.')

        self._fp = h5py.File(self.stats_fn, 'r')
        self.stat_type = self._fp.attrs.get('stat_type')
        self.region_size = self._fp.attrs.get('block_size')
        self.stat_blocks = self._fp[STAT_BLOCKS_H5_NAME]
        self.num_blocks = 0
        blocks_index = defaultdict(dict)
        for block_name, block_data in self.stat_blocks.items():
            blocks_index[
                (block_data.attrs.get('chrm'), block_data.attrs.get('strand'))][
                    block_data.attrs.get('start')] = block_name
            self.num_blocks += 1
        self.blocks_index = dict(blocks_index)

        self.cov_thresh = self._fp.attrs.get(COV_THRESH_H5_NAME)
        most_signif_grp = self._fp[MOST_SIGNIF_H5_NAME]
        # read full most significant array into memory
        self.most_signif_stats = most_signif_grp[MOST_SIGNIF_H5_NAME][:]
        self.most_signif_chrm_map = dict(
            (v,k) for k,v in most_signif_grp['chrm_ids'].attrs.items())
        # LevelStats doesn't have damp counts
        try:
            self.cov_damp_counts = dict(self._fp[
                COV_DAMP_COUNTS_H5_NAME].attrs.items())
        except:
            self.cov_damp_counts = None

        return

    def _create_new_stats_file(self):
        # try to remove file for overwriting old results
        try:
            os.remove(self.stats_fn)
        except:
            pass
        # open file for writing
        self._fp = h5py.File(self.stats_fn, 'w')

        # save attributes to file and open stats blocks group
        self._fp.attrs['stat_type'] = self.stat_type
        self._fp.attrs['block_size'] = self.region_size
        self.stat_blocks = self._fp.create_group(STAT_BLOCKS_H5_NAME)

        # save coverage damp counts and threshold attributes
        self._fp.attrs[COV_THRESH_H5_NAME] = self.cov_thresh
        self.cov_damp_counts_grp = self._fp.create_group(
            COV_DAMP_COUNTS_H5_NAME)
        self.cov_damp_counts_grp.attrs[
            'unmod'] = self.cov_damp_counts['unmod']
        self.cov_damp_counts_grp.attrs[
            'mod'] = self.cov_damp_counts['mod']

        # storage for most significant stats
        self.most_signif_sites = self._fp.create_group(MOST_SIGNIF_H5_NAME)
        self.running_most_signif_sites = np.empty(
            shape=(self.num_most_signif,),
            dtype=[(str('damp_frac'), 'f8'), (str('frac'), 'f8'),
                   (str('pos'), 'u4'), (str('cov'), 'u4'),
                   (str('control_cov'), 'u4'), (str('valid_cov'), 'u4'),
                   (str('chrm'), 'u4'), (str('strand'), 'S1')])
        self.running_most_signif_sites[:] = np.NAN
        # store a queue of completed stat batches to be concatenated and stored
        # as a group to avoid too many array copy and sorting ops
        self.queued_stat_batches = []
        # store chromosomes names in dict for storing most signif array
        self.curr_chrm_id = 0
        self.chrm_names = {}
        self.chrm_id_grp = self.most_signif_sites.create_group('chrm_ids')

        self.is_empty = True

        return

    def __init__(self, stats_fn, stat_type=None, region_size=None,
                 cov_damp_counts=None, cov_thresh=None, num_most_signif=None,
                 most_signif_num_batches=MOST_SIGNIF_NUM_BATCHES_DEFAULT):
        """Parse or open for writing a standard (per-genomic base) Tombo
        statistics file.

        Example::

            stats = tombo_stats.TomboStats('path/to/stats.file')
            for chrm, strand, pos, frac, damp_frac, valid_cov in \
                    stats.iter_most_signif_sites():
                # do stuff

        Args:
            stats_fn (str): filename for previously saved tombo stats
            stat_type (str): type of statistic (model_compare, de_novo, or
                sample_compare); only applicable for new file writing
            region_size (int): size of chunked storage blocks; only applicable
                for new file writing
            cov_damp_counts (tuple): pseudo-counts for modified and un-modified
                reads to compute ``damp_frac``
            cov_thresh (int): only sites with coverage greater than or equal to
                this value will be stored
            num_most_signif (int): number of most significant sites to be stored
                for faster access
            most_signif_num_batches (int): number of region batches to store
                before re-computing the most significant array (default: 10)

        Warning:

            If all arguments are provided the current file's contents will be
               deleted.

            Intended to open a fresh ``ModelStats`` file for writing.
        """
        self.stats_fn = stats_fn

        if any(arg is None for arg in (stat_type, region_size, cov_damp_counts,
                                       cov_thresh, num_most_signif)):
            self.open_for_writing = False
            # open file for reading
            try:
                self._parse_stats()
            except:
                raise th.TomboError(
                    'Invalid statistics file provided. Try running ' +
                    'tombo/scripts/convert_stats.py if this stats file ' +
                    'was created before Tombo v1.3.1')
        else:
            self.open_for_writing = True
            # set class attributes
            self.stat_type = stat_type
            self.region_size = region_size
            self.curr_block_num = 0
            self.cov_damp_counts = dict(zip(('unmod', 'mod'), cov_damp_counts))
            self.cov_thresh = cov_thresh
            self.num_most_signif = num_most_signif
            self.most_signif_num_batches = most_signif_num_batches
            # open file for writing
            self._create_new_stats_file()

        if self.stat_type not in PER_READ_STATS:
            if self.stat_type in LEVEL_STATS_TXTS:
                raise th.TomboError(
                    'This appears to be a group-comparison stats file. Open ' +
                    'with tombo_stats.LevelStats.')
            raise th.TomboError(
                'This file is not a valid ModelStats file. `stat_type` ' +
                'listed as "' + self.stat_type + '".')
        self.is_model_stats = True
        self._stat_slot = str('damp_frac')
        self._stat_text = 'Est. Frac. Alternate: {0:.2g}'
        self._stat_transform = lambda pos_stat: 1 - pos_stat[self._stat_slot]

        return

    def _update_most_signif(self):
        tmp_most_signif = np.concatenate(
            [self.running_most_signif_sites,] + self.queued_stat_batches)
        tmp_most_signif.sort(kind='mergesort', order=self._stat_slot)
        self.running_most_signif_sites = tmp_most_signif[:self.num_most_signif]
        self.queued_stat_batches = []
        return

    def _add_to_most_signif(self, reg_stats, chrm, strand):
        if chrm not in self.chrm_names:
            self.chrm_names[chrm] = self.curr_chrm_id
            self.curr_chrm_id += 1

        self.queued_stat_batches.append(append_fields(
            base=reg_stats, names=(str('chrm'), str('strand')),
            data=(list(repeat(self.chrm_names[chrm], reg_stats.shape[0])),
                  list(repeat(strand, reg_stats.shape[0]))),
            dtypes=('u4', 'S1')))
        if len(self.queued_stat_batches) >= self.most_signif_num_batches:
            self._update_most_signif()

        return

    def _write_stat_block(self, reg_stats):
        """Write region statistics block to file.
        """
        try:
            block_data = self.stat_blocks.create_group(
                'Block_' + unicode(self.curr_block_num))
            self.curr_block_num += 1
        except:
            th.warning_message('Statistics file not opened for writing.')
            return

        block_data.attrs['chrm'] = reg_stats.chrm
        block_data.attrs['strand'] = reg_stats.strand
        block_data.attrs['start'] = reg_stats.start

        damp_frac = calc_damp_fraction(
            self.cov_damp_counts, reg_stats.reg_frac_standard_base,
            reg_stats.valid_cov)
        reg_stats_arr = np.array(
            [pos_stats for pos_stats in zip(
                damp_frac, reg_stats.reg_frac_standard_base,
                reg_stats.reg_poss, reg_stats.reg_cov,
                reg_stats.ctrl_cov, reg_stats.valid_cov)
             if not np.isnan(pos_stats[0])],
            dtype=[
                (str('damp_frac'), 'f8'), (str('frac'), 'f8'),
                (str('pos'), 'u4'), (str('cov'), 'u4'),
                (str('control_cov'), 'u4'), (str('valid_cov'), 'u4')])
        block_data.create_dataset(
            'block_stats', data=reg_stats_arr, compression="gzip")

        self._add_to_most_signif(reg_stats_arr, reg_stats.chrm, reg_stats.strand)

        #self._fp.flush()
        self.is_empty = False

        return

    def _close_write(self):
        # process any remaining batches
        if len(self.queued_stat_batches) >= 1:
            self._update_most_signif()
        # trim the array if necessary
        if np.isnan(self.running_most_signif_sites[self._stat_slot][-1]):
            # not as many signif sites were stored as requested so trim array
            first_nan = np.where(np.isnan(
                self.running_most_signif_sites[self._stat_slot]))[0][0]
            self.running_most_signif_sites = self.running_most_signif_sites[
                :first_nan,]
        # add dataset to file
        self.most_signif_sites.create_dataset(
            MOST_SIGNIF_H5_NAME, data=self.running_most_signif_sites,
            compression="gzip")
        # and add chrm ids map to file (store in reverse order of useful dict,
        # since int's can't be hdf5 keys
        for chrm_name, chrm_id in self.chrm_names.items():
            self.chrm_id_grp.attrs[chrm_name] = chrm_id

        return

    def close(self):
        """Close open HDF5 file and write most significant sites if open for
        writing
        """
        if self.open_for_writing:
            self._close_write()
        self._fp.close()
        return


    # Reading functions
    def _get_chrm_name(self, pos_stat):
        return self.most_signif_chrm_map[pos_stat['chrm']]

    def iter_stat_seqs(
            self, genome_index, before_bases, after_bases, include_pos=True):
        """Iterate through most significant genomic sites returning the genomic
        sequence surrounding each position.

        Args:
            genome_index (:class:`tombo.tombo_helper.Fasta`): genome index
                object
            before_bases (int): number of sequence bases before positions to
                include
            after_bases (int): number of sequence bases after positions to
                include
            include_pos (bool): yeild (pos_seq, chrm, strand, start, end) for
                each site (default: True)
        """
        for pos_stat in self.most_signif_stats:
            chrm, strand, pos = (self._get_chrm_name(pos_stat),
                                 pos_stat['strand'].decode(),
                                 pos_stat['pos'])
            if strand == '+':
                start, end = pos - before_bases, pos + after_bases + 1
                try:
                    stat_seq = genome_index.get_seq(chrm, start, end)
                except th.TomboError:
                    continue
            else:
                start, end = pos - after_bases, pos + before_bases + 1
                try:
                    stat_seq = th.rev_comp(genome_index.get_seq(
                        chrm, start, end))
                except th.TomboError:
                    continue
            if include_pos:
                yield stat_seq, chrm, strand, start, end
            else:
                yield stat_seq

        return

    def iter_most_signif_sites(self):
        """Iterate through statistics table yeilding (chrm, strand, pos, stat).
        """
        for pos_stat in self.most_signif_stats:
            yield (
                self._get_chrm_name(pos_stat), pos_stat['strand'].decode(),
                pos_stat['pos'],
                self._stat_transform(pos_stat[self._stat_slot]))

        return

    def get_most_signif_regions(
            self, num_bases, num_regions, unique_pos=True,
            prepend_loc_to_text=False):
        """Select regions centered on locations with the largest fraction of
        modified bases

        Args:
            num_bases (int): number of bases to output
            num_regions (int): number of regions to output
            unique_pos (bool): get only unique positions (optional; default
                True) intervals may overlap, but identified significant
                position is outside other intervals
            prepend_loc_to_text (bool): pre-prend most significant location to
               the region text (can be off for interval near start/end of
               sequence records)

        Returns:
            A list of :class:`tombo.tombo_helper.intervalData` objects
        """
        selected_regs = []
        used_intervals = defaultdict(set)
        for i, pos_stat in enumerate(self.most_signif_stats):
            int_start = max(0, pos_stat['pos'] - int(num_bases / 2.0))
            chrm = self._get_chrm_name(pos_stat)
            strand = pos_stat['strand'].decode()
            if (not unique_pos or
                pos_stat['pos'] not in used_intervals[(chrm, strand)]):
                used_intervals[(chrm, strand)].update(
                    range(int_start, int_start + num_bases))
                int_text = self._stat_text.format(
                    self._stat_transform(pos_stat))
                if prepend_loc_to_text:
                    int_text = '{0}:{1:d}:{2}'.format(
                        chrm, pos_stat['pos'] + 1, strand) + " " + int_text
                selected_regs.append(th.intervalData(
                    chrm=chrm, start=int_start, end=int_start + num_bases,
                    strand=strand, reg_id='{:03d}'.format(i),
                    reg_text=int_text))
                if len(selected_regs) >= num_regions: break

        if len(selected_regs) == 0:
            th.error_message_and_exit(
                'No locations identified. Most likely an empty ' +
                'statistics file.')
        if len(selected_regs) < num_regions:
            th.warning_message(
                'Fewer unique significant locations more than ' +
                '[--num-bases]/2 apart were identified. Continuing with ' +
                str(len(selected_regs)) + ' unique locations. Must raise ' +
                '--num-most-significant-stored in order to see more most ' +
                'significant stats.')

        return selected_regs

    def compute_motif_stats(
            self, motif_descs, genome_index,
            stats_per_block=None, total_stats_limit=None):
        """Compute lists of statistic values and whether this site represents a
        match to the provided motifs

        Args:
            motif_descs (list; see
                :class:`tombo.tombo_helper.parse_motif_descs`): containing
                tuples with :class:`tombo.tombo_helper.TomboMotif` and
                motif/modification names
            genome_index (:class:`tombo.tombo_helper.Fasta`): genome index
            stats_per_block (int): statistics to include in calculations
                per-block (`--multiprocess-region-size`)
            total_stats_limit (int): maximum total statistics to include in
               computation (Default: include all stats)

        Returns:
            Dictionary with (key) motif/modification name and (value) list of
            tuples containing statistic value and boolean motif match
        """
        return _compute_motif_stats(
            self, motif_descs, genome_index, stats_per_block=stats_per_block,
            total_stats_limit=total_stats_limit)

    def compute_ground_truth_stats(self, ground_truth_locs):
        """Compute lists of statistic values and ground truth modification
        status (boolean)

        Args:
            ground_truth_locs: list containing tuples of 1) modified locations
                (from :class:`tombo.tombo_helper.parse_locs_file`), 2)
                unmodified locations and mod names.

        Returns:
            Dictionary with (key) motif/modification name and (value) list of
            tuples containing statistic value and boolean motif match
        """
        return _compute_ground_truth_stats(self, ground_truth_locs)

    def compute_ctrl_motif_stats(
            self, ctrl_stats, motif_descs, genome_index,
            stats_per_block=None, total_stats_limit=None):
        """Compute lists of statistic values and whether this site represents a
        match to the provided motifs

        Args:
            ctrl_stats (:class:`tombo.tombo_stats.ModelStats`): control
                statistics
            motif_descs (list; see
                :class:`tombo.tombo_helper.parse_motif_descs`): containing
                tuples with :class:`tombo.tombo_helper.TomboMotif` and
                motif/modification names
            genome_index (:class:`tombo.tombo_helper.Fasta`): genome index
            stats_per_block (int): statistics to include in calculations
                per-block (`--multiprocess-region-size`)
            total_stats_limit (int): maximum total statistics to include in
                computation (Default: include all stats)

        Returns:
            Dictionary with (key) motif/modification name and (value) list of
            tuples containing statistic value and boolean native or control
            sample stats
        """
        return _compute_ctrl_motif_stats(
            self, ctrl_stats, motif_descs, genome_index,
            stats_per_block=stats_per_block,
            total_stats_limit=total_stats_limit)

    def __iter__(self):
        """Iterator over all statistics blocks, yeilding chrm, strand, start,
        end, block_stats
        """
        self.iter_all_cs = iter(sorted(self.blocks_index))
        try:
            self.iter_curr_cs = next(self.iter_all_cs)
        except StopIteration:
            self.iter_curr_cs = None
            self.iter_curr_cs_blocks = iter([])
        else:
            self.iter_curr_cs_blocks = iter(
                self.blocks_index[self.iter_curr_cs].items())

        return self

    def __next__(self):
        try:
            next_start, next_block_name = next(self.iter_curr_cs_blocks)
        except StopIteration:
            # move to next chromosome and strand
            # this will raise a second StopIteration
            # when the end of the blocks is hit
            self.iter_curr_cs = next(self.iter_all_cs)
            self.iter_curr_cs_blocks = iter(sorted(
                self.blocks_index[self.iter_curr_cs].items()))
            next_start, next_block_name = next(self.iter_curr_cs_blocks)

        chrm, strand = self.iter_curr_cs
        return (chrm, strand, next_start, next_start + self.region_size,
                self.stat_blocks[next_block_name]['block_stats'][:])

    # for python2 compatibility
    def next(self):
        """Return next statistics block from file including (chrm, strand,
        block start, block end and statistics table ``numpy structured array``)
        """
        return self.__next__()

    def get_pos_stat(self, chrm, strand, pos, missing_value=None):
        """Extract statistic value from the requested genomic position.
        """
        try:
            pos_block_start = np.floor_divide(
                pos, self.region_size) * self.region_size
            # TODO: blocks may have missing data (consider full sized blocks
            # for random disk access to single or range of elements)
            #block_pos = np.remainder(pos, self.region_size)
            block_name = self.blocks_index[(chrm, strand)][pos_block_start]
            block_data = self.stat_blocks[block_name]['block_stats'][:]
            pos_index = np.where(block_data['pos'] == pos)[0]
            if len(pos_index) != 1: raise KeyError
            pos_stat = self._stat_transform(block_data[pos_index[0]])
        except KeyError:
            pos_stat = missing_value

        return pos_stat

    def get_reg_stats(self, chrm, strand, start, end):
        if (chrm, strand) not in self.blocks_index: return
        reg_stats = []
        for reg_start, block_name in sorted(
                self.blocks_index[(chrm, strand)].items()):
            # if this is an overlapping block
            if reg_start < end and reg_start + self.region_size > start:
                # extract stats overlapping this region
                reg_stats_block = self.stat_blocks[block_name]['block_stats'][:]
                reg_stats.append(reg_stats_block[
                    np.logical_and(
                        np.greater_equal(reg_stats_block['pos'], start),
                        np.less(reg_stats_block['pos'], end))])

        if len(reg_stats) == 0:
            return
        elif len(reg_stats) == 1:
            return reg_stats[0]
        return np.vstack(reg_stats)


class LevelStats(ModelStats):
    def _create_new_stats_file(self):
        # try to remove file for overwriting old results
        try:
            os.remove(self.stats_fn)
        except:
            pass
        # open file for writing
        self._fp = h5py.File(self.stats_fn, 'w')

        # save attributes to file and open stats blocks group
        self._fp.attrs['stat_type'] = self.stat_type
        self._fp.attrs['block_size'] = self.region_size
        self.stat_blocks = self._fp.create_group(STAT_BLOCKS_H5_NAME)

        # save coverage damp counts and threshold attributes
        self._fp.attrs[COV_THRESH_H5_NAME] = self.cov_thresh

        # storage for most significant stats
        self.most_signif_sites = self._fp.create_group(MOST_SIGNIF_H5_NAME)
        self.running_most_signif_sites = np.empty(
            shape=(self.num_most_signif,),
            dtype=[(str('stat'), 'f8'), (str('pos'), 'u4'),
                   (str('cov'), 'u4'), (str('control_cov'), 'u4'),
                   (str('chrm'), 'u4'), (str('strand'), 'S1')])
        self.running_most_signif_sites[:] = np.NAN
        # store a queue of completed stat batches to be concatenated and stored
        # as a group to avoid too many array copy and sorting ops
        self.queued_stat_batches = []
        # store chromosomes names in dict for storing most signif array
        self.curr_chrm_id = 0
        self.chrm_names = {}
        self.chrm_id_grp = self.most_signif_sites.create_group('chrm_ids')

        self.is_empty = True

        return

    def __init__(self, stats_fn, stat_type=None, region_size=None,
                 cov_thresh=None, num_most_signif=None,
                 most_signif_num_batches=MOST_SIGNIF_NUM_BATCHES_DEFAULT):
        """Parse or open for writing a standard (per-genomic base) Tombo
        statistics file.

        Example::

            stats = tombo_stats.TomboStats('path/to/stats.file')
            for chrm, strand, pos, nl_pval in \
                    stats.iter_most_signif_sites():
                # do stuff

        Args:
            stats_fn (str): filename for previously saved tombo stats
            stat_type (str): type of statistic (ks, ttest, utest); only
                applicable for new file writing
            region_size (int): size of chunked storage blocks; only applicable
                for new file writing
            cov_thresh (int): only sites with coverage greater than or equal to
                this value will be stored
            num_most_signif (int): number of most significant sites to be stored
                for faster access
            most_signif_num_batches (int): number of region batches to store
                before re-computing the most significant array (default: 10)

        Warning:

            If all arguments are provided the current file's contents will be
               deleted.

            Intended to open a fresh ``ModelStats`` file for writing.
        """
        self.stats_fn = stats_fn

        if any(arg is None for arg in (
                stat_type, region_size, cov_thresh, num_most_signif)):
            self.open_for_writing = False
            # open file for reading
            try:
                self._parse_stats()
            except:
                raise th.TomboError(
                    'Invalid statistics file provided. Try running ' +
                    'tombo/scripts/convert_stats.py if this stats file ' +
                    'was created before Tombo v1.3.1')
        else:
            self.open_for_writing = True
            # set class attributes
            self.stat_type = stat_type
            self.region_size = region_size
            self.curr_block_num = 0
            self.cov_thresh = cov_thresh
            self.num_most_signif = num_most_signif
            self.most_signif_num_batches = most_signif_num_batches
            # open file for writing
            self._create_new_stats_file()

        if self.stat_type not in LEVEL_STATS_TXTS:
            if self.stat_type in PER_READ_STATS:
                raise th.TomboError(
                    'This appears to be a model-based comparison stats ' +
                    'file. Open with tombo_stats.ModelStats.')
            raise th.TomboError(
                'This file is not a valid LevelStats file. `stat_type` ' +
                'listed as "' + self.stat_type + '".')

        # set values to access statistics consistently
        self.is_model_stats = False
        self._stat_slot = str('stat')
        if self.stat_type in (KS_TEST_TXT, U_TEST_TXT, T_TEST_TXT):
            self._stat_text = '-log10(p-value): {0:.2g}'
            self._stat_transform = lambda pos_stat: -np.log10(
                pos_stat[self._stat_slot])
        elif self.stat_type == KS_STAT_TEST_TXT:
            self._stat_text = 'D Statistic: {0:.2g}'
            self._stat_transform = lambda pos_stat: (
                1 - pos_stat[self._stat_slot])
        elif self.stat_type == U_STAT_TEST_TXT:
            self._stat_text = 'Common Language Marginal Effect: {0:.2g}'
            self._stat_transform = lambda pos_stat: -pos_stat[self._stat_slot]
        elif self.stat_type == T_STAT_TEST_TXT:
            self._stat_text = "Cohen's D: {0:.2g}"
            self._stat_transform = lambda pos_stat: -pos_stat[self._stat_slot]
        else:
            raise th.TomboError('Unknown statistic type.')

        return

    def _write_stat_block(self, grp_stats):
        """Write region group statistics block to file.
        """
        try:
            block_data = self.stat_blocks.create_group(
                'Block_' + unicode(self.curr_block_num))
            self.curr_block_num += 1
        except:
            th.warning_message('Statistics file not opened for writing.')
            return

        block_data.attrs['chrm'] = grp_stats.chrm
        block_data.attrs['strand'] = grp_stats.strand
        block_data.attrs['start'] = grp_stats.start

        grp_stats_arr = np.array(
            [pos_stats for pos_stats in zip(
                grp_stats.reg_stats, grp_stats.reg_poss,
                grp_stats.reg_cov, grp_stats.ctrl_cov)
             if not np.isnan(pos_stats[0])],
            dtype=[(str('stat'), 'f8'), (str('pos'), 'u4'),
                   (str('cov'), 'u4'), (str('control_cov'), 'u4')])
        block_data.create_dataset(
            'block_stats', data=grp_stats_arr, compression="gzip")

        self._add_to_most_signif(
            grp_stats_arr, grp_stats.chrm, grp_stats.strand)

        self.is_empty = False

        return

def TomboStats(stat_fn):
    """Load per-reference location Tombo statistics. Safe method to load statistics for read-only purposes.

    Returns:
        Either :class:`tombo.tombo_stats.ModelStats` or :class:`tombo.tombo_stats.LevelStats` depending on the type of Tombo statistics file provided.
    """
    try:
        stats = ModelStats(stat_fn)
    except th.TomboError:
        stats = LevelStats(stat_fn)

    return stats

class PerReadStats(object):
    """Store and accses per-read modified base testing statistics

    .. automethod:: __init__
    """
    # TODO add attributes
    def _parse_per_read_stats(self):
        self._fp = h5py.File(self.per_read_stats_fn, 'r')
        self.stat_type = self._fp.attrs.get('stat_type')
        self.region_size = self._fp.attrs.get('block_size')
        self.per_read_blocks = self._fp[STAT_BLOCKS_H5_NAME]
        self.num_blocks = 0
        blocks_index = defaultdict(dict)
        for block_name, block_data in self.per_read_blocks.items():
            blocks_index[
                (block_data.attrs.get('chrm'), block_data.attrs.get('strand'))][
                    block_data.attrs.get('start')] = block_name
            self.num_blocks += 1
        self.blocks_index = dict(blocks_index)

        return

    def _create_new_per_read_stats_file(self):
        # try to remove file for overwriting old results
        try:
            os.remove(self.per_read_stats_fn)
        except:
            pass
        # open file for writing
        self._fp = h5py.File(self.per_read_stats_fn, 'w')

        # save attributes to file and open stats blocks group
        self.curr_block_num = 0
        self._fp.attrs['stat_type'] = self.stat_type
        self._fp.attrs['block_size'] = self.region_size
        self.per_read_blocks = self._fp.create_group(STAT_BLOCKS_H5_NAME)

        return

    def __init__(self, per_read_stats_fn, stat_type=None, region_size=None):
        """Open per-read statistics file.

        Examples::

            per_read_stats = tombo_stats.PerReadStats(\
                'path/to/sample.tombo.per_read_stats')
            int_data = tombo_helper.intervalData(
                chrm='chr20', start=10000, end=10100, strand='+')
            reg_per_read_stats = per_read_stats.get_region_per_read_stats(
                int_data, num_reads=10)

        Args:

            per_read_stats_fn (str): filename containing (or to write) per-read
                Tombo statistics
            stat_type (str): type of statistic (model_compare, de_novo, or
                sample_compare); only applicable for new file writing
            region_size (int): size of chunked storage blocks; only applicable
                for new file writing

        Warning:

            If ``stat_type`` and ``region_size`` are provided the current
            file's contents will be deleted.

            Intended to open a fresh ``PerReadStats`` file for writing.
        """
        self.per_read_stats_fn = per_read_stats_fn
        if stat_type is None or region_size is None:
            # open file for reading
            try:
                self._parse_per_read_stats()
            except:
                th.error_message_and_exit(
                    'Non-existent or invalid per-read statistics file ' +
                    'provided.')
        else:
            # set class attributes
            self.stat_type = stat_type
            self.region_size = region_size

            self._create_new_per_read_stats_file()

        self.are_pvals = self.stat_type != ALT_MODEL_TXT

        self._stat_slot = str('stat')
        self._stat_text = '-log10(p-value): {0:.2g}'
        self._stat_transform = lambda pos_stat: -np.log10(
            pos_stat[self._stat_slot])

        return

    def _write_per_read_block(
            self, per_read_block, read_id_lookup, chrm, strand, start):
        """Write region statistics block to file.
        """
        try:
            block_data = self.per_read_blocks.create_group(
                'Block_' + unicode(self.curr_block_num))
            self.curr_block_num += 1
        except:
            th.warning_message(
                'Per-read statistics file not opened for writing.')
            return

        block_data.attrs['chrm'] = chrm
        block_data.attrs['strand'] = strand
        block_data.attrs['start'] = start
        block_data.create_dataset(
            'block_stats', data=per_read_block, compression="gzip")
        # add lookup dict for read_id slot stored in table to save space and
        # avoid memory leak due to vlen slots in h5py datasets
        dt = h5py.special_dtype(vlen=str)
        read_ids = np.array(list(read_id_lookup.keys()), dtype=dt)
        read_ids_ds = block_data.create_dataset(
            'read_ids', read_ids.shape, dtype=dt, compression="gzip")
        read_ids_ds[...] = read_ids
        block_data.create_dataset(
            'read_id_vals', data=np.array(list(read_id_lookup.values())),
            compression="gzip")

        self._fp.flush()

        return

    def get_region_per_read_stats(self, interval_data, num_reads=None):
        """Extract per-read statistics over the specifed interval.

        Args:

            interval_data (:class:`tombo.tombo_helper.intervalData`):
                genomic interval
            num_reads (int): randomly select this many reads (default: inlcude
                all reads)

        Returns:
            `np.array` structured array containing ``pos``, ``stat`` and
            ``read_id`` for per-read stats over requested interval
        """
        try:
            cs_blocks = self.blocks_index[(
                interval_data.chrm, interval_data.strand)]
        except KeyError:
            return

        int_block_stats = []
        for block_start, block_name in cs_blocks.items():
            if (interval_data.end < block_start or
                interval_data.start > block_start + self.region_size): continue
            # extract stats from FAST5
            block_stats = self.per_read_blocks[block_name]['block_stats'][:]
            reg_poss = block_stats['pos']
            reg_read_stats = block_stats[self._stat_slot]
            # extract and convert read_ids back into strings
            if 'read_id_vals' in self.per_read_blocks[block_name]:
                block_read_id_lookup = dict([
                    (read_id_val, read_id) for read_id, read_id_val in
                    zip(self.per_read_blocks[block_name]['read_ids'].value,
                        self.per_read_blocks[block_name]['read_id_vals'].value)])
            else:
                # read_ids previously stored (inefficiently) as attributes
                # so parse read_ids attributes for backwards compatibility
                block_read_id_lookup = dict([
                    (read_id_val, read_id) for read_id, read_id_val in
                    self.per_read_blocks[block_name]['read_ids'].attrs.items()])
            reg_read_ids = [
                block_read_id_lookup[r_id] for r_id in block_stats['read_id']]
            int_block_stats.append(np.array(
                list(zip(reg_poss, reg_read_stats, reg_read_ids)),
                dtype=[(str('pos'), 'u4'), (str('stat'), 'f8'),
                       (str('read_id'), object)]))

        if len(int_block_stats) == 0:
            return

        if len(int_block_stats) == 1:
            int_block_stats = int_block_stats[0]
        else:
            int_block_stats = np.concatenate(int_block_stats)

        all_int_stats = int_block_stats[
            (int_block_stats['pos'] >= interval_data.start) &
            (int_block_stats['pos'] < interval_data.end)]

        read_ids = set(all_int_stats['read_id'])
        if num_reads is not None and num_reads < len(read_ids):
            int_plot_reads = set(random.sample(read_ids, num_reads))
            all_int_stats = all_int_stats[
                np.array([read_id in int_plot_reads
                          for read_id in all_int_stats['read_id']])]

        return all_int_stats

    def compute_motif_stats(
            self, motif_descs, genome_index,
            stats_per_block=None, total_stats_limit=None):
        """Compute lists of statistic values and whether this site represents a
        match to the provided motifs

        Args:
            motif_descs (list; see
                :class:`tombo.tombo_helper.parse_motif_descs`): containing
                tuples with :class:`tombo.tombo_helper.TomboMotif` and
                motif/modification names
            genome_index (:class:`tombo.tombo_helper.Fasta`): genome index
            stats_per_block (int): statistics to include in calculations
                per-block (`--multiprocess-region-size`)
            total_stats_limit (int): maximum total statistics to include in
                computation (Default: include all stats)

        Returns:
            Dictionary with (key) motif/modification name and (value) list of
            tuples containing statistic value and boolean motif match
        """
        return _compute_motif_stats(
            self, motif_descs, genome_index, stats_per_block=stats_per_block,
            total_stats_limit=total_stats_limit)

    def compute_ground_truth_stats(self, ground_truth_locs):
        """Compute lists of statistic values and ground truth modification
        status (boolean)

        Args:
            ground_truth_locs: list containing tuples of 1) modified locations
                (from :class:`tombo.tombo_helper.parse_locs_file`), 2)
                unmodified locations and mod names.

        Returns:
            Dictionary with (key) motif/modification name and (value) list of
            tuples containing statistic value and boolean motif match
        """
        return _compute_ground_truth_stats(self, ground_truth_locs)

    def compute_ctrl_motif_stats(
            self, pr_ctrl_stats, motif_descs, genome_index,
            stats_per_block=None, total_stats_limit=None):
        """Compute lists of statistic values and whether this site represents a
        match to the provided motifs

        Args:
            pr_ctrl_stats (:class:`tombo.tombo_stats.PerReadStats`): control
                per-read statistics
            motif_descs (list; see
                :class:`tombo.tombo_helper.parse_motif_descs`): containing
                tuples with :class:`tombo.tombo_helper.TomboMotif` and
                motif/modification names
            genome_index (:class:`tombo.tombo_helper.Fasta`): genome index
            stats_per_block (int): statistics to include in calculations
                per-block (`--multiprocess-region-size`)
            total_stats_limit (int): maximum total statistics to include in
                computation (Default: include all stats)

        Returns:
            Dictionary with (key) motif/modification name and (value) list of
            tuples containing statistic value and boolean native or control
            sample stats
        """
        return _compute_ctrl_motif_stats(
            self, pr_ctrl_stats, motif_descs, genome_index,
            stats_per_block=stats_per_block,
            total_stats_limit=total_stats_limit)

    def __iter__(self):
        """
        Iterator over all statistics blocks, yeilding chrm, strand,
        start, end, block_stats
        """
        self.iter_all_cs = iter(list(self.blocks_index))
        self.iter_curr_cs = next(self.iter_all_cs)
        self.iter_curr_cs_blocks = iter(
            self.blocks_index[self.iter_curr_cs].items())
        return self

    def __next__(self):
        try:
            next_start, next_block_name = next(self.iter_curr_cs_blocks)
        except StopIteration:
            # move to next chromosome and strand
            # this will raise a second StopIteration
            # when the end of the blocks is hit
            self.iter_curr_cs = next(self.iter_all_cs)
            self.iter_curr_cs_blocks = iter(sorted(
                self.blocks_index[self.iter_curr_cs].items()))
            next_start, next_block_name = next(self.iter_curr_cs_blocks)

        chrm, strand = self.iter_curr_cs
        return (chrm, strand, next_start, next_start + self.region_size,
                self.per_read_blocks[next_block_name]['block_stats'][:])

    # for python2 compatibility
    def next(self):
        """Return next per-read statistics block from file including (chrm,
        strand, block start, block end and per-read statistics table
        ``numpy structured array``)
        """
        return self.__next__()

    def close(self):
        """Close HDF5 file
        """
        self._fp.close()
        return

    def get_reg_stats(self, chrm, strand, start, end):
        if (chrm, strand) not in self.blocks_index: return
        reg_stats = []
        for reg_start, block_name in sorted(
                self.blocks_index[(chrm, strand)].items()):
            # if this is an overlapping block
            if reg_start < end and reg_start + self.region_size > start:
                # extract stats overlapping this region
                reg_stats_block = self.per_read_blocks[
                    block_name]['block_stats'][:]
                reg_stats.append(reg_stats_block[
                    np.logical_and(
                        np.greater_equal(reg_stats_block['pos'], start),
                        np.less(reg_stats_block['pos'], end))])

        if len(reg_stats) == 0:
            return
        elif len(reg_stats) == 1:
            return reg_stats[0]
        return np.vstack(reg_stats)


################################
##### Base-by-base Testing #####
################################

def compute_posterior_samp_dists(
        ctrl_means, ctrl_sds, ctrl_cov, ctrl_reg_data, std_ref,
        prior_weights, min_test_reads, fm_offset):
    dnstrm_bases = std_ref.kmer_width - std_ref.central_pos - 1
    gnm_begin_lag = (
        std_ref.central_pos if ctrl_reg_data.strand == '+' else dnstrm_bases)
    gnm_end_lag = (
        dnstrm_bases if ctrl_reg_data.strand == '+' else std_ref.central_pos)
    reg_seq = ctrl_reg_data.copy().update(
        start=ctrl_reg_data.start - gnm_begin_lag - fm_offset,
        end=ctrl_reg_data.end + gnm_end_lag + fm_offset).add_seq().seq
    if ctrl_reg_data.strand == '-':
        reg_seq = th.rev_comp(reg_seq)

    reg_ref_means, reg_ref_sds = std_ref.get_exp_levels_from_seq_with_gaps(
        reg_seq, ctrl_reg_data.strand == '-')

    # compute vectorized weighted means for new mean and sd estimates
    post_ref_means = ((
        (prior_weights[0] * reg_ref_means) + (ctrl_cov * ctrl_means)) /
                      (prior_weights[0] + ctrl_cov))
    post_ref_sds = ((
        (prior_weights[1] * reg_ref_sds) + (ctrl_cov * ctrl_sds)) /
                    (prior_weights[1] + ctrl_cov))

    # This bit should work, but the SD estimates seem to be incorrect
    # and the computation is likely far too much for likely very similar
    # results from a weighted average with the prior SD.
    """
    # compute posterior sds; see example formula here:
    # https://www.statlect.com/fundamentals-of-statistics/\
    #     normal-distribution-Bayesian-estimation
    # optimizations to matrix ops applied here
    def compute_mean_diff_factor(event_means, ref_mean):
        valid_indices = ~np.isnan(event_means)
        num_valid_indices = sum(valid_indices)
        if num_valid_indices < min_test_reads:
            return np.NAN
        n = float(num_valid_indices + prior_weights[0])
        c1 = (n - 1) / n
        c2 = -1 / n
        mean_diffs = event_means[valid_indices] - ref_mean
        mds_sum = mean_diffs.sum()
        return sum(((i_md * c1) + (mds_sum - i_md) * c2) * i_md
                   for i, i_md in enumerate(mean_diffs))
    mean_diff_factors = np.array([
        compute_mean_diff_factor(b_events, ref_mean)
        for b_events, ref_mean in zip(ctrl_base_events, reg_ref_means)])
    post_ref_sds = np.sqrt((
        mean_diff_factors + (prior_weights[1] * np.square(reg_ref_sds))) / (
            ctrl_cov + prior_weights[1]))
    """

    return post_ref_means, post_ref_sds

def get_reads_ref(
        reg_data, min_test_reads, fm_offset, std_ref=None,
        prior_weights=None, est_mean=False):
    """Get mean and standard deviation of levels from a sample across the genome
    """
    central_func = np.mean if est_mean else np.median
    reg_size = reg_data.end - reg_data.start + (fm_offset * 2)
    level_means, level_sds = np.empty(reg_size), np.empty(reg_size)
    level_means[:] = np.NAN
    level_sds[:] = np.NAN

    # expand region to include fm_offset
    bases_levels = reg_data.copy().update(
        start=reg_data.start - fm_offset,
        end=reg_data.end + fm_offset).get_base_levels()
    valid_indices = np.logical_not(np.isnan(bases_levels))
    cov = valid_indices.sum(axis=1)
    cov_regs = np.where(np.diff(np.concatenate([
        [False,], np.greater_equal(cov, min_test_reads), [False,]])))[0]
    if len(cov_regs) == 0:
        return level_means, level_sds, {}

    for cov_start, cov_end in zip(cov_regs[:-1:2], cov_regs[1::2]):
        level_means[cov_start:cov_end] = np.array([
            central_func(b_levels[valid_indices[cov_start + i]])
            for i, b_levels in enumerate(bases_levels[cov_start:cov_end])])
        level_sds[cov_start:cov_end] = np.array([
            np.std(b_levels[valid_indices[cov_start + i]])
            for i, b_levels in enumerate(bases_levels[cov_start:cov_end])])

    if std_ref is not None:
        if prior_weights is None:
            prior_weights = (MEAN_PRIOR_CONST, SD_PRIOR_CONST)
        level_means, level_sds = compute_posterior_samp_dists(
            level_means, level_sds, cov, reg_data, std_ref,
            prior_weights, min_test_reads, fm_offset)

    # convert coverage to a dict for later lookup
    cov = dict(zip(range(reg_data.start - fm_offset,
                         reg_data.end + fm_offset), cov))

    # mask regions with 0 estimated SD (likely very low coverage)
    zero_sd_sites = level_sds == 0
    level_means[zero_sd_sites] = np.NAN
    level_sds[zero_sd_sites] = np.NAN

    return level_means, level_sds, cov

def compute_sample_compare_read_stats(
        r_data, ctrl_means, ctrl_sds, fm_offset=FM_OFFSET_DEFAULT,
        reg_data=None):
    """Compute signficance statistics using comparison of two sequenceing
    samples method for a single read within a specified genomic region.

    Args:

        r_data (:class:`tombo.tombo_helper.readData`): read data
        ctrl_means (`np.array::np.float64`): mean level values from control set
            of reads
        ctrl_sds (`np.array::np.float64`): level SD values from control set of
            reads
        fm_offset (int): Fisher's Method offset for computing locally combined
            p-values (optional; default: 1)
        reg_data (:class:`tombo.tombo_helper.intervalData`): region to test
            (default: test whole read)

    Returns:
        Read testing results, positions tested and the read_id

        1) r_pvals (`np.array::np.float64`): p-values for testing over specified region
        2) r_poss (`np.array::np.int64`): genomic positions for returned p-values
        3) read_id (str): read identifier
    """
    reg_start = reg_data.start if reg_data is not None else r_data.start
    reg_size = (reg_data.end - reg_data.start if reg_data is not None
                else r_data.end - r_data.start)

    def comp_clip_and_flip():
        with h5py.File(r_data.fn, 'r') as fast5_data:
            r_means = th.get_single_slot_read_centric(
                fast5_data, 'norm_mean', r_data.corr_group)
            read_id = th.get_raw_read_slot(fast5_data).attrs.get('read_id')
        if r_means is None:
            raise th.TomboError(
                'Read does not contain re-squiggled level means.')

        read_start, read_end = r_data.start, r_data.end
        if read_start + fm_offset < reg_start:
            num_start_clip = reg_start - (read_start + fm_offset)
            read_start = reg_start - fm_offset
            if r_data.strand == '+':
                r_means = r_means[num_start_clip:]
            else:
                r_means = r_means[:-num_start_clip]
        if read_end - fm_offset > reg_start + reg_size:
            num_end_clip = (read_end - fm_offset) - (reg_start + reg_size)
            read_end = reg_start + reg_size + fm_offset
            if r_data.strand == '+':
                r_means = r_means[:-num_end_clip]
            else:
                r_means = r_means[num_end_clip:]

        # flip means to match genomic positions
        if r_data.strand == '-':
            r_means = r_means[::-1]

        return r_means, read_start, read_end, read_id

    def get_read_comp_z_score(r_means, read_start, read_end):
        r_z_scores = np.abs(
            r_means - ctrl_means[read_start - reg_start + fm_offset:
                                 read_end - reg_start + fm_offset]) / ctrl_sds[
                                     read_start - reg_start + fm_offset:
                                     read_end - reg_start + fm_offset]

        return r_z_scores

    def get_pvals(r_z_scores):
        # mask out nan z-scores for efficient CDF computation
        r_poss = np.where(np.logical_not(np.isnan(r_z_scores)))[0]
        r_pvals = np.empty(r_z_scores.shape)
        r_pvals[:] = np.NAN
        valid_r_pvals = stats.norm.cdf(-r_z_scores[r_poss]) * 2.0
        r_pvals[r_poss] = valid_r_pvals

        return r_pvals, r_poss

    r_means, read_start, read_end, read_id = comp_clip_and_flip()
    r_z_scores = get_read_comp_z_score(r_means, read_start, read_end)

    if np.sum(np.logical_not(np.isnan(r_z_scores))) == 0:
        raise th.TomboError('No valid z-scores in read.')
    r_pvals, r_poss = get_pvals(r_z_scores)
    if fm_offset > 0:
        r_pvals = calc_window_fishers_method(r_pvals, fm_offset)
        r_poss = np.where(np.logical_not(np.isnan(r_pvals)))[0]
        r_pvals = r_pvals[r_poss]
    else:
        r_pvals = r_pvals[r_poss]

    r_poss += read_start

    return {SAMP_COMP_TXT:r_pvals}, {SAMP_COMP_TXT:r_poss}, read_id

def compute_de_novo_read_stats(
        r_data, std_ref, fm_offset=FM_OFFSET_DEFAULT, reg_data=None):
    """Compute signficance statistics using de novo comparison to a canonical
    model method for a single read within a specified genomic region.

    Args:

        r_data (:class:`tombo.tombo_helper.readData`): read data
        std_ref (:class:`tombo.tombo_stats.TomboModel`): canonical expected
           signal level model
        fm_offset (int): Fisher's Method offset for computing locally combined
            p-values (optional; default: 1)
        reg_data (:class:`tombo.tombo_helper.intervalData`): region to test
            (default: test whole read)

    Returns:
        Read testing results, positions tested and the read_id

        1) r_pvals (`np.array::np.float64`): p-values for testing over specified region
        2) r_poss (`np.array::np.int64`): genomic positions for returned p-values
        3) read_id (str): read identifier
    """
    reg_start = reg_data.start if reg_data is not None else r_data.start
    reg_size = (reg_data.end - reg_data.start if reg_data is not None
                else r_data.end - r_data.start)
    dnstrm_bases = std_ref.kmer_width - std_ref.central_pos - 1
    gnm_begin_lag = (std_ref.central_pos if r_data.strand == '+' else
                     dnstrm_bases)
    gnm_end_lag = (dnstrm_bases if r_data.strand == '+' else
                   std_ref.central_pos)

    def de_novo_clip_and_flip():
        with h5py.File(r_data.fn, 'r') as fast5_data:
            r_means, r_seq = th.get_multiple_slots_read_centric(
                fast5_data, ['norm_mean', 'base'], r_data.corr_group)
            read_id = th.get_raw_read_slot(fast5_data).attrs.get('read_id')

        if r_means is None or r_seq is None:
            raise th.TomboError(
                'Read does not contain valid re-squiggled data.')
        r_seq = b''.join(r_seq).decode()

        read_start, read_end = r_data.start, r_data.end
        # clip read if it extends outside the current genomic region, so
        # stats are only computed within this region
        if read_start + gnm_begin_lag + fm_offset < reg_start:
            num_start_clip = reg_start - (
                read_start + gnm_begin_lag + fm_offset)
            read_start = reg_start - gnm_begin_lag - fm_offset
            if r_data.strand == '+':
                r_means = r_means[num_start_clip:]
                r_seq = r_seq[num_start_clip:]
            else:
                r_means = r_means[:-num_start_clip]
                r_seq = r_seq[:-num_start_clip]
        if read_end - gnm_end_lag - fm_offset > reg_start + reg_size:
            num_end_clip = (read_end - gnm_end_lag - fm_offset) - (
                reg_start + reg_size)
            read_end = reg_start + reg_size + gnm_end_lag + fm_offset
            if r_data.strand == '+':
                r_means = r_means[:-num_end_clip]
                r_seq = r_seq[:-num_end_clip]
            else:
                r_means = r_means[num_end_clip:]
                r_seq = r_seq[num_end_clip:]

        # if this read does not cover enough of this region for stat
        # computation raise an error to be handled below
        if len(r_seq) < std_ref.kmer_width:
            raise th.TomboError(
                'Read does not contain information in this region.')

        r_ref_means, r_ref_sds = std_ref.get_exp_levels_from_seq(
            r_seq, r_data.strand == '-')

        if r_data.strand == '-':
            # reverse means to match genomic order
            r_means = r_means[::-1]
        # clip means that don't have testing model data
        # note that this is still extended by fm_offset
        r_means = r_means[gnm_begin_lag:-gnm_end_lag]
        read_start += gnm_begin_lag
        read_end -= gnm_end_lag

        return (r_means, r_ref_means, r_ref_sds, read_start,
                read_end, read_id)


    (r_means, r_ref_means, r_ref_sds,
     read_start, read_end, read_id) = de_novo_clip_and_flip()

    z_scores = np.abs(r_means - r_ref_means) / r_ref_sds
    r_pvals = stats.norm.cdf(-z_scores) * 2.0
    if fm_offset > 0:
        r_pvals = calc_window_fishers_method(r_pvals, fm_offset)

    # ignore errors in max over NAN values if fisher's method was used
    with np.errstate(invalid='ignore'):
        r_pvals = np.maximum(r_pvals, SMALLEST_PVAL)

    r_poss = np.arange(read_start, read_end)

    return {DE_NOVO_TXT:r_pvals}, {DE_NOVO_TXT:r_poss}, read_id

def calc_llh_ratio(
        reg_means, reg_ref_means, reg_ref_vars, reg_alt_means, reg_alt_vars):
    """Compute log likelihood ratio. This is about 10X slower than the cython
    version in tombo._c_helper, but has been kept for debugging purposes.
    """
    # compute log likelihood ratio
    # positive value means standard base fits data better
    # negative value means alternative base fits data better
    return (np.sum(np.square(reg_means - reg_alt_means) / reg_alt_vars) +
            np.sum(np.log(reg_alt_vars))) - (
                np.sum(np.square(reg_means - reg_ref_means) / reg_ref_vars) +
                np.sum(np.log(reg_ref_vars)))

def trim_seq_and_means(
        seq, means, r_start, reg_start, reg_end, strand,
        kmer_width, central_pos, max_motif_bb, max_motif_ab):
    """Return trimmed arrays to the genomic region specified.

    Args:
       seq: read genomic sequence (read-centric for rev strand)
       means: read base level means
       r_start: genomic start position
       reg_start: genomic region start
       reg_end: genomic region end
       strand: read mapped strand
       kmer_width: standard and alt k-mer width
       central_pos: central position within standard k-mer
       max_motif_bb: maximum (over alt models) bases required upstream for
           motif search
       max_motif_ab: maximum (over alt models) bases required downstream for
           motif search

    Arrays are:
        1) read centric k-mers (for expected level lookup)
        2) read-centric k-mer model-able means
        3) genome-centric alt test-able start
        4) motif search seq
    """
    r_end = r_start + means.shape[0]
    # save un-clipped motif search seq
    motif_search_seq = seq

    # clip read if it extends outside the current genomic region, so
    # stats are only computed within this region
    num_start_clip, num_end_clip = 0, 0
    if r_start + kmer_width - 1 < reg_start:
        if strand == '+':
            num_start_clip = reg_start - (r_start + kmer_width - 1)
        else:
            num_end_clip = reg_start - (r_start + kmer_width - 1)
        r_start = reg_start - (kmer_width - 1)
    if r_end - kmer_width + 1 > reg_end:
        if strand == '+':
            num_end_clip = r_end - kmer_width + 1 - reg_end
        else:
            num_start_clip = r_end - kmer_width + 1 - reg_end

    # clip sequence to that required for expected level lookups
    seq = seq[num_start_clip:]
    if num_end_clip > 0:
        seq = seq[:-num_end_clip]

    # clip extra bits from means to test-able means
    means = means[num_start_clip + central_pos:]
    means = means[:-(num_end_clip + kmer_width - central_pos - 1)]

    # if this read does not cover enough of this region for stat
    # alternate stat computation raise an error to be handled below
    if means.shape[0] < kmer_width:
        raise th.TomboError('Read sequence too short in this region.')

    kmers = th.get_seq_kmers(seq, kmer_width)
    if len(kmers) != means.shape[0]:
        raise th.TomboError('Mismatching k-mer and mean levels.')

    # shift r_start to first alt test-able position
    r_start += kmer_width - 1

    # clip and/or pad sequence for motif searching
    if num_start_clip + kmer_width - 1 - max_motif_bb >= 0:
        motif_search_seq = motif_search_seq[
            num_start_clip + kmer_width - 1 - max_motif_bb:]
    else:
        motif_search_seq = 'N' * -(
            num_start_clip + kmer_width - 1 - max_motif_bb) + motif_search_seq
    if num_end_clip + kmer_width - 1 - max_motif_ab >= 0:
        motif_search_seq = motif_search_seq[
            :-(num_end_clip + kmer_width - 1 - max_motif_ab)]
    else:
        motif_search_seq = motif_search_seq + 'N' * -(
            num_end_clip + kmer_width - 1 - max_motif_ab)
    if (len(motif_search_seq) - max_motif_bb - max_motif_bb !=
        means.shape[0] - ((kmer_width - 1) * 2)):
        th.TomboError('Motif search sequence not correct length.')

    return kmers, means, r_start, motif_search_seq

def compute_alt_model_read_stats(
        r_data, std_ref, alt_refs, use_standard_llhr=False, reg_data=None):
    """Compute signficance statistics using comparison of read signal to
    canonical and alternative models method for a single read within a
    specified genomic region.

    Args:
        r_data (:class:`tombo.tombo_helper.readData`): read data
        std_ref (:class:`tombo.tombo_stats.TomboModel`): canonical expected
            signal level model
        alt_refs (list of :class:`tombo.tombo_stats.AltModel`): alternative
            expected signal level models
        use_standard_llhr (bool): compute standard likelihood ratio; for details
            see https://nanoporetech.github.io/tombo/modified_base_detection.html#alternative-model-method
            (optional; default: False)
        reg_data (:class:`tombo.tombo_helper.intervalData`): region to test
            (default: test whole read)

    Returns:
        Read testing results, positions tested and the read_id

        1) r_llhrs (`np.array::np.float64`): log-likelihood ratios (or psuedo-llhrs) for testing over specified region
        2) r_poss (`np.array::np.int64`): genomic positions for returned p-values
        3) read_id (str): read identifier
    """
    reg_start = reg_data.start if reg_data is not None else r_data.start
    reg_end = reg_data.end if reg_data is not None else r_data.end

    # number of sequence positions to keep in order to look up both motif hits
    # and expected levels
    max_motif_bb = max([
        alt_ref.motif.mod_pos - 1 for _, alt_ref in alt_refs])
    max_motif_ab = max([
        alt_ref.motif.motif_len - alt_ref.motif.mod_pos
        for _, alt_ref in alt_refs])

    def alt_clip_and_flip():
        """Get read information in order to test requested region of this read

        base means, seq, k-mers, genomic start and read id
        """
        with h5py.File(r_data.fn, 'r') as fast5_data:
            r_means, r_seq = th.get_multiple_slots_read_centric(
                fast5_data, ['norm_mean', 'base'], r_data.corr_group)
            read_id = th.get_raw_read_slot(fast5_data).attrs.get('read_id')

        if r_means is None or r_seq is None:
            raise th.TomboError(
                'Read does not contain valid re-squiggled data.')
        r_seq = b''.join(r_seq).decode()

        r_kmers, r_means, r_start, motif_search_seq = trim_seq_and_means(
            r_seq, r_means, r_data.start, reg_start, reg_end, r_data.strand,
            std_ref.kmer_width, std_ref.central_pos, max_motif_bb, max_motif_ab)

        return r_kmers, r_means, motif_search_seq, r_start, read_id


    r_kmers, r_means, motif_search_seq, r_start, read_id = alt_clip_and_flip()
    testable_len = r_means.shape[0] - std_ref.kmer_width + 1
    r_ref_means, r_ref_sds = std_ref.get_exp_levels_from_kmers(r_kmers)
    r_ref_vars = np.square(r_ref_sds)

    all_poss = {}
    all_llhrs = {}
    for alt_name, alt_ref in alt_refs:
        alt_base_poss = []
        log_lh_ratios = []
        # note search space is clipped since all k-mers covering the position
        # of interest must be valid
        alt_i_motif_search_seq = motif_search_seq[
            max_motif_bb - (alt_ref.motif.mod_pos - 1):]
        if max_motif_ab - (alt_ref.motif.motif_len - alt_ref.motif.mod_pos) > 0:
            alt_i_motif_search_seq = alt_i_motif_search_seq[
                :-(max_motif_ab - (
                    alt_ref.motif.motif_len - alt_ref.motif.mod_pos))]
        for alt_m in alt_ref.motif.motif_pat.finditer(alt_i_motif_search_seq):
            alt_pos = alt_m.start()
            if r_data.strand == '+':
                alt_base_poss.append(r_start + alt_pos)
            else:
                alt_base_poss.append(r_start + testable_len - alt_pos - 1)
            r_pos_alt_means, r_pos_alt_sds = alt_ref.get_exp_levels_from_kmers(
                r_kmers[alt_pos:alt_pos + alt_ref.kmer_width])
            pos_args = [r_means[alt_pos:alt_pos + std_ref.kmer_width],
                        r_ref_means[alt_pos:alt_pos + std_ref.kmer_width],
                        r_pos_alt_means]
            if CONST_SD_MODEL:
                const_var = r_ref_vars[alt_pos]
                if use_standard_llhr:
                    pos_lh_ratio = c_calc_llh_ratio_const_var(
                        *(pos_args + const_var))
                else:
                    pos_lh_ratio = c_calc_scaled_llh_ratio_const_var(
                        *(pos_args + [const_var, OCLLHR_SCALE,
                                      OCLLHR_HEIGHT, OCLLHR_POWER]))
            else:
                if use_standard_llhr:
                    pos_lh_ratio = c_calc_llh_ratio(
                        *(pos_args + [
                            r_ref_vars[alt_pos:alt_pos + std_ref.kmer_width],
                            np.square(r_pos_alt_sds)]))
                else:
                    raise th.TomboError(
                        'Variable SD scaled likelihood ratio not implemented.')
            log_lh_ratios.append(pos_lh_ratio)

        all_llhrs[alt_name] = np.array(log_lh_ratios)
        all_poss[alt_name] = np.array(alt_base_poss)

    return all_llhrs, all_poss, read_id

def apply_per_read_thresh(
        reg_base_stats, single_read_thresh, lower_thresh, stat_type,
        stat_locs, ctrl_cov=None):
    reg_cov = np.array([base_stats.shape[0] for base_stats in reg_base_stats])

    if lower_thresh is not None:
        # filter base statistics that fall between the upper and lower
        # stat threshold for the log likelihood statistic
        reg_base_stats = [
            base_stats[np.logical_or(base_stats <= lower_thresh,
                                     base_stats >= single_read_thresh)]
            for base_stats in reg_base_stats]
        valid_cov = np.array([base_stats.shape[0]
                              for base_stats in reg_base_stats])
    elif stat_type == ALT_MODEL_TXT:
        # filter base statistics that fall between the upper and lower
        # stat threshold for the log likelihood statistic
        reg_base_stats = [base_stats[np.abs(base_stats) >= single_read_thresh]
                          for base_stats in reg_base_stats]
        valid_cov = np.array([base_stats.shape[0]
                              for base_stats in reg_base_stats])
    else:
        valid_cov = reg_cov

    if stat_type == SAMP_COMP_TXT:
        ctrl_cov = [ctrl_cov[pos] if pos in ctrl_cov else 0
                    for pos in stat_locs]
    else:
        # convert to list since python2 repeat objects can't be pickled
        ctrl_cov = list(repeat(0, stat_locs.shape[0]))

    reg_frac_std_base = np.array([
        np.greater_equal(
            base_stats, single_read_thresh).sum() / base_stats.shape[0]
        if base_stats.shape[0] > 0 else np.NAN
        for base_stats in reg_base_stats])

    return reg_frac_std_base, reg_cov, ctrl_cov, valid_cov

def collate_reg_stats(
        stats, stat_locs, read_ids, per_read_q, reg_data, single_read_thresh,
        lower_thresh, stat_type, stat_name, ctrl_cov):
    stats = np.concatenate(stats)
    stat_locs = np.concatenate(stat_locs)
    # remove nans possibly introduced by fisher's method calculcations
    valid_poss = ~np.isnan(stats)
    stat_locs = stat_locs[valid_poss]
    stats = stats[valid_poss]
    assert stat_locs.shape[0] == stats.shape[0], '\t'.join(map(str, (
        stat_locs.shape[0], stats.shape[0])))

    if per_read_q is not None:
        valid_read_ids = [
            rep_r_id for rep_r_id, is_valid in zip(
                [rep_r_id for r_id, r_len in read_ids
                 for rep_r_id in repeat(r_id, r_len)], valid_poss) if is_valid]
        read_id_lookup = dict((
            (read_id, read_id_val)
            for read_id_val, read_id in enumerate(set(valid_read_ids))))
        conv_read_ids = np.array([
            read_id_lookup[r_id] for r_id in valid_read_ids])
        assert conv_read_ids.shape[0] == stat_locs.shape[0]
        per_read_block = np.array(
            list(zip(stat_locs, stats, conv_read_ids)),
            dtype=[(str('pos'), 'u4'), (str('stat'), 'f8'),
                   (str('read_id'), 'u4')])
        per_read_q.put((
            stat_name, (per_read_block, read_id_lookup, reg_data.chrm,
                        reg_data.strand, reg_data.start)))

    # get order of all bases from position array
    as_stat_locs = np.argsort(stat_locs)
    # sort all positions from all reads
    stat_locs = stat_locs[as_stat_locs]
    # get unique tested genomic positions across all reads
    us_stat_locs = np.unique(stat_locs)

    if stat_locs.shape[0] == 0:
        raise th.TomboError('No valid positions in this region.')

    # then sort the stats array by genomic position and
    # split into stats by genomic base position
    reg_base_stats = np.split(
        stats[as_stat_locs],
        np.where(np.concatenate([[0,], np.diff(stat_locs)]) > 0)[0])

    reg_frac_std_base, reg_cov, ctrl_cov, valid_cov = apply_per_read_thresh(
        reg_base_stats, single_read_thresh, lower_thresh,
        stat_type, stat_locs, ctrl_cov)

    return th.regionStats(
        reg_frac_std_base, us_stat_locs, reg_data.chrm, reg_data.strand,
        reg_data.start, reg_cov, ctrl_cov, valid_cov)

def compute_reg_stats(
        reg_data, fm_offset, min_test_reads,
        single_read_thresh, lower_thresh, ctrl_reg_data, std_ref,
        alt_refs, use_standard_llhr, per_read_q, stat_type, prior_weights):
    if stat_type == SAMP_COMP_TXT:
        ctrl_means, ctrl_sds, ctrl_cov = get_reads_ref(
            ctrl_reg_data, min_test_reads, fm_offset, std_ref, prior_weights)
    else:
        # TODO get region sequence and expected levels/sds here
        # instead of for each read
        # after that add per-read stat computation to API
        ctrl_cov = None

    # store multiple alt references lists (default to stat_type for de novo or
    # sample comp testing)
    stat_names = [stat_type,] if stat_type != ALT_MODEL_TXT else list(
        zip(*alt_refs))[0]
    reg_read_stats, stat_locs, reg_ids = [
        dict((stat_name, []) for stat_name in stat_names) for _ in range(3)]
    for r_data in reg_data.reads:
        try:
            if stat_type == SAMP_COMP_TXT:
                r_stats, r_poss, read_id = compute_sample_compare_read_stats(
                    r_data, ctrl_means, ctrl_sds, fm_offset, reg_data)
            elif stat_type == DE_NOVO_TXT:
                r_stats, r_poss, read_id = compute_de_novo_read_stats(
                    r_data, std_ref, fm_offset, reg_data)
            else:
                r_stats, r_poss, read_id = compute_alt_model_read_stats(
                    r_data, std_ref, alt_refs, use_standard_llhr, reg_data)
        except th.TomboError:
            continue
        if r_stats is None: continue
        for stat_name, stat_r_stats in r_stats.items():
            reg_read_stats[stat_name].append(stat_r_stats)
            reg_ids[stat_name].append((read_id, stat_r_stats.shape[0]))
            stat_locs[stat_name].append(r_poss[stat_name])

    if sum(len(stat_reg_read_stats)
           for stat_reg_read_stats in reg_read_stats.values()) == 0:
        raise th.TomboError('Reads contains no statistics in this region.')

    reg_stats = [
        (stat_name, collate_reg_stats(
            stat_name_stats, stat_locs[stat_name], reg_ids[stat_name],
            per_read_q, reg_data, single_read_thresh, lower_thresh, stat_type,
            stat_name,  ctrl_cov))
        for stat_name, stat_name_stats in reg_read_stats.items()]

    return reg_stats


###################################
########## Level Testing ##########
###################################

def compute_ks_tests(samp_base_levels, ctrl_base_levels, return_stat):
    def compute_pos_ks_test(pos_samp_levels, pos_ctrl_levels):
        """Compute effect size statistic or p-value of two-sample
        Kolmogorov-Smirnov test

        Using definition from
        https://github.com/scipy/scipy/blob/v0.14.0/scipy/stats/stats.py#L3886
        """
        samp_n, ctrl_n = pos_samp_levels.shape[0], pos_ctrl_levels.shape[0]
        pos_all_levels = np.concatenate([pos_samp_levels, pos_ctrl_levels])
        samp_cdf = np.searchsorted(
            pos_samp_levels, pos_all_levels, side='right') / samp_n
        ctrl_cdf = np.searchsorted(
            pos_ctrl_levels, pos_all_levels, side='right') / ctrl_n
        d = np.max(np.absolute(samp_cdf - ctrl_cdf))
        if return_stat:
            # subtract 1 so most significant are smallest values
            return 1 - d
        en = np.sqrt(samp_n * ctrl_n / float(samp_n + ctrl_n))
        return stats.distributions.kstwobign.sf((en + 0.12 + 0.11 / en) * d)


    samp_valid_indices = np.logical_not(np.isnan(samp_base_levels))
    ctrl_valid_indices = np.logical_not(np.isnan(ctrl_base_levels))
    return np.array([compute_pos_ks_test(
        np.sort(pos_samp_levels[samp_valid_indices[i]]),
        np.sort(pos_ctrl_levels[ctrl_valid_indices[i]]))
                     for i, (pos_samp_levels, pos_ctrl_levels) in enumerate(zip(
                             samp_base_levels, ctrl_base_levels))])

def compute_u_tests(samp_base_levels, ctrl_base_levels, return_stat):
    def compute_pos_u_test(pos_samp_levels, pos_ctrl_levels):
        """Compute effect size statistic or p-value of two-sample u-test
        """
        samp_n, ctrl_n = pos_samp_levels.shape[0], pos_ctrl_levels.shape[0]
        tot_comps = samp_n * ctrl_n
        pos_all_levels = np.concatenate([pos_samp_levels, pos_ctrl_levels])

        ranks = np.empty(samp_n + ctrl_n, int)
        ranks[pos_all_levels.argsort()] = np.arange(1, samp_n + ctrl_n + 1)
        samp_ranks_sum = ranks[:samp_n].sum()
        #ctrl_ranks_sum = ranks[samp_n:].sum()

        u1 = samp_ranks_sum - (samp_n * (samp_n + 1)) / 2
        #u2 = ctrl_ranks_sum - (ctrl_n * (ctrl_n + 1)) / 2
        u2 = tot_comps - u1
        u = min(u1, u2)

        mu = tot_comps / 2
        if return_stat:
            # this will be negative (flip sign for stat transform)
            return (u - mu) / mu

        rhou = np.sqrt(tot_comps * (tot_comps + 1) / 12)
        z = (u - mu) / rhou
        return stats.norm.cdf(z) * 2.0


    samp_valid_indices = np.logical_not(np.isnan(samp_base_levels))
    ctrl_valid_indices = np.logical_not(np.isnan(ctrl_base_levels))
    return np.array([compute_pos_u_test(
        np.sort(pos_samp_levels[samp_valid_indices[i]]),
        np.sort(pos_ctrl_levels[ctrl_valid_indices[i]]))
                     for i, (pos_samp_levels, pos_ctrl_levels) in enumerate(zip(
                             samp_base_levels, ctrl_base_levels))])

def compute_t_tests(samp_base_levels, ctrl_base_levels, return_stat):
    def compute_pos_t_test(pos_samp_levels, pos_ctrl_levels):
        """Compute effect size statistic or p-value of two-sample t-test
        """
        samp_n, ctrl_n = pos_samp_levels.shape[0], pos_ctrl_levels.shape[0]
        tot_comps = samp_n * ctrl_n

        samp_mean, samp_sd = c_mean_std(pos_samp_levels)
        ctrl_mean, ctrl_sd = c_mean_std(pos_ctrl_levels)

        if return_stat:
            # this will be negative (flip sign for stat transform)
            return -np.abs(samp_mean - ctrl_mean) / np.sqrt(
                ((samp_sd ** 2) + (ctrl_sd ** 2)) / 2)

        sp = np.sqrt((((samp_n - 1) * (samp_sd ** 2)) +
                      (ctrl_n - 1) * (ctrl_sd ** 2)) /
                     (samp_n + ctrl_n - 2))
        t = -np.abs(samp_mean - ctrl_mean) / (
            sp * np.sqrt((1 / samp_n) + (1 / ctrl_n)))

        # t dist with samp_n + ctrl_n - 2 d.o.f.
        return stats.t.cdf(t, samp_n + ctrl_n - 2) * 2.0


    samp_valid_indices = np.logical_not(np.isnan(samp_base_levels))
    ctrl_valid_indices = np.logical_not(np.isnan(ctrl_base_levels))
    return np.array([compute_pos_t_test(
        np.sort(pos_samp_levels[samp_valid_indices[i]]),
        np.sort(pos_ctrl_levels[ctrl_valid_indices[i]]))
                     for i, (pos_samp_levels, pos_ctrl_levels) in enumerate(zip(
                             samp_base_levels, ctrl_base_levels))])

def compute_group_reg_stats(
        reg_data, ctrl_reg_data, fm_offset, min_test_reads, stat_type):
    samp_base_levels = reg_data.copy().update(
        start=reg_data.start - fm_offset,
        end=reg_data.end + fm_offset).get_base_levels()
    ctrl_base_levels = ctrl_reg_data.copy().update(
        start=ctrl_reg_data.start - fm_offset,
        end=ctrl_reg_data.end + fm_offset).get_base_levels()

    # get regions with coverage greater than min_test_reads
    samp_cov = np.logical_not(np.isnan(samp_base_levels)).sum(axis=1)
    ctrl_cov = np.logical_not(np.isnan(ctrl_base_levels)).sum(axis=1)
    cov_regs = np.where(np.diff(np.concatenate([[False,], np.logical_and(
        np.greater_equal(samp_cov, min_test_reads),
        np.greater_equal(ctrl_cov, min_test_reads)), [False,]])))[0]
    if len(cov_regs) == 0:
        return []

    reg_stats, reg_poss, reg_cov, reg_ctrl_cov = [], [], [], []
    for cov_start, cov_end in zip(cov_regs[:-1:2], cov_regs[1::2]):
        if cov_end - cov_start < (fm_offset * 2) + 1: continue
        if stat_type in (KS_TEST_TXT, KS_STAT_TEST_TXT):
            cov_reg_stats = compute_ks_tests(
                samp_base_levels[cov_start:cov_end],
                ctrl_base_levels[cov_start:cov_end],
                stat_type == KS_STAT_TEST_TXT)
        elif stat_type in (U_TEST_TXT, U_STAT_TEST_TXT):
            cov_reg_stats = compute_u_tests(
                samp_base_levels[cov_start:cov_end],
                ctrl_base_levels[cov_start:cov_end],
                stat_type == U_STAT_TEST_TXT)
        elif stat_type in (T_TEST_TXT, T_STAT_TEST_TXT):
            cov_reg_stats = compute_t_tests(
                samp_base_levels[cov_start:cov_end],
                ctrl_base_levels[cov_start:cov_end],
                stat_type == T_STAT_TEST_TXT)
        else:
            raise NotImplementedError('Unrecognized test type.')
        if fm_offset > 0:
            if stat_type in (KS_TEST_TXT, U_TEST_TXT, T_TEST_TXT):
                cov_reg_stats = calc_window_fishers_method(
                    cov_reg_stats, fm_offset)
            else:
                cov_reg_stats = calc_window_means(cov_reg_stats, fm_offset)
        reg_stats.append(cov_reg_stats)
        reg_poss.append(np.arange(reg_data.start - fm_offset + cov_start,
                                  reg_data.start - fm_offset + cov_end))
        reg_cov.append(samp_cov[cov_start:cov_end])
        reg_ctrl_cov.append(ctrl_cov[cov_start:cov_end])

    return [(stat_type, th.groupStats(
        np.concatenate(reg_stats), np.concatenate(reg_poss),
        reg_data.chrm, reg_data.strand, reg_data.start,
        np.concatenate(reg_cov), np.concatenate(reg_ctrl_cov))),]


##############################################
########## Testing Multi-processing ##########
##############################################

def _test_signif_worker(
        region_q, stats_q, progress_q, per_read_q, reads_index, fm_offset,
        min_test_reads, single_read_thresh, lower_thresh, ctrl_reads_index,
        std_ref, alt_refs, use_standard_llhr, stat_type, prior_weights):
    ctrl_reg_data = None
    while True:
        try:
            reg_data = region_q.get(block=False)
        except queue.Empty:
            # sometimes throws false empty error with get(block=False)
            sleep(0.01)
            if not region_q.empty(): continue
            break

        if ctrl_reads_index is not None:
            ctrl_reg_data = reg_data.copy().add_reads(ctrl_reads_index)
        reg_data.add_reads(reads_index)
        if len(reg_data.reads) == 0:
            progress_q.put(1)
            continue

        try:
            if stat_type in (ALT_MODEL_TXT, DE_NOVO_TXT, SAMP_COMP_TXT):
                stat_type_reg_stats = compute_reg_stats(
                    reg_data, fm_offset, min_test_reads, single_read_thresh,
                    lower_thresh, ctrl_reg_data, std_ref, alt_refs,
                    use_standard_llhr, per_read_q, stat_type, prior_weights)
            else:
                stat_type_reg_stats = compute_group_reg_stats(
                    reg_data, ctrl_reg_data, fm_offset, min_test_reads,
                    stat_type)
        except th.TomboError:
            progress_q.put(1)
            continue
        for stat_type_i_reg_stats in stat_type_reg_stats:
            stats_q.put(stat_type_i_reg_stats)
        progress_q.put(1)

    return

if _PROFILE_SIGNIF:
    _test_signif_wrapper = _test_signif_worker
    def _test_signif_worker(*args):
        import cProfile
        cProfile.runctx('_test_signif_wrapper(*args)', globals(), locals(),
                        filename='test_signif.prof')
        return

def _get_stats_queue(
        stats_q, stats_conn, min_test_reads, stats_file_bn,
        stat_type, alt_names, reg_size, cov_damp_counts, num_most_signif):
    # multiple files with alt_names
    all_stats = {}
    if stat_type == ALT_MODEL_TXT:
        for alt_name in alt_names:
            all_stats[alt_name] = ModelStats(
                stats_file_bn + '.' + alt_name + '.tombo.stats',
                stat_type=stat_type, region_size=reg_size,
                cov_damp_counts=cov_damp_counts, cov_thresh=min_test_reads,
                num_most_signif=num_most_signif)
    elif stat_type in (DE_NOVO_TXT, SAMP_COMP_TXT):
        all_stats[stat_type] = ModelStats(
            stats_file_bn + '.tombo.stats',
            stat_type=stat_type, region_size=reg_size,
            cov_damp_counts=cov_damp_counts, cov_thresh=min_test_reads,
            num_most_signif=num_most_signif)
    else:
        all_stats[stat_type] = LevelStats(
            stats_file_bn + '.tombo.stats',
            stat_type=stat_type, region_size=reg_size,
            cov_thresh=min_test_reads, num_most_signif=num_most_signif)

    while True:
        try:
            stat_name, reg_stats = stats_q.get(block=False)
            all_stats[stat_name]._write_stat_block(reg_stats)
        except queue.Empty:
            # wait for main process to send indicator that all regions
            # have been processed
            if stats_conn.poll():
                sleep(0.1)
                break
            sleep(0.1)
            continue

    # Clear leftover values from queues
    while not stats_q.empty():
        stat_name, reg_stats = stats_q.get(block=False)
        all_stats[stat_name]._write_stat_block(reg_stats)

    for stat_name_all_stats in all_stats.values():
        stat_name_all_stats.close()
    stats_conn.send(True)

    return

if _PROFILE_SIGNIF_STATS_OUT:
    _get_stats_queue_wrapper = _get_stats_queue
    def _get_stats_queue(*args):
        import cProfile
        cProfile.runctx('_get_stats_queue_wrapper(*args)', globals(), locals(),
                        filename='test_signif_stats_out.prof')
        return

def _get_per_read_queue(
        per_read_q, per_read_conn, per_read_bn, stat_type, alt_names,
        region_size):
    per_read_stats = {}
    if stat_type == ALT_MODEL_TXT:
        for alt_name in alt_names:
            per_read_stats[alt_name] = PerReadStats(
                per_read_bn + '.' + alt_name + '.tombo.per_read_stats',
                stat_type, region_size)
    else:
        per_read_stats[stat_type] = PerReadStats(
                per_read_bn + '.tombo.per_read_stats',
                stat_type, region_size)

    while True:
        try:
            stat_name, per_read_block = per_read_q.get(block=False)
            per_read_stats[stat_name]._write_per_read_block(*per_read_block)
            del per_read_block
        except queue.Empty:
            if per_read_conn.poll():
                sleep(0.1)
                break
            sleep(0.1)
            continue

    # Clear leftover values from queues
    while not per_read_q.empty():
        stat_name, per_read_block = per_read_q.get(block=False)
        per_read_stats[stat_name]._write_per_read_block(*per_read_block)
        del per_read_block
    for stat_name_per_read_stats in per_read_stats.values():
        stat_name_per_read_stats.close()

    # indicate that the process has closed
    per_read_conn.send(True)

    return

if _PROFILE_SIGNIF_PER_READ:
    _get_per_read_queue_wrapper = _get_per_read_queue
    def _get_per_read_queue(*args):
        import cProfile
        cProfile.runctx(
            '_get_per_read_queue_wrapper(*args)', globals(), locals(),
            filename='test_signif_per_read.prof')
        return

def _get_progress_queue(progress_q, prog_conn, num_regions):
    th.status_message(
        'Performing modified base detection across genomic regions.')
    bar = tqdm(total=num_regions, smoothing=0)

    tot_num_rec_proc = 0
    while True:
        try:
            iter_val = progress_q.get(block=False)
            tot_num_rec_proc += iter_val
            bar.update(iter_val)
        except queue.Empty:
            if prog_conn.poll():
                break
            sleep(0.1)
            continue

    bar.close()
    prog_conn.send(tot_num_rec_proc)

    return

def test_significance(
        reads_index, stat_type, stats_file_bn, region_size, num_processes,
        min_test_reads, num_most_signif, per_read_bn=None,
        single_read_thresh=None, lower_thresh=None, cov_damp_counts=None,
        fm_offset=None, ctrl_reads_index=None, std_ref=None, alt_refs=None,
        use_standard_llhr=False, prior_weights=None):
    """Test for significant shifted signal in mutliprocessed batches
    """
    region_q = Queue()
    stats_q = Queue(STAT_BLOCKS_QUEUE_LIMIT)
    progress_q = Queue()
    per_read_q = Queue(STAT_BLOCKS_QUEUE_LIMIT) \
                 if per_read_bn else None
    # split chromosomes into separate regions to process independently
    num_regions = 0
    # TODO add min coverage when adding ctrl coverage (instead of adding
    # as is done currently)
    for chrm, strand, reg_start in reads_index.iter_cov_regs(
            1, region_size, ctrl_reads_index):
        region_q.put(th.intervalData(
            chrm=chrm, start=reg_start, end=reg_start + region_size,
            strand=strand))
        num_regions += 1
    # wait for queue items to register in queue and avoid queue appearing empty
    sleep(0.1)

    test_args = (
        region_q, stats_q, progress_q, per_read_q, reads_index, fm_offset,
        min_test_reads, single_read_thresh, lower_thresh, ctrl_reads_index,
        std_ref, alt_refs, use_standard_llhr, stat_type, prior_weights)
    test_ps = []
    for p_id in range(num_processes):
        p = Process(target=_test_signif_worker, args=test_args)
        p.start()
        test_ps.append(p)

    # start queue getter processes
    if VERBOSE:
        main_prog_conn, prog_conn = Pipe()
        prog_p = Process(target=_get_progress_queue,
                            args=(progress_q, prog_conn, num_regions))
        prog_p.daemon = True
        prog_p.start()

    # main region stats queue getter
    alt_names = None if alt_refs is None else list(zip(*alt_refs))[0]
    main_stats_conn, stats_conn = Pipe()
    stats_p = Process(target=_get_stats_queue, args=(
        stats_q, stats_conn, min_test_reads, stats_file_bn,
        stat_type, alt_names, region_size, cov_damp_counts,
        num_most_signif))
    stats_p.daemon = True
    stats_p.start()

    # per-read stats queue getter
    if per_read_bn is not None:
        main_per_read_conn, per_read_conn = Pipe()
        per_read_p = Process(
            target=_get_per_read_queue,
            args=(per_read_q, per_read_conn, per_read_bn, stat_type,
                  alt_names, region_size))
        per_read_p.daemon = True
        per_read_p.start()

    # wait for test processes to finish
    for test_p in test_ps:
        test_p.join()

    # in a very unlikely case the progress queue could die while the
    # main process remains active and thus we would have a deadlock here
    if VERBOSE and prog_p.is_alive():
        # send signal to getter queue to finish and return results
        main_prog_conn.send(True)
        # returns total number of processed reads if that is needed
        main_prog_conn.recv()

    if per_read_bn is not None:
        main_per_read_conn.send(True)
        main_per_read_conn.recv()

    main_stats_conn.send(True)
    main_stats_conn.recv()

    return


################################################
########## Aggregate Multi-processing ##########
################################################

def _write_stats(
        stats_q, stats_fn, stat_type, region_size, cov_damp_counts,
        min_test_reads, num_most_signif, num_blocks, num_processes):
    all_stats = ModelStats(
        stats_fn, stat_type=stat_type, region_size=region_size,
        cov_damp_counts=cov_damp_counts, cov_thresh=min_test_reads,
        num_most_signif=num_most_signif)
    if VERBOSE:
        bar = tqdm(total=num_blocks, smoothing=0)
    num_agg_ps_finished = 0
    while True:
        try:
            agg_stats = stats_q.get(block=False)
            if agg_stats is None:
                num_agg_ps_finished += 1
                if num_agg_ps_finished >= num_processes: break
                continue

            if VERBOSE: bar.update(1)
            ((reg_frac_std_base, reg_cov, ctrl_cov, valid_cov),
             chrm, strand, start, us_reg_poss) = agg_stats
            all_stats._write_stat_block(
                th.regionStats(reg_frac_std_base, us_reg_poss, chrm, strand,
                               start, reg_cov, ctrl_cov, valid_cov))
        except queue.Empty:
            sleep(0.1)

    if VERBOSE: bar.close()
    all_stats.close()
    if all_stats.is_empty:
        th.error_message_and_exit(
            'No genomic positions contain --minimum-test-reads.')

    return

def _agg_stats_worker(
        pr_stats_q, stats_q, stat_type, single_read_thresh, lower_thresh):
    while True:
        try:
            block_pr_stats = pr_stats_q.get(block=False)
            # None value indicates that per-reads blocks have been exhausted
            if block_pr_stats is None:
                stats_q.put(None)
                break
            chrm, strand, start, end, block_stats = block_pr_stats

            block_stats.sort(order=str('pos'))
            reg_poss = block_stats['pos']
            us_reg_poss = np.unique(reg_poss)

            reg_base_stats = np.split(
                block_stats['stat'], np.where(np.concatenate(
                    [[0,], np.diff(reg_poss)]) > 0)[0])

            reg_stats = apply_per_read_thresh(
                 reg_base_stats, single_read_thresh, lower_thresh,
                 stat_type, reg_poss)
            stats_q.put((reg_stats, chrm, strand, start, us_reg_poss))
        except queue.Empty:
            sleep(0.1)

    return

def _load_stats_batches(pr_stats_fn, pr_stats_q, num_processes):
    pr_stats = PerReadStats(pr_stats_fn)
    for pr_block in pr_stats:
        pr_stats_q.put(pr_block)
    for _ in range(num_processes):
        pr_stats_q.put(None)

    return

def aggregate_per_read_stats(
        pr_stats_fn, single_read_thresh, lower_thresh, stats_fn,
        cov_damp_counts, min_test_reads, num_most_signif, num_processes):
    if VERBOSE: th.status_message(
            'Loading and aggregating per-read statistics.')

    # pre-load per-read stats queue
    pr_stats = PerReadStats(pr_stats_fn)
    stat_type, num_blocks, region_size = (
        pr_stats.stat_type, pr_stats.num_blocks, pr_stats.region_size)
    pr_stats.close()

    pr_stats_q = Queue(STAT_BLOCKS_QUEUE_LIMIT)
    stats_q = Queue(STAT_BLOCKS_QUEUE_LIMIT)
    write_stats_q = Queue(STAT_BLOCKS_QUEUE_LIMIT)

    load_stats_p = Process(target=_load_stats_batches, args=(
        pr_stats_fn, pr_stats_q, num_processes))
    load_stats_p.daemon = True
    load_stats_p.start()

    agg_stats_ps = []
    for p_id in range(num_processes):
        agg_p = Process(target=_agg_stats_worker, args=(
            pr_stats_q, stats_q, stat_type, single_read_thresh, lower_thresh))
        agg_p.daemon = True
        agg_p.start()
        agg_stats_ps.append(agg_p)

    write_stats_p = Process(target=_write_stats, args=(
        stats_q, stats_fn, stat_type, region_size, cov_damp_counts,
        min_test_reads, num_most_signif, num_blocks, num_processes))
    write_stats_p.daemon = True
    write_stats_p.start()

    # wait for processes to complete
    load_stats_p.join()
    for agg_p in agg_stats_ps:
        agg_p.join()
    write_stats_p.join()

    return


##########################
##### Main Functions #####
##########################

def _est_ref_main(args):
    global VERBOSE
    VERBOSE = not args.quiet
    th.VERBOSE = VERBOSE

    if min(args.upstream_bases, args.downstream_bases) == 0:
        th.error_message_and_exit(
            'Context upstream and downstream must be greater ' +
            'than 0 for model estimation.')

    std_ref = estimate_kmer_model(
        args.fast5_basedirs, args.corrected_group, args.basecall_subgroups,
        args.minimum_test_reads, args.upstream_bases, args.downstream_bases,
        args.minimum_kmer_observations, args.kmer_specific_sd,
        args.coverage_threshold, args.estimate_mean,
        args.multiprocess_region_size, args.processes)
    std_ref.write_model(args.tombo_model_filename)

    return

def _est_alt_ref_main(args):
    global VERBOSE
    VERBOSE = not args.quiet
    th.VERBOSE = VERBOSE

    alt_ref = estimate_alt_model(
        args.fast5_basedirs, args.control_fast5_basedirs,
        args.corrected_group, args.basecall_subgroups,
        args.tombo_model_filename, args.seq_sample_type,
        args.alternate_model_base, args.alt_fraction_percentile,
        args.minimum_kmer_observations, args.save_density_basename,
        args.kernel_density_bandwidth, args.alternate_density_filename,
        args.control_density_filename, args.processes)
    # returns None when profiling method
    if alt_ref is None: return
    alt_ref.name = args.alternate_model_name
    alt_ref.write_model(args.alternate_model_filename)

    return

def _est_motif_alt_ref_main(args):
    global VERBOSE
    VERBOSE = not args.quiet
    th.VERBOSE = VERBOSE

    alt_ref = estimate_motif_alt_model(
        args.fast5_basedirs, args.corrected_group, args.basecall_subgroups,
        args.motif_description, args.upstream_bases, args.downstream_bases,
        args.valid_locations_filename, args.minimum_kmer_observations,
        args.minimum_test_reads, args.coverage_threshold,
        args.multiprocess_region_size, args.processes)
    alt_ref.name = args.alternate_model_name
    alt_ref.write_model(args.alternate_model_filename)

    return

def _estimate_scale_main(args):
    global VERBOSE
    VERBOSE = not args.quiet
    th.VERBOSE = VERBOSE

    if VERBOSE: th.status_message('Getting files list.')
    try:
        if not os.path.isdir(args.fast5s_basedir):
            th.error_message_and_exit(
                'Provided [fast5-basedir] is not a directory.')
        fast5s_basedir = (
            args.fast5s_basedir if args.fast5s_basedir.endswith('/') else
            args.fast5s_basedir + '/')
        fast5_fns = th.get_files_list(fast5s_basedir)
    except OSError:
        th.error_message_and_exit(
            'Reads base directory, a sub-directory or an old (hidden) ' +
            'index file does not appear to be accessible. Check ' +
            'directory permissions.')
    if len(fast5_fns) < 1:
        th.error_message_and_exit(
            'No files identified in the specified ' +
            'directory or within immediate subdirectories.')

    th.status_message('Global scaling estimate: ' +
                       unicode(estimate_global_scale(fast5_fns)))

    return

def _de_novo_main(
        args, lower_thresh, single_read_thresh, seq_samp_type, reads_index):
    if seq_samp_type is None:
        seq_samp_type = th.get_seq_sample_type(reads_index=reads_index)
    std_ref = TomboModel(
        ref_fn=args.tombo_model_filename, seq_samp_type=seq_samp_type,
        reads_index=reads_index)

    stat_type = DE_NOVO_TXT
    lower_thresh, single_read_thresh = (
        (lower_thresh, single_read_thresh) if single_read_thresh
        is not None else DE_NOVO_THRESH[seq_samp_type.name])
    if VERBOSE: th.status_message(
            'Performing de novo model testing against canonical model.')
    test_significance(
        reads_index, stat_type, args.statistics_file_basename,
        args.multiprocess_region_size, args.processes, args.minimum_test_reads,
        args.num_most_significant_stored, args.per_read_statistics_basename,
        single_read_thresh, lower_thresh, args.coverage_dampen_counts,
        fm_offset=args.fishers_method_context, std_ref=std_ref)

    return

def _alt_main(
        args, lower_thresh, single_read_thresh, seq_samp_type, reads_index):
    if seq_samp_type is None:
        seq_samp_type = th.get_seq_sample_type(reads_index=reads_index)
    std_ref = TomboModel(
        ref_fn=args.tombo_model_filename, seq_samp_type=seq_samp_type,
        reads_index=reads_index)

    stat_type = ALT_MODEL_TXT
    lower_thresh, single_read_thresh = (
        (lower_thresh, single_read_thresh) if single_read_thresh
        is not None else LLR_THRESH[seq_samp_type.name])
    if VERBOSE: th.status_message('Performing alternative model testing.')
    alt_refs = load_alt_refs(
        args.alternate_model_filenames, args.alternate_bases,
        reads_index, std_ref, seq_samp_type)
    if len(alt_refs) == 0:
        th.error_message_and_exit('No alternative models successfully loaded.')

    if VERBOSE: th.status_message(
            'Performing specific alternate base(s) testing.')
    test_significance(
        reads_index, stat_type, args.statistics_file_basename,
        args.multiprocess_region_size, args.processes,
        args.minimum_test_reads, args.num_most_significant_stored,
        args.per_read_statistics_basename, single_read_thresh, lower_thresh,
        args.coverage_dampen_counts,
        std_ref=std_ref, alt_refs=list(alt_refs.items()),
        use_standard_llhr=args.standard_log_likelihood_ratio)

    return

def _model_samp_comp_main(
        args, lower_thresh, single_read_thresh, seq_samp_type, reads_index):
    stat_type = SAMP_COMP_TXT
    if single_read_thresh is None:
        if seq_samp_type is None:
            seq_samp_type = th.get_seq_sample_type(reads_index=reads_index)
        lower_thresh, single_read_thresh = SAMP_COMP_THRESH[seq_samp_type.name]
    if VERBOSE: th.status_message(
            'Performing two-sample comparison significance testing.')
    ctrl_reads_index = th.TomboReads(
        args.control_fast5_basedirs, args.corrected_group,
        args.basecall_subgroups)

    # load expected levels ref for posterior computation
    std_ref = None if args.sample_only_estimates else TomboModel(
        ref_fn=args.tombo_model_filename, seq_samp_type=seq_samp_type,
        reads_index=reads_index)

    test_significance(
        reads_index, stat_type, args.statistics_file_basename,
        args.multiprocess_region_size, args.processes, args.minimum_test_reads,
        args.num_most_significant_stored, args.per_read_statistics_basename,
        single_read_thresh, lower_thresh, args.coverage_dampen_counts,
        fm_offset=args.fishers_method_context,
        ctrl_reads_index=ctrl_reads_index, std_ref=std_ref,
        prior_weights=args.model_prior_weights)

    return

def _level_samp_comp_main(args, reads_index):
    if args.statistic_type == 'ks':
        stat_type = KS_TEST_TXT if args.store_p_value else KS_STAT_TEST_TXT
    elif args.statistic_type == 'u':
        stat_type = U_TEST_TXT if args.store_p_value else U_STAT_TEST_TXT
    elif args.statistic_type == 't':
        stat_type = T_TEST_TXT if args.store_p_value else T_STAT_TEST_TXT
    else:
        raise th.TomboError('Invalid statistic type.')

    if VERBOSE: th.status_message(
            'Performing two-sample group comparison significance testing.')
    ctrl_reads_index = th.TomboReads(
        args.alternate_fast5_basedirs, args.corrected_group,
        args.basecall_subgroups)

    test_significance(
        reads_index, stat_type, args.statistics_file_basename,
        args.multiprocess_region_size, args.processes,
        args.minimum_test_reads, args.num_most_significant_stored,
        fm_offset=args.fishers_method_context,
        ctrl_reads_index=ctrl_reads_index)

    return

def _test_shifts_main(args):
    global VERBOSE
    VERBOSE = not args.quiet
    th.VERBOSE = VERBOSE

    # Check extra requirements for alternative model testing
    if args.action_command == 'alternative_model':
        if args.print_available_models:
            _print_alt_models()
            sys.exit()
        if args.fast5_basedirs is None or args.statistics_file_basename is None:
            th.error_message_and_exit(
                'Must provide both a set of FAST5 read files ' +
                '(--fast5-basedirs) and an output file basename ' +
                '(--statistics-file-basename).')
        if (args.alternate_model_filenames is None and
            args.alternate_bases is None):
            th.error_message_and_exit(
                'Must provide an alterntive model against which to test.\n\t' +
                'Run with --print-available-models option to see possible ' +
                'values for the --alternate-bases option.')

    if 'single_read_threshold' in args:
        if args.single_read_threshold is None:
            lower_thresh = None
            single_read_thresh = None
        elif len(args.single_read_threshold) == 1:
            single_read_thresh = args.single_read_threshold[0]
            lower_thresh = None
        else:
            if len(args.single_read_threshold) > 2:
                th.warning_message(
                    'Only 1 or 2 values may be passed as single-read ' +
                    'thresholds. Only using the first 2 options provided.')
            lower_thresh = args.single_read_threshold[0]
            single_read_thresh = args.single_read_threshold[1]

    try:
        if args.seq_sample_type is None:
            seq_samp_type = None
        else:
            # sample compare does not have seq_sample_type in the namespace
            seq_samp_type = th.seqSampleType(DNA_SAMP_TYPE, False) \
                            if args.seq_sample_type == DNA_SAMP_TYPE else \
                               th.seqSampleType(RNA_SAMP_TYPE, True)
    except AttributeError:
        seq_samp_type = None

    reads_index = th.TomboReads(
        args.fast5_basedirs, args.corrected_group, args.basecall_subgroups)

    if args.action_command == 'de_novo':
        _de_novo_main(
            args, lower_thresh, single_read_thresh, seq_samp_type, reads_index)
    elif args.action_command == 'alternative_model':
        _alt_main(
            args, lower_thresh, single_read_thresh, seq_samp_type, reads_index)
    elif args.action_command == 'model_sample_compare':
        _model_samp_comp_main(
            args, lower_thresh, single_read_thresh, seq_samp_type, reads_index)
    elif args.action_command == 'level_sample_compare':
        _level_samp_comp_main(args, reads_index)
    else:
        th.error_message_and_exit('Invalid Tombo detect_modifications command.')

    return

def _aggregate_per_read_main(args):
    global VERBOSE
    VERBOSE = not args.quiet
    th.VERBOSE = VERBOSE

    if len(args.single_read_threshold) == 1:
        lower_thresh = None
        single_read_thresh = args.single_read_threshold[0]
    else:
        if len(args.single_read_threshold) > 2:
            th.warning_message(
                'Only 1 or 2 values may be passed as single-read ' +
                'thresholds. Only using the first 2 options provided.')
        lower_thresh = args.single_read_threshold[0]
        single_read_thresh = args.single_read_threshold[1]

    aggregate_per_read_stats(
        args.per_read_statistics_filename, single_read_thresh, lower_thresh,
        args.statistics_filename, args.coverage_dampen_counts,
        args.minimum_test_reads, args.num_most_significant_stored,
        args.processes)

    return


if __name__ == '__main__':
    sys.stderr.write('This is a module. See commands with `tombo -h`')
    sys.exit(1)
