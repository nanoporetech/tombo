from __future__ import division, unicode_literals, absolute_import

from builtins import int, range, dict, map, zip

import sys

# Future warning from cython in h5py
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
import h5py

import numpy as np

from operator import itemgetter

if sys.version_info[0] > 2:
    unicode = str

from . import tombo_helper as th


VERBOSE = False


##############################
###### Filter Functions ######
##############################

def clear_filters(fast5s_dir, corr_grp):
    """Clear filters applied to this directories index files
    """
    reads_index = th.TomboReads([fast5s_dir,], corr_grp, remove_filtered=False)

    th.status_message('Clearing all filters.')
    reads_index.replace_index(dict(
        (chrm_strand, [rd._replace(filtered=False) for rd in cs_reads_index])
        for chrm_strand, cs_reads_index in reads_index))

    reads_index.write_index_file()
    th.status_message('All filters successfully cleared!')

    return

def print_filter_mess(
        num_filt_reads, prev_unfilt_reads, total_reads, fast5s_dir, filter_text):
    if prev_unfilt_reads == 0:
        th.error_message_and_exit(
            'No unfiltered reads present in current Tombo index.')

    th.status_message(
        'Filtered {:d} reads ({:.1%} of previously filtered and '.format(
            num_filt_reads, float(num_filt_reads) / prev_unfilt_reads) +
        '{:.1%} of all valid reads)'.format(
            float(num_filt_reads) / total_reads) +
        ' reads due to ' + filter_text + ' filter from ' + fast5s_dir + '.')
    return

def filter_reads_for_stuck(fast5s_dir, corr_grp, obs_filter):
    """Filter reads based on some observation per base threshold criteria
    """
    def read_is_stuck(fast5_fn, s_grp):
        try:
            with h5py.File(fast5_fn, 'r') as fast5_data:
                base_lens = th.get_single_slot_read_centric(
                    fast5_data, 'length', s_grp)
                if base_lens is None: return True
                return any(np.percentile(base_lens, pctl) > thresh
                           for pctl, thresh in obs_filter)
        except:
            return True


    reads_index = th.TomboReads([fast5s_dir,], corr_grp, remove_filtered=False)

    th.status_message('Filtering stuck reads.')
    filt_reads_index = {}
    prev_unfilt_reads, num_filt_reads, total_reads = 0, 0, 0
    for chrm_strand, cs_reads in reads_index:
        cs_filt_reads = []
        for rd in cs_reads:
            total_reads += 1
            if not rd.filtered:
                prev_unfilt_reads += 1
                if read_is_stuck(rd.fn, rd.corr_group):
                    num_filt_reads += 1
                    rd = rd._replace(filtered=True)
            cs_filt_reads.append(rd)
        filt_reads_index[chrm_strand] = cs_filt_reads

    print_filter_mess(num_filt_reads, prev_unfilt_reads, total_reads,
                      fast5s_dir, 'observations per base')

    reads_index.replace_index(filt_reads_index)
    reads_index.write_index_file()

    return

def filter_reads_for_coverage(fast5s_dir, corr_grp, frac_to_filter):
    """Filter reads at higher coverage regions
    """
    reads_index = th.TomboReads([fast5s_dir,], corr_grp, remove_filtered=False)

    th.status_message('Filtering reads to obtain more uniform coverage.')
    filt_reads_index = dict((cs, []) for cs in reads_index.get_all_cs())
    unfilt_reads = []
    unfilt_reads_cov = []
    total_reads = 0
    for chrm_strand, cs_reads in reads_index:
        # TODO: perform coverage computation iteratively for larger fractions
        # of reads requested to avoid canceling out high coverage locations

        # compute coverage
        max_end = max(rd.end for rd in cs_reads)
        cs_coverage = np.zeros(max_end, dtype=np.int64)
        for rd in cs_reads:
            total_reads += 1
            if rd.filtered:
                # add previously filtered reads straight back to new index
                filt_reads_index[chrm_strand].append(rd)
            cs_coverage[rd.start:rd.end] += 1
        # assign coverage value to each read
        for rd in (rd for rd in cs_reads if not rd.filtered):
            # add approximate coverage from middle of read
            # faster than mean over the whole read
            unfilt_reads_cov.append(cs_coverage[
                rd.start + ((rd.end - rd.start) // 2)])
            unfilt_reads.append((chrm_strand, rd))

    prev_unfilt_reads = len(unfilt_reads)
    if prev_unfilt_reads == 0:
        th.error_message_and_exit(
            'No unfiltered reads present in current Tombo index.')
    num_filt_reads = int(frac_to_filter * prev_unfilt_reads)
    print_filter_mess(num_filt_reads, prev_unfilt_reads, total_reads,
                      fast5s_dir, 'even coverage')

    # create probabilities array with coverage values normalized to sum to 1
    unfilt_reads_cov = np.array(unfilt_reads_cov, dtype=np.float)
    unfilt_reads_p = unfilt_reads_cov / unfilt_reads_cov.sum()
    # randomly chose reads to filter
    filt_indices = set(np.random.choice(
        prev_unfilt_reads, size=num_filt_reads, replace=False, p=unfilt_reads_p))
    for i, (chrm_strand, rd) in enumerate(unfilt_reads):
        if i in filt_indices:
            rd = rd._replace(filtered=True)
        filt_reads_index[chrm_strand].append(rd)

    reads_index.replace_index(filt_reads_index)
    reads_index.write_index_file()

    return

def filter_reads_for_qscore(fast5s_dir, bc_grp, corr_grp, q_score_thresh):
    """Filter reads based on mean q-score
    """
    def read_fails_q_score(fast5_fn, s_grp):
        try:
            with h5py.File(fast5_fn, 'r') as fast5_data:
                r_q_scores = fast5_data['/Analyses/' + bc_grp + '/' + s_grp +
                                        '/Fastq'][:].decode().split('\n')[3]
                return th.get_mean_q_score(r_q_scores) < q_score_thresh
        except:
            return True


    reads_index = th.TomboReads([fast5s_dir,], corr_grp, remove_filtered=False)

    th.status_message('Filtering reads below a mean q-score cutoff.')
    filt_reads_index = dict((cs, []) for cs in reads_index.get_all_cs())
    num_filt_reads, prev_unfilt_reads, total_reads = 0, 0, 0
    for chrm_strand, cs_reads in reads_index:
        for rd in cs_reads:
            total_reads += 1
            if rd.filtered:
                filt_reads_index[chrm_strand].append(rd)
                continue
            prev_unfilt_reads += 1
            if rd.mean_q_score is None:
                filter_read = read_fails_q_score(
                    rd.fn, rd.corr_group.split('/')[-1])
            else:
                filter_read = rd.mean_q_score < q_score_thresh
            if filter_read:
                num_filt_reads += 1
                rd = rd._replace(filtered=True)
            filt_reads_index[chrm_strand].append(rd)

    print_filter_mess(num_filt_reads, prev_unfilt_reads, total_reads,
                      fast5s_dir, 'q-score')

    reads_index.replace_index(filt_reads_index)
    reads_index.write_index_file()

    return

def filter_reads_for_signal_matching(fast5s_dir, corr_grp, sig_match_thresh):
    """Filter reads based on observed to expected signal matching score
    """
    def read_fails_matching_score(fast5_fn, corr_group):
        try:
            with h5py.File(fast5_fn, 'r') as fast5_data:
                return fast5_data['/Analyses/' + corr_group].attrs.get(
                    'signal_match_score') > sig_match_thresh
        except:
            return True


    reads_index = th.TomboReads([fast5s_dir,], corr_grp, remove_filtered=False)

    th.status_message('Filtering reads above a signal matching score threshold.')
    filt_reads_index = dict((cs, []) for cs in reads_index.get_all_cs())
    num_filt_reads, prev_unfilt_reads, total_reads = 0, 0, 0
    for chrm_strand, cs_reads in reads_index:
        for rd in cs_reads:
            total_reads += 1
            if rd.filtered:
                filt_reads_index[chrm_strand].append(rd)
                continue
            prev_unfilt_reads += 1
            if rd.sig_match_score is None:
                filter_read = read_fails_matching_score(rd.fn, rd.corr_group)
            else:
                filter_read = rd.sig_match_score > sig_match_thresh
            if filter_read:
                num_filt_reads += 1
                rd = rd._replace(filtered=True)
            filt_reads_index[chrm_strand].append(rd)

    print_filter_mess(num_filt_reads, prev_unfilt_reads, total_reads,
                      fast5s_dir, 'signal matching')

    reads_index.replace_index(filt_reads_index)
    reads_index.write_index_file()

    return

def filter_reads_for_genome_pos(
        fast5s_dir, corr_grp, include_regs, include_partial=False):
    """Filter reads to include or exclude genomic regions
    """
    def read_included(start, end, chrm_include_regs):
        if chrm_include_regs is None:
            return True
        if include_partial:
            # include all reads partially overlapping regions
            return any(not (start > i_end or end < i_start)
                       for i_start, i_end in chrm_include_regs)

        # test if read is completely contained within the interval
        return any((start >= i_start and end <= i_end)
                   for i_start, i_end in chrm_include_regs)


    reads_index = th.TomboReads([fast5s_dir,], corr_grp, remove_filtered=False)

    th.status_message('Filtering reads outside of the specified ' +
                      'genomic location.')
    filt_reads_index = dict((cs, []) for cs in reads_index.get_all_cs())
    num_filt_reads, prev_unfilt_reads, total_reads = 0, 0, 0
    for (chrm, strand), cs_reads in reads_index:
        do_filter_cs_reads = chrm not in include_regs
        for rd in cs_reads:
            total_reads += 1
            if rd.filtered:
                filt_reads_index[(chrm, strand)].append(rd)
                continue
            prev_unfilt_reads += 1
            if do_filter_cs_reads or not read_included(
                    rd.start, rd.end, include_regs[chrm]):
                num_filt_reads += 1
                rd = rd._replace(filtered=True)
            filt_reads_index[(chrm, strand)].append(rd)

    print_filter_mess(num_filt_reads, prev_unfilt_reads, total_reads,
                      fast5s_dir, 'genomic position')

    reads_index.replace_index(filt_reads_index)
    reads_index.write_index_file()

    return


###################################
###### Filter Main Functions ######
###################################

def _clear_filters_main(args):
    for fast5s_dir in args.fast5_basedirs:
        clear_filters(fast5s_dir, args.corrected_group)

    return

def _filter_stuck_main(args):
    obs_filter = th.parse_obs_filter(args.obs_per_base_filter)
    for fast5s_dir in args.fast5_basedirs:
        filter_reads_for_stuck(fast5s_dir, args.corrected_group, obs_filter)

    return

def _filter_coverage_main(args):
    if not 0 < args.percent_to_filter < 100:
        th.error_message_and_exit(
            '--percent-to-filter must be between 0 and 100.')

    for fast5s_dir in args.fast5_basedirs:
        filter_reads_for_coverage(
            fast5s_dir, args.corrected_group, args.percent_to_filter / 100.0)

    return

def _filter_q_score_main(args):
    if not 0 < args.q_score < 40:
        th.error_message_and_exit('--q-score must be between 0 and 40.')

    for fast5s_dir in args.fast5_basedirs:
        filter_reads_for_qscore(
            fast5s_dir, args.basecall_group, args.corrected_group, args.q_score)

    return

def _filter_signal_matching_main(args):
    if not 0 < args.signal_matching_score < 10:
        th.error_message_and_exit(
            '--signal-matching-score must be between 0 and 10.')

    for fast5s_dir in args.fast5_basedirs:
        filter_reads_for_signal_matching(
            fast5s_dir, args.corrected_group, args.signal_matching_score)

    return

def _filter_genome_pos_main(args):
    include_regs = th.parse_genome_regions(args.include_regions)

    for fast5s_dir in args.fast5_basedirs:
        filter_reads_for_genome_pos(
            fast5s_dir, args.corrected_group, include_regs,
            args.include_partial_overlap)

    return

def filter_main(args):
    global VERBOSE
    VERBOSE = not args.quiet
    th.VERBOSE = VERBOSE

    if args.action_command == 'clear_filters':
        _clear_filters_main(args)
    elif args.action_command == 'genome_locations':
        _filter_genome_pos_main(args)
    elif args.action_command == 'stuck':
        _filter_stuck_main(args)
    elif args.action_command == 'level_coverage':
        _filter_coverage_main(args)
    elif args.action_command == 'q_score':
        _filter_q_score_main(args)
    elif args.action_command == 'raw_signal_matching':
        _filter_signal_matching_main(args)
    else:
        th.error_message_and_exit('Invalid Tombo filter command.')

    return


if __name__ == '__main__':
    sys.stderr.write('This is a module. See commands with `tombo -h`')
    sys.exit(1)
