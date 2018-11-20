from __future__ import division, unicode_literals, absolute_import

from builtins import int, range, dict, map, zip

import io
import sys

import numpy as np

from collections import defaultdict

if sys.version_info[0] > 2:
    unicode = str

# import tombo functions
from . import tombo_stats as ts
from . import tombo_helper as th

from ._default_parameters import SMALLEST_PVAL

VERBOSE = False

OUT_HEADER='track type={0} name="{1}_{2}_{3}{4}" ' + \
    'description="{1} {2} {3}{5}"\n'
BG_TYPE = 'bedgraph'
WIG_TYPE = 'wig'
OUT_TYPES = {WIG_TYPE:'wiggle_0', BG_TYPE:'bedGraph'}
GROUP_NAME='sample'
CTRL_NAME='control'

COV_WIG_TYPE = 'coverage'

# event table slot values
SIG_SLOT = 'norm_mean'
SD_SLOT = 'norm_stdev'
DWELL_SLOT = 'length'
SIG_WIG_TYPE = 'signal'
DIFF_WIG_TYPE = 'difference'
SD_WIG_TYPE = 'signal_sd'
DWELL_WIG_TYPE = 'dwell'

# stat table slot values
POS_SLOT = 'pos'
FRAC_SLOT = 'frac'
DFRAC_SLOT = 'damp_frac'
STAT_SLOT = 'stat'
VCOV_SLOT = 'valid_cov'

FRAC_WIG_TYPE = 'fraction'
DFRAC_WIG_TYPE = 'dampened_fraction'
STAT_WIG_TYPE = 'statistic'
VCOV_WIG_TYPE = 'valid_coverage'

FRAC_WIG_NAME = 'fraction_modified_reads'
DFRAC_WIG_NAME = 'dampened_fraction_modified_reads'
STAT_WIG_NAME = 'statistic'
VCOV_WIG_NAME = 'valid_coverage'


########################
###### WIG Output ######
########################

def open_browser_files(wig_base, group_text, type_name, out_type=WIG_TYPE):
    group_w_dot = '' if group_text == '' else '.' + group_text
    group_w_us = '' if group_text == '' else '_' + group_text
    group_w_space = '' if group_text == '' else ' ' + group_text
    plus_wig_fp = io.open(
        wig_base + '.' + type_name + group_w_dot + '.plus.' + out_type, 'wt')
    minus_wig_fp = io.open(
        wig_base + '.' + type_name + group_w_dot + '.minus.' + out_type, 'wt')
    plus_wig_fp.write(OUT_HEADER.format(
        OUT_TYPES[out_type], wig_base, type_name, 'fwd_strand',
        group_w_us, group_w_space))
    minus_wig_fp.write(OUT_HEADER.format(
        OUT_TYPES[out_type], wig_base, type_name, 'rev_strand',
        group_w_us, group_w_space))

    return plus_wig_fp, minus_wig_fp

def _write_cs_data(wig_fp, chrm, cs_poss, cs_vals):
    wig_fp.write("variableStep chrom={} span=1\n".format(chrm))
    wig_fp.write('\n'.join(['{:d} {:.4f}'.format(x[0] + 1, x[1])
                            for x in zip(cs_poss, cs_vals)]) + '\n')

    return

def _write_cs_int_data(wig_fp, chrm, cs_poss, cs_vals):
    wig_fp.write("variableStep chrom={} span=1\n".format(chrm))
    wig_fp.write('\n'.join(['{:d} {:d}'.format(x[0] + 1, x[1])
                            for x in zip(cs_poss, cs_vals)]) + '\n')

    return

def write_frac_wigs(
        all_stats, wig_base, do_frac, do_damp, do_stats, do_vcov, fasta_fn,
        motif_descs):
    if VERBOSE: th.status_message('Parsing and outputting statistics wiggles.')
    filter_motifs = not(fasta_fn is None or motif_descs is None)
    if filter_motifs:
        genome_index = th.Fasta(fasta_fn)
        motifs = [
            motif for motif in th.parse_motif_descs('::'.join(motif_descs))]
    else:
        motifs = [(None, ''),]

    if do_frac:
        frac_fps = dict(
            kv for _, mod_name in motifs
            for kv in zip(
                    (('+', mod_name), ('-', mod_name)),
                    open_browser_files(wig_base, mod_name, FRAC_WIG_NAME)))
    if do_damp:
        damp_fps = dict(
            kv for _, mod_name in motifs
            for kv in zip(
                    (('+', mod_name), ('-', mod_name)),
                    open_browser_files(wig_base, mod_name, DFRAC_WIG_NAME)))
    if do_stats:
        stat_fps = dict(
            kv for _, mod_name in motifs
            for kv in zip(
                    (('+', mod_name), ('-', mod_name)),
                    open_browser_files(wig_base, mod_name, STAT_WIG_NAME)))
    if do_vcov:
        vcov_fps = dict(
            kv for _, mod_name in motifs
            for kv in zip(
                    (('+', mod_name), ('-', mod_name)),
                    open_browser_files(wig_base, mod_name, VCOV_WIG_NAME)))


    def write_cs_stats(
            curr_chrm, curr_strand, curr_poss,
            curr_fracs, curr_dampf, curr_stats, curr_vcov):
        curr_poss = np.concatenate(curr_poss)
        if do_frac: curr_fracs = np.concatenate(curr_fracs)
        if do_damp: curr_dampf = np.concatenate(curr_dampf)
        if do_stats: curr_stats = np.concatenate(curr_stats)
        if do_vcov: curr_vcov = np.concatenate(curr_vcov)

        if filter_motifs:
            chrm_seq = genome_index.get_seq(curr_chrm)
        for motif, mod_name in motifs:
            if do_frac:
                frac_fp = frac_fps[(curr_strand, mod_name)]
                motif_fracs = curr_fracs
            if do_damp:
                damp_fp = damp_fps[(curr_strand, mod_name)]
                motif_damp = curr_dampf
            if do_stats:
                stat_fp = stat_fps[(curr_strand, mod_name)]
                motif_stats = curr_stats
            if do_vcov:
                vcov_fp = vcov_fps[(curr_strand, mod_name)]
                motif_vcov = curr_vcov
            m_curr_poss = curr_poss
            if filter_motifs:
                if curr_strand == '-':
                    chrm_motif_poss = np.array([
                        m.start() + motif.motif_len - motif.mod_pos
                        for m in motif.rev_comp_pat.finditer(chrm_seq)])
                else:
                    chrm_motif_poss = np.array([
                        m.start() + motif.mod_pos - 1
                        for m in motif.motif_pat.finditer(chrm_seq)])
                valid_poss = np.isin(
                    curr_poss, chrm_motif_poss, assume_unique=True)
                m_curr_poss = curr_poss[valid_poss]
                if do_frac: motif_fracs = curr_fracs[valid_poss]
                if do_damp: motif_damp = curr_dampf[valid_poss]
                if do_stats: motif_stats = curr_stats[valid_poss]
                if do_vcov: motif_vcov = curr_vcov[valid_poss]

            # write current chrm/strand data
            if do_frac:
                _write_cs_data(frac_fp, curr_chrm, m_curr_poss, motif_fracs)
            if do_damp:
                _write_cs_data(damp_fp, curr_chrm, m_curr_poss, motif_damp)
            if do_stats:
                _write_cs_data(stat_fp, curr_chrm, m_curr_poss, motif_stats)
            if do_vcov:
                _write_cs_int_data(vcov_fp, curr_chrm, m_curr_poss, motif_vcov)

        return


    (curr_chrm, curr_strand, curr_poss, curr_fracs, curr_dampf, curr_stats,
     curr_vcov) = (None, None, [], [], [], [], [])
    for chrm, strand, start, end, block_stats in all_stats:
        if chrm != curr_chrm or strand != curr_strand:
            if len(curr_poss) > 0:
                write_cs_stats(
                    curr_chrm, curr_strand, curr_poss, curr_fracs,
                    curr_dampf, curr_stats, curr_vcov)

            # set new chrm and strand and empty lists
            curr_chrm, curr_strand = chrm, strand
            curr_poss, curr_fracs, curr_dampf, curr_stats, curr_vcov = (
                [], [], [], [], [])

        # store block statistics
        curr_poss.append(block_stats[POS_SLOT])
        if do_frac: curr_fracs.append(1 - block_stats[FRAC_SLOT])
        if do_damp: curr_dampf.append(1 - block_stats[DFRAC_SLOT])
        if do_stats: curr_stats.append(all_stats._stat_transform(block_stats))
        if do_vcov: curr_vcov.append(block_stats[VCOV_SLOT])

    # write last chrm/strand data
    if len(curr_poss) > 0:
        write_cs_stats(
            curr_chrm, curr_strand, curr_poss,
            curr_fracs, curr_dampf, curr_stats, curr_vcov)

    if do_frac:
        for wig_fp in frac_fps.values():
            wig_fp.close()
    if do_damp:
        for wig_fp in damp_fps.values():
            wig_fp.close()
    if do_stats:
        for wig_fp in stat_fps.values():
            wig_fp.close()
    if do_vcov:
        for wig_fp in vcov_fps.values():
            wig_fp.close()

    return

def filter_cs_nans(cs_vals):
    valid_poss = np.where(~np.isnan(cs_vals))[0]
    valid_vals = cs_vals[valid_poss]
    return valid_poss, valid_vals

def write_slot_mean_wig(
        reads_index, chrm_sizes, wig_base, group_name, wig_type, slot_name):
    if VERBOSE: th.status_message(
            'Parsing and outputting ' + group_name + ' ' + wig_type + '.')
    plus_sd_fp, minus_sd_fp = open_browser_files(wig_base, group_name, wig_type)
    for chrm, strand, cs_vals, _ in th.iter_mean_slot_values(
            reads_index, chrm_sizes, slot_name):
        sd_fp = plus_sd_fp if strand == '+' else minus_sd_fp
        cs_poss, cs_vals = filter_cs_nans(cs_vals)
        _write_cs_data(sd_fp, chrm, cs_poss, cs_vals)

    plus_sd_fp.close()
    minus_sd_fp.close()

    return

def write_signal_and_diff_wigs(
        reads_index, ctrl_reads_index, chrm_sizes,
        wig_base, group_name, write_sig, write_diff):
    if VERBOSE: th.status_message(
            'Parsing and outputting signal means and differences.')
    # open all file pointers
    if write_sig:
        plus_sig1_fp, minus_sig1_fp = open_browser_files(
            wig_base, group_name, SIG_WIG_TYPE)
        if ctrl_reads_index is not None:
            plus_sig2_fp, minus_sig2_fp = open_browser_files(
                wig_base, CTRL_NAME, SIG_WIG_TYPE)
    if write_diff:
        plus_diff_fp, minus_diff_fp = open_browser_files(
            wig_base, '', DIFF_WIG_TYPE)

    # iterate over mean signal values for all chrm/strand combinations with
    # coverage in either sample. None returned if one sample is not covered
    for chrm, strand, cs_sig_means1, cs_sig_means2 in th.iter_mean_slot_values(
            reads_index, chrm_sizes, SIG_SLOT, ctrl_reads_index):
        # compute valid positions since it will either be used here for signal
        # output or for diff below
        # note small wasted effort for diff only output when second sample
        # does not have coverage
        if cs_sig_means1 is not None:
            cs_poss1, cs_means1 = filter_cs_nans(cs_sig_means1)
            if write_sig:
                sig1_fp = plus_sig1_fp if strand == '+' else minus_sig1_fp
                _write_cs_data(sig1_fp, chrm, cs_poss1, cs_means1)

        if cs_sig_means2 is not None:
            # ocmpute filtered poss since it will be used for either signal
            # diff (or both) outputs
            cs_poss2, cs_means2 = filter_cs_nans(cs_sig_means2)
            if write_sig:
                sig2_fp = plus_sig2_fp if strand == '+' else minus_sig2_fp
                _write_cs_data(sig2_fp, chrm, cs_poss2, cs_means2)

            # write diff values if both samples have coverage
            if cs_sig_means1 is not None and write_diff:
                diff_fp = plus_diff_fp if strand == '+' else minus_diff_fp
                valid_diff_poss = np.intersect1d(
                    cs_poss1, cs_poss2, assume_unique=True)
                cs_diffs = (cs_sig_means1[valid_diff_poss] -
                            cs_sig_means2[valid_diff_poss])
                _write_cs_data(diff_fp, chrm, valid_diff_poss, cs_diffs)

    return

def write_cov_wig(reads_index, out_base, group_text):
    if VERBOSE: th.status_message('Getting and writing ' + group_text +
                                  ' coverage bedgraphs.')
    plus_bg_fp, minus_bg_fp = open_browser_files(
        out_base, group_text, COV_WIG_TYPE, BG_TYPE)
    for (chrm, strand, cs_cov,
         cs_cov_starts) in reads_index.iter_coverage_regions():
        # extract only values from each region and convert to str
        cs_cov = np.char.mod('%d', cs_cov)
        cs_cov_starts = np.char.mod('%d', cs_cov_starts)

        bg_fp = plus_bg_fp if strand == '+' else minus_bg_fp
        bg_fp.write(
            '\n'.join('\t'.join((
                chrm, cs_cov_starts[i], cs_cov_starts[i + 1],
                cs_cov[i])) for i in range(cs_cov.shape[0])) + '\n')

    plus_bg_fp.close()
    minus_bg_fp.close()

    return

def write_all_browser_files(
        fast5s_dirs, ctrl_fast5s_dirs, corr_grp, bc_subgrps,
        stats_fn, wig_base, wig_types, motif_descs, fasta_fn):
    if fast5s_dirs is not None:
        reads_index = th.TomboReads(fast5s_dirs, corr_grp, bc_subgrps)
        if reads_index.is_empty():
            th.error_message_and_exit('No reads present in --fast5-basedirs.')

    group_name = '' if ctrl_fast5s_dirs is None else GROUP_NAME
    if ctrl_fast5s_dirs is not None:
        ctrl_reads_index = th.TomboReads(ctrl_fast5s_dirs, corr_grp, bc_subgrps)
        chrm_sizes = th.get_chrm_sizes(reads_index, ctrl_reads_index)

        if COV_WIG_TYPE in wig_types:
            write_cov_wig(ctrl_reads_index, wig_base, CTRL_NAME)
        if SD_WIG_TYPE in wig_types:
            write_slot_mean_wig(
                ctrl_reads_index, chrm_sizes, wig_base, CTRL_NAME,
                SD_WIG_TYPE, SD_SLOT)
        if DWELL_WIG_TYPE in wig_types:
            write_slot_mean_wig(
                ctrl_reads_index, chrm_sizes, wig_base, CTRL_NAME,
                DWELL_WIG_TYPE, DWELL_SLOT)

        # need to do signal and difference call once either with or
        # w/o second set of files (unlike coverage, sds and length
        if SIG_WIG_TYPE in wig_types or DIFF_WIG_TYPE in wig_types:
            write_signal_and_diff_wigs(
                reads_index, ctrl_reads_index, chrm_sizes,
                wig_base, group_name, SIG_WIG_TYPE in wig_types,
                DIFF_WIG_TYPE in wig_types)
    elif fast5s_dirs is not None:
        chrm_sizes = th.get_chrm_sizes(reads_index)
        if SIG_WIG_TYPE in wig_types:
            write_signal_and_diff_wigs(
                reads_index, None, chrm_sizes, wig_base,
                group_name, SIG_WIG_TYPE in wig_types, False)

    if COV_WIG_TYPE in wig_types:
        write_cov_wig(reads_index, wig_base, group_name)
    if SD_WIG_TYPE in wig_types:
        write_slot_mean_wig(
            reads_index, chrm_sizes, wig_base, CTRL_NAME, SD_WIG_TYPE, SD_SLOT)
    if DWELL_WIG_TYPE in wig_types:
        write_slot_mean_wig(
            reads_index, chrm_sizes, wig_base, CTRL_NAME,
            DWELL_WIG_TYPE, DWELL_SLOT)
    if any(wig_type in wig_types for wig_type in (
            FRAC_WIG_TYPE, DFRAC_WIG_TYPE, STAT_WIG_TYPE, VCOV_WIG_TYPE)):
        if VERBOSE: th.status_message('Loading statistics from file.')
        all_stats = ts.TomboStats(stats_fn)
        if all_stats.is_model_stats and any((
                FRAC_WIG_TYPE in wig_types, DFRAC_WIG_TYPE in wig_types,
                VCOV_WIG_TYPE in wig_types)):
            th.TomboError('Cannot output fraction, dampened_fraction or ' +
                          'valid_coverage for LevelStats statistics.')
        if not all_stats.is_model_stats and STAT_WIG_TYPE in wig_types:
            th.TomboError('Cannot output stat for ModelStats statistics.')
        write_frac_wigs(
            all_stats, wig_base, FRAC_WIG_TYPE in wig_types,
            DFRAC_WIG_TYPE in wig_types, STAT_WIG_TYPE in wig_types,
            VCOV_WIG_TYPE in wig_types, motif_descs, fasta_fn)

    return


##########################
###### FASTA Output ######
##########################

def write_most_signif(
        fast5s_dirs, fasta_fn, num_regions, corr_grp, bc_subgrps, seqs_fn,
        num_bases, stats_fn):
    if VERBOSE: th.status_message('Loading statistics from file.')
    plot_intervals = ts.TomboStats(stats_fn).get_most_signif_regions(
        num_bases, num_regions, prepend_loc_to_text=True)

    # get each regions sequence either from reads or fasta index
    if fasta_fn is None:
        reads_index = th.TomboReads(fast5s_dirs, corr_grp, bc_subgrps)
        for p_int in plot_intervals:
            p_int.add_reads(reads_index).add_seq()
    else:
        genome_index = th.Fasta(fasta_fn)
        for p_int in plot_intervals:
            p_int.add_seq(genome_index)

    if VERBOSE: th.status_message('Outputting region seqeuences.')
    with io.open(seqs_fn, 'wt') as seqs_fp:
        for p_int in plot_intervals:
            reg_seq = (p_int.seq if p_int.strand == '+' else
                       th.rev_comp(p_int.seq))
            seqs_fp.write('>{0}\n{1}\n'.format(p_int.reg_text, ''.join(reg_seq)))

    return


############################
###### Main functions ######
############################

def _browser_files_main(args):
    global VERBOSE
    VERBOSE = not args.quiet
    th.VERBOSE = VERBOSE
    ts.VERBOSE = VERBOSE

    if (any(data_type in args.file_types
            for data_type in [SIG_WIG_TYPE, DIFF_WIG_TYPE, COV_WIG_TYPE,
                              SD_WIG_TYPE, DWELL_WIG_TYPE]) and
        args.fast5_basedirs is None):
        th.error_message_and_exit(
            'Must provide a fast5 basedir to output signal, difference, ' +
            'coverage, signal_sd and/or length browser files.')
    if (any(wig_type in args.file_types for wig_type in (
            FRAC_WIG_TYPE, DFRAC_WIG_TYPE, STAT_WIG_TYPE, VCOV_WIG_TYPE)) and
        args.statistics_filename is None):
        th.error_message_and_exit(
            'Must provide a statistics filename to output ' +
            'fraction or valid coverage browser files.')
    if (DIFF_WIG_TYPE in args.file_types and
        args.control_fast5_basedirs is None):
        th.error_message_and_exit(
            'Must provide two sets of FAST5s ' + \
            'to output difference wiggle files.')
    if (args.control_fast5_basedirs is not None and
        args.fast5_basedirs is None):
        th.error_message_and_exit(
            'Cannot provide a control FAST5 set of directories ' +
            'without a sample set of FAST5 directories.')

    write_all_browser_files(
        args.fast5_basedirs, args.control_fast5_basedirs, args.corrected_group,
        args.basecall_subgroups, args.statistics_filename,
        args.browser_file_basename, args.file_types, args.genome_fasta,
        args.motif_descriptions)

    return

def _write_signif_diff_main(args):
    global VERBOSE
    VERBOSE = not args.quiet
    th.VERBOSE = VERBOSE
    ts.VERBOSE = VERBOSE

    if args.fast5_basedirs is None and args.genome_fasta is None:
        th.error_message_and_exit(
            'Must provide either FAST5 directory(ies) or a fasta file.')

    write_most_signif(
        args.fast5_basedirs, args.genome_fasta, args.num_regions,
        args.corrected_group, args.basecall_subgroups, args.sequences_filename,
        args.num_bases, args.statistics_filename)

    return


if __name__ == '__main__':
    sys.stderr.write('This is a module. See commands with `tombo -h`')
    sys.exit(1)
