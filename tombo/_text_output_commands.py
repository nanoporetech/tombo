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
OUT_TYPES = {'wig':'wiggle_0', 'bedgraph':'bedGraph'}
GROUP1_NAME='sample'
GROUP2_NAME='control'


########################
###### WIG Output ######
########################

def open_browser_files(wig_base, group_text, type_name, out_type='wig'):
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

def write_frac_wigs(all_stats, wig_base, do_frac, do_damp, do_valid_cov):
    if VERBOSE: th._status_message(
            'Parsing and outputting statistics wiggles.')
    if do_frac:
        plus_frac_fp, minus_frac_fp = open_browser_files(
            wig_base, '', 'fraction_modified_reads')
    if do_damp:
        plus_damp_fp, minus_damp_fp = open_browser_files(
            wig_base, '', 'dampened_fraction_modified_reads')
    if do_valid_cov:
        plus_vcov_fp, minus_vcov_fp = open_browser_files(
            wig_base, '', 'valid_coverage')

    (curr_chrm, curr_strand, curr_poss, curr_fracs, curr_damp_fracs,
     curr_valid_cov) = (None, None, [], [], [], [])
    all_stats.order_by_pos()
    for chrm, strand, pos, frac, damp_frac, valid_cov in all_stats.iter_fracs():
        if chrm != curr_chrm or strand != curr_strand:
            if len(curr_poss) > 0:
                # write current chrm/strand data
                if do_frac:
                    wig_fp = plus_frac_fp if curr_strand == '+' else minus_frac_fp
                    _write_cs_data(wig_fp, curr_chrm, curr_poss, curr_fracs)
                if do_damp:
                    wig_fp = plus_damp_fp if curr_strand == '+' else minus_damp_fp
                    _write_cs_data(wig_fp, curr_chrm, curr_poss, curr_damp_fracs)
                if do_valid_cov:
                    wig_fp = plus_vcov_fp if curr_strand == '+' else minus_vcov_fp
                    _write_cs_int_data(
                        wig_fp, curr_chrm, curr_poss, curr_valid_cov)

            # set new chrm and strand and empty lists
            curr_chrm, curr_strand = chrm, strand
            curr_poss, curr_fracs, curr_damp_fracs, curr_valid_cov = (
                [], [], [], [])

        # store position statistics
        curr_poss.append(pos)
        if do_frac:
            curr_fracs.append(1 - frac)
        if do_damp:
            curr_damp_fracs.append(1 - damp_frac)
        if do_valid_cov:
            curr_valid_cov.append(valid_cov)

    # write last chrm/strand data
    if len(curr_poss) > 0:
        if do_frac:
            wig_fp = plus_frac_fp if curr_strand == '+' else minus_frac_fp
            _write_cs_data(wig_fp, curr_chrm, curr_poss, curr_fracs)
        if do_damp:
            wig_fp = plus_damp_fp if curr_strand == '+' else minus_damp_fp
            _write_cs_data(wig_fp, curr_chrm, curr_poss, curr_damp_fracs)
        if do_valid_cov:
            wig_fp = plus_vcov_fp if curr_strand == '+' else minus_vcov_fp
            _write_cs_int_data(wig_fp, curr_chrm, curr_poss, curr_valid_cov)

    if do_frac:
        plus_frac_fp.close()
        minus_frac_fp.close()
    if do_damp:
        plus_damp_fp.close()
        minus_damp_fp.close()
    if do_valid_cov:
        plus_vcov_fp.close()
        minus_vcov_fp.close()

    return

def filter_cs_nans(cs_vals):
    valid_poss = np.where(~np.isnan(cs_vals))[0]
    valid_vals = cs_vals[valid_poss]
    return valid_poss, valid_vals

def write_length_wig(
        raw_read_coverage, chrm_sizes, wig_base, group_name):
    if VERBOSE: th._status_message('Parsing and outputting ' + group_name +
                                   ' dwell times.')
    plus_dwell_fp, minus_dwell_fp = open_browser_files(
        wig_base, group_name, 'dwell')
    for chrm, strand, cs_vals in th.iter_mean_slot_values(
            raw_read_coverage, chrm_sizes, 'length'):
        dwell_fp = plus_dwell_fp if strand == '+' else minus_dwell_fp
        cs_poss, cs_vals = filter_cs_nans(cs_vals)
        _write_cs_data(dwell_fp, chrm, cs_poss, cs_vals)

    plus_dwell_fp.close()
    minus_dwell_fp.close()

    return

def write_signal_sd_wig(
        raw_read_coverage, chrm_sizes, wig_base, group_name):
    if VERBOSE: th._status_message('Parsing and outputting ' + group_name +
                                   ' signal SDs.')
    plus_sd_fp, minus_sd_fp = open_browser_files(
        wig_base, group_name, 'signal_sd')
    for chrm, strand, cs_vals in th.iter_mean_slot_values(
            raw_read_coverage, chrm_sizes, 'norm_stdev'):
        sd_fp = plus_sd_fp if strand == '+' else minus_sd_fp
        cs_poss, cs_vals = filter_cs_nans(cs_vals)
        _write_cs_data(sd_fp, chrm, cs_poss, cs_vals)

    plus_sd_fp.close()
    minus_sd_fp.close()

    return

def write_signal_and_diff_wigs(
        raw_read_coverage1, raw_read_coverage2, chrm_sizes,
        wig_base, group1_name, write_sig, write_diff):
    if VERBOSE: th._status_message(
            'Parsing and outputting signal means and differences.')
    # open all file pointers
    if write_sig:
        plus_sig1_fp, minus_sig1_fp = open_browser_files(
            wig_base, group1_name, 'signal')
        if raw_read_coverage2 is not None:
            plus_sig2_fp, minus_sig2_fp = open_browser_files(
                wig_base, GROUP2_NAME, 'signal')
    if write_diff:
        plus_diff_fp, minus_diff_fp = open_browser_files(
            wig_base, '', 'difference')

    # iterate over mean signal values for all chrm/strand combinations with
    # coverage in either sample. None returned if one sample is not covered
    for chrm, strand, cs_sig_means1, cs_sig_means2 in th.iter_mean_slot_values(
            raw_read_coverage1, chrm_sizes, 'norm_mean', raw_read_coverage2):
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

def write_cov_wig(raw_read_coverage, out_base, group_text):
    if VERBOSE: th._status_message('Getting and writing ' + group_text +
                                   ' coverage bedgraphs.')
    plus_bg_fp, minus_bg_fp = open_browser_files(
        out_base, group_text, 'coverage', 'bedgraph')
    for chrm, strand, cs_cov, cs_cov_starts in th.get_coverage_regions(
            raw_read_coverage):
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
        f5_dirs1, f5_dirs2, corr_grp, bc_subgrps,
        stats_fn, wig_base, wig_types, cov_damp_counts):
    if f5_dirs1 is not None:
        raw_read_coverage1 = th.parse_fast5s(
            f5_dirs1, corr_grp, bc_subgrps, sample_name='sample')
        if len(raw_read_coverage1) == 0:
            th._error_message_and_exit(
                'No reads present in --fast5-basedirs.')

    group1_name = '' if f5_dirs2 is None else GROUP1_NAME
    if f5_dirs2 is not None:
        raw_read_coverage2 = th.parse_fast5s(
            f5_dirs2, corr_grp, bc_subgrps, sample_name='control')
        chrm_sizes = th.get_chrm_sizes(
            raw_read_coverage1, raw_read_coverage2)

        if 'coverage' in wig_types:
            write_cov_wig(raw_read_coverage2, wig_base, GROUP2_NAME)
        if 'signal_sd' in wig_types:
            write_signal_sd_wig(
                raw_read_coverage2, chrm_sizes, wig_base, GROUP2_NAME)
        if 'dwell' in wig_types:
            write_length_wig(raw_read_coverage2, chrm_sizes,
                             wig_base, GROUP2_NAME)

        # need to do signal and difference call once either with or
        # w/o second set of files (unlike coverage, sds and length
        if 'signal' in wig_types or 'difference' in wig_types:
            write_signal_and_diff_wigs(
                raw_read_coverage1, raw_read_coverage2, chrm_sizes,
                wig_base, group1_name, 'signal' in wig_types,
                'difference' in wig_types)
    elif f5_dirs1 is not None:
        chrm_sizes = th.get_chrm_sizes(raw_read_coverage1)
        if 'signal' in wig_types:
            write_signal_and_diff_wigs(
                raw_read_coverage1, None, chrm_sizes, wig_base,
                group1_name, 'signal' in wig_types, False)

    if 'coverage' in wig_types:
        write_cov_wig(raw_read_coverage1, wig_base, group1_name)
    if 'signal_sd' in wig_types:
        write_signal_sd_wig(
            raw_read_coverage1, chrm_sizes, wig_base, group1_name)
    if 'dwell' in wig_types:
        write_length_wig(raw_read_coverage1, chrm_sizes, wig_base, group1_name)
    if any(wig_type in wig_types for wig_type in (
            'fraction', 'dampened_fraction', 'valid_coverge')):
        if VERBOSE: th._status_message('Loading statistics from file.')
        all_stats = ts.TomboStats(stats_fn)
        if 'dampened_fraction' in wig_types:
            all_stats.calc_damp_fraction(cov_damp_counts)
        write_frac_wigs(all_stats, wig_base,
                        'fraction' in wig_types,
                        'dampened_fraction' in wig_types,
                        'valid_coverage' in wig_types)

    return


##########################
###### FASTA Output ######
##########################

def write_most_signif(
        f5_dirs, fasta_fn, num_regions, corr_grp, bc_subgrps, seqs_fn,
        num_bases, stats_fn, cov_damp_counts):
    if VERBOSE: th._status_message('Loading statistics from file.')
    plot_intervals = ts.TomboStats(stats_fn).get_most_signif_regions(
        num_bases, num_regions, cov_damp_counts=cov_damp_counts)

    # get each regions sequence either from reads or fasta index
    if fasta_fn is None:
        raw_read_coverage = th.parse_fast5s(f5_dirs, corr_grp, bc_subgrps)
        all_reg_data = th.get_region_sequences(
            plot_intervals, raw_read_coverage)
    else:
        genome_index = th.Fasta(fasta_fn)
        all_reg_data = [
            int_i._replace(
                seq=genome_index.get_seq(int_i.chrm, int_i.start, int_i.end))
            for int_i in plot_intervals if int_i.chrm in genome_index]

    if VERBOSE: th._status_message('Outputting region seqeuences.')
    with io.open(seqs_fn, 'wt') as seqs_fp:
        for int_i in all_reg_data:
            reg_seq = int_i.seq
            if int_i.strand == '-':
                reg_seq = th.rev_comp(reg_seq)
            seqs_fp.write('>{0}:{1:d}:{2} {3}\n{4}\n'.format(
                int_i.chrm, int(int_i.start + (num_bases // 2)),
                int_i.strand, int_i.reg_text, ''.join(reg_seq)))

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
            for data_type in ['signal', 'difference', 'coverage',
                              'signal_sd', 'dwell']) and
        args.fast5_basedirs is None):
        th._error_message_and_exit(
            'Must provide a fast5 basedir to output signal, difference, ' +
            'coverage, signal_sd and/or length browser files.')
    if (any(wig_type in args.file_types for wig_type in (
            'fraction', 'dampened_fraction', 'valid_coverage')) and
        args.statistics_filename is None):
        th._error_message_and_exit(
            'Must provide a statistics filename to output ' +
            'fraction or valid coverage browser files.')
    if ('difference' in args.file_types and
        args.control_fast5_basedirs is None):
        th._error_message_and_exit(
            'Must provide two sets of FAST5s ' + \
            'to output difference wiggle files.')
    if (args.control_fast5_basedirs is not None and
        args.fast5_basedirs is None):
        th._error_message_and_exit(
            'Cannot provide a control FAST5 set of directories ' +
            'without a sample set of FAST5 directories.')
    if (args.coverage_dampen_counts is None and
        'dampened_fraction' in args.file_types):
        th._error_message_and_exit(
            'Cannot compute dampened fractions without ' +
            '--coverage-dampened-counts values.')

    write_all_browser_files(
        args.fast5_basedirs, args.control_fast5_basedirs, args.corrected_group,
        args.basecall_subgroups, args.statistics_filename,
        args.browser_file_basename, args.file_types,
        args.coverage_dampen_counts)

    return

def _write_signif_diff_main(args):
    global VERBOSE
    VERBOSE = not args.quiet
    th.VERBOSE = VERBOSE
    ts.VERBOSE = VERBOSE

    if args.fast5_basedirs is None and args.genome_fasta is None:
        th._error_message_and_exit(
            'Must provide either FAST5 directory(ies) or a fasta file.')

    write_most_signif(
        args.fast5_basedirs, args.genome_fasta, args.num_regions,
        args.corrected_group, args.basecall_subgroups, args.sequences_filename,
        args.num_bases, args.statistics_filename, args.coverage_dampen_counts)

    return


if __name__ == '__main__':
    raise NotImplementedError(
        'This is a module. See commands with `tombo -h`')
