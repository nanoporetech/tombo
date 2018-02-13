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

WIG_HEADER='track type=wiggle_0 name="{0}_{1}_{2}{3}" ' + \
    'description="{0} {1} {2}{4}"\n'
GROUP1_NAME='sample'
GROUP2_NAME='control'

def _write_wiggle(wig_base, group_text, data_values, type_name,
                 filter_zeros=False):
    group_w_dot = '' if group_text == '' else '.' + group_text
    group_w_us = '' if group_text == '' else '_' + group_text
    group_w_space = '' if group_text == '' else ' ' + group_text
    plus_wig_fp = io.open(
        wig_base + '.' + type_name + group_w_dot + '.plus.wig', 'wt')
    minus_wig_fp = io.open(
        wig_base + '.' + type_name + group_w_dot + '.minus.wig', 'wt')
    plus_wig_fp.write(WIG_HEADER.format(
        wig_base, type_name, 'fwd_strand', group_w_us, group_w_space))
    minus_wig_fp.write(WIG_HEADER.format(
        wig_base, type_name, 'rev_strand', group_w_us, group_w_space))
    for (chrm, strand), cs_values in data_values.items():
        wig_fp = plus_wig_fp if strand == '+' else minus_wig_fp
        wig_fp.write("variableStep chrom={} span=1\n".format(chrm))
        wig_fp.write('\n'.join([
            unicode(int(pos) + 1) + " " + unicode(round(val, 4))
            for pos, val in enumerate(cs_values)
            if not (np.isnan(val) or (
                    filter_zeros and np.equal(val, 0.0)))]) + '\n')

    plus_wig_fp.close()
    minus_wig_fp.close()

    return

def write_frac_wigs(all_stats, wig_base, do_frac, do_damp, cov_damp_counts):
    if VERBOSE: sys.stderr.write('Parsing statistics.\n')
    raw_chrm_strand_stats = defaultdict(list)
    for stat in all_stats:
        raw_chrm_strand_stats[(
            stat['chrm'].decode(), stat['strand'].decode())].append(
                (stat['pos'], stat['frac'], stat['valid_cov']))

    if do_frac:
        all_frac = {}
    if do_damp:
        all_damp_frac = {}
    for chrm_strand, stats in raw_chrm_strand_stats.items():
        cs_poss, raw_cs_frac, raw_cs_cov = map(np.array, zip(*stats))
        max_pos = max(cs_poss)

        if do_frac:
            cs_frac = np.empty(max_pos + 1)
            cs_frac[:] = np.nan
            # fraction is stored as fraction of unmodified bases, but
            # higher values show better in a wig file, so flip the fractions
            np.put(cs_frac, cs_poss, 1 - raw_cs_frac)
            all_frac[chrm_strand] = cs_frac
        if do_damp:
            cs_damp_frac = np.empty(max_pos + 1)
            cs_damp_frac[:] = np.nan
            non_mod_counts = np.round(raw_cs_frac * raw_cs_cov)
            # compute dampened fraction of modified reads by adding psuedo-counts
            # to the modified and un-modified counts (equivalent to a beta prior
            # on the fraction estimation as a binomial variable)
            raw_cs_damp_frac = (non_mod_counts + cov_damp_counts[0]) / (
                raw_cs_cov + sum(cov_damp_counts))
            # fraction is stored as fraction of unmodified bases, but
            # higher values show better in a wig file, so flip the fractions
            np.put(cs_damp_frac, cs_poss, 1 - raw_cs_damp_frac)
            all_damp_frac[chrm_strand] = cs_damp_frac

    if VERBOSE: sys.stderr.write('Writing fraction wigs.\n')
    if do_frac:
        _write_wiggle(wig_base, '', all_frac, 'fraction_modified_reads')
    if do_damp:
        _write_wiggle(wig_base, '', all_damp_frac,
                      'dampened_fraction_modified_reads')

    return

def write_length_wig(
        raw_read_coverage, chrm_sizes, wig_base, group_name):
    if VERBOSE: sys.stderr.write('Parsing events lengths.\n')
    base_lens = th.get_all_mean_lengths(raw_read_coverage, chrm_sizes)

    if VERBOSE: sys.stderr.write('Writing length wig.\n')
    _write_wiggle(wig_base, group_name, base_lens, 'dwell')

    return

def write_signal_sd_wig(
        raw_read_coverage, chrm_sizes, wig_base, group_name):
    if VERBOSE: sys.stderr.write('Parsing signal SDs.\n')
    base_sds = th.get_all_mean_stdev(raw_read_coverage, chrm_sizes)

    if VERBOSE: sys.stderr.write('Writing signal SD wig.\n')
    _write_wiggle(wig_base, group_name, base_sds, 'signalSd')

    return

def write_signal_and_diff_wigs(
        raw_read_coverage1, raw_read_coverage2, chrm_sizes,
        wig_base, group1_name, write_sig, write_diff):
    if VERBOSE: sys.stderr.write('Parsing mean base signals.\n')
    base_means1 = th.get_all_mean_levels(raw_read_coverage1, chrm_sizes)
    if raw_read_coverage2 is not None:
        base_means2 = th.get_all_mean_levels(raw_read_coverage2, chrm_sizes)

        if write_diff:
            if VERBOSE: sys.stderr.write(
                    'Calculating signal differences.\n')
            sig_diffs = {}
            for chrm, strand in [(c, s) for c in chrm_sizes for s in ('+', '-')]:
                # calculate difference and set no coverage
                # (nan) values to zero
                sig_diffs[(chrm, strand)] \
                    = base_means1[(chrm, strand)] - \
                    base_means2[(chrm, strand)]
            if VERBOSE: sys.stderr.write('Writing differnce wig.\n')
            _write_wiggle(wig_base, '', sig_diffs, 'difference')
        if write_sig:
            if VERBOSE: sys.stderr.write('Writing signal wigs.\n')
            _write_wiggle(wig_base, GROUP2_NAME, base_means2, 'signal')

    if write_sig:
        _write_wiggle(wig_base, group1_name, base_means1, 'signal')

    return

def write_cov_wig(raw_read_coverage, wig_base, group_text):
    read_coverage = th.get_coverage(raw_read_coverage)

    if VERBOSE: sys.stderr.write('Writing coverage wig.\n')
    _write_wiggle(wig_base, group_text, read_coverage, 'coverage', True)

    return

def write_all_wiggles(
        f5_dirs1, f5_dirs2, corrected_group, basecall_subgroups,
        stats_fn, wig_base, wig_types, cov_damp_counts):
    if f5_dirs1 is not None:
        raw_read_coverage1 = th.parse_fast5s(
            f5_dirs1, corrected_group, basecall_subgroups)
        if len(raw_read_coverage1) == 0:
            th._error_message_and_exit(
                'No reads present in --fast5-basedirs.')

    group1_name = '' if f5_dirs2 is None else GROUP1_NAME
    if f5_dirs2 is not None:
        raw_read_coverage2 = th.parse_fast5s(
            f5_dirs2, corrected_group, basecall_subgroups)
        chrm_sizes = th.get_chrm_sizes(
            raw_read_coverage1, raw_read_coverage2)

        if VERBOSE: sys.stderr.write('Writing wiggles.\n')
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
        if VERBOSE: sys.stderr.write('Writing wiggles.\n')
        if 'signal' in wig_types:
            write_signal_and_diff_wigs(
                raw_read_coverage1, None, chrm_sizes, wig_base,
                group1_name, 'signal' in wig_types, False)

    if 'coverage' in wig_types:
        write_cov_wig(raw_read_coverage1, wig_base, group1_name)
    if 'signal_sd' in wig_types:
        write_signal_sd_wig(raw_read_coverage1, chrm_sizes,
                            wig_base, group1_name)
    if 'dwell' in wig_types:
        write_length_wig(raw_read_coverage1, chrm_sizes,
                         wig_base, group1_name)
    if any(wig_type in wig_types for wig_type in (
            'fraction', 'dampened_fraction')):
        if VERBOSE: sys.stderr.write('Loading statistics from file.\n')
        all_stats, stat_type = ts.parse_stats(stats_fn)
        write_frac_wigs(all_stats, wig_base, 'fraction' in wig_types,
                        'dampened_fraction' in wig_types, cov_damp_counts)

    return

def write_most_signif(
        f5_dirs, fasta_fn, num_regions, corrected_group,
        basecall_subgroups, seqs_fn, num_bases, stats_fn):
    if VERBOSE: sys.stderr.write('Loading statistics from file.\n')
    all_stats, stat_type = ts.parse_stats(stats_fn)
    plot_intervals = ts.get_most_signif_regions(
        all_stats, num_bases, num_regions)

    # get each regions sequence either from reads or fasta index
    if fasta_fn is None:
        raw_read_coverage = th.parse_fast5s(
            f5_dirs, corrected_group, basecall_subgroups)
        all_reg_data = th.get_region_sequences(plot_intervals, raw_read_coverage)
    else:
        genome_index = th.Fasta(fasta_fn)
        all_reg_data = [
            int_i._replace(
                seq=genome_index.get_seq(int_i.chrm, int_i.start, int_i.end))
            for int_i in plot_intervals if int_i.chrm in genome_index.index]

    if VERBOSE: sys.stderr.write('Outputting region seqeuences.\n')
    with io.open(seqs_fn, 'wt') as seqs_fp:
        for int_i in all_reg_data:
            reg_seq = int_i.seq
            if int_i.strand == '-':
                reg_seq = th.rev_comp(reg_seq)
            seqs_fp.write('>{0}:{1:d}:{2} {3}\n{4}\n'.format(
                int_i.chrm, int(int_i.start + (num_bases // 2)),
                int_i.strand, int_i.reg_text, ''.join(reg_seq)))

    return


def wiggle_main(args):
    global VERBOSE
    VERBOSE = not args.quiet
    th.VERBOSE = VERBOSE
    ts.VERBOSE = VERBOSE

    if (any(data_type in args.wiggle_types
            for data_type in ['signal', 'difference', 'coverage',
                              'signal_sd', 'dwell']) and
        args.fast5_basedirs is None):
        th._error_message_and_exit(
            'Must provide a fast5 basedir to output signal, difference, ' +
            'coverage, signal_sd and/or length wiggle files.')
    if (any(wig_type in args.wiggle_types for wig_type in (
            'fraction', 'dampened_fraction')) and
        args.statistics_filename is None):
        th._error_message_and_exit(
            'Must provide a statistics filename to output ' +
            'fraction wiggle files.')
    if ('difference' in args.wiggle_types and
        args.control_fast5_basedirs is None):
        th._error_message_and_exit(
            'Must provide two sets of FAST5s ' + \
            'to output difference wiggle files.')
    if (args.control_fast5_basedirs is not None and
        args.fast5_basedirs is None):
        th._error_message_and_exit(
            'Cannot provide a control FAST5 set of directories ' +
            'without a sample set of FAST5 directories.')

    write_all_wiggles(
        args.fast5_basedirs, args.control_fast5_basedirs, args.corrected_group,
        args.basecall_subgroups, args.statistics_filename, args.wiggle_basename,
        args.wiggle_types, args.coverage_dampen_counts)

    return

def write_signif_diff_main(args):
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
        args.num_bases, args.statistics_filename)

    return


if __name__ == '__main__':
    raise NotImplementedError(
        'This is a module. See commands with `tombo -h`')
