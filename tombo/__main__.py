from __future__ import unicode_literals, absolute_import

import sys

from . import _option_parsers
from ._version import TOMBO_VERSION

def main(args=None):
    """The main routine."""
    if args is None:
        args = sys.argv[1:]

    commands = [
        ('Pre-processing:', [
            ('annotate_raw_with_fastqs','Add basecalled sequence ' +
             'from FASTQs to raw FAST5s.',
             _option_parsers.get_add_fastqs_parser()),
        ]),
        ('Re-squiggle:', [
            ('resquiggle','Re-annotate raw signal with ' +
             'genomic alignment from existing basecalls.',
             _option_parsers.get_resquiggle_parser()),
        ]),
        ('Modified Base Detection:',[
            ('test_significance','Test for shifts in signal ' +
             'indicative of non-canonical bases.',
             _option_parsers.get_test_signif_parser()),
            ('aggregate_per_read_stats','Aggregate per-read statistics ' +
             'to produce a genomic base statistics file.',
             _option_parsers.get_aggregate_per_read_parser()),
        ]),
        ('Text Output Commands:', [
            ('write_wiggles','Write text outputs for genome browser ' +
             'visualization and bioinformatic processing (wiggle file format).',
             _option_parsers.get_wiggle_parser()),
            ('write_most_significant_fasta',
             'Write sequence centered on most modified genomic locations.',
             _option_parsers.get_write_signif_diff_parser()),
        ]),
        ('Genome Anchored Plotting Commands:', [
            ('plot_max_coverage',
             'Plot raw signal in regions with maximal coverage.',
             _option_parsers.get_max_cov_parser()),
            ('plot_genome_location',
             'Plot raw signal at defined genomic locations.',
             _option_parsers.get_genome_loc_parser()),
            ('plot_motif_centered',
             'Plot raw signal at a specific motif.',
             _option_parsers.get_motif_loc_parser()),
            ('plot_max_difference',
             'Plot raw signal where signal differs most between two read groups.',
             _option_parsers.get_max_diff_parser()),
            ('plot_most_significant',
             'Plot raw signal at most modified locations.',
             _option_parsers.get_signif_diff_parser()),
            ('plot_motif_with_stats',
             'Plot example signal and statistic distributions around a ' +
             'motif of interst.',
             _option_parsers.get_signif_motif_parser()),
            ('plot_per_read',
             'Plot per read modified base probabilities.',
             _option_parsers.get_per_read_parser()),
        ]),
        ('Other Plotting Commands:', [
            ('plot_roc','Plot ROC curve from known motif(s).',
             _option_parsers.get_roc_parser()),
            ('plot_per_read_roc','Plot per-read ROC curve from known motif(s).',
             _option_parsers.get_per_read_roc_parser()),
            ('plot_kmer','Plot signal distributions acorss kmers.',
             _option_parsers.get_kmer_dist_parser()),
            ('cluster_most_significant',
             'Clustering traces at bases with most significant stats.',
             _option_parsers.get_cluster_signif_diff_parser()),
        ]),
        ('Read Filtering (Only effects tombo index file):', [
            ('clear_filters',
             'Clear filters to process all successfully re-squiggled reads.',
             _option_parsers.get_clear_filters_parser()),
            ('filter_stuck',
             'Apply filter based on observations per base thresholds.',
             _option_parsers.get_filter_stuck_parser()),
            ('filter_coverage',
             'Apply filter to downsample for more even coverage.',
             _option_parsers.get_filter_coverage_parser()),
        ]),
        ('Model Estimation and Event-based Re-squiggle:', [
            ('estimate_reference',
             'Estimate reference tombo model derived from the provided reads.',
             _option_parsers.get_est_ref_parser()),
            ('estimate_alt_reference', 'Estimate alternative tombo model from ' +
             'a sample containing canonical bases spiked with a single ' +
             'non-standard base.',
             _option_parsers.get_est_alt_ref_parser()),
            ('estimate_scale', 'Estimate a global scaling parameter from a ' +
             'sub-set of reads.',
             _option_parsers.get_estimate_scale_parser()),
            ('event_resquiggle', 'Re-annotate raw signal with genomic ' +
             'alignment from existing basecalls using event table.',
             _option_parsers.get_event_resquiggle_parser()),
        ]),
    ]
    desc = '\n\n'.join([
        grp + '\n' + '\n'.join([
            '\t{0: <30}{1}'.format(cmd, cmd_help)
            for cmd, cmd_help, cmd_parser in cmds])
        for grp, cmds in commands])

    import argparse
    parser = argparse.ArgumentParser(
        prog='tombo',
        description='********** TOMBO *********\n\nTombo is a suite of tools ' +
        'primarily for the identification of modified nucleotides from ' +
        'nanopore sequencing data.\n\nTombo also provides tools for the ' +
        'analysis and visualization of raw nanopore signal.',
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument(
        '-v', '--version', action='version',
        version='tombo version: {}'.format(TOMBO_VERSION),
        help='show tombo version and exit.')
    subparsers = parser.add_subparsers(
        title='commands', description=desc,
        help='Additional help available for subcommands.')

    # fill subparser with parsers and linked main functions
    for grp, cmds in commands:
        for cmd, cmd_help, cmd_parser in cmds:
            subparser_cmd = subparsers.add_parser(
                cmd, parents=[cmd_parser,], add_help=False)
            subparser_cmd.set_defaults(subcmd=cmd, group=grp)

    args = parser.parse_args(args)

    if args.subcmd == 'resquiggle':
        from . import resquiggle
        resquiggle.resquiggle_main(args)
    elif args.subcmd == 'event_resquiggle':
        from . import _event_resquiggle
        _event_resquiggle.event_resquiggle_main(args)
    elif args.subcmd == 'test_significance':
        from . import tombo_stats
        tombo_stats.test_shifts_main(args)
    elif args.subcmd == 'aggregate_per_read_stats':
        from . import tombo_stats
        tombo_stats.aggregate_per_read_main(args)
    elif args.subcmd == 'estimate_reference':
        from . import tombo_stats
        tombo_stats.est_ref_main(args)
    elif args.subcmd == 'estimate_alt_reference':
        from . import tombo_stats
        tombo_stats.est_alt_ref_main(args)
    elif args.subcmd == 'estimate_scale':
        from . import tombo_stats
        tombo_stats.estimate_scale_main(args)
    elif args.subcmd == 'annotate_raw_with_fastqs':
        from . import tombo_helper
        tombo_helper.annotate_reads_with_fastq_main(args)
    elif args.subcmd == 'clear_filters':
        from . import tombo_helper
        tombo_helper.clear_filters_main(args)
    elif args.subcmd == 'filter_stuck':
        from . import tombo_helper
        tombo_helper.filter_stuck_main(args)
    elif args.subcmd == 'filter_coverage':
        from . import tombo_helper
        tombo_helper.filter_coverage_main(args)
    elif args.group == 'Text Output Commands:':
        from . import text_output_commands
        if args.subcmd == 'write_wiggles':
            text_output_commands.wiggle_main(args)
        else:
            text_output_commands.write_signif_diff_main(args)
    else:
        from . import plot_commands
        plot_commands.plot_main(args)

    return

if __name__ == "__main__":
    main()
