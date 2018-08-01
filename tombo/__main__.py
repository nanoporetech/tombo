from __future__ import unicode_literals, absolute_import

import sys

from . import _option_parsers
from ._version import TOMBO_VERSION

import argparse
class SubcommandHelpFormatter(argparse.RawDescriptionHelpFormatter):
    def _format_action(self, action):
        parts = super(SubcommandHelpFormatter, self)._format_action(action)
        if action.nargs == argparse.PARSER:
            parts = "\n".join(parts.split("\n")[1:])
        return parts

def main(args=None):
    """The main routine."""
    if args is None:
        args = sys.argv[1:]

    # seperate re-squiggle command since all others are nested
    rsqgl_help = [
        ('resquiggle',
         'Re-annotate raw signal with genomic alignment from ' +
         'existing basecalls.', _option_parsers.get_resquiggle_parser()),]
    nested_commands = [
        ('preprocess', 'Pre-process nanopore reads for Tombo processing.', [
            ('annotate_raw_with_fastqs','Add basecalled sequence ' +
             'from FASTQs to raw FAST5s.',
             _option_parsers.get_add_fastqs_parser()),
        ]),
        ('filter', 'Apply filter to Tombo index file for specified criterion.', [
            ('clear_filters',
             'Clear filters to process all successfully re-squiggled reads.',
             _option_parsers.get_clear_filters_parser()),
            ('genome_locations',
             'Filter reads based on mapping location.',
             _option_parsers.get_filter_genome_pos_parser()),
            ('raw_signal_matching',
             'Filter reads with poor raw to expected signal matching.',
             _option_parsers.get_filter_signal_matching_parser()),
            ('q_score',
             'Filter reads with poor mean basecalling quality.',
             _option_parsers.get_filter_qscore_parser()),
            ('level_coverage',
             'Filter reads for more even coverage.',
             _option_parsers.get_filter_coverage_parser()),
            ('stuck',
             'Filter reads with more "stuck" bases.',
             _option_parsers.get_filter_stuck_parser()),
        ]),
        ('detect_modifications', 'Perform statistical testing to detect ' +
         'non-standard nucleotides.', [
             ('de_novo', 'Test for shifts in raw signal against a ' +
              'canonical base model.',
              _option_parsers.get_de_novo_test_signif_parser()),
             ('alternative_model', 'Test for shifts in raw signal which match ' +
              'those of a specific known non-canonical base.',
              _option_parsers.get_alt_test_signif_parser()),
             ('sample_compare', 'Test for shifts in raw signal against signal ' +
              'levels derived from a canonical base only sample (PCR/IVT).',
              _option_parsers.get_samp_comp_test_signif_parser()),
             ('aggregate_per_read_stats','Aggregate Tombo per-read statistics ' +
              'to produce a genomic base statistics file.',
              _option_parsers.get_aggregate_per_read_parser()),
         ]),
        ('text_output', 'Output Tombo results in text files.', [
            ('browser_files', 'Write text outputs for genome browser ' +
             'visualization and bioinformatic processing (wiggle or ' +
             'bedGraph file format).',
             _option_parsers.get_browser_files_parser()),
            ('signif_sequence_context',
             'Write genomic/transcriptomic sequence centered on most ' +
             'modified genomic locations.',
             _option_parsers.get_write_signif_diff_parser()),
        ]),
        ('plot', 'Save plots to visualize raw nanopore signal or ' +
         'testing results.', [
            ('max_coverage',
             'Plot raw signal in regions with maximal coverage.',
             _option_parsers.get_max_cov_parser()),
            ('genome_locations',
             'Plot raw signal at defined genomic locations.',
             _option_parsers.get_genome_loc_parser()),
            ('motif_centered',
             'Plot raw signal at a specific motif.',
             _option_parsers.get_motif_loc_parser()),
            ('max_difference',
             'Plot raw signal where signal differs most between two ' +
             'read groups.', _option_parsers.get_max_diff_parser()),
            ('most_significant',
             'Plot raw signal at most modified locations.',
             _option_parsers.get_signif_diff_parser()),
            ('motif_with_stats',
             'Plot example signal and statistic distributions around a ' +
             'motif of interst.', _option_parsers.get_signif_motif_parser()),
            ('per_read',
             'Plot per-read modified base probabilities.',
             _option_parsers.get_per_read_parser()),
            ('roc','Plot ROC curve from known motif(s).',
             _option_parsers.get_roc_parser()),
            ('per_read_roc','Plot per-read ROC curve from known motif(s).',
             _option_parsers.get_per_read_roc_parser()),
            ('kmer','Plot signal distributions acorss kmers.',
             _option_parsers.get_kmer_dist_parser()),
            ('cluster_most_significant',
             'Clustering traces at bases with most significant stats.',
             _option_parsers.get_cluster_signif_diff_parser()),
        ]),
        ('build_model', 'Create canonical and alternative base Tombo models.', [
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

    desc = ('Tombo command groups (additional help available ' +
            'within each command group):\n' + '\n'.join([
                '\t{0: <25}{1}'.format(grp_name, grp_help)
                for grp_name, grp_help, _ in rsqgl_help + nested_commands]))
    parser = argparse.ArgumentParser(
        prog='tombo',
        description='********** Tombo *********\n\nTombo is a suite of tools ' +
        'primarily for the identification of modified nucleotides from ' +
        'nanopore sequencing data.\n\nTombo also provides tools for the ' +
        'analysis and visualization of raw nanopore signal.\n\n' + desc,
        formatter_class=SubcommandHelpFormatter)
    parser.add_argument(
        '-v', '--version', action='version',
        version='Tombo version: {}'.format(TOMBO_VERSION),
        help='show Tombo version and exit.')

    # Tombo command groups
    service_subparsers = parser.add_subparsers(dest="service_command")

    # seperate re-squiggle command since all others are nested
    rsqgl_parser = service_subparsers.add_parser(
        rsqgl_help[0][0], parents=[rsqgl_help[0][2],],
        add_help=False)
    # resquiggle is both the service parser and action parser
    rsqgl_parser.set_defaults(action_command=rsqgl_help[0][0])

    for grp_name, grp_help, grp_sub_cmds in nested_commands:
        grp_desc = '\n'.join([
            '\t{0: <30}{1}'.format(cmd_name, cmd_help)
            for cmd_name, cmd_help, _ in grp_sub_cmds])
        grp_parser = service_subparsers.add_parser(
            grp_name, formatter_class=SubcommandHelpFormatter,
            description=grp_desc)
        grp_subparser = grp_parser.add_subparsers(
            title=grp_name, dest="action_command")
        for cmd_name, cmd_help, cmd_parser in grp_sub_cmds:
            subparser_cmd = grp_subparser.add_parser(
                cmd_name, parents=[cmd_parser,], add_help=False)

    try:
        save_args = args
        args = parser.parse_args(args)
    except:
        # catch for re-squiggle advanced args printing
        import re
        if any(re.match(rsqgl_help[0][0], val) for val in args) and any(
                re.match(_option_parsers.printadv_opt[0], val)
                for val in args):
            args.extend(['foo', 'foo'])
            args = parser.parse_args(args)
        else:
            raise

    if args.service_command is None:
        parser.print_help()
        sys.stderr.write('\nTombo error: Must provide a tombo command group.\n')
        sys.exit(2)

    # if no second level parser is provided print that command groups help
    if args.action_command is None:
        save_args.append('-h')
        parser.parse_args(save_args)

    if args.action_command == 'resquiggle':
        from . import resquiggle
        resquiggle._resquiggle_main(args)

    elif args.action_command == 'annotate_raw_with_fastqs':
        from . import _preprocess
        _preprocess.annotate_reads_with_fastq_main(args)

    elif args.service_command == 'detect_modifications':
        from . import tombo_stats
        if args.action_command == 'aggregate_per_read_stats':
            tombo_stats._aggregate_per_read_main(args)
        else:
            tombo_stats._test_shifts_main(args)

    elif args.action_command == 'event_resquiggle':
        from . import _event_resquiggle
        _event_resquiggle._event_resquiggle_main(args)
    elif args.service_command == 'build_model':
        from . import tombo_stats
        if args.action_command == 'estimate_reference':
            tombo_stats._est_ref_main(args)
        elif args.action_command == 'estimate_alt_reference':
            tombo_stats._est_alt_ref_main(args)
        elif args.action_command == 'estimate_scale':
            tombo_stats._estimate_scale_main(args)
        else:
            from . import tombo_helper
            tombo_helper._error_message_and_exit(
                'Invalid Tombo build_model command.')

    elif args.service_command == 'filter':
        from . import _filter_reads
        _filter_reads.filter_main(args)

    elif args.service_command == 'text_output':
        from . import _text_output_commands
        if args.action_command == 'browser_files':
            _text_output_commands._browser_files_main(args)
        elif args.action_command == 'signif_sequence_context':
            _text_output_commands._write_signif_diff_main(args)
        else:
            from . import tombo_helper
            tombo_helper._error_message_and_exit(
                'Invalid Tombo text_output command.')

    elif args.service_command == 'plot':
        from . import _plot_commands
        _plot_commands.plot_main(args)

    else:
        from . import tombo_helper
        tombo_helper._error_message_and_exit('Invalid Tombo command.')

    return

if __name__ == "__main__":
    main()
