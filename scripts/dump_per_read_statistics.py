#! /usr/bin/env python
"""Dump per-read statistics

This script takes a per-read statistics file and dumps its contents out into a tab-separated values file. The columns in the file are 'chrm', 'pos', 'strand', 'read_id' and 'stat'.
"""
from tombo import tombo_stats
from tombo.tombo_helper import intervalData
import argparse
import sys
import os

def extract_per_read_stats(input_file, output_file):
    """Dump per-read statistics to tab-separated values"""
    if not os.path.isfile(input_file):
        sys.exit('"{}" is not a valid file'.format(input_file))

    pr_stats = tombo_stats.PerReadStats(input_file)

    with open(output_file, 'w') as out_fp:
        out_fp.write('{}\t{}\t{}\t{}\t{}\n'.format(
            'chrm', 'pos', 'strand', 'read_id', 'stat'))
        for (chrm, strand), cs_blocks in pr_stats.blocks_index.items():
            for start, block_name in cs_blocks.items():
                for pos, stat, read_id in pr_stats.get_region_per_read_stats(
                        intervalData(chrm, start, start + pr_stats.region_size,
                                    strand)):
                    out_fp.write('{}\t{}\t{}\t{}\t{}\n'.format(
                        chrm, pos, strand, read_id, stat))


def main(input_file, output_file):
    extract_per_read_stats(input_file, output_file)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        'input_file',
        type=str,
        help="the per-read statistics file produced by tombo"
    )
    parser.add_argument(
        '-o',
        '--output_file',
        default='per_read_stats.txt',
        type=str,
        help="the name of the output tsv file (default: 'per_read_stats.txt')"
    )
    args = parser.parse_args()
    main(args.input_file, args.output_file)
