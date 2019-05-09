#! /usr/bin/env python
from tombo import tombo_stats
from tombo.tombo_helper import intervalData
import sys
import os

def main():
    if len(sys.argv) < 2:
        sys.exit('missing argument for filename')
    fname = sys.argv[1]
    if not os.path.isfile(fname):
        sys.exit('"{}" is not a valid file'.format(fname))

    pr_stats = tombo_stats.PerReadStats(fname)

    with open('per_read_stats.txt', 'w') as out_fp:
        out_fp.write('{}\t{}\t{}\t{}\t{}\n'.format(
            'chrm', 'pos', 'strand', 'read_id', 'stat'))
        for (chrm, strand), cs_blocks in pr_stats.blocks_index.items():
            for start, block_name in cs_blocks.items():
                for pos, stat, read_id in pr_stats.get_region_per_read_stats(
                        intervalData(chrm, start, start + pr_stats.region_size,
                                    strand)):
                    out_fp.write('{}\t{}\t{}\t{}\t{}\n'.format(
                        chrm, pos, strand, read_id, stat))

if __name__ == '__main__':
    main()
