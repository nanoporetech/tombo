from tombo import tombo_stats
from tombo.tombo_helper import intervalData

pr_stats = tombo_stats.PerReadStats('path/to/stats.per_read_stats')

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
