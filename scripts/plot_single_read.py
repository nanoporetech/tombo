from tombo import tombo_helper as th
from tombo._plot_commands import plot_single_read, test_r_imports

test_r_imports()

# plot single read
read_fn = 'path/to/read.fast5'
plot_single_read(fast5_fn=read_fn)


# plot first 10k raw obs from several reads
num_reads = 10
reads_index = th.TomboReads(['path/to/fast5/basedir/',])
for i, r_data in enumerate(reads_index.iter_reads()):
    plot_single_read(
        fast5_fn=r_data.fn, num_obs=10000,
        png_fn='single_read_raw_signal.num_' + str(i) + '.png')
    if i > num_reads: break


# plot best reads
best_reads = sorted(
    (r_data.sig_match_score, r_data)
    for r_data in reads_index.iter_reads())[:num_reads]
for i, (r_score, r_data) in enumerate(best_reads):
    plot_single_read(
        fast5_fn=r_data.fn, num_obs=10000,
        png_fn='best_reads.num_' + str(i) + '.png')


# plot worst scoring reads
reads_index = th.TomboReads(['path/to/fast5/basedir/',], remove_filtered=False)
worst_reads = sorted(
    (r_data.sig_match_score, r_data)
    for r_data in reads_index.iter_reads())[::-1][:num_reads]
for i, (r_score, r_data) in enumerate(worst_reads):
    plot_single_read(
        fast5_fn=r_data.fn, num_obs=10000,
        png_fn='worst_reads.num_' + str(i) + '.png')
