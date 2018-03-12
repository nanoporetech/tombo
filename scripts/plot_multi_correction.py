from tombo.plot_commands import plot_multi_corrections

from rpy2.robjects.packages import importr
ggplot = importr(str('ggplot2'))

f5_dirs = ['test_data/native_reads',]
corrected_group = 'RawEventCorrected'
basecall_subgroups = ['BaseCalled_template',]
pdf_fn = 'tombo_results.multi_corrected.pdf'

# either provide list of genome locations or num_regions
genome_locations = None
num_regions = 10

num_reads_per_plot = 10
num_obs = 500
include_orig_bcs = False

plot_multi_corrections(
    f5_dirs, corrected_group, basecall_subgroups, pdf_fn,
    num_reads_per_plot, num_regions, num_obs, include_orig_bcs,
    genome_locations)
