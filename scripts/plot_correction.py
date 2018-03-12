from tombo.plot_commands import plot_corrections

from rpy2.robjects.packages import importr
ggplot = importr(str('ggplot2'))

f5_dirs = ['test_data/native_reads',]
corrected_group = 'RawEventCorrected'
basecall_subgroups = ['BaseCalled_template',]
pdf_fn = 'tombo_results.corrected.pdf'

reg_type = 'random' # or start or end
num_obs = 500
num_reads = 10

plot_corrections(
    f5_dirs, corrected_group, basecall_subgroups, pdf_fn,
    reg_type, num_obs, num_reads)
