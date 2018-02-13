from __future__ import unicode_literals, absolute_import

from builtins import map

import sys
import argparse

if sys.version_info[0] > 2:
    unicode = str

from ._default_parameters import SEG_PARAMS_TABLE, ALGN_PARAMS_TABLE, \
    LLR_THRESH, HYPO_THRESH, ALTERNATE_MODELS

ALT_BASES = tuple(set(alt_name.split('_')[1] for alt_name in ALTERNATE_MODELS))


##################################
###### Positional arguments ######
##################################

basedir_opt=('fast5_basedir', {
    'type':unicode,
    'help':'Directory containing fast5 files. All files ending in "fast5" ' +
    'found recursively within this base directory will be processed.'})
fasta_pos_opt=(
    'genome_fasta', {'type':unicode, 'help':'Path to fasta file for mapping.'})


############################
###### Text arguments ######
############################

minimap2_opt=('--minimap2-executable', {
    'type':unicode, 'help':'Path to minimap2 executable.'})
minindx_opt=('--minimap2-index', {
    'type':unicode, 'help':'Path to minimap2 index (with map-ont preset) '
    'file corresponding to the [genome_fasta] provided.'})
bwamem_opt=('--bwa-mem-executable', {
    'type':unicode, 'help':'Path to bwa-mem executable.'})
graphmap_opt=('--graphmap-executable', {
    'type':unicode, 'help':'Path to graphmap executable.'})

poremod_opt=('--pore-model-filename', {
    'type':unicode,
    'help':'File containing kmer model parameters (level_mean ' +
    'and level_stdv) used in order to compute kmer-based corrected pA ' +
    'values. E.g. https://github.com/jts/nanopolish/blob/master/etc/' +
    'r9-models/template_median68pA.5mers.model'})
tbmod_opt=('--tombo-model-filename', {
    'type':unicode, 'help':'Tombo model filename. If no file is provided, ' +
    'the default DNA or RNA Tombo model will be used.'})
tbmod_w_opt=('--tombo-model-filename', {
    'type':unicode, 'help':'Filename to save Tombo model.'})
atbmod_opt=('--alternate-model-filename', {
    'type':unicode,
    'help':'Tombo model for alternative likelihood ratio significance testing.'})
hidden_tbmod_opt=('--tombo-model-filename', {
    'type':unicode, 'help':argparse.SUPPRESS})
hidden_atbmod_opt=('--alternate-model-filename', {
    'type':unicode, 'help':argparse.SUPPRESS})
hidden_atbmods_opt=('--alternate-model-filenames', {
    'type':unicode, 'nargs':'+', 'help':argparse.SUPPRESS})
altname_opt=('--alternate-model-name', {
    'type':unicode,
    'help':'A short name to associate with this alternate model (e.g. 5mC, ' +
    '6mA, etc.). This text will be included in output filenames when this ' +
    'model is used for testing.'})

failed_opt=('--failed-reads-filename', {
    'help':'Output failed read filenames with assoicated error. Default: ' +
    'Do not store failed reads.'})

sfast5dir_opt = ('--fast5-basedir', {
    'type':unicode, 'help':'Directory containing fast5 files.'})
fast5dir_opt = ('--fast5-basedirs', {
    'type':unicode, 'nargs':'+', 'help':'Directories containing fast5 files.'})
ctrlfast5dir_opt=('--control-fast5-basedirs', {
    'type':unicode, 'nargs':'+',
    'help':'Set of directories containing fast5 files for control reads, ' +
    'containing only canonical nucleotides.'})

corrgrp_opt=('--corrected-group', {
    'type':unicode, 'default':'RawGenomeCorrected_000',
    'help':'FAST5 group created by resquiggle command. Default: %(default)s'})
correvntgrp_opt=('--corrected-group', {
    'type':unicode, 'default':'RawGenomeCorrected_000',
    'help':'FAST5 group created by resquiggle command. Default: %(default)s'})
newcorrgrp_opt=('--new-corrected-group', {
    'type':unicode, 'default':'RawModelCorrected_000',
    'help':'FAST5 group created by resquiggle command. Default: %(default)s'})
bcgrp_opt=('--basecall-group', {
    'type':unicode, 'default':'Basecall_1D_000',
    'help':'FAST5 group obtain original basecalls (under Analyses group). ' +
    'Default: %(default)s'})
bcsubgrps_opt=('--basecall-subgroups', {
    'type':unicode, 'default':['BaseCalled_template',], 'nargs':'+',
    'help':'FAST5 subgroup(s) (under /Analyses/[--basecall-group]/) containing ' +
    'basecalls and created within [--corrected-group] containing re-squiggle ' +
    'results. Default: %(default)s'})
bcsubgrp_opt=('--basecall-subgroup', {
    'type':unicode, 'default':'BaseCalled_template',
    'help':'FAST5 subgroup (under /Analyses/[--basecall-group]/) under which ' +
    'to store basecalls from FASTQs. Default: %(default)s'})

gnmloc_opt=('--genome-locations', {
    'type':unicode, 'nargs':'+',
    'help':'Genomic locations at which to plot signal. Format locations ' +
    'as "chrm:position[:strand] [chrm2:position2[:strand2] ...]" ' +
    '(strand not applicable for all applications)'})
fasta_opt=('--genome-fasta', {
    'type':unicode,
    'help':'FASTA file used to re-squiggle. For faster sequence access.'})
motif_opt=('--motif', {
    'type':unicode,
    'help':'Motif of interest at which to plot signal and statsitics. ' +
    'Supports IUPAC single letter codes (use T for RNA).'})

obsfilt_opt=('--obs-per-base-filter', {
    'type':unicode, 'nargs':'+', 'default':[],
    'help':'Filter reads baseed on observations per base percentile ' +
    'thresholds. Format thresholds as "percentile:thresh ' +
    '[pctl2:thresh2 ...]". For example to filter reads with 99th ' +
    'pctl > 200 obs/base or max > 5k obs/base use "99:200 100:5000".'})

fastqs_opt = ('--fastq-filenames', {
    'type':unicode, 'nargs':'+',
    'help':'FASTQ filenames containing basecalls to be added to ' +
    'raw FAST5 files.'})

wigfn_opt=('--wiggle-basename', {
    'type':unicode, 'default':'tombo_results',
    'help':'Basename for output wiggle files. Two files (plus and minus ' +
    'strand) will be produced for each --wiggle-types supplied. ' +
    'Filenames formatted as "[wiggle-basename].[wiggle-type].' +
    '[sample|control]?.[plus|minus].wig". Default: %(default)s'})
pdf_opt=('--pdf-filename', {
    'type':unicode, 'help':'PDF filename to store plot(s). Default: %(default)s'})
statfn_opt=('--statistics-filename', {
    'type':unicode, 'help':"File to save/load base by base statistics."})
statbsnm_opt=('--statistics-file-basename', {
    'type':unicode,
    'help':"File base name to save base by base statistics from testing. " +
    "Filenames will be [--statistics-file-basename]." +
    "[--alternate-bases]?.tombo.stats"})
rdata_opt=('--r-data-filename', {
    'type':unicode,
    'help':"Filename to save R data structure. Default: Don't save"})
seqs_opt=('--sequences-filename', {
    'type':unicode,
    'help':'File for sequences from selected regions. Sequences will be ' +
    'stored in FASTA format. Default: %(default)s.'})
densbn_opt=('--save-density-basename', {
    'type':unicode,
    'help':"Basename to save alternative model estimation density " +
    "estimation information. See scripts/debug_est_alt.R for info use " +
    "example. Default: Don't save."})
altden_opt=('--alternate-density-filename', {
    'type':unicode,
    'help':'File containing k-mer level kernel density estimates for the ' +
    'alternative sample saved using --save-density-basename.'})
ctrlden_opt=('--control-density-filename', {
    'type':unicode,
    'help':'File containing k-mer level kernel density estimates for the ' +
    'control sample saved using --save-density-basename.'})
prstatbn_opt=('--per-read-statistics-basename', {
    'type':unicode,
    'help':'Base for binary files containing per-read statistics from ' +
    'statistical testing. Filenames will be [--per-read-statistics-basename].' +
    '[--alternate-bases]?.tombo.per_read_stats'})
prstat_opt=('--per-read-statistics-filename', {
    'type':unicode,
    'help':'Binary file containing per-read statistics from ' +
    'statistical testing.'})

statfns_opt=('--statistics-filenames', {
    'type':unicode, 'nargs':'+',
    'help':"Files to load base by base statistics."})
motifdesc_opt=('--motif-descriptions', {
    'type':unicode, 'nargs':'+',
    'help':'Ground truth, motif centered, modified base descriptions for ' +
    'computing ROC and PR curves. Each statistics file is associated with ' +
    'a set of motif descriptions. Format descriptions as: "motif:mod_pos:name' +
    '[::motif2:mod_pos2:name2...]". The mod_pos indicated the modified base ' +
    'within the motif (1-based index). Example: CCWGG:2:"dcm 5mC"::GATC:2:' +
    '"dam 6mA" would assess the performance of a single Tombo statistics ' +
    'file for identification of E. coli dam and dcm methylation.'})


############################
###### Int arguments ######
############################

proc_opt=('--processes', {
    'type':int, 'help':'Number of processes. Default: %(default)d'})
alignproc_opt=('--align-processes', {
    'type':int, 'default':1,
    'help':'Number of processes to use for parsing and aligning ' +
    'original basecalls. Each process will independently load the ' +
    'genome into memory, so use caution with larger genomes ' +
    '(e.g. human). Default: %(default)d'})
alignthrds_opt=('--align-threads-per-process', {
    'type':int,
    'help':'Number of threads to use for aligner system call. ' +
    'Default: [--processes] / (2 * [--align-processes)]'})
rsqglproc_opt=('--resquiggle-processes', {
    'type':int,
    'help':'Number of processes to use for resquiggle algorithm. ' +
    'Default: [--processes] / 2'})
batchsize_opt=('--alignment-batch-size', {
    'default':1000, 'type':int,
    'help':'Number of reads included in each alignment call. ' +
    'Note: A new system mapping call is made for each batch ' +
    '(including loading of the genome), so it is advised to use ' +
    'larger values for larger genomes. Default: %(default)d'})
mpreg_opt=('--multiprocess-region-size', {
    'default':10000, 'type':int,
    'help':'Size of regions over which to multiprocesses statistic ' +
    'computation. For very deep samples a smaller value is recommmended ' +
    'in order to control memory consumption. Default: %(default)d'})

timeout_opt=('--timeout', {
    'type':int, 'help':'Timeout in seconds for processing a single read. ' +
    'Default: No timeout.'})
cpt_opt=('--cpts-limit', {
    'type':int, 'help':'Maximum number of changepoints to find within a ' +
    'single indel group. Default: No limit.'})

kmerthresh_opt=('--num-kmer-threshold', {
    'default':1, 'type':int,
    'help':'Observations of each k-mer required to include a read in ' +
    'read level averages. Default: %(default)d'})

minobs_opt=('--min-obs-per-base', {
    'type':int,
    'help':'Minimum raw observations to assign to a genomic base. ' +
    'Default: %(default)d'})
covthresh_opt=('--coverage-threshold', {
    'type':int,
    'help':'Maximum mean coverage per region when estimating k-mer model ' +
    '(limits compute time for deep samples). Default: %(default)d'})
maxbase_opt=('--max-bases-shift', {
    'type':int, 'default':3,
    'help':'Maximum bases to shift raw signal from event_resquiggle ' +
    'assignment. Default: %(default)d'})
bmaxbase_opt=('--base-score-max-bases-shift', {
    'type':int, 'default':4,
    'help':'Maximum bases to shift raw signal from first round of ' +\
    'model re-squiggle. Default: %(default)d'})
minkmer_opt=('--minimum-kmer-observations', {
    'type':int,
    'help':'Number of each k-mer observations required in order to produce ' +
    'a reference (genomic locations for standard reference and per-read ' +
    'for alternative reference). Default: %(default)d'})

numbases_opt=('--num-bases', {
    'type':int,
    'help':'Number of bases to plot/output. Default: %(default)d'})
numreads_opt=('--num-reads', {
    'type':int, 'help':'Number of reads to plot. Default: %(default)d'})
numobs_opt=('--num-obs', {
    'type':int, 'default':500,
    'help':'Number of observations to plot. Default: %(default)d'})
numreg_opt=('--num-regions', {
    'type':int, 'help':'Number of regions to plot. Default: %(default)d'})

slides_opt=('--slide-span', {
    'type':int, 'default':0,
    'help':'Number of bases offset over which to search when computing ' +
    'distances for signal cluster plotting. Default: 0 (exact position)'})
cntxt_opt=('--num-context', {
    'type':int, 'default':5,
    'help':'Number of context bases around motif. Default: %(default)d'})
numstat_opt=('--num-statistics', {
    'type':int, 'default':200,
    'help':'Number of motif centered regions to include in ' +
    'statistic distributions. Default: %(default)d'})
ovpltthresh_opt=('--overplot-threshold', {
    'type':int, 'default':50,
    'help':'Coverage level to trigger alternative plot type ' +
    'instead of raw signal. Default: %(default)d'})

fmo_opt=('--fishers-method-context', {
    'type':int, 'default':1,
    'help':'Number of context bases up and downstream over which to compute ' +
    "Fisher's method combined p-values. Note: Not applicable " +
    "for alternative model likelihood ratio tests. Default: %(default)d."})
szo_opt=('--stouffer-z-context', {
    'type':int, 'default':1,
    'help':'Number of context bases up and downstream over which to compute ' +
    "Stouffer's Z combined z-scores. Default: %(default)d."})
regcntxt_opt=('--region-context', {
    'type':int, 'default':2,
    'help':'Number of context bases up and downstream of poorly fit ' +
    'regions to perform model re-squiggle. Default: %(default)d.'})
brcntxt_opt=('--base-score-region-context', {
    'type':int, 'default':4,
    'help':'Number of context bases up and downstream of poorly fit regions ' +
    'to perform iterative base-score model ' +
    're-squiggle. Default: %(default)d.'})
bsiters_opt=('--base-score-iterations', {
    'type':int, 'default':2,
    'help':'Number of iterations through each read to perform ' +
    '(computationally expensive) base space model re-squiggle algorithm. ' +
    'Default: %(default)d.'})
minreads_opt=('--minimum-test-reads', {
    'type':int,
    'help':'Number of reads required at a position to perform significance ' +
    'testing or contribute to model estimation. Default: %(default)d'})

segpars_opt=('--segmentation-parameters', {
    'type':int, 'nargs':3,
    'help':'Specify the 3 parameters for segmentation 1) running neighboring ' +
    'windows width 2) minimum raw observations per genomic base 3) mean raw ' +
    'observations per event. Sample type defaults: ' +
    ' || '.join((bst + ' : ' + ' '.join(map(str, params)))
                for bst, params in SEG_PARAMS_TABLE.items())})
segpars2_opt=('--segmentation-parameters', {
    'type':int, 'nargs':2,
    'help':'Specify the 2 parameters for segmentation 1) running neighboring ' +
    'windows width 2) minimum raw observations per genomic base. Sample type ' +
    'defaults:\n' +
    ' || '.join((bst + ' : ' + ' '.join(map(str, params[:2])))
                for bst, params in SEG_PARAMS_TABLE.items())})


###############################
###### Boolean arguments ######
###############################

skpidx_opt=('--skip-index', {
    'default':False, 'action':'store_true',
    'help':'Skip creation of tombo index. This drastically slows downstream '+
    'tombo commands. Default stores tombo index named ".[--fast5-basedir].' +
    '[--corrected-group].tombo.index" to be loaded automatically for ' +
    'downstream commands.'})
ovrwrt_opt=('--overwrite', {
    'default':False, 'action':'store_true',
    'help':'Overwrite previous corrected group in FAST5 files. Note: ' +
    'only effects --corrected-group or --new-corrected-group.'})

estmean_opt=('--estimate-mean', {
    'default':False, 'action':'store_true',
    'help':"Use the mean instead of median for model level estimation. Note:" +
    " This can cause poor fits due to outliers"})
kmspec_opt=('--kmer-specific-sd', {
    'default':False, 'action':'store_true',
    'help':"Estimate standard deviation for each k-mers individually."})
incldsd_opt=('--include-event-stdev', {
    'default':False, 'action':'store_true',
    'help':'Include corrected event standard deviation in output FAST5 data.'})
fitscl_opt=('--fit-scale-per-read', {
    'default':False, 'action':'store_true',
    'help':'Fit the scaling parameter for each read. If not set then a ' +
    'global scaling parameter is estimated from a random subset of reads, ' +
    'which has could provide more robust results.'})

readmean_opt=('--read-mean', {
    'default':False, 'action':'store_true',
    'help':'Plot k-mer means across whole reads as opposed to ' +
    'individual k-mer event levels.'})
boxc_opt=('--box-center', {
    'default':False, 'action':'store_true',
    'help':"Plot a box around the central base."})
deepcov_opt=('--deepest-coverage', {
    'default':False, 'action':'store_true',
    'help':'Plot the deepest coverage regions.'})
noplot_opt=('--dont-plot', {
    'default':False, 'action':'store_true',
    'help':"Don't plot result. Useful to produce only R data file."})
origbcs_opt=('--include-original-basecalls', {
    'default':False, 'action':'store_true',
    'help':"Include original basecalls in plots."})
pstdmod_opt=('--plot-standard-model', {
    'default':False, 'action':'store_true',
    'help':"Add default standard model distribution to the plot."})

quiet_opt=(('--quiet', '-q'), {
    'default':False, 'action':'store_true',
    'help':"Don't print status information."})


##############################
###### Float arguments ######
##############################

otlthresh_opt=('--outlier-threshold', {
    'default':5, 'type':float,
    'help':'Windosrize the signal at this number of scale values. ' +
    'Negative value disables outlier clipping. Default: %(default)f'})
pvalthrsh_opt=('--p-value-threshold', {
    'default':0.1, 'type':float,
    'help':'P-value threshold to identify regions to apply model ' +
    're-squiggle algorithm. Default: %(default)f'})
snglrdthrsh_opt=('--single-read-threshold', {
    'type':float,
    'help':'P-value or log likelihood ratio threshold when computing ' +
    'fraction of significant reads at each genomic position. Default: ' +
    'p-value:{0:.2g}; likelihood ratio:{1:.2g}'.format(HYPO_THRESH, LLR_THRESH)})
altfrac_opt=('--alt-fraction-percentile', {
    'default':1, 'type':float,
    'help':'When esitmating the alternative base incorporation rate, this ' +
    'percent of k-mers are assumed to have significantly shifted signal so ' +
    'the alternative distribution minimally overlaps the standard base ' +
    'distribution. Default: %(default)f'})
kernden_opt=('--kernel-density-bandwidth', {
    'default':0.05, 'type':float,
    'help':'Bandwidth applied when performing Gaussian kernal density ' +
    'esitmation on standard and alternative base signal distributions. ' +
    'Default: %(default)f'})
pctfilt_opt=('--percent-to-filter', {
    'type':float, 'default':10,
    'help':'Percentage of all reads to filter. Reads are randomly selected ' +
    'weighted according to the approximate coverage at the mapped genomic ' +
    'location. This can be useful in modeling and testing.'})
fxdscl_opt=('--fixed-scale', {
    'type':float,
    'help':'Fixed scaling parameter to use for raw signal normalization.'})
cvgdmp_opt=('--coverage-dampen-counts', {
    'type':float, 'nargs':2, 'default':[2,0.5],
    'help':'Dampen fraction modified estimates for low coverage sites. Two ' +
    'parameters are psuedo unmodified and modified read counts. This is ' +
    'equivalent to a beta prior on the fraction estimate. Default: %(default)s'})

sigapars_opt=('--signal-align-parameters', {
    'type':float, 'nargs':4,
    'help':'Specify the 4 parameters for signal to genome sequence alignment ' +
    'algorithm 1) match expected value 2) skip penalty 3) bandwidth 4) mean ' +
    'signal segmentation half-normal score threshold. Sample type defaults: ' +
    ' || '.join((bst + ' : ' + ' '.join(map(str, params)))
                for bst, params in ALGN_PARAMS_TABLE.items())})


##############################
###### Choice arguments ######
##############################

normtype_opt=('--normalization-type', {
    'type':unicode,
    'default':'median', 'choices':('median', 'pA', 'pA_raw', 'none'),
    'help':'Choices: "none": raw 16-bit DAQ values, "pA_raw": pA as in the ' +
    'ONT events (using offset, range and digitization), "pA": k-mer-based ' +
    'correction for pA drift as in nanopolish (requires ' +
    '[--pore-model-filename]), "median": median and MAD from raw signal. ' +
    'Default: %(default)s'})

upstrmbs_opt=('--upstream-bases', {
    'default':1, 'type':int, 'choices':(0,1,2,3,4),
    'help':'Upstream bases in k-mer. Default: %(default)d'})
dnstrmbs_opt=('--downstream-bases', {
    'default':2, 'type':int, 'choices':(0,1,2,3,4),
    'help':'Downstream bases in k-mer. Default: %(default)d'})
altbs_opt=('--alternate-model-base', {
    'type':unicode, 'choices':('A','C','G','T'),
    'help':'Non-standard base is an alternative to this base.'})
modbs_opt=('--alternate-bases', {
    'type':unicode, 'choices':ALT_BASES, 'nargs':'+',
    'help':'Default non-standard base model for testing.'})
paltmod_opt=('--plot-alternate-model', {
    'type':unicode, 'choices':ALT_BASES,
    'help':'Add alternative model distribution to the plot.'})
regtype_opt=('--region-type', {
    'type':unicode, 'default':'random', 'choices':['random', 'start', 'end'],
    'help':'Region to plot within each read. Default: random'})
ovplttype_opt=('--overplot-type', {
    'type':unicode, 'default':'Downsample',
    'choices':['Downsample', 'Boxplot', 'Quantile', 'Density'],
    'help':'Plot type for regions with higher coverage. Default: Downsample'})
wigtypes_opt=('--wiggle-types', {
    'type':unicode, 'default':['coverage', 'fraction'], 'nargs':'+',
    'choices':['coverage', 'fraction', 'dampened_fraction', 'signal',
               'signal_sd', 'dwell', 'difference'],
    'help':'Data types of wiggles to produce. Default: "coverage fraction"'})

dna_opt=('--dna', {
    'dest':'bio_sample_type', 'action':'store_const', 'const':'DNA',
    'help':'Explicitly select default DNA model. Default: Automatically ' +
    'determine from FAST5s'})
rna_opt=('--rna', {
    'dest':'bio_sample_type', 'action':'store_const', 'const':'RNA',
    'help':'Explicitly select default RNA model. Default: Automatically ' +
    'determine from FAST5s'})


###########################
###### Help argument ######
###########################

help_opt=(('--help', '-h'), {
    'action':'help',
    'help':"Print this help message and exit"})


##############################
###### Helper functions ######
##############################

OUTPUT_BASE = 'tombo_results'

def add_misc_args(parser):
    misc_args = parser.add_argument_group('Miscellaneous Arguments')
    misc_args.add_argument(*quiet_opt[0], **quiet_opt[1])
    misc_args.add_argument(*help_opt[0], **help_opt[1])

    return misc_args, parser

def add_default_args(parser):
    fast5_args = parser.add_argument_group('FAST5 Data Arguments')
    fast5_args.add_argument(corrgrp_opt[0], **corrgrp_opt[1])
    fast5_args.add_argument(bcsubgrps_opt[0], **bcsubgrps_opt[1])

    misc_args = add_misc_args(parser)

    return fast5_args, misc_args, parser

def add_comp_dist_args(parser):
    alt_args = parser.add_argument_group('Comparison Arguments')
    alt_args.add_argument(ctrlfast5dir_opt[0], **ctrlfast5dir_opt[1])
    alt_args.add_argument(pstdmod_opt[0], **pstdmod_opt[1])
    alt_args.add_argument(paltmod_opt[0], **paltmod_opt[1])
    alt_args.add_argument(hidden_tbmod_opt[0], **hidden_tbmod_opt[1])
    alt_args.add_argument(hidden_atbmod_opt[0], **hidden_atbmod_opt[1])

    return alt_args, parser


#####################################
###### Main re-squiggle parser ######
#####################################

def get_eventless_resquiggle_parser():
    parser = argparse.ArgumentParser(
        description='Re-segment raw nanopore signal to match with mapped ' +
        'portion of a known genomic sequence guided by a k-mer model.',
        add_help=False)
    req_args = parser.add_argument_group('Required Arguments')
    req_args.add_argument(basedir_opt[0], **basedir_opt[1])
    req_args.add_argument(fasta_pos_opt[0], **fasta_pos_opt[1])

    mod_args = parser.add_argument_group('Model Parameters')
    mod_args.add_argument(dna_opt[0], **dna_opt[1])
    mod_args.add_argument(rna_opt[0], **rna_opt[1])
    mod_args.add_argument(hidden_tbmod_opt[0], **hidden_tbmod_opt[1])

    alg_args = parser.add_argument_group(
        'Event to Sequence Assignment Parameters')
    alg_args.add_argument(segpars_opt[0], **segpars_opt[1])
    alg_args.add_argument(sigapars_opt[0], **sigapars_opt[1])

    sig_args = parser.add_argument_group( 'Signal Scaling Parameters')
    sig_args.add_argument(fitscl_opt[0], **fitscl_opt[1])
    sig_args.add_argument(fxdscl_opt[0], **fxdscl_opt[1])
    sig_args.add_argument(otlthresh_opt[0], **otlthresh_opt[1])

    io_args = parser.add_argument_group('Input/Output Arguments')
    io_args.add_argument(minindx_opt[0], **minindx_opt[1])
    io_args.add_argument(skpidx_opt[0], **skpidx_opt[1])
    io_args.add_argument(failed_opt[0], **failed_opt[1])
    io_args.add_argument(incldsd_opt[0], **incldsd_opt[1])

    filt_args = parser.add_argument_group('Read Filtering Argument')
    filt_args.add_argument(obsfilt_opt[0], **obsfilt_opt[1])

    multi_args = parser.add_argument_group('Multiprocessing Arguments')
    multi_args.add_argument(proc_opt[0], default=1, **proc_opt[1])

    fast5_args = parser.add_argument_group('FAST5 Data Arguments')
    fast5_args.add_argument(corrgrp_opt[0], **corrgrp_opt[1])
    fast5_args.add_argument(bcgrp_opt[0], **bcgrp_opt[1])
    fast5_args.add_argument(bcsubgrps_opt[0], **bcsubgrps_opt[1])
    fast5_args.add_argument(ovrwrt_opt[0], **ovrwrt_opt[1])

    misc_args, parser = add_misc_args(parser)

    return parser


#############################################
###### Alternaitve re-squiggle parsers ######
#############################################

def get_event_resquiggle_parser():
    parser = argparse.ArgumentParser(
        description='Re-segment raw nanopore signal to match with mapped ' +
        'portion of a known genomic sequence guided by stored events ' +
        'from basecaller.', add_help=False)
    req_args = parser.add_argument_group('Required Arguments')
    req_args.add_argument(basedir_opt[0], **basedir_opt[1])
    req_args.add_argument(fasta_pos_opt[0], **fasta_pos_opt[1])

    mapper_args = parser.add_argument_group(
        'Mapper Arguments (One mapper is required)')
    mapper_args.add_argument(minimap2_opt[0], **minimap2_opt[1])
    mapper_args.add_argument(minindx_opt[0], **minindx_opt[1])
    mapper_args.add_argument(bwamem_opt[0], **bwamem_opt[1])
    mapper_args.add_argument(graphmap_opt[0], **graphmap_opt[1])
    mapper_args.add_argument(batchsize_opt[0], **batchsize_opt[1])

    norm_args = parser.add_argument_group('Signal Processing Arguments')
    norm_args.add_argument(normtype_opt[0], **normtype_opt[1])
    norm_args.add_argument(poremod_opt[0], **poremod_opt[1])
    norm_args.add_argument(otlthresh_opt[0], **otlthresh_opt[1])
    norm_args.add_argument(segpars2_opt[0], **segpars2_opt[1])

    filt_args = parser.add_argument_group('Read Filtering Arguments')
    filt_args.add_argument(obsfilt_opt[0], **obsfilt_opt[1])
    filt_args.add_argument(timeout_opt[0], **timeout_opt[1])
    filt_args.add_argument(cpt_opt[0], **cpt_opt[1])

    io_args = parser.add_argument_group('Input/Output Arguments')
    io_args.add_argument(skpidx_opt[0], **skpidx_opt[1])
    io_args.add_argument(ovrwrt_opt[0], **ovrwrt_opt[1])
    io_args.add_argument(failed_opt[0], **failed_opt[1])
    io_args.add_argument(incldsd_opt[0], **incldsd_opt[1])

    fast5_args = parser.add_argument_group('FAST5 Data Arguments')
    fast5_args.add_argument(correvntgrp_opt[0], **correvntgrp_opt[1])
    fast5_args.add_argument(bcgrp_opt[0], **bcgrp_opt[1])
    fast5_args.add_argument(bcsubgrps_opt[0], **bcsubgrps_opt[1])

    multi_args = parser.add_argument_group('Multiprocessing Arguments')
    multi_args.add_argument(proc_opt[0], default=2, **proc_opt[1])
    multi_args.add_argument(alignproc_opt[0], **alignproc_opt[1])
    multi_args.add_argument(alignthrds_opt[0], **alignthrds_opt[1])
    multi_args.add_argument(rsqglproc_opt[0], **rsqglproc_opt[1])

    misc_args, parser = add_misc_args(parser)

    return parser

def get_model_resquiggle_parser():
    parser = argparse.ArgumentParser(
        description='Re-segment raw nanopore signal at local regions that do ' +
        'not match a tombo model.', add_help=False)
    req_args = parser.add_argument_group('Required Argument')
    req_args.add_argument(fast5dir_opt[0], required=True, **fast5dir_opt[1])

    reg_args = parser.add_argument_group('Region Selection Arguments')
    reg_args.add_argument(szo_opt[0], **szo_opt[1])
    reg_args.add_argument(regcntxt_opt[0], **regcntxt_opt[1])
    reg_args.add_argument(pvalthrsh_opt[0], **pvalthrsh_opt[1])

    modr_args = parser.add_argument_group('Model Re-squiggle Arguments')
    modr_args.add_argument(dna_opt[0], **dna_opt[1])
    modr_args.add_argument(rna_opt[0], **rna_opt[1])
    modr_args.add_argument(maxbase_opt[0], **maxbase_opt[1])
    modr_args.add_argument(minobs_opt[0], default=3, **minobs_opt[1])
    modr_args.add_argument(hidden_tbmod_opt[0], **hidden_tbmod_opt[1])

    brsqgl_args = parser.add_argument_group('Base Scoring Arguments')
    brsqgl_args.add_argument(bsiters_opt[0], **bsiters_opt[1])
    brsqgl_args.add_argument(brcntxt_opt[0], **brcntxt_opt[1])
    brsqgl_args.add_argument(bmaxbase_opt[0], **bmaxbase_opt[1])

    fast5_args = parser.add_argument_group('FAST5 Data Arguments')
    fast5_args.add_argument(corrgrp_opt[0], **corrgrp_opt[1])
    fast5_args.add_argument(newcorrgrp_opt[0], **newcorrgrp_opt[1])
    fast5_args.add_argument(bcsubgrps_opt[0], **bcsubgrps_opt[1])

    out_args = parser.add_argument_group('Output Arguments')
    out_args.add_argument(ovrwrt_opt[0], **ovrwrt_opt[1])
    out_args.add_argument(failed_opt[0], **failed_opt[1])
    out_args.add_argument(incldsd_opt[0], **incldsd_opt[1])

    multi_args = parser.add_argument_group('Multiprocessing Argument')
    multi_args.add_argument(proc_opt[0], default=1, **proc_opt[1])

    misc_args, parser = add_misc_args(parser)

    return parser


###################################
###### Pre-processing parser ######
###################################

def get_add_fastqs_parser():
    parser = argparse.ArgumentParser(
        description='Annotate raw FAST5 files with basecalls from a FASTQ.',
        add_help=False)
    req_args = parser.add_argument_group('Required Arguments')
    req_args.add_argument(sfast5dir_opt[0], required=True, **sfast5dir_opt[1])
    req_args.add_argument(fastqs_opt[0], required=True, **fastqs_opt[1])

    fast5_args = parser.add_argument_group('FAST5 Data Arguments')
    fast5_args.add_argument(bcgrp_opt[0], **bcgrp_opt[1])
    fast5_args.add_argument(bcsubgrp_opt[0], **bcsubgrp_opt[1])
    fast5_args.add_argument(ovrwrt_opt[0], **ovrwrt_opt[1])

    misc_args, parser = add_misc_args(parser)

    return parser


######################################
###### Model estimation parsers ######
######################################

def get_est_ref_parser():
    parser = argparse.ArgumentParser(
        description='Estimate standard tombo model for use in re-squiggle ' +
        'and  testing without an amplified (un-modified) sample.', add_help=False)
    req_args = parser.add_argument_group('Required Arguments')
    req_args.add_argument(fast5dir_opt[0], required=True, **fast5dir_opt[1])
    req_args.add_argument(tbmod_w_opt[0], required=True, **tbmod_w_opt[1])

    stat_args = parser.add_argument_group('Modeling Arguments')
    stat_args.add_argument(estmean_opt[0], **estmean_opt[1])
    stat_args.add_argument(kmspec_opt[0], **kmspec_opt[1])
    stat_args.add_argument(upstrmbs_opt[0], **upstrmbs_opt[1])
    stat_args.add_argument(dnstrmbs_opt[0], **dnstrmbs_opt[1])

    filt_args = parser.add_argument_group('Filtering Arguments')
    filt_args.add_argument(minreads_opt[0], default=10, **minreads_opt[1])
    filt_args.add_argument(covthresh_opt[0],  default=100, **covthresh_opt[1])
    filt_args.add_argument(minkmer_opt[0], default=5, **minkmer_opt[1])

    multi_args = parser.add_argument_group('Multiprocessing Arguments')
    multi_args.add_argument(mpreg_opt[0], **mpreg_opt[1])
    multi_args.add_argument(proc_opt[0], default=1, **proc_opt[1])

    fast5_args, misc_args, parser = add_default_args(parser)

    return parser

def get_est_alt_ref_parser():
    parser = argparse.ArgumentParser(
        description='Estimate alternative k-mer reference model for use ' +
        'in testing for specific modification types. [--fast5-basedirs] ' +
        'should contain a sample spiked with a single known randomly ' +
        'incorporated base.', add_help=False)
    req_args = parser.add_argument_group('Required Arguments')
    req_args.add_argument(atbmod_opt[0], required=True, **atbmod_opt[1])
    req_args.add_argument(altname_opt[0], required=True, **altname_opt[1])
    req_args.add_argument(altbs_opt[0], required=True, **altbs_opt[1])

    dens_args = parser.add_argument_group(
        'Signal Data Arguments (Must provide either FAST5 dirs or ' +
        'previous density estimates)')
    dens_args.add_argument(fast5dir_opt[0], **fast5dir_opt[1])
    dens_args.add_argument(ctrlfast5dir_opt[0], **ctrlfast5dir_opt[1])
    dens_args.add_argument(altden_opt[0], **altden_opt[1])
    dens_args.add_argument(ctrlden_opt[0], **ctrlden_opt[1])

    mod_args = parser.add_argument_group('Standard Model Arguments')
    mod_args.add_argument(dna_opt[0], **dna_opt[1])
    mod_args.add_argument(rna_opt[0], **rna_opt[1])
    mod_args.add_argument(tbmod_opt[0], **tbmod_opt[1])

    stat_args = parser.add_argument_group('Model Fitting Arguments')
    stat_args.add_argument(altfrac_opt[0], **altfrac_opt[1])
    stat_args.add_argument(kernden_opt[0], **kernden_opt[1])

    filt_args = parser.add_argument_group('Filtering Argument')
    filt_args.add_argument(minkmer_opt[0], default=1000, **minkmer_opt[1])

    io_args = parser.add_argument_group('Output Argument')
    io_args.add_argument(densbn_opt[0], **densbn_opt[1])

    multi_args = parser.add_argument_group('Multiprocessing Arguments')
    multi_args.add_argument(proc_opt[0], default=1, **proc_opt[1])

    fast5_args, misc_args, parser = add_default_args(parser)

    return parser

def get_estimate_scale_parser():
    parser = argparse.ArgumentParser(
        description='Estimate a global scaling parameter from a ' +
        'sub-set of reads.', add_help=False)
    req_args = parser.add_argument_group('Required Arguments')
    req_args.add_argument(basedir_opt[0], **basedir_opt[1])

    misc_args = parser.add_argument_group('Miscellaneous Arguments')
    misc_args.add_argument(*quiet_opt[0], **quiet_opt[1])
    misc_args.add_argument(*help_opt[0], **help_opt[1])

    return parser


#########################################
###### Significance testing parser ######
#########################################

def get_test_signif_parser():
    parser = argparse.ArgumentParser(
        description='Test for significant shifts in raw nanopore signal ' +
        'against either a model, a set of two models or another sequencing ' +
        'sample.', add_help=False)
    req_args = parser.add_argument_group('Required Argument')
    req_args.add_argument(fast5dir_opt[0], required=True, **fast5dir_opt[1])
    req_args.add_argument(statbsnm_opt[0], required=True, **statbsnm_opt[1])

    alt_args = parser.add_argument_group(
        'Comparison Arguments (Default: De novo testing against default ' +
        'standard model)')
    alt_args.add_argument(modbs_opt[0], **modbs_opt[1])
    alt_args.add_argument(ctrlfast5dir_opt[0], **ctrlfast5dir_opt[1])
    alt_args.add_argument(dna_opt[0], **dna_opt[1])
    alt_args.add_argument(rna_opt[0], **rna_opt[1])
    alt_args.add_argument(hidden_tbmod_opt[0], **hidden_tbmod_opt[1])
    alt_args.add_argument(hidden_atbmods_opt[0], **hidden_atbmods_opt[1])

    test_args = parser.add_argument_group('Significance Test Arguments')
    test_args.add_argument(fmo_opt[0], **fmo_opt[1])
    test_args.add_argument(minreads_opt[0], default=1, **minreads_opt[1])
    test_args.add_argument(snglrdthrsh_opt[0], **snglrdthrsh_opt[1])

    io_args = parser.add_argument_group('Output Argument')
    io_args.add_argument(prstatbn_opt[0], **prstatbn_opt[1])

    multi_args = parser.add_argument_group('Multiprocessing Arguments')
    multi_args.add_argument(mpreg_opt[0], **mpreg_opt[1])
    multi_args.add_argument(proc_opt[0], default=1, **proc_opt[1])

    fast5_args, misc_args, parser = add_default_args(parser)

    return parser


############################
###### Filter parsers ######
############################

def get_clear_filters_parser():
    parser = argparse.ArgumentParser(
        description='Clear all filters applied to re-squiggled reads.',
        add_help=False)
    req_args = parser.add_argument_group('Required Argument')
    req_args.add_argument(fast5dir_opt[0], required=True, **fast5dir_opt[1])

    fast5_args = parser.add_argument_group('FAST5 Data Argument')
    fast5_args.add_argument(corrgrp_opt[0], **corrgrp_opt[1])

    misc_args, parser = add_misc_args(parser)

    return parser

def get_filter_stuck_parser():
    parser = argparse.ArgumentParser(
        description='Filter reads based on observations ' +
        'per base thresholds.', add_help=False)
    req_args = parser.add_argument_group('Required Argument')
    req_args.add_argument(fast5dir_opt[0], required=True, **fast5dir_opt[1])

    filter_args = parser.add_argument_group('Read Filtering Argument')
    filter_args.add_argument(obsfilt_opt[0], **obsfilt_opt[1])

    fast5_args = parser.add_argument_group('FAST5 Data Argument')
    fast5_args.add_argument(corrgrp_opt[0], **corrgrp_opt[1])

    misc_args, parser = add_misc_args(parser)

    return parser

def get_filter_coverage_parser():
    parser = argparse.ArgumentParser(
        description='Filter reads by downsampling for more even coverage.',
        add_help=False)
    req_args = parser.add_argument_group('Required Arguments')
    req_args.add_argument(fast5dir_opt[0], required=True, **fast5dir_opt[1])

    filter_args = parser.add_argument_group('Read Filtering Argument')
    filter_args.add_argument(pctfilt_opt[0], **pctfilt_opt[1])

    fast5_args = parser.add_argument_group('FAST5 Data Arguments')
    fast5_args.add_argument(corrgrp_opt[0], **corrgrp_opt[1])

    misc_args, parser = add_misc_args(parser)

    return parser


##############################################
###### Genome-anchored plotting parsers ######
##############################################

def get_max_cov_parser():
    parser = argparse.ArgumentParser(
        description='Plot raw signal in regions with the deepest coverage.',
        add_help=False)
    req_args = parser.add_argument_group('Required Argument')
    req_args.add_argument(fast5dir_opt[0], required=True, **fast5dir_opt[1])

    alt_args, parser = add_comp_dist_args(parser)

    ovplt_args = parser.add_argument_group('Overplotting Arguments')
    ovplt_args.add_argument(ovpltthresh_opt[0], **ovpltthresh_opt[1])
    ovplt_args.add_argument(ovplttype_opt[0], **ovplttype_opt[1])

    reg_args = parser.add_argument_group('Plotting Region Arguments')
    reg_args.add_argument(numreg_opt[0], default=10, **numreg_opt[1])
    reg_args.add_argument(numbases_opt[0], default=21, **numbases_opt[1])

    out_args = parser.add_argument_group('Output Argument')
    out_args.add_argument(pdf_opt[0], default=OUTPUT_BASE + '.max_coverage.pdf',
                          **pdf_opt[1])

    fast5_args, misc_args, parser = add_default_args(parser)

    return parser

def get_genome_loc_parser():
    parser = argparse.ArgumentParser(
        description='Plot raw signal at defined genomic locations.',
        add_help=False)
    req_args = parser.add_argument_group('Required Arguments')
    req_args.add_argument(fast5dir_opt[0], required=True, **fast5dir_opt[1])
    req_args.add_argument(gnmloc_opt[0], required=True, **gnmloc_opt[1])

    alt_args, parser = add_comp_dist_args(parser)

    ovplt_args = parser.add_argument_group('Overplotting Arguments')
    ovplt_args.add_argument(ovpltthresh_opt[0], **ovpltthresh_opt[1])
    ovplt_args.add_argument(ovplttype_opt[0], **ovplttype_opt[1])

    reg_args = parser.add_argument_group('Plotting Region Argument')
    reg_args.add_argument(numbases_opt[0], default=21, **numbases_opt[1])

    out_args = parser.add_argument_group('Output Argument')
    out_args.add_argument(pdf_opt[0],
                          default=OUTPUT_BASE + '.genome_locations.pdf',
                          **pdf_opt[1])

    fast5_args, misc_args, parser = add_default_args(parser)

    return parser

def get_motif_loc_parser():
    parser = argparse.ArgumentParser(
        description='Plot raw signal centered on a motif of interest.',
        add_help=False)
    req_args = parser.add_argument_group('Required Arguments')
    req_args.add_argument(fast5dir_opt[0], required=True, **fast5dir_opt[1])
    req_args.add_argument(motif_opt[0], required=True, **motif_opt[1])
    req_args.add_argument(fasta_opt[0], required=True, **fasta_opt[1])

    alt_args, parser = add_comp_dist_args(parser)

    ovplt_args = parser.add_argument_group('Overplotting Arguments')
    ovplt_args.add_argument(ovpltthresh_opt[0], **ovpltthresh_opt[1])
    ovplt_args.add_argument(ovplttype_opt[0], **ovplttype_opt[1])

    reg_args = parser.add_argument_group('Plotting Region Arguments')
    reg_args.add_argument(numreg_opt[0], default=10, **numreg_opt[1])
    reg_args.add_argument(numbases_opt[0], default=21, **numbases_opt[1])
    reg_args.add_argument(deepcov_opt[0], **deepcov_opt[1])

    out_args = parser.add_argument_group('Output Argument')
    out_args.add_argument(pdf_opt[0], default=OUTPUT_BASE + '.motif_centered.pdf',
                          **pdf_opt[1])

    fast5_args, misc_args, parser = add_default_args(parser)

    return parser

def get_max_diff_parser():
    parser = argparse.ArgumentParser(
        description='Plot raw signal where the average signal differs ' +
        'most between two samples.', add_help=False)
    req_args = parser.add_argument_group('Required Arguments')
    req_args.add_argument(fast5dir_opt[0], required=True, **fast5dir_opt[1])
    req_args.add_argument(ctrlfast5dir_opt[0], required=True,
                          **ctrlfast5dir_opt[1])

    ovplt_args = parser.add_argument_group('Overplotting Arguments')
    ovplt_args.add_argument(ovpltthresh_opt[0], **ovpltthresh_opt[1])
    ovplt_args.add_argument(ovplttype_opt[0], **ovplttype_opt[1])

    reg_args = parser.add_argument_group('Plotting Region Arguments')
    reg_args.add_argument(numreg_opt[0], default=10, **numreg_opt[1])
    reg_args.add_argument(numbases_opt[0], default=21, **numbases_opt[1])

    out_args = parser.add_argument_group('Output Arguments')
    out_args.add_argument(pdf_opt[0], default=OUTPUT_BASE + '.max_difference.pdf',
                          **pdf_opt[1])
    out_args.add_argument(seqs_opt[0], **seqs_opt[1])

    fast5_args, misc_args, parser = add_default_args(parser)

    return parser

def get_signif_diff_parser():
    parser = argparse.ArgumentParser(
        description='Plot raw signal at most significant genomic locations ' +
        'from previous test_significance results.', add_help=False)
    req_args = parser.add_argument_group('Required Arguments')
    req_args.add_argument(fast5dir_opt[0], required=True, **fast5dir_opt[1])
    req_args.add_argument(statfn_opt[0], required=True, **statfn_opt[1])

    alt_args, parser = add_comp_dist_args(parser)

    ovplt_args = parser.add_argument_group('Overplotting Arguments')
    ovplt_args.add_argument(ovpltthresh_opt[0], **ovpltthresh_opt[1])
    ovplt_args.add_argument(ovplttype_opt[0], **ovplttype_opt[1])

    reg_args = parser.add_argument_group('Plotting Region Arguments')
    reg_args.add_argument(numreg_opt[0], default=10, **numreg_opt[1])
    reg_args.add_argument(numbases_opt[0], default=21, **numbases_opt[1])

    out_args = parser.add_argument_group('Output Arguments')
    out_args.add_argument(pdf_opt[0],
                          default=OUTPUT_BASE + '.significant_difference.pdf',
                          **pdf_opt[1])
    out_args.add_argument(seqs_opt[0], **seqs_opt[1])

    fast5_args, misc_args, parser = add_default_args(parser)

    return parser

def get_signif_motif_parser():
    parser = argparse.ArgumentParser(
        description='Plot raw signal centered on motif of interest along ' +
        'with test statistic distributions at relative genomic position.',
        add_help=False)
    req_args = parser.add_argument_group('Required Arguments')
    req_args.add_argument(fast5dir_opt[0], required=True, **fast5dir_opt[1])
    req_args.add_argument(motif_opt[0], required=True, **motif_opt[1])
    req_args.add_argument(statfn_opt[0], required=True, **statfn_opt[1])

    alt_args, parser = add_comp_dist_args(parser)

    ovplt_args = parser.add_argument_group('Overplotting Argument')
    ovplt_args.add_argument(ovpltthresh_opt[0], **ovpltthresh_opt[1])

    reg_args = parser.add_argument_group('Plotting Region Arguments')
    reg_args.add_argument(numreg_opt[0], default=3, **numreg_opt[1])
    reg_args.add_argument(cntxt_opt[0], **cntxt_opt[1])
    reg_args.add_argument(numstat_opt[0], **numstat_opt[1])

    seq_args = parser.add_argument_group(
        'Sequence Argument (for faster sequence access)')
    seq_args.add_argument(fasta_opt[0], **fasta_opt[1])

    out_args = parser.add_argument_group('Output Argument')
    out_args.add_argument(pdf_opt[0],
                          default=OUTPUT_BASE + '.motif_statistics.pdf',
                          **pdf_opt[1])

    fast5_args, misc_args, parser = add_default_args(parser)

    return parser

def get_per_read_parser():
    parser = argparse.ArgumentParser(
        description='Plot non-standard base statistic per read at specified ' +
        'genomic locations.', add_help=False)
    req_args = parser.add_argument_group('Required Arguments')
    req_args.add_argument(gnmloc_opt[0], required=True, **gnmloc_opt[1])
    req_args.add_argument(prstat_opt[0], required=True, **prstat_opt[1])

    seq_args = parser.add_argument_group(
        'Sequence Arguments (Provide either FAST5s dir or genome FASTA)')
    seq_args.add_argument(fasta_opt[0], **fasta_opt[1])
    seq_args.add_argument(fast5dir_opt[0], **fast5dir_opt[1])

    reg_args = parser.add_argument_group('Plotting Region Arguments')
    reg_args.add_argument(numreads_opt[0], default=100, **numreads_opt[1])
    reg_args.add_argument(numbases_opt[0], default=51, **numbases_opt[1])
    reg_args.add_argument(boxc_opt[0], **boxc_opt[1])

    out_args = parser.add_argument_group('Output Argument')
    out_args.add_argument(
        pdf_opt[0], default=OUTPUT_BASE + '.per_read_stats.pdf', **pdf_opt[1])

    fast5_args, misc_args, parser = add_default_args(parser)

    return parser


####################################
###### Other plotting parsers ######
####################################

def get_correction_parser():
    parser = argparse.ArgumentParser(
        description='Plot segments before and after event-based re-resquiggle.',
        add_help=False)
    req_args = parser.add_argument_group('Required Argument')
    req_args.add_argument(fast5dir_opt[0], required=True, **fast5dir_opt[1])

    type_args = parser.add_argument_group('Region Type Argument')
    type_args.add_argument(regtype_opt[0], **regtype_opt[1])

    reg_args = parser.add_argument_group('Plotting Region Arguments')
    reg_args.add_argument(numreads_opt[0], default=10, **numreads_opt[1])
    reg_args.add_argument(numobs_opt[0], **numobs_opt[1])

    out_args = parser.add_argument_group('Output Argument')
    out_args.add_argument(pdf_opt[0], default=OUTPUT_BASE + '.corrected.pdf',
                          **pdf_opt[1])

    fast5_args, misc_args, parser = add_default_args(parser)

    return parser

def get_multi_correction_parser():
    parser = argparse.ArgumentParser(
        description='Plot segments before and after event-based ' +
        're-squiggle for multiple reads anchored at a genomic location.',
        add_help=False)
    req_args = parser.add_argument_group('Required Argument')
    req_args.add_argument(fast5dir_opt[0], required=True, **fast5dir_opt[1])

    reg_args = parser.add_argument_group('Plotting Region Arguments')
    reg_args.add_argument(gnmloc_opt[0], **gnmloc_opt[1])
    reg_args.add_argument(numreg_opt[0], default=10, **numreg_opt[1])
    reg_args.add_argument(numreads_opt[0], default=5, **numreads_opt[1])
    reg_args.add_argument(numobs_opt[0], **numobs_opt[1])

    out_args = parser.add_argument_group('Output Arguments')
    out_args.add_argument(pdf_opt[0],
                          default=OUTPUT_BASE + '.multi_corrected.pdf',
                          **pdf_opt[1])
    out_args.add_argument(origbcs_opt[0], **origbcs_opt[1])

    fast5_args, misc_args, parser = add_default_args(parser)

    return parser

def get_kmer_dist_parser():
    parser = argparse.ArgumentParser(
        description='Plot signal level distribution across kmers.',
        add_help=False)
    req_args = parser.add_argument_group('Required Argument')
    req_args.add_argument(fast5dir_opt[0], required=True, **fast5dir_opt[1])

    proc_args = parser.add_argument_group('Data Processing Arguments')
    proc_args.add_argument(upstrmbs_opt[0], **upstrmbs_opt[1])
    proc_args.add_argument(dnstrmbs_opt[0], **dnstrmbs_opt[1])
    proc_args.add_argument(readmean_opt[0], **readmean_opt[1])
    proc_args.add_argument(kmerthresh_opt[0], **kmerthresh_opt[1])

    reg_args = parser.add_argument_group('Plotting Region Arguments')
    reg_args.add_argument(numreads_opt[0], default=100, **numreads_opt[1])

    out_args = parser.add_argument_group('Output Arguments')
    out_args.add_argument(pdf_opt[0],
                         default=OUTPUT_BASE + '.kmer_distribution.pdf',
                         **pdf_opt[1])
    out_args.add_argument(rdata_opt[0], **rdata_opt[1])
    out_args.add_argument(noplot_opt[0], **noplot_opt[1])

    fast5_args, misc_args, parser = add_default_args(parser)

    return parser

def get_roc_parser():
    parser = argparse.ArgumentParser(
        description='Plot ROC curve given known motif(s).',
        add_help=False)
    req_args = parser.add_argument_group('Required Argument')
    req_args.add_argument(statfns_opt[0], required=True, **statfns_opt[1])
    req_args.add_argument(motifdesc_opt[0], required=True, **motifdesc_opt[1])
    req_args.add_argument(fasta_opt[0], required=True, **fasta_opt[1])

    out_args = parser.add_argument_group('Output Arguments')
    out_args.add_argument(pdf_opt[0],
                         default=OUTPUT_BASE + '.roc.pdf',
                         **pdf_opt[1])

    filt_args = parser.add_argument_group('Filtering Arguments')
    filt_args.add_argument(minreads_opt[0], default=1, **minreads_opt[1])

    misc_args, parser = add_misc_args(parser)

    return parser

def get_cluster_signif_diff_parser():
    parser = argparse.ArgumentParser(
        description='Cluster signal trace differences at most significant ' +
        'signal shifts from previous test_significance results.',
        add_help=False)
    req_args = parser.add_argument_group('Required Arguments')
    req_args.add_argument(fast5dir_opt[0], required=True, **fast5dir_opt[1])
    req_args.add_argument(ctrlfast5dir_opt[0], required=True,
                          **ctrlfast5dir_opt[1])
    req_args.add_argument(statfn_opt[0], required=True, **statfn_opt[1])

    fasta_args = parser.add_argument_group('FASTA Sequence Argument')
    fasta_args.add_argument(fasta_opt[0], **fasta_opt[1])

    multi_args = parser.add_argument_group('Multiprocessing Argument')
    multi_args.add_argument(proc_opt[0], default=1, **proc_opt[1])

    reg_args = parser.add_argument_group('Plotting Region Arguments')
    reg_args.add_argument(numreg_opt[0], default=10, **numreg_opt[1])
    reg_args.add_argument(numbases_opt[0], default=21, **numbases_opt[1])
    reg_args.add_argument(slides_opt[0], **slides_opt[1])

    out_args = parser.add_argument_group('Output Arguments')
    out_args.add_argument(pdf_opt[0],
                          default=OUTPUT_BASE + '.signal_clusters.pdf',
                          **pdf_opt[1])
    out_args.add_argument(rdata_opt[0], **rdata_opt[1])

    fast5_args, misc_args, parser = add_default_args(parser)

    return parser


#################################
###### Text output parsers ######
#################################

def get_wiggle_parser():
    parser = argparse.ArgumentParser(
        description='Write wiggle files for specified data types.',
        add_help=False)
    data_args = parser.add_argument_group('Data Arguments')
    data_args.add_argument(fast5dir_opt[0], **fast5dir_opt[1])
    data_args.add_argument(ctrlfast5dir_opt[0], **ctrlfast5dir_opt[1])
    data_args.add_argument(statfn_opt[0], **statfn_opt[1])

    out_args = parser.add_argument_group('Output Arguments')
    out_args.add_argument(wigfn_opt[0], **wigfn_opt[1])
    out_args.add_argument(wigtypes_opt[0], **wigtypes_opt[1])

    stat_args = parser.add_argument_group('Statistical Argument')
    stat_args.add_argument(cvgdmp_opt[0], **cvgdmp_opt[1])

    fast5_args, misc_args, parser = add_default_args(parser)

    return parser

def get_write_signif_diff_parser():
    parser = argparse.ArgumentParser(
        description='Write sequence at genomic locations with most ' +
        'significant difference from previous test_significance results.',
        add_help=False)
    req_args = parser.add_argument_group('Required Argument')
    req_args.add_argument(statfn_opt[0], required=True, **statfn_opt[1])

    seq_args = parser.add_argument_group(
        'Sequence Arguments (Provide either FAST5s dir or genome FASTA)')
    seq_args.add_argument(fasta_opt[0], **fasta_opt[1])
    seq_args.add_argument(fast5dir_opt[0], **fast5dir_opt[1])

    reg_args = parser.add_argument_group('Region Selection Arguments')
    reg_args.add_argument(numreg_opt[0], default=100, **numreg_opt[1])
    reg_args.add_argument(numbases_opt[0], default=15, **numbases_opt[1])

    out_args = parser.add_argument_group('Output Arguments')
    out_args.add_argument(
        seqs_opt[0], default=OUTPUT_BASE + '.significant_regions.fasta',
        **seqs_opt[1])

    fast5_args, misc_args, parser = add_default_args(parser)

    return parser


if __name__ == '__main__':
    raise NotImplementedError(
        'This is a module. See commands with `tombo -h`')
