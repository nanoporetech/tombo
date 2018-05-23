from __future__ import unicode_literals, absolute_import

from builtins import map

import sys
import argparse

if sys.version_info[0] > 2:
    unicode = str

from ._default_parameters import (
    SEG_PARAMS_TABLE, ALGN_PARAMS_TABLE, LLR_THRESH, SAMP_COMP_THRESH,
    DE_NOVO_THRESH, ALTERNATE_MODELS, MAX_SCALING_ITERS, ALT_EST_PCTL,
    COV_DAMP_COUNTS, SIG_MATCH_THRESH)

ALT_BASES = tuple(set(alt_name.split('_')[1] for alt_name in ALTERNATE_MODELS))


##################################
###### Positional arguments ######
##################################

basedir_opt=('fast5_basedir', {
    'type':unicode,
    'help':'Directory containing fast5 files. All files ending in "fast5" ' +
    'found recursively within this base directory will be processed.'})
fasta_pos_opt=('reference', {
    'type':unicode, 'help':'Reference genome/transcriptome FASTA file ' +
    'or minimap2 index (with "map-ont" preset) for mapping.'})
fasta_event_opt=('reference_fasta', {
    'type':unicode, 'help':'Reference genome/transcriptome FASTA file ' +
    'for mapping.'})
# put re-squiggle positional arguments in one argument to allow printing
# hidden arguments help
rsqgl_pos_opt=('fast5s_and_reference', {
    'type':unicode, 'nargs':'*',
    'help':'Directory containing fast5 files and a genome/transcriptome ' +
    'reference. Directory will be searched recursively for files ending ' +
    'in ".fast5". Reference may be a FASTA file or minimap2 index file.'})


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
bcgrp_opt=('--basecall-group', {
    'type':unicode, 'default':'Basecall_1D_000',
    'help':'FAST5 group obtain original basecalls (under Analyses group). ' +
    'Default: %(default)s'})
bcsubgrps_opt=('--basecall-subgroups', {
    'type':unicode, 'default':['BaseCalled_template',], 'nargs':'+',
    'help':'FAST5 subgroup(s) (under /Analyses/[--basecall-group]/) ' +
    'containing basecalls and created within [--corrected-group] ' +
    'containing re-squiggle results. Default: %(default)s'})
bcsubgrp_opt=('--basecall-subgroup', {
    'type':unicode, 'default':'BaseCalled_template',
    'help':'FAST5 subgroup (under /Analyses/[--basecall-group]/) under which ' +
    'to store basecalls from FASTQs. Default: %(default)s'})

gnmloc_opt=('--genome-locations', {
    'type':unicode, 'nargs':'+',
    'help':'Genomic locations at which to plot signal. Format locations ' +
    'as "chrm:position[:strand] [chrm2:position2[:strand2] ...]" ' +
    '(strand not applicable for all applications)'})
incldreg_opt=('--include-regions', {
    'type':unicode, 'nargs':'+',
    'help':'Filter out reads not falling completely within include regions. ' +
    'Omit start and end coordinates to include an entire chromosome/sequence ' +
    'record. Format regions as "chrm[:start-end] [chrm2[:start2-end2] ...]".'})

fasta_opt=('--genome-fasta', {
    'type':unicode,
    'help':'FASTA file used to re-squiggle. For faster sequence access.'})
motif_opt=('--motif', {
    'type':unicode,
    'help':'Motif of interest at which to plot signal and statsitics. ' +
    'Supports IUPAC single letter codes (use T for RNA).'})

obsfilt_opt=('--obs-per-base-filter', {
    'type':unicode, 'nargs':'+', 'default':[],
    'help':'Filter reads based on observations per base percentile ' +
    'thresholds. Format thresholds as "percentile:thresh ' +
    '[pctl2:thresh2 ...]". For example to filter reads with 99th ' +
    'pctl > 200 obs/base or max > 5k obs/base use "99:200 100:5000".'})

fastqs_opt = ('--fastq-filenames', {
    'type':unicode, 'nargs':'+',
    'help':'FASTQ filenames containing basecalls to be added to ' +
    'raw FAST5 files.'})
seqsum_opt = ('--sequencing-summary-filenames', {
    'type':unicode, 'nargs':'+',
    'help':'Sequencing summary filenames produced by albacore. These can ' +
    'make annotation of raw FAST5 files with FASTQ sequence much faster.'})

brsrfn_opt=('--browser-file-basename', {
    'type':unicode, 'default':'tombo_results',
    'help':'Basename for output browser files. Two files (plus and minus ' +
    'strand) will be produced for each --file-types supplied. ' +
    'Filenames formatted as "[browser-file-basename].[file-type].' +
    '[sample|control]?.[plus|minus].[wig|bedgraph]". Default: %(default)s'})
pdf_opt=('--pdf-filename', {
    'type':unicode, 'help':'PDF filename to store plot(s). ' +
    'Default: %(default)s'})
statfn_opt=('--statistics-filename', {
    'type':unicode, 'help':'File to save/load genomic base anchored ' +
    'statistics.'})
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
prstats_opt=('--per-read-statistics-filenames', {
    'type':unicode, 'nargs':'+',
    'help':'Binary files containing per-read statistics from ' +
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
thrpp_opt=('--threads-per-process', {
    'type':int,
    'help':'Number of file input/output and mapping threads per compute ' +
    'process [--processes]. This should likely be left at 1, but may ' +
    'improve performance on some systems. Default: %(default)d'})

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

covthresh_opt=('--coverage-threshold', {
    'type':int,
    'help':'Maximum mean coverage per region when estimating k-mer model ' +
    '(limits compute time for deep samples). Default: %(default)d'})
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
minreads_opt=('--minimum-test-reads', {
    'type':int,
    'help':'Number of reads required at a position to perform significance ' +
    'testing or contribute to model estimation. Default: %(default)d'})

spb_opt=('--statistics-per-block', {
    'type':int, 'default':100000,
    'help':'Number of randomly selected per-read, per-base statistics to ' +
    'extract from each genomic block for plotting. Default: %(default)d'})
tsl_opt=('--total-statistics-limit', {
    'type':int, 'default':5000000,
    'help':'Total per-read statistics to be extracted for plotting. ' +
    'Avoids memory overflow for large runs. Default: %(default)d'})

segpars_opt=('--segmentation-parameters', {
    'type':int, 'nargs':3,
    'help':'Specify the 3 parameters for segmentation 1) running neighboring ' +
    'windows width 2) minimum raw observations per genomic base 3) mean raw ' +
    'observations per event. Sample type defaults: ' +
    ' || '.join((bst + ' : ' + ' '.join(map(str, params)))
                for bst, params in SEG_PARAMS_TABLE.items())})
hidsegpars_opt=('--segmentation-parameters', {
    'type':int, 'nargs':3, 'help':argparse.SUPPRESS})
segpars2_opt=('--segmentation-parameters', {
    'type':int, 'nargs':2,
    'help':'Specify the 2 parameters for segmentation 1) running neighboring ' +
    'windows width 2) minimum raw observations per genomic base. Sample type ' +
    'defaults:\n' +
    ' || '.join((bst + ' : ' + ' '.join(map(str, params[:2])))
                for bst, params in SEG_PARAMS_TABLE.items())})
msi_opt=('--max-scaling-iterations', {
    'type':int, 'default':MAX_SCALING_ITERS,
    'help':'Maximum re-squiggle iterations to perform. At each iteration ' +
    'the signal normalization parameters are re-fit. Higher values ' +
    'recommended for highly modified reads. Default: %(default)d'})
hidmsi_opt=('--max-scaling-iterations', {
    'type':int, 'default':MAX_SCALING_ITERS, 'help':argparse.SUPPRESS})


###############################
###### Boolean arguments ######
###############################

skpidx_opt=('--skip-index', {
    'default':False, 'action':'store_true',
    'help':'Skip creation of tombo index. This drastically slows downstream '+
    'tombo commands. Default stores tombo index named ".[--fast5-basedir].' +
    '[--corrected-group].tombo.index" to be loaded automatically for ' +
    'downstream commands.'})
hidskpidx_opt=('--skip-index', {
    'default':False, 'action':'store_true', 'help':argparse.SUPPRESS})
ovrwrt_opt=('--overwrite', {
    'default':False, 'action':'store_true',
    'help':'Overwrite previous corrected group in FAST5 files. Note: ' +
    'only effects --corrected-group or --new-corrected-group.'})
ignrlock_opt=('--ignore-read-locks', {
    'default':False, 'action':'store_true',
    'help':'Ignore Tombo locks, used to ensure that reads are only accessed ' +
    'from a single resquiggle processes avoiding potential file corruption.'})
hidignrlock_opt=('--ignore-read-locks', {
    'default':False, 'action':'store_true', 'help':argparse.SUPPRESS})
printadv_opt=('--print-advanced-arguments', {
    'default':False, 'action':'store_true',
    'help':'Print advanced re-squiggle arguments and exit.'})
printalt_opt=('--print-available-models', {
    'default':False, 'action':'store_true',
    'help':'Print available alternative models and exit.'})

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
hidincldsd_opt=('--include-event-stdev', {
    'default':False, 'action':'store_true', 'help':argparse.SUPPRESS})
fitscl_opt=('--fit-global-scale', {
    'default':False, 'action':'store_true',
    'help':'Fit a global scaling parameter for all reads. Otherwise fit ' +
    'the scaling parameter for each read. Global parameter estimated from ' +
    'a random subset of reads, which may produce more robust results for ' +
    'some samples.'})
hidfitscl_opt=('--fit-global-scale', {
    'default':False, 'action':'store_true', 'help':argparse.SUPPRESS})
sss_opt=('--skip-sequence-rescaling', {
    'default':False, 'action':'store_true',
    'help':'Skip sequence-based re-scaling. Otherwise, after re-squiggle, ' +
    'signal normalization parameters are re-fit (using Theil-Sen estimator).'})
hidsss_opt=('--skip-sequence-rescaling', {
    'default':False, 'action':'store_true', 'help':argparse.SUPPRESS})
stdllhr_opt=('--standard-log-likelihood-ratio', {
    'default':False, 'action':'store_true',
    'help':'Use a standard log likelihood ratio (LLR) statistic. Default ' +
    'is to use an outlier-robust LLR-like statistic. Detail in full ' +
    'online documentation.'})

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
hidotlthresh_opt=('--outlier-threshold', {
    'default':5, 'type':float, 'help':argparse.SUPPRESS})

snglrdthrsh_opt=('--single-read-threshold', {
    'type':float, 'nargs':'+',
    'help':(
        'P-value or log likelihood ratio threshold when computing ' +
        'fraction of significant reads at each genomic position. If two ' +
        'values are provided, statistics between these values are not ' +
        'considered.')})
dnthresh_opt=('--single-read-threshold', {
    'type':float, 'nargs':'+',
    'help':(
        'P-value threshold when computing fraction of significant reads at ' +
        'each genomic position. If two values are provided, statistics ' +
        'between these values are not considered. ' +
        'Default thresholds: ' +
        ', '.join(bst + ':' + (str(thresh[1]) if thresh[0] is None else
                               str(thresh[0]) + '-' + str(thresh[1])) + ' '
                  for bst, thresh in DE_NOVO_THRESH.items()))})
scompthresh_opt=('--single-read-threshold', {
    'type':float, 'nargs':'+',
    'help':(
        'P-value threshold when computing fraction of significant reads at ' +
        'each genomic position. If two values are provided, statistics ' +
        'between these values are not considered. ' +
        'Default thresholds: ' +
        ', '.join(bst + ':' + (str(thresh[1]) if thresh[0] is None else
                               str(thresh[0]) + '-' + str(thresh[1])) + ' '
                  for bst, thresh in SAMP_COMP_THRESH.items()))})
altthresh_opt=('--single-read-threshold', {
    'type':float, 'nargs':'+',
    'help':(
        'Log likelihood ratio threshold when computing fraction of ' +
        'significant reads at each genomic position. If two values ' +
        'are provided, statistics between these values are not considered. ' +
        'Default thresholds: ' +
        ', '.join(bst + ':' + (str(thresh[1]) if thresh[0] is None else
                               str(thresh[0]) + '-' + str(thresh[1])) + ' '
                  for bst, thresh in LLR_THRESH.items()))})

altfrac_opt=('--alt-fraction-percentile', {
    'default':ALT_EST_PCTL, 'type':float,
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
    'location. This can be useful in modeling and testing. ' +
    'Default: %(default)f'})
qscr_opt=('--q-score', {
    'type':float,
    'help':'Q-score threshold for filtering low quality reads. ' +
    'Default: %(default)f'})
sms_opt=('--signal-matching-score', {
    'type':float,
    'help':'Mean half normal z-score threshold for filtering reads with ' +
    'poor raw to expected signal matching. Signal type defaults: ' +
    ' || '.join(bst + ' : ' + str(params)
                for bst, params in SIG_MATCH_THRESH.items())})
fxdscl_opt=('--fixed-scale', {
    'type':float,
    'help':'Fixed scaling parameter to use for raw signal normalization.'})
hidfxdscl_opt=('--fixed-scale', {
    'type':float, 'help':argparse.SUPPRESS})
cvgdmp_opt=('--coverage-dampen-counts', {
    'type':float, 'nargs':2, 'default':COV_DAMP_COUNTS,
    'help':'Dampen fraction modified estimates for low coverage sites. Two ' +
    'parameters are unmodified and modified psuedo read counts. This is ' +
    'equivalent to a beta prior on the fraction estimate. Set to "0 0" to ' +
    'disable dampened fraction estimation. Default: %(default)s'})

sigapars_opt=('--signal-align-parameters', {
    'type':float, 'nargs':5,
    'help':'Specify the 4 parameters for signal to genome sequence alignment ' +
    'algorithm 1) match expected value 2) skip penalty 3) bandwidth 4) save ' +
    'bandwidth (if read fails with bandwith) 5) z-score winsorizing ' +
    'threshold. Sample type defaults: ' + ' || '.join(
        (bst + ' : ' + ' '.join(map(str, params)))
        for bst, params in ALGN_PARAMS_TABLE.items())})
hidsigapars_opt=('--signal-align-parameters', {
    'type':float, 'nargs':5, 'help':argparse.SUPPRESS})


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
ovplttype_opt=('--overplot-type', {
    'type':unicode, 'default':'Downsample',
    'choices':['Downsample', 'Boxplot', 'Quantile', 'Density'],
    'help':'Plot type for regions with higher coverage. Default: Downsample'})
ftypes_opt=('--file-types', {
    'type':unicode, 'default':['coverage', ], 'nargs':'+',
    'choices':['coverage', 'valid_coverage', 'fraction', 'dampened_fraction',
               'signal', 'signal_sd', 'dwell', 'difference'],
    'help':'Data types of genome browser files to produce. Produced coverage ' +
    'files are in bedGraph format, while all other file types  will be in ' +
    'wiggle format (https://genome.ucsc.edu/goldenpath/help/wiggle.html). ' +
    'Default: "coverage"'})

dna_opt=('--dna', {
    'dest':'bio_sample_type', 'action':'store_const', 'const':'DNA',
    'help':'Explicitly select canonical DNA model. Default: Automatically ' +
    'determine from FAST5s'})
rna_opt=('--rna', {
    'dest':'bio_sample_type', 'action':'store_const', 'const':'RNA',
    'help':'Explicitly select canonical RNA model. Default: Automatically ' +
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

def add_common_testing_args(parser):
    io_args = parser.add_argument_group('Output Argument')
    io_args.add_argument(prstatbn_opt[0], **prstatbn_opt[1])

    multi_args = parser.add_argument_group('Multiprocessing Arguments')
    multi_args.add_argument(mpreg_opt[0], **mpreg_opt[1])
    multi_args.add_argument(proc_opt[0], default=1, **proc_opt[1])

    return io_args, multi_args

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

def get_resquiggle_parser():
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

    filt_args = parser.add_argument_group('Read Filtering Argument')
    filt_args.add_argument(obsfilt_opt[0], **obsfilt_opt[1])
    filt_args.add_argument(qscr_opt[0], default=0, **qscr_opt[1])
    filt_args.add_argument(sms_opt[0], **sms_opt[1])

    multi_args = parser.add_argument_group('Multiprocessing Arguments')
    multi_args.add_argument(proc_opt[0], default=1, **proc_opt[1])
    multi_args.add_argument(thrpp_opt[0], default=1, **thrpp_opt[1])

    fast5_args = parser.add_argument_group('FAST5 Data Arguments')
    fast5_args.add_argument(corrgrp_opt[0], **corrgrp_opt[1])
    fast5_args.add_argument(bcgrp_opt[0], **bcgrp_opt[1])
    fast5_args.add_argument(bcsubgrps_opt[0], **bcsubgrps_opt[1])
    fast5_args.add_argument(ovrwrt_opt[0], **ovrwrt_opt[1])

    io_args = parser.add_argument_group('Input/Output Arguments')
    io_args.add_argument(failed_opt[0], **failed_opt[1])

    hid_args = parser.add_argument_group('Advanced Arguments')
    hid_args.add_argument(printadv_opt[0], **printadv_opt[1])
    hid_args.add_argument(hidden_tbmod_opt[0], **hidden_tbmod_opt[1])
    hid_args.add_argument(hidsegpars_opt[0], **hidsegpars_opt[1])
    hid_args.add_argument(hidsigapars_opt[0], **hidsigapars_opt[1])
    hid_args.add_argument(hidsss_opt[0], **hidsss_opt[1])
    hid_args.add_argument(hidmsi_opt[0], **hidmsi_opt[1])
    hid_args.add_argument(hidfitscl_opt[0], **hidfitscl_opt[1])
    hid_args.add_argument(hidfxdscl_opt[0], **hidfxdscl_opt[1])
    hid_args.add_argument(hidotlthresh_opt[0], **hidotlthresh_opt[1])
    hid_args.add_argument(hidskpidx_opt[0], **hidskpidx_opt[1])
    hid_args.add_argument(hidincldsd_opt[0], **hidincldsd_opt[1])
    hid_args.add_argument(hidignrlock_opt[0], **hidignrlock_opt[1])

    misc_args, parser = add_misc_args(parser)

    return parser

def print_advanced_resquiggle():
    parser = argparse.ArgumentParser(
        description='Hidden parameters to the resquiggle command.',
        add_help=False, usage='')
    hid_args = parser.add_argument_group('Hidden Arguments')
    hid_args.add_argument(tbmod_opt[0], **tbmod_opt[1])
    hid_args.add_argument(segpars_opt[0], **segpars_opt[1])
    hid_args.add_argument(sigapars_opt[0], **sigapars_opt[1])
    hid_args.add_argument(sss_opt[0], **sss_opt[1])
    hid_args.add_argument(msi_opt[0], **msi_opt[1])
    hid_args.add_argument(fitscl_opt[0], **fitscl_opt[1])
    hid_args.add_argument(fxdscl_opt[0], **fxdscl_opt[1])
    hid_args.add_argument(otlthresh_opt[0], **otlthresh_opt[1])
    hid_args.add_argument(skpidx_opt[0], **skpidx_opt[1])
    hid_args.add_argument(incldsd_opt[0], **incldsd_opt[1])
    hid_args.add_argument(ignrlock_opt[0], **ignrlock_opt[1])

    hid_args.add_argument(*help_opt[0], **help_opt[1])

    return parser.parse_args(['-h',])


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
    req_args.add_argument(fasta_event_opt[0], **fasta_event_opt[1])

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

    seqsum_args = parser.add_argument_group('Sequencing Summary Argument')
    seqsum_args.add_argument(seqsum_opt[0], **seqsum_opt[1])

    multi_args = parser.add_argument_group('Multiprocessing Argument')
    multi_args.add_argument(proc_opt[0], default=1, **proc_opt[1])

    misc_args, parser = add_misc_args(parser)

    return parser


######################################
###### Model estimation parsers ######
######################################

def get_est_ref_parser():
    parser = argparse.ArgumentParser(
        description='Estimate standard tombo model for use in re-squiggle ' +
        'and  testing without an amplified (un-modified) sample.',
        add_help=False)
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

def get_de_novo_test_signif_parser():
    parser = argparse.ArgumentParser(
        description='Test for significant shifts in raw nanopore signal ' +
        'against either a canonical model.', add_help=False)
    req_args = parser.add_argument_group('Required Argument')
    req_args.add_argument(fast5dir_opt[0], required=True, **fast5dir_opt[1])
    req_args.add_argument(statbsnm_opt[0], required=True, **statbsnm_opt[1])

    alt_args = parser.add_argument_group('Comparison Arguments')
    alt_args.add_argument(dna_opt[0], **dna_opt[1])
    alt_args.add_argument(rna_opt[0], **rna_opt[1])
    alt_args.add_argument(hidden_tbmod_opt[0], **hidden_tbmod_opt[1])

    test_args = parser.add_argument_group('Significance Test Arguments')
    test_args.add_argument(fmo_opt[0], **fmo_opt[1])
    test_args.add_argument(minreads_opt[0], default=1, **minreads_opt[1])
    test_args.add_argument(dnthresh_opt[0], **dnthresh_opt[1])

    io_args, multi_args = add_common_testing_args(parser)
    fast5_args, misc_args, parser = add_default_args(parser)

    return parser

def get_alt_test_signif_parser():
    parser = argparse.ArgumentParser(
        description='Test for significant shifts in raw nanopore signal ' +
        'which match a specific non-canonical base model.', add_help=False)
    req_args = parser.add_argument_group('Required Argument')
    req_args.add_argument(fast5dir_opt[0], **fast5dir_opt[1])
    req_args.add_argument(statbsnm_opt[0], **statbsnm_opt[1])

    alt_args = parser.add_argument_group('Comparison Arguments')
    alt_args.add_argument(modbs_opt[0], **modbs_opt[1])
    alt_args.add_argument(printalt_opt[0], **printalt_opt[1])
    alt_args.add_argument(dna_opt[0], **dna_opt[1])
    alt_args.add_argument(rna_opt[0], **rna_opt[1])
    alt_args.add_argument(hidden_tbmod_opt[0], **hidden_tbmod_opt[1])
    alt_args.add_argument(hidden_atbmods_opt[0], **hidden_atbmods_opt[1])

    test_args = parser.add_argument_group('Significance Test Arguments')
    test_args.add_argument(minreads_opt[0], default=1, **minreads_opt[1])
    test_args.add_argument(altthresh_opt[0], **altthresh_opt[1])
    test_args.add_argument(stdllhr_opt[0], **stdllhr_opt[1])

    io_args, multi_args = add_common_testing_args(parser)
    fast5_args, misc_args, parser = add_default_args(parser)

    return parser

def get_samp_comp_test_signif_parser():
    parser = argparse.ArgumentParser(
        description='Test for significant shifts in raw nanopore signal ' +
        'against either a model, a set of two models or another sequencing ' +
        'sample.', add_help=False)
    req_args = parser.add_argument_group('Required Argument')
    req_args.add_argument(fast5dir_opt[0], required=True, **fast5dir_opt[1])
    req_args.add_argument(statbsnm_opt[0], required=True, **statbsnm_opt[1])

    alt_args = parser.add_argument_group('Comparison Arguments')
    alt_args.add_argument(ctrlfast5dir_opt[0], **ctrlfast5dir_opt[1])

    test_args = parser.add_argument_group('Significance Test Arguments')
    test_args.add_argument(fmo_opt[0], **fmo_opt[1])
    test_args.add_argument(minreads_opt[0], default=1, **minreads_opt[1])
    test_args.add_argument(scompthresh_opt[0], **scompthresh_opt[1])

    io_args, multi_args = add_common_testing_args(parser)
    fast5_args, misc_args, parser = add_default_args(parser)

    return parser

def get_aggregate_per_read_parser():
    parser = argparse.ArgumentParser(
        description='Aggregate per-read statistics to produce a standard ' +
        '(genomic base) statistics file.', add_help=False)
    req_args = parser.add_argument_group('Required Argument')
    req_args.add_argument(prstat_opt[0], required=True, **prstat_opt[1])
    req_args.add_argument(statbsnm_opt[0], required=True, **statbsnm_opt[1])
    req_args.add_argument(snglrdthrsh_opt[0], required=True, **snglrdthrsh_opt[1])

    test_args = parser.add_argument_group('Significance Test Arguments')
    test_args.add_argument(minreads_opt[0], default=1, **minreads_opt[1])

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

def get_filter_qscore_parser():
    parser = argparse.ArgumentParser(
        description='Filter reads to remove low quality reads.',
        add_help=False)
    req_args = parser.add_argument_group('Required Arguments')
    req_args.add_argument(fast5dir_opt[0], required=True, **fast5dir_opt[1])

    filter_args = parser.add_argument_group('Read Filtering Argument')
    filter_args.add_argument(qscr_opt[0], default=7, **qscr_opt[1])

    fast5_args = parser.add_argument_group('FAST5 Data Arguments')
    fast5_args.add_argument(corrgrp_opt[0], **corrgrp_opt[1])
    fast5_args.add_argument(bcgrp_opt[0], **bcgrp_opt[1])

    misc_args, parser = add_misc_args(parser)

    return parser

def get_filter_signal_matching_parser():
    parser = argparse.ArgumentParser(
        description='Filter reads with poor raw to expected signal matching.',
        add_help=False)
    req_args = parser.add_argument_group('Required Arguments')
    req_args.add_argument(fast5dir_opt[0], required=True, **fast5dir_opt[1])
    req_args.add_argument(sms_opt[0], required=True, **sms_opt[1])

    fast5_args = parser.add_argument_group('FAST5 Data Arguments')
    fast5_args.add_argument(corrgrp_opt[0], **corrgrp_opt[1])

    misc_args, parser = add_misc_args(parser)

    return parser

def get_filter_genome_pos_parser():
    parser = argparse.ArgumentParser(
        description='Filter reads based on genome mapping location.',
        add_help=False)
    req_args = parser.add_argument_group('Required Arguments')
    req_args.add_argument(fast5dir_opt[0], required=True, **fast5dir_opt[1])
    req_args.add_argument(incldreg_opt[0], **incldreg_opt[1])

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
    out_args.add_argument(pdf_opt[0],
                          default=OUTPUT_BASE + '.motif_centered.pdf',
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

    stat_args = parser.add_argument_group('Statistical Argument')
    stat_args.add_argument(cvgdmp_opt[0], **cvgdmp_opt[1])

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
    req_args.add_argument(fasta_opt[0], required=True,  **fasta_opt[1])

    alt_args, parser = add_comp_dist_args(parser)

    ovplt_args = parser.add_argument_group('Overplotting Argument')
    ovplt_args.add_argument(ovpltthresh_opt[0], **ovpltthresh_opt[1])

    reg_args = parser.add_argument_group('Plotting Region Arguments')
    reg_args.add_argument(numreg_opt[0], default=3, **numreg_opt[1])
    reg_args.add_argument(cntxt_opt[0], **cntxt_opt[1])
    reg_args.add_argument(numstat_opt[0], **numstat_opt[1])

    stat_args = parser.add_argument_group('Statistical Argument')
    stat_args.add_argument(cvgdmp_opt[0], **cvgdmp_opt[1])

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

    stat_args = parser.add_argument_group('Statistical Argument')
    stat_args.add_argument(cvgdmp_opt[0], **cvgdmp_opt[1])

    misc_args, parser = add_misc_args(parser)

    return parser

def get_per_read_roc_parser():
    parser = argparse.ArgumentParser(
        description='Plot per-read ROC curve given known motif(s).',
        add_help=False)
    req_args = parser.add_argument_group('Required Argument')
    req_args.add_argument(prstats_opt[0], required=True, **prstats_opt[1])
    req_args.add_argument(motifdesc_opt[0], required=True, **motifdesc_opt[1])
    req_args.add_argument(fasta_opt[0], required=True, **fasta_opt[1])

    limit_args = parser.add_argument_group('Down-sampling Arguments')
    limit_args.add_argument(spb_opt[0], **spb_opt[1])
    limit_args.add_argument(tsl_opt[0], **tsl_opt[1])

    out_args = parser.add_argument_group('Output Arguments')
    out_args.add_argument(pdf_opt[0],
                         default=OUTPUT_BASE + '.per_reads_roc.pdf',
                         **pdf_opt[1])

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

def get_browser_files_parser():
    parser = argparse.ArgumentParser(
        description='Write wiggle files for specified data types.',
        add_help=False)
    data_args = parser.add_argument_group('Data Arguments')
    data_args.add_argument(fast5dir_opt[0], **fast5dir_opt[1])
    data_args.add_argument(ctrlfast5dir_opt[0], **ctrlfast5dir_opt[1])
    data_args.add_argument(statfn_opt[0], **statfn_opt[1])

    out_args = parser.add_argument_group('Output Arguments')
    out_args.add_argument(brsrfn_opt[0], **brsrfn_opt[1])
    out_args.add_argument(ftypes_opt[0], **ftypes_opt[1])

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

    stat_args = parser.add_argument_group('Statistical Argument')
    stat_args.add_argument(cvgdmp_opt[0], **cvgdmp_opt[1])

    fast5_args, misc_args, parser = add_default_args(parser)

    return parser


if __name__ == '__main__':
    raise NotImplementedError(
        'This is a module. See commands with `tombo -h`')
