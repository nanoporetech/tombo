from __future__ import unicode_literals

###############################
##### Model Name Defaults #####
###############################

# default model names
STANDARD_MODELS = {
    'DNA':'tombo.DNA.model',
    'RNA':'tombo.RNA.200mV.model',
}
ALTERNATE_MODELS = {
    'DNA_5mC':'tombo.DNA.5mC.model',
    'DNA_6mA':'tombo.DNA.6mA.model',
}


################################
##### Re-squiggle Defaults #####
################################

# table containing default segmentation parameters for different sample types
#   1) running neighboring window width for segmentation scoring
#   2) minimum observations per genomic base
#   3) mean number of observations per event during segmentation
SEG_PARAMS_TABLE = {
    'RNA':(8, 4, 10),
    'DNA':(5, 3, 5),
}

# table containing default signal to sequence assignment parameters
# for different sample types
#   1) expected value for matching event to sequence
#   2) penalty for skipped sequence position
#   3) adaptive bandwidth
#   4) signal segmentation mean half-normal score threshold
ALGN_PARAMS_TABLE = {
    'RNA':(4, 10, 1400, 2.0),
    'DNA':(4.2, 4.2, 1200, 1.75),
}

# factor of extra raw signal above minimum to add around skipped bases for
# raw signal segment detection
EXTRA_SIG_FACTOR = 1.1

MASK_FILL_Z_SCORE = -10
MASK_BASES = 50

START_BANDWIDTH = 5000
START_SEQ_WINDOW = 500
BAND_BOUNDARY_THRESH = 5

DEL_FIX_WINDOW = 2
MAX_DEL_FIX_WINDOW = 8
MAX_RAW_CPTS = 200
MIN_EVENT_TO_SEQ_RATIO = 1.1


############################
##### Testing Defaults #####
############################

LLR_THRESH = 0.0
SAMP_COMP_THRESH = 0.1
DE_NOVO_THRESH = 0.5


#####################################
##### Model Estimation Defaults #####
#####################################

ALT_EST_BATCH = 1000
MAX_KMER_OBS = 10000
MIN_KMER_OBS_TO_EST = 50
KERNEL_DENSITY_RANGE = (-5,5)


##########################
##### Misc. Defaults #####
##########################

SMALLEST_PVAL = 1e-50

# got quantiles from analysis of stability after shift-only normalization
ROBUST_QUANTS = (46.5, 53.5)

# minimum standard deviation for a genomic position when estimating spread
# from a sample
MIN_POSITION_SD = 0.01

# number of points at which to estimate the k-mer siganl densities
NUM_DENS_POINTS = 500

# number of reads to estimate global scale parameter
NUM_READS_FOR_SCALE = 1000

# number of points to plot in the ROC curve plotting command
ROC_PLOT_POINTS = 1000
