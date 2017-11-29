GENOME_FASTA='e_coli.fasta'
ALT_MODEL_STATS_FN='e_coli.5mC.tombo.stats'
MODEL_OR_PCR_STATS_FN='e_coli.5mC.pcr_or_model.tombo.stats'
IS_MODEL_STATS=True
COV_THRESH=10
MOTIF='CCWGG'
# second position in motif is the ground truth modified base
MOD_POS=2
NUM_ROC_PLOT_POINTS=1000


# load modules
import h5py

import numpy as np
from numpy.lib import recfunctions

from tombo.tombo_helper import parse_fasta, rev_comp, parse_motif

# parse genome and motif
genome_seq = parse_fasta(GENOME_FASTA)
motif_pat = parse_motif(MOTIF)
motif_len = len(MOTIF)
mod_base = MOTIF[MOD_POS - 1]

# parse valid stat positions from alternative model stats
with h5py.File(ALT_MODEL_STATS_FN) as dat:
    stats = dat['stats'].value
stats = stats[stats['cov'] >= COV_THRESH]
stats = stats[np.logical_or(stats['frac'] > 0, stats['alt_frac'] > 0)]

frac_tot = stats['frac'] + stats['alt_frac']
valid_cov = (frac_tot * stats['cov']).astype(np.int32)
valid_frac = stats['frac'] / frac_tot

new_stats = recfunctions.append_fields(
    stats, ('valid_frac', 'valid_cov'), (valid_frac, valid_cov),
    (np.float64, np.int32), usemask=False)
new_stats = new_stats[new_stats['valid_cov'] >= COV_THRESH]
new_stats.sort(order='valid_frac')

# determine locations that match the motif
vs_alt_stat_has_mod = np.zeros(new_stats.shape[0], dtype=np.int32)
for i, pos_stat in enumerate(new_stats):
    if pos_stat['strand'] == '+':
        stat_seq = genome_seq[pos_stat['chrm']][
            pos_stat['pos'] - MOD_POS + 1:
            pos_stat['pos'] + motif_len - MOD_POS + 1]
    else:
        stat_seq = rev_comp(genome_seq[pos_stat['chrm']][
            pos_stat['pos'] - motif_len + MOD_POS:
            pos_stat['pos'] + MOD_POS])
    if motif_pat.match(stat_seq) is not None:
        vs_alt_stat_has_mod[i] = 1

# compute true and false positive rates
vs_alt_tp_rate = np.cumsum(vs_alt_stat_has_mod)
vs_alt_tp_rate = vs_alt_tp_rate / float(vs_alt_tp_rate[-1])
vs_alt_fp_rate = np.cumsum(np.logical_not(vs_alt_stat_has_mod))
vs_alt_fp_rate = vs_alt_fp_rate / float(vs_alt_fp_rate[-1])

vs_alt_tp_rate = vs_alt_tp_rate[::vs_alt_tp_rate.shape[0] / NUM_ROC_PLOT_POINTS]
vs_alt_fp_rate = vs_alt_fp_rate[::vs_alt_fp_rate.shape[0] / NUM_ROC_PLOT_POINTS]




# now compute true and false positive rates for either a versus standard model or versus PCR sample stats
with h5py.File(MODEL_OR_PCR_STATS_FN) as dat:
   stats = dat['stats'].value
stats = stats[np.logical_and(
    stats['cov'] >= 10,
    np.logical_or(IS_MODEL_STATS, stats['control_cov'] >= 10))]
stats.sort(order='frac')

vs_pcr_or_model_stat_has_mod = []
for i, pos_stat in enumerate(stats):
    if pos_stat['strand'] == '+':
        stat_seq = genome_seq[pos_stat['chrm']][
            pos_stat['pos'] - MOD_POS + 1:
            pos_stat['pos'] + motif_len - MOD_POS + 1]
    else:
        stat_seq = rev_comp(genome_seq[pos_stat['chrm']][
            pos_stat['pos'] - motif_len + MOD_POS:
            pos_stat['pos'] + MOD_POS])
    if motif_pat.match(stat_seq) is not None:
        vs_pcr_or_model_stat_has_mod.append(True)
    elif stat_seq[MOD_POS - 1] == mod_base:
        vs_pcr_or_model_stat_has_mod.append(False)

vs_pcr_or_model_tp_rate = np.cumsum(vs_pcr_or_model_stat_has_mod)
vs_pcr_or_model_tp_rate = vs_pcr_or_model_tp_rate / float(vs_pcr_or_model_tp_rate[-1])
vs_pcr_or_model_fp_rate = np.cumsum(np.logical_not(vs_pcr_or_model_stat_has_mod))
vs_pcr_or_model_fp_rate = vs_pcr_or_model_fp_rate / float(vs_pcr_or_model_fp_rate[-1])

vs_pcr_or_model_tp_rate = vs_pcr_or_model_tp_rate[::vs_pcr_or_model_tp_rate.shape[0] / 1000]
vs_pcr_or_model_fp_rate = vs_pcr_or_model_fp_rate[::vs_pcr_or_model_fp_rate.shape[0] / 1000]

with open('roc_output.txt', 'w') as roc_fp:
    roc_fp.write('TP\tFP\tComparison\n')
    roc_fp.write('\n'.join(
        '{:0.4}\t{:0.4g}\t{}'.format(tpr, fpr, 'Alternative 5mC Model')
        for tpr, fpr in zip(vs_alt_tp_rate, vs_alt_fp_rate)) + '\n')
    roc_fp.write('\n'.join(
        '{:0.4}\t{:0.4g}\t{}'.format(tpr, fpr, 'PCR or Standard Model Comparison')
        for tpr, fpr in zip(vs_pcr_or_model_tp_rate, vs_pcr_or_model_fp_rate)) + '\n')
