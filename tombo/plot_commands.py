from __future__ import division, unicode_literals, absolute_import

from future.utils import native

from builtins import int, range, dict, map, zip

import os
import io
import re
import sys
import queue

import numpy as np
import multiprocessing as mp

from tqdm import tqdm
from time import sleep
from operator import itemgetter
from collections import defaultdict
from itertools import repeat, groupby
from pkg_resources import resource_string

if sys.version_info[0] > 2:
    unicode = str

# import tombo functions
from . import tombo_stats as ts
from . import tombo_helper as th

from ._default_parameters import SMALLEST_PVAL

VERBOSE = False

# quantiles and especially density plots at leat 3 values
# to be meaningful
QUANT_MIN = 3

# plotting names for strands
FWD_STRAND = 'Forward Strand'
REV_STRAND = 'Reverse Strand'

# try global rpy2 imports
try:
    from rpy2 import robjects as r
    from rpy2.robjects.packages import importr
except:
    # pass here and raise error when main functions are actually called
    # in order to give specific error message
    pass

_PROFILE_PLOT_MAX = False


####################
#### ROC Curves ####
####################

def parse_motif_descs(stat_motif_descs):
    parsed_motif_descs = []
    try:
        for motif_desc in stat_motif_descs.split('::'):
            raw_motif, mod_pos, mod_name = motif_desc.split(':')
            motif = th.TomboMotif(raw_motif, int(mod_pos))
            parsed_motif_descs.append((motif, mod_name))
    except:
        th._error_message_and_exit(
            'Invalid motif decriptions format. Format descriptions as: ' +
            '"motif:mod_pos:name[::motif2:mod_pos2:name2...]".')

    return parsed_motif_descs

def plot_roc(stats_fns, motif_descs, fasta_fn, min_reads, pdf_fn,
             cov_damp_counts):
    if len(motif_descs) != len(stats_fns):
        th._error_message_and_exit(
            'Must provide exactly one set of motif descriptions for ' +
            'each statistics file.')

    if VERBOSE: th._status_message('Parsing motifs.')
    motif_descs = [parse_motif_descs(stat_motif_descs)
                   for stat_motif_descs in motif_descs]
    mod_names = [mod_name for stat_mds in motif_descs
                 for _, mod_name in stat_mds]
    if len(mod_names) != len(set(mod_names)):
        th._error_message_and_exit('Modified base names are not unique.')

    if VERBOSE: th._status_message('Parsing genome.')
    genome_index = th.Fasta(fasta_fn)

    if VERBOSE: th._status_message('Computing accuracy statistics.')
    tp_rates, fp_rates, precisions, mod_names_for_r = [], [], [], []
    for stats_fn, stat_motif_descs in zip(stats_fns, motif_descs):
        if not os.path.isfile(stats_fn):
            th._warning_message('Statistics file does not exist. Skipping: ' +
                                stats_fn)
            continue
        stats = ts.TomboStats(stats_fn)
        stats.filter_coverage(min_reads)
        if stats.is_empty():
            th._warning_message(
                'No locations pass coverage threshold. Skipping: ' + stats_fn)
            continue
        stats.order_by_frac(cov_damp_counts)

        for motif, mod_name in stat_motif_descs:
            if (stats.stat_type == ts.ALT_MODEL_TXT and
                next(stats.iter_stat_seqs(
                    genome_index, motif.mod_pos - 1,
                    motif.motif_len - motif.mod_pos))[
                        motif.mod_pos - 1] != motif.mod_base):
                th._warning_message(
                    'Cannot assess modified base accuracy with alternative ' +
                    'model testing to another canonical base. Skipping: ' +
                    mod_name)
                continue
            mod_tp_rate, mod_fp_rate, mod_precision = ts.get_motif_stats(
                motif, stats, genome_index)
            # print auc and average precision
            auc = np.sum(mod_tp_rate[:-1] *
                         (mod_fp_rate[1:] - mod_fp_rate[:-1]))
            # TODO compute precision recall summary stat
            if VERBOSE: sys.stderr.write('\t'.join((
                    '', mod_name.ljust(30), 'AUC:',
                    '{:.4f}'.format(auc))) + '\n')
            tp_rates.extend(mod_tp_rate)
            fp_rates.extend(mod_fp_rate)
            precisions.extend(mod_precision)
            mod_names_for_r.extend(repeat(mod_name, len(mod_tp_rate)))

    if VERBOSE: th._status_message('Plotting.')
    rocDat = r.DataFrame({
        'TP':r.FloatVector(tp_rates),
        'FP':r.FloatVector(fp_rates),
        'Precision':r.FloatVector(precisions),
        'Comparison':r.StrVector(mod_names_for_r)})
    r.r(resource_string(__name__, 'R_scripts/plotROC.R').decode())
    r.r('pdf("' + pdf_fn + '", height=4, width=6)')
    r.globalenv[str('plotROC')](rocDat)
    r.r('dev.off()')

    return

def plot_per_read_roc(
        pr_stats_fns, motif_descs, fasta_fn, pdf_fn,
        stats_per_block, total_stats_limit):
    if len(motif_descs) != len(pr_stats_fns):
        th._error_message_and_exit(
            'Must provide exactly one set of motif descriptions for ' +
            'each statistics file.')

    if VERBOSE: th._status_message('Parsing motifs.')
    motif_descs = [parse_motif_descs(stat_motif_descs)
                   for stat_motif_descs in motif_descs]
    mod_names = [mod_name for stat_mds in motif_descs
                 for _, mod_name in stat_mds]
    if len(mod_names) != len(set(mod_names)):
        th._error_message_and_exit('Modified base names are not unique.')

    if VERBOSE: th._status_message('Parsing genome.')
    genome_index = th.Fasta(fasta_fn)

    if VERBOSE: th._status_message('Extracting per-read statistics.')
    all_motif_stats = {}
    all_motif_stats_for_r = {}
    for pr_stats_fn, stat_motif_descs in zip(pr_stats_fns, motif_descs):
        if not os.path.isfile(pr_stats_fn):
            th._warning_message('Statistics file does not exist. Skipping: ' +
                                pr_stats_fn)
            continue
        pr_stats = ts.PerReadStats(pr_stats_fn)
        for motif, mod_name in stat_motif_descs:
            all_motif_stats[mod_name] = []
        before_bases = max((motif.mod_pos for motif, _ in stat_motif_descs)) - 1
        after_bases = max((motif.motif_len - motif.mod_pos
                           for motif, _ in stat_motif_descs))
        total_num_stats = 0
        for chrm, strand, start, end, block_stats in pr_stats:
            if strand == '+':
                seq_start = max(start - before_bases, 0)
                seq_end = end + after_bases
            else:
                seq_start = max(start - after_bases, 0)
                seq_end = end + before_bases

            reg_seq = genome_index.get_seq(chrm, seq_start, seq_end)
            # randomly sub-sample per-read stats here
            if block_stats.shape[0] > stats_per_block:
                block_stats = block_stats[np.random.choice(
                    block_stats.shape[0], stats_per_block, replace=False)]
            total_num_stats += block_stats.shape[0]
            for r_pos_stat in block_stats:
                # extract position sequence
                if strand == '+':
                    r_pos_seq = reg_seq[
                        r_pos_stat['pos'] - seq_start - before_bases:
                            r_pos_stat['pos'] - seq_start + after_bases + 1]
                else:
                    r_pos_seq = th.rev_comp(reg_seq[
                            r_pos_stat['pos'] - seq_start - after_bases:
                            r_pos_stat['pos'] - seq_start + before_bases + 1])

                # add statistic and whether the sequence matches each motif
                for motif, mod_name in stat_motif_descs:
                    if r_pos_seq[before_bases] != motif.mod_base: continue
                    all_motif_stats[mod_name].append((
                        r_pos_stat['stat'],
                        bool(motif.motif_pat.match(
                            r_pos_seq[before_bases - motif.mod_pos + 1:]))))

            if total_num_stats >= total_stats_limit:
                break

        for mod_name, stats in all_motif_stats.items():
            unzip_stats = list(zip(*stats))
            all_motif_stats_for_r[mod_name] = r.DataFrame({
                'stat':r.FloatVector(unzip_stats[0]),
                'motif_match':r.BoolVector(unzip_stats[1])})

    # python2 rpy2 ListVector can't take unicode keys
    if sys.version_info[0] < 3:
        conv_all_motif_stats_for_r = {}
        for k, v in all_motif_stats_for_r.items():
            conv_all_motif_stats_for_r[k.encode()] = v
        all_motif_stats_for_r = conv_all_motif_stats_for_r
    all_motif_stats_for_r = r.ListVector(all_motif_stats_for_r)

    if VERBOSE: th._status_message('Computing accuracy statistics.')
    tp_rates, fp_rates, precisions, mod_names_for_r = [], [], [], []
    for mod_name, mod_stats in all_motif_stats.items():
        ordered_mod_tf = list(zip(*sorted(mod_stats)))[1]
        mod_tp_rate, mod_fp_rate, mod_precision = ts.compute_accuracy_rates(
            ordered_mod_tf)
        auc = np.sum(mod_tp_rate[:-1] * (mod_fp_rate[1:] - mod_fp_rate[:-1]))
        # TODO compute precision recall summary stat
        if VERBOSE: sys.stderr.write('\t'.join((
                '', mod_name.ljust(30), 'AUC:',
                '{:.4f}'.format(auc))) + '\n')
        tp_rates.extend(mod_tp_rate)
        fp_rates.extend(mod_fp_rate)
        precisions.extend(mod_precision)
        mod_names_for_r.extend(repeat(mod_name, len(mod_tp_rate)))

    rocDat = r.DataFrame({
        'TP':r.FloatVector(tp_rates),
        'FP':r.FloatVector(fp_rates),
        'Precision':r.FloatVector(precisions),
        'Comparison':r.StrVector(mod_names_for_r)})

    if VERBOSE: th._status_message('Plotting.')
    r.r(resource_string(__name__, 'R_scripts/plotROCPerRead.R').decode())
    r.r('pdf("' + pdf_fn + '", height=4, width=6)')
    r.globalenv[str('plotROCPerRead')](rocDat, all_motif_stats_for_r)
    r.r('dev.off()')

    return


###################################
#### K-mer Signal Distribution ####
###################################

def plot_kmer_dist(
        f5_dirs1, corrected_group, basecall_subgroups, pdf_fn,
        read_mean, upstrm_bases, dnstrm_bases, kmer_thresh, num_reads,
        r_struct_fn, dont_plot):
    kmer_width = upstrm_bases + dnstrm_bases + 1
    reads_added = 0
    all_kmers = defaultdict(list)

    raw_read_coverage = th.parse_fast5s(
        f5_dirs1, corrected_group, basecall_subgroups)
    if VERBOSE: th._status_message('Extracting read levels.')
    files = [r_data for cs_r_data in raw_read_coverage.values()
             for r_data in cs_r_data]
    np.random.shuffle(files)
    for r_data in files:
        r_means, r_seq = th.get_multiple_slots_read_centric(
            r_data, ['norm_mean', 'base'])
        if r_means is None: continue
        r_seq = b''.join(r_seq).decode()

        read_kmers = defaultdict(list)
        for kmer, event_mean in zip(
                [r_seq[i:i+kmer_width]
                 for i in range(len(r_seq)-kmer_width+1)],
                r_means[upstrm_bases:]):
            read_kmers[kmer].append(event_mean)
        # if every k-mer is present (unless kmer is greater than 4) and
        # each k-mer has the requested number of occurences
        if kmer_thresh == 0 or (
                len(read_kmers) == 4 ** kmer_width and
                (kmer_thresh == 1 or
                 min(len(x) for x in read_kmers.values()) > kmer_thresh)):
            reads_added += 1
            for kmer, kmer_means in read_kmers.items():
                if read_mean:
                    all_kmers[kmer].append((
                        np.mean(kmer_means), reads_added))
                else:
                    all_kmers[kmer].extend(
                        zip(kmer_means, repeat(reads_added)))

        if reads_added >= num_reads:
            break

    if reads_added in (0,1):
        th._error_message_and_exit(
            'No valid reads present.\n\t\tCheck that [--corrected-group] ' +
            'matches value used in resquiggle.\n\t\tAlso consider lowering ' +
            '[--num-kmer-threshold] especially for k-mer lengths greater than 4.')
    if reads_added < num_reads:
        th._warning_message(
            'Fewer valid reads present than requested.\n\tConsider ' +
            'lowering [--num-kmer-threshold] especially for k-mer lengths ' +
            'greater than 4.')

    if VERBOSE: th._status_message('Preparing plot data.')
    kmer_levels = [kmer for means, kmer in sorted([
        (np.mean(list(map(itemgetter(0), means))), kmer)
        for kmer, means in all_kmers.items()])]

    plot_kmers, plot_bases, plot_means, plot_r_ids = zip(*[
        (kmer, kmer[upstrm_bases], sig_mean, read_i)
        for kmer in kmer_levels
        for sig_mean, read_i in all_kmers[kmer]])

    kmerDat = r.DataFrame({
        'Kmer':r.FactorVector(
            r.StrVector(plot_kmers),
            ordered=True, levels=r.StrVector(kmer_levels)),
        'Base':r.StrVector(plot_bases),
        'Signal':r.FloatVector(plot_means),
        'Read':r.StrVector(plot_r_ids)})
    # df to plot kmers as tile of colors requires gridExtra R package
    try:
        importr(str('gridExtra'))
        baseDat = r.DataFrame({
            'Kmer':r.FactorVector(
                r.StrVector([kmer for kmer in kmer_levels
                             for _ in range(kmer_width)]),
                ordered=True, levels=r.StrVector(kmer_levels)),
            'Base':r.StrVector([kmer[i] for kmer in kmer_levels
                                for i in range(kmer_width)]),
            'Position':r.IntVector([
                i - upstrm_bases for kmer in kmer_levels
                for i in range(kmer_width)])})
    except:
        th._warning_message(
            'Install R package `gridExtra` for ' +
            'visual kmer display. Using text kmer display.')
        baseDat = r.NA_Character

    if r_struct_fn is None:
        r_struct_fn = r.NA_Character
    else:
        r_struct_fn = r.StrVector([r_struct_fn,])
    dont_plot_r = r.BoolVector([dont_plot,])

    if VERBOSE: th._status_message('Plotting.')
    r.r(resource_string(__name__, 'R_scripts/plotKmerDist.R').decode())
    if not dont_plot: r.r('pdf("' + pdf_fn + '", height=7, width=10)')
    if read_mean:
        r.globalenv[str('plotKmerDistWReadPath')](
            kmerDat, baseDat, r_struct_fn, dont_plot_r)
    else:
        r.globalenv[str('plotKmerDist')](
            kmerDat, baseDat, r_struct_fn, dont_plot_r)
    if not dont_plot: r.r('dev.off()')

    return



########################################
#### General data parsing functions ####
########################################

def get_read_correction_data(
        r_data, reg_type, num_obs, region_name=None, start_at_zero=False,
        r_start=None, r_strand=None):
    read_corr_data = th.parse_read_correction_data(r_data)
    if read_corr_data is None:
        return None, None, None, None
    (read_id, signal_data, raw_offset, scale_values, old_segs, old_align_vals,
     new_align_vals, events_end, new_segs) = read_corr_data

    if np.issubdtype(type(native(reg_type)), np.integer):
        if r_strand == '+':
            reg_start = int(new_segs[reg_type - r_start] - (num_obs / 2) - 1)
        else:
            reg_start = int((new_segs[len(new_segs) - (reg_type - r_start) - 1]
                             - num_obs) + (num_obs / 2))
    elif reg_type == 'start':
        reg_start = 0
    elif reg_type == 'end':
        reg_start = events_end - num_obs
    elif reg_type == 'random':
        reg_start = np.random.randint(0, events_end - num_obs)
    else:
        raise NotImplementedError(
            'Invalid reg_type (int or str) to extract read correction data')

    norm_reg_signal, _ = ts.normalize_raw_signal(
        signal_data, raw_offset + reg_start, num_obs, scale_values=scale_values)

    # calculate running difference
    min_seg_len = 4
    sig_cs = np.cumsum(np.insert(norm_reg_signal, 0, 0))
    running_diffs = np.abs((2 * sig_cs[min_seg_len:-min_seg_len]) -
                           sig_cs[:-2*min_seg_len] -
                           sig_cs[2*min_seg_len:])

    # note that I need to check that both new and old segments are
    # in the region as the start of the same genomic position can
    # shift in raw space (i.e. only the old or new position could be
    # in the region of interest)
    old_segs_in_reg = np.where(np.logical_and(
            reg_start <= old_segs, old_segs < reg_start + num_obs))[0]
    old_reg_segs = old_segs[old_segs_in_reg]
    new_segs_in_reg = np.where(np.logical_and(
            reg_start <= new_segs, new_segs < reg_start + num_obs))[0]
    new_reg_segs = new_segs[new_segs_in_reg]

    i_old_segs = iter(old_segs)
    i_new_segs = iter(new_segs)
    align_vals = [((old_b, next(i_old_segs) if old_b != '-' else -1),
                   (new_b, next(i_new_segs) if new_b != '-' else -1))
                  for old_b, new_b in zip(old_align_vals, new_align_vals)]
    reg_align_vals = [
        ((old_b, old_pos, old_pos in old_reg_segs),
         (new_b, new_pos, new_pos in new_reg_segs))
        for (old_b, old_pos), (new_b, new_pos) in align_vals
        if old_pos in old_reg_segs or new_pos in new_reg_segs]

    # summarize alignment for old and new segments
    old_is_del, old_is_mismatch, new_is_ins = [], [], []
    last_was_del = False
    for (old_b, old_pos, old_in_reg), (
            new_b, new_pos, new_in_reg) in reg_align_vals:
        if old_b == '-' and new_in_reg:
            new_is_ins.append(True)
        elif new_b == '-' and old_in_reg:
            old_is_del.append(True)
            old_is_mismatch.append(False)
            last_was_del = True
        else:
            if new_in_reg:
                new_is_ins.append(False)
            if old_in_reg:
                if last_was_del:
                    old_is_del.append(True)
                    last_was_del = False
                else:
                    old_is_del.append(False)
                old_is_mismatch.append(old_b != new_b)

    old_bases, old_reg_segs = [], []
    if (len(reg_align_vals) > 0 and
        sum(map(itemgetter(2), map(itemgetter(0), reg_align_vals))) > 0):
        old_bases, old_reg_segs = zip(*[
            (b, pos) for b, pos, in_reg in map(itemgetter(0), reg_align_vals)
            if in_reg])
    new_bases, new_reg_segs = [], []
    if (len(reg_align_vals) > 0 and
        sum(map(itemgetter(2), map(itemgetter(1), reg_align_vals))) > 0):
        new_bases, new_reg_segs = zip(*[
            (b, pos) for b, pos, in_reg in map(itemgetter(1), reg_align_vals)
            if in_reg])

    # bring positions to zero start if aligning multiple sequences
    sig_range = list(range(reg_start, reg_start + num_obs))
    if start_at_zero:
        old_reg_segs = [
            old_seg_pos - reg_start for old_seg_pos in old_reg_segs]
        new_reg_segs = [
            new_seg_pos - reg_start for new_seg_pos in new_reg_segs]
        sig_range = list(range(0, num_obs))

    old_dat = {
        'Position':r.FloatVector(old_reg_segs),
        'Base':r.StrVector(old_bases),
        'IsDel':r.BoolVector(old_is_del),
        'IsMismatch':r.BoolVector(old_is_mismatch),
        'Read':r.StrVector([read_id for _ in range(len(old_bases))])}
    new_dat = {
        'Position':r.FloatVector(new_reg_segs),
        'Base':r.StrVector(new_bases),
        'IsIns':r.BoolVector(new_is_ins),
        'Read':r.StrVector([read_id for _ in range(len(new_bases))])}
    sig_dat = {
        'Signal':r.FloatVector(norm_reg_signal),
        'Position':r.FloatVector(sig_range),
        'Read':r.StrVector([
            read_id for _ in range(len(norm_reg_signal))])}
    diff_dat = {
        'Signal':r.FloatVector(running_diffs),
        'Position':r.FloatVector(sig_range[
            min_seg_len:len(running_diffs) + min_seg_len]),
        'Read':r.StrVector([
            read_id for _ in range(len(running_diffs))])}
    # add region is applicable
    if region_name is not None:
        old_dat['Region'] = r.StrVector([
            region_name for _ in range(len(old_bases))])
        new_dat['Region'] = r.StrVector([
            region_name for _ in range(len(new_bases))])
        sig_dat['Region'] = r.StrVector([
            region_name for _ in range(len(norm_reg_signal))])
        diff_dat['Region'] = r.StrVector([
            region_name for _ in range(len(running_diffs))])

    old_dat = r.DataFrame(old_dat)
    new_dat = r.DataFrame(new_dat)
    sig_dat = r.DataFrame(sig_dat)
    diff_dat = r.DataFrame(diff_dat)

    return old_dat, new_dat, sig_dat, diff_dat

def get_read_reg_events(r_data, int_start, int_end):
    r_means = th.get_single_slot_genome_centric(r_data, 'norm_mean')
    if r_means is None: return None
    if r_data.start > int_start and r_data.end < int_end:
        # handle reads that are contained in a region
        start_overlap = int_end - r_data.start
        end_overlap = r_data.end - int_start
        # create region with nan values
        region_means = np.empty(int_end - int_start)
        region_means[:] = np.NAN
        region_means[-start_overlap:end_overlap] = r_means[
            -end_overlap:start_overlap]
    elif r_data.start > int_start:
        # handle reads that start in middle of region
        start_overlap = int_end - r_data.start
        # create region with nan values
        region_means = np.empty(int_end - int_start)
        region_means[:] = np.NAN
        region_means[-start_overlap:] = r_means[:start_overlap]
    elif r_data.end < int_end:
        # handle reads that end inside region
        end_overlap = r_data.end - int_start
        # create region with nan values
        region_means = np.empty(int_end - int_start)
        region_means[:] = np.NAN
        region_means[:end_overlap] = r_means[-end_overlap:]
    else:
        region_means = r_means[
            int_start - r_data.start:int_end - r_data.start]

    return region_means

def get_reg_events(reg_reads, int_start, int_end, strand,
                   read_rows=False, num_reads=None):
    reg_events = [
        get_read_reg_events(r_data, int_start, int_end)
        for r_data in reg_reads if strand is None or r_data.strand == strand]
    reg_events = [r_means for r_means in reg_events
                  if r_means is not None]
    if num_reads is not None:
        reg_events = reg_events[:num_reads]

    if read_rows:
        return np.row_stack(reg_events)
    return np.column_stack(reg_events)

def get_r_event_data(
        all_reg_data, plot_types, overplot_thresh, group_num='Group1'):
    Position, Signal, Strand, Region = [], [], [], []
    for reg_plot_sig, reg_data in zip(plot_types, all_reg_data):
        if reg_plot_sig != 'Density': continue

        for strand in ('+', '-'):
            if sum(r_data.strand == strand
                   for r_data in reg_data.reads) == 0:
                continue
            reg_events = get_reg_events(
                reg_data.reads, reg_data.start, reg_data.end, strand)
            for pos, base_read_means in enumerate(reg_events):
                # skip bases with zero or 1 read as ggplot won't
                # be able to estimate the density
                if sum(~np.isnan(base_read_means)) < 2:
                    continue
                # remove nan  regions of reads from partial overlaps
                base_read_means = base_read_means[
                    ~np.isnan(base_read_means)]
                Position.extend(repeat(
                    pos + reg_data.start, base_read_means.shape[0]))
                Signal.extend(base_read_means)
                Strand.extend(repeat(
                    FWD_STRAND if strand == '+' else REV_STRAND,
                    base_read_means.shape[0]))
                Region.extend(repeat(
                    reg_data.reg_id, base_read_means.shape[0]))

    return r.DataFrame({
        'Position':r.IntVector(Position),
        'Signal':r.FloatVector(Signal),
        'Strand':r.FactorVector(
            r.StrVector(Strand),
            ordered=True, levels=r.StrVector((FWD_STRAND, REV_STRAND))),
        'Region':r.StrVector(Region),
        'Group':r.StrVector(list(repeat(group_num, len(Position))))})

def get_r_boxplot_data(
        all_reg_data, plot_types, overplot_thresh, group_num='Group1'):
    (Position, SigMin, Sig25, SigMed, Sig75, SigMax, Strand, Region) = (
        [], [], [], [], [], [], [], [])
    for reg_plot_sig, reg_data in zip(plot_types, all_reg_data):
        if reg_plot_sig != 'Boxplot': continue

        for strand in ('+', '-'):
            if sum(r_data.strand == strand
                   for r_data in reg_data.reads) == 0:
                continue
            reg_events = get_reg_events(
                reg_data.reads, reg_data.start, reg_data.end, strand)
            for pos, base_read_means in enumerate(reg_events):
                # skip regions with no coverage
                if sum(~np.isnan(base_read_means)) == 0:
                    continue
                # remove nan  regions of reads from partial overlaps
                base_read_means = base_read_means[
                    ~np.isnan(base_read_means)]
                Position.append(pos + reg_data.start)
                SigMin.append(np.percentile(base_read_means, 0))
                Sig25.append(np.percentile(base_read_means, 25))
                SigMed.append(np.percentile(base_read_means, 50))
                Sig75.append(np.percentile(base_read_means, 75))
                SigMax.append(np.percentile(base_read_means, 100))
                Strand.append(
                    FWD_STRAND if strand == '+' else REV_STRAND)
                Region.append(reg_data.reg_id)

    return r.DataFrame({
        'Position':r.IntVector(Position),
        'SigMin':r.FloatVector(SigMin),
        'Sig25':r.FloatVector(Sig25),
        'SigMed':r.FloatVector(SigMed),
        'Sig75':r.FloatVector(Sig75),
        'SigMax':r.FloatVector(SigMax),
        'Strand':r.FactorVector(
            r.StrVector(Strand),
            ordered=True, levels=r.StrVector((FWD_STRAND, REV_STRAND))),
        'Region':r.StrVector(Region),
        'Group':r.StrVector(list(repeat(group_num, len(Position))))})

def get_r_quant_data(
        all_reg_data, plot_types, overplot_thresh, group_num='Group1',
        pos_offest=0, pcntls=[1,10,20,30,40,49]):
    upper_pcntls = [100 - pcntl for pcntl in pcntls]
    Position, Lower, Upper, Strand, Region = [], [], [], [], []
    for reg_plot_sig, reg_data in zip(plot_types, all_reg_data):
        if reg_plot_sig != 'Quantile': continue

        for strand in ('+', '-'):
            if sum(r_data.strand == strand
                   for r_data in reg_data.reads) == 0:
                continue
            reg_events = get_reg_events(
                reg_data.reads, reg_data.start, reg_data.end, strand)
            for pos, base_read_means in enumerate(reg_events):
                # skip regions with no coverage
                if sum(~np.isnan(base_read_means)) == 0:
                    continue
                # remove nan regions of reads from partial overlaps
                base_read_means = base_read_means[
                    ~np.isnan(base_read_means)]
                Position.extend(list(repeat(
                    pos + reg_data.start + pos_offest, len(pcntls))))
                Lower.extend(np.percentile(
                    base_read_means, pcntls, interpolation='nearest'))
                Upper.extend(np.percentile(
                    base_read_means, upper_pcntls,
                    interpolation='nearest'))
                Strand.extend(
                    list(repeat(FWD_STRAND if strand == '+' else
                                REV_STRAND, len(pcntls))))
                Region.extend(list(repeat(reg_data.reg_id, len(pcntls))))

    return r.DataFrame({
        'Position':r.FloatVector(Position),
        'Lower':r.FloatVector(Lower),
        'Upper':r.FloatVector(Upper),
        'Strand':r.FactorVector(
            r.StrVector(Strand),
            ordered=True, levels=r.StrVector((FWD_STRAND, REV_STRAND))),
        'Region':r.StrVector(Region),
        'Group':r.StrVector(list(repeat(group_num, len(Position))))})

def get_r_raw_signal_data(
        all_reg_data, plot_types, overplot_thresh, group_num='Group1'):
    not_warned = True
    Position, Signal, Read, Strand, Region = [], [], [], [], []
    for reg_plot_sig, reg_data in zip(plot_types, all_reg_data):
        if not reg_plot_sig in ('Signal', 'Downsample'): continue

        reg_reads = reg_data.reads
        if reg_plot_sig == 'Downsample':
            plus_reads = [r_data for r_data in reg_data.reads
                          if r_data.strand == '+']
            minus_reads = [r_data for r_data in reg_data.reads
                           if r_data.strand == '-']
            # randomly select reads to plot if too many
            if len(plus_reads) > overplot_thresh:
                np.random.shuffle(plus_reads)
                plus_reads = plus_reads[:overplot_thresh]
            if len(minus_reads) > overplot_thresh:
                np.random.shuffle(minus_reads)
                minus_reads = minus_reads[:overplot_thresh]
            reg_reads = plus_reads + minus_reads
        for r_num, r_data in enumerate(reg_reads):
            try:
                (r_sig, overlap_seg_data, start_offset,
                 scale_vals) = th.get_raw_signal(
                    r_data, reg_data.start, reg_data.end)
                r_sig, _ = ts.normalize_raw_signal(
                    r_sig, 0, r_sig.shape[0], scale_values=scale_vals)
            except:
                if not_warned:
                    not_warned = False
                    th._warning_message(
                        'Genome resolved raw signal could not be retrieved ' +
                        'for some reads. Ensure that reads have been ' +
                        're-squiggled and that all data slot corresponding ' +
                        'accordingly.')
                continue
            for base_i, (b_start, b_end) in enumerate(zip(
                    overlap_seg_data[:-1], overlap_seg_data[1:])):
                Position.extend(
                    reg_data.start + base_i + start_offset +
                    np.linspace(0, 1, b_end - b_start, endpoint=False))
                Signal.extend(r_sig[b_start-overlap_seg_data[0]:
                                    b_end-overlap_seg_data[0]])
                Read.extend(list(repeat(
                    unicode(r_num) + '_' + group_num, b_end - b_start)))
                Strand.extend(list(repeat(
                    FWD_STRAND if r_data.strand == '+' else
                    REV_STRAND, b_end - b_start)))
                Region.extend(list(repeat(reg_data.reg_id, b_end - b_start)))

    return r.DataFrame({
        'Position':r.FloatVector(Position),
        'Signal':r.FloatVector(Signal),
        'Read':r.StrVector(Read),
        'Strand':r.FactorVector(
            r.StrVector(Strand),
            ordered=True, levels=r.StrVector((FWD_STRAND, REV_STRAND))),
        'Region':r.StrVector(Region),
        'Group':r.StrVector(list(repeat(group_num, len(Position))))})

def get_plot_types_data(plot_args, quant_offset=0):
    SignalData = get_r_raw_signal_data(*plot_args)
    QuantData = get_r_quant_data(*plot_args, pos_offest=quant_offset)
    BoxData = get_r_boxplot_data(*plot_args)
    EventData = get_r_event_data(*plot_args)

    return SignalData, QuantData, BoxData, EventData

def get_base_r_data(all_reg_data, zero_start=False, is_rna=False):
    BaseStart, Bases, BaseRegion, BaseStrand = [], [], [], []
    for reg_data in all_reg_data:
        # skip regions without sequence data
        if reg_data.seq is None:
            continue
        if reg_data.strand == '+' or reg_data.strand is None:
            for i, base in enumerate(reg_data.seq):
                if is_rna and base == 'T':
                    base = 'U'
                if zero_start:
                    BaseStart.append(unicode(i))
                else:
                    BaseStart.append(unicode(i + reg_data.start))
                Bases.append(base)
                BaseRegion.append(reg_data.reg_id)
                BaseStrand.append(FWD_STRAND)
        if reg_data.strand == '-' or reg_data.strand is None:
            for i, base in enumerate(reg_data.seq):
                base = base.translate(th.COMP_BASES)
                if is_rna and base == 'T':
                    base = 'U'
                if zero_start:
                    BaseStart.append(unicode(i))
                else:
                    BaseStart.append(unicode(i + reg_data.start))
                Bases.append(base)
                BaseRegion.append(reg_data.reg_id)
                BaseStrand.append(REV_STRAND)

    return r.DataFrame({
        'Position':r.FloatVector(BaseStart),
        'Base':r.StrVector(Bases),
        'Region':r.StrVector(BaseRegion),
        'Strand':r.FactorVector(
            r.StrVector(BaseStrand),
            ordered=True, levels=r.StrVector((FWD_STRAND, REV_STRAND)))})


def get_model_r_data(all_reg_model_data):
    Position, Strand, Mean, SD, Region = [], [], [], [], []
    for reg_id, strand, fwd_model_data, rev_model_data in all_reg_model_data:
        if strand == '+' or strand is None:
            for pos, base_model_mean, base_model_sd in fwd_model_data:
                Position.append(pos)
                Strand.append(FWD_STRAND)
                Mean.append(base_model_mean)
                SD.append(base_model_sd)
                Region.append(reg_id)
        if strand == '-' or strand is None:
            for pos, base_model_mean, base_model_sd in rev_model_data:
                Position.append(pos)
                Strand.append(REV_STRAND)
                Mean.append(base_model_mean)
                SD.append(base_model_sd)
                Region.append(reg_id)

    return r.DataFrame({
        'Position':r.FloatVector(Position),
        'Strand':r.FactorVector(
            r.StrVector(Strand),
            ordered=True, levels=r.StrVector((FWD_STRAND, REV_STRAND))),
        'Mean':r.FloatVector(Mean),
        'SD':r.FloatVector(SD),
        'Region':r.StrVector(Region)})

def get_reg_r_stats(all_reg_stats, are_pvals=True):
    Stats, Position, Read, Region = [], [], [], []
    OrdRead, OrdRegion = [], []
    for reg_id, reg_stats in all_reg_stats:
        if are_pvals:
            reg_stats = -np.log10(reg_stats)
        OrdRead.extend(ts.order_reads(reg_stats))
        OrdRegion.extend(repeat(reg_id, reg_stats.shape[0]))
        for read_i, read_stats in enumerate(reg_stats):
            for pos, pos_stat in enumerate(read_stats):
                Stats.append(pos_stat)
                Position.append(pos)
                Read.append(read_i)
                Region.append(reg_id)

    StatsData = r.DataFrame({
        'Stats':r.FloatVector(Stats),
        'Position':r.FloatVector(Position),
        'Read':r.FloatVector(Read),
        'Region':r.StrVector(Region)})
    OrdData = r.DataFrame({
        'Read':r.FloatVector(OrdRead),
        'Region':r.StrVector(OrdRegion)})
    return StatsData, OrdData



###########################################################################
#### Correction plotting (deprecated, but still accessible via script) ####
###########################################################################

def plot_corrections(
        f5_dirs1, corrected_group, basecall_subgroups, pdf_fn,
        reg_type, num_obs, num_reads):
    th._warning_message('The plot_correction command may be deprecated in ' +
                        'future versions of Tombo.')
    if VERBOSE: th._status_message('Preparing plot data.')
    OldSegDat, NewSegDat, SigDat, DiffDat = [], [], [], []
    raw_read_coverage = th.parse_fast5s(
        f5_dirs1, corrected_group, basecall_subgroups)
    files = [r_data for cs_r_data in raw_read_coverage.values()
             for r_data in cs_r_data]
    for r_data in files:
        old_dat, new_dat, signal_dat, diff_dat = get_read_correction_data(
            r_data, reg_type, num_obs)
        if old_dat is None:
            # skip reads that don't have correction slots b/c they
            # couldn't be corrected
            continue
        OldSegDat.append(old_dat)
        NewSegDat.append(new_dat)
        SigDat.append(signal_dat)
        DiffDat.append(diff_dat)
        if len(OldSegDat) >= num_reads:
            break
    if len(OldSegDat) == 0:
        th._error_message_and_exit(
            'No reads were able to be processed. This command is ' +
            'only applicable to reads processed with event_resquiggle. ' +
            'Also check that --corrected-group and --basecall-subgroup ' +
            'match the event_resquiggle command.')
    if VERBOSE and len(OldSegDat) < num_reads:
        th._warning_message(
            'Fewer reads than requested were able to ' +
            'be processed. Likely too few reads provided or ' +
            'those provided were not corrected.')
    OldSegDat = r.DataFrame.rbind(*OldSegDat)
    NewSegDat = r.DataFrame.rbind(*NewSegDat)
    SigDat = r.DataFrame.rbind(*SigDat)
    DiffDat = r.DataFrame.rbind(*DiffDat)

    if VERBOSE: th._status_message('Plotting.')
    r.r(resource_string(__name__, 'R_scripts/plotReadCorr.R').decode())
    r.r('pdf("' + pdf_fn + '", height=7, width=11)')
    r.globalenv[str('plotReadCorr')](OldSegDat, NewSegDat, SigDat, DiffDat)
    r.r('dev.off()')

    return

def plot_multi_corrections(
        f5_dirs1, corrected_group, basecall_subgroups, pdf_fn,
        num_reads_per_plot, num_regions, num_obs, include_orig_bcs, genome_locs):
    th._warning_message('The plot_multi_correction command may be deprecated ' +
                        'in future versions of Tombo.')
    num_regions = num_regions if num_regions % 2 == 0 else \
                  num_regions + 1
    raw_read_coverage = th.parse_fast5s(
        f5_dirs1, corrected_group, basecall_subgroups)
    read_coverage = th.get_coverage(raw_read_coverage)

    if genome_locs is None:
        coverage_regions = []
        for (chrm, strand), cs_coverage in read_coverage.items():
            reg_covs, reg_lens = zip(*[
                (x, len(list(y))) for x, y in groupby(cs_coverage)])
            coverage_regions.extend(zip(
                reg_covs, [start + (reg_len // 2) for start, reg_len in
                           zip(np.cumsum(np.insert(reg_lens, 0, 0)), reg_lens)],
                repeat(chrm), repeat(strand)))

        # randomly select regions with at least num_reads_to_plot regions
        coverage_regions = [
            (chrm, reg_center, strand)
            for stat, reg_center, chrm, strand in
            coverage_regions if stat >= num_reads_per_plot]
        np.random.shuffle(coverage_regions)
        plot_locs = list(zip(
            ['{:03d}'.format(rn) for rn in range(num_regions)],
            coverage_regions[:num_regions]))
        if len(plot_locs) < num_regions:
            th._warning_message(
                'Fewer regions contain minimum ' +
                'number of reads than requested.')
    else:
        if VERBOSE: th._status_message('Parsing genome locations.')
        parsed_locs = th.parse_genome_locations(genome_locs, default_strand='+')
        plot_locs = [
            ('{:03d}'.format(i), (chrm, int(pos) - 1, strand))
            for i, (chrm, pos, strand) in enumerate(parsed_locs)]
        # filter regions with no coverage
        plot_locs = [
            (reg_i, (chrm, start, strand))
            for (reg_i, (chrm, start, strand)) in plot_locs
            if (chrm, strand) in read_coverage and
            read_coverage[(chrm, strand)][start] > 0]
        if len(plot_locs) < len(parsed_locs):
            th._warning_message(
                'Some regions did not contain read coverage.')

    if len(plot_locs) == 0:
        th._error_message_and_exit(
            'No regions contain minimum number of reads.')

    if VERBOSE: th._status_message('Preparing plot data.')
    OldSegDat, NewSegDat, SigDat = [], [], []
    for reg_i, (chrm, reg_center, strand) in plot_locs:
        reg_num_reads = 0
        ## get num_reads_per_region reads from this region
        reg_reads = [
            r_data for r_data in raw_read_coverage[(chrm, strand)]
            if r_data.start <= reg_center - (num_obs / 2.0) and
            r_data.end > reg_center + (num_obs / 2.0) and
            r_data.strand == strand]
        for r_data in reg_reads:
            try:
                old_dat, new_dat, signal_dat, diff_dat \
                    = get_read_correction_data(
                        r_data, reg_center, num_obs, reg_i, True,
                        r_data.start, r_data.strand)
            # some FAST5 files give an error:
            #     "Can't read data (Inflate() failed)"
            except (IOError, KeyError) as e:
                continue
            if old_dat is None:
                # skip reads that don't have correction slots b/c they
                # couldn't be corrected
                continue
            OldSegDat.append(old_dat)
            NewSegDat.append(new_dat)
            SigDat.append(signal_dat)
            reg_num_reads += 1
            if reg_num_reads >= num_reads_per_plot:
                break
        if reg_num_reads < num_reads_per_plot:
            # TODO: figure out if we should warn here
            pass
    try:
        OldSegDat = r.DataFrame.rbind(*OldSegDat)
    except:
        OldSegDat = None
    NewSegDat = r.DataFrame.rbind(*NewSegDat)
    SigDat = r.DataFrame.rbind(*SigDat)

    if VERBOSE: th._status_message('Plotting.')
    r.r(resource_string(__name__, 'R_scripts/plotMultiReadCorr.R').decode())
    r.r('pdf("' + pdf_fn + '", height=5, width=11)')
    if include_orig_bcs and OldSegDat is not None:
        r.globalenv[str('plotMultiReadCorr')](OldSegDat, NewSegDat, SigDat)
    else:
        r.globalenv[str('plotMultiReadCorrNoOrig')](NewSegDat, SigDat)
    r.r('dev.off()')

    return


########################################
#### Base plotting linker functions ####
########################################

def get_plots_titles(all_reg_data, all_reg_data2, overplot_type,
                     overplot_thresh, model_plot=False):
    strand_cov = []
    for reg_i in range(len(all_reg_data)):
        reg_cov1 = [
            sum(r_data.strand == '+' for r_data in all_reg_data[reg_i].reads),
            sum(r_data.strand == '-' for r_data in all_reg_data[reg_i].reads)]
        reg_cov2 = [] if all_reg_data2 is None else [
            sum(r_data.strand == '+' for r_data in all_reg_data2[reg_i].reads),
            sum(r_data.strand == '-' for r_data in all_reg_data2[reg_i].reads)]
        strand_cov.append(reg_cov1 + reg_cov2)
    # downsample can be plotted with any number of reads on either strand
    if overplot_type == "Downsample":
        plot_types = ["Downsample" for _ in strand_cov]
        dnspl_stars = [
            ['*' if grp_r_ovr_cov > overplot_thresh else ''
             for grp_r_ovr_cov in r_cov] for r_cov in strand_cov]
    else:
        # need to remove 0 cov as might be plotting on just one strand
        gt0_strand_cov = [[x for x in covs if x > 0]
                          for covs in strand_cov]
        plot_types = [
            'Signal' if (max(covs) < overplot_thresh or
                         min(covs) < QUANT_MIN)
            else overplot_type for covs in gt0_strand_cov]
        dnspl_stars = [['' for _ in r_cov] for r_cov in strand_cov]

    titles = []
    for int_i, r_cov, r_ovp in zip(all_reg_data, strand_cov, dnspl_stars):
        if all_reg_data2 is None:
            if int_i.strand is None:
                reg_title = int_i.chrm + ' ' + int_i.reg_text + \
                        " ::: Coverage: " + unicode(r_cov[0]) + r_ovp[0] + \
                        " + " + unicode(r_cov[1]) + r_ovp[1] + " -"
            else:
                cov_str = unicode(r_cov[0]) + r_ovp[0] if int_i.strand == '+' \
                          else unicode(r_cov[1]) + r_ovp[1]
                reg_title = int_i.chrm + (
                    ":" + int_i.strand if int_i.strand else '') + \
                    ' ' + int_i.reg_text + " ::: Coverage: " + cov_str
            if model_plot and overplot_type in (
                    'Density', 'Quantile', 'Boxplot'):
                reg_title += ' (Model in Black)'
            titles.append(reg_title)
        else:
            if int_i.strand is None:
                titles.append(
                    int_i.chrm + ' ' + int_i.reg_text +
                    " ::: Coverage: Sample (Red): " +
                    unicode(r_cov[0]) + r_ovp[0] + " + " +
                    unicode(r_cov[1]) + r_ovp[1] + " -; Control (Black): " +
                    unicode(r_cov[2]) + r_ovp[2] + " + " +
                    unicode(r_cov[3]) + r_ovp[3] + " -")
            else:
                cov_str = (
                    'Sample (Red): ' + unicode(r_cov[0]) + r_ovp[0] +
                    '; Control (Black): ' + unicode(r_cov[2]) + r_ovp[2]
                ) if int_i.strand == '+' else (
                    'Sample (Red): ' + unicode(r_cov[1]) + r_ovp[1] +
                    '; Control (Black): ' + unicode(r_cov[3]) + r_ovp[3])
                titles.append(
                    int_i.chrm + ":" + int_i.strand + ' ' + int_i.reg_text +
                    " ::: Coverage: " + cov_str)

    Titles = r.DataFrame({
        'Title':r.StrVector(titles),
        'Region':r.StrVector([int_i.reg_id for int_i in all_reg_data])})

    return Titles, plot_types

def plot_single_sample(
        plot_intervals, raw_read_coverage, overplot_thresh,
        overplot_type, pdf_fn):
    if VERBOSE: th._status_message('Preparing plot data.')
    all_reg_data = th.get_region_reads(plot_intervals, raw_read_coverage)
    if len(all_reg_data) == 0:
        th._error_message_and_exit('No reads in any selected regions.')
    if len(all_reg_data) < len(plot_intervals):
        th._warning_message('Some selected regions contain no reads.')
    rna = th.is_rna(raw_read_coverage)

    Titles, plot_types = get_plots_titles(
        all_reg_data, None, overplot_type, overplot_thresh)

    BasesData = get_base_r_data(all_reg_data, is_rna=rna)
    SignalData, QuantData, BoxData, EventData = get_plot_types_data(
        (all_reg_data, plot_types, overplot_thresh, 'Group1'))

    if VERBOSE: th._status_message('Plotting.')
    r.r(resource_string(__name__, 'R_scripts/plotSingleRun.R').decode())
    r.r('pdf("' + pdf_fn + '", height=5, width=11)')
    r.globalenv[str('plotSingleRun')](SignalData, QuantData, BoxData,
                                      EventData, BasesData, Titles)
    r.r('dev.off()')

    return

def filter_and_merge_group_regs(g1_data, g2_data):
    filt_g1, filt_g2, merged_reg_data, both_no_cov = [], [], [], []
    for r1, r2 in zip(g1_data, g2_data):
        both_reads = r1.reads + r2.reads
        if len(both_reads) > 0:
            merged_reg_data.append(r1._replace(reads=both_reads))
            filt_g1.append(r1)
            filt_g2.append(r2)
        else:
            both_no_cov.append(':'.join(map(str, (
                r1.chrm, unicode(r1.start) + '-' + unicode(r1.end),
                r1.strand))))

    if len(both_no_cov) > 0 and VERBOSE:
        th._warning_message(
            'Some regions include no reads: ' + '\t'.join(both_no_cov))

    if len(merged_reg_data) == 0:
        th._error_message_and_exit('No reads in any selected regions.')

    return merged_reg_data, filt_g1, filt_g2

def plot_two_samples(
        plot_intervals, raw_read_coverage1, raw_read_coverage2,
        overplot_thresh, overplot_type, pdf_fn, seqs_fn=None):
    if VERBOSE: th._status_message('Preparing plot data.')
    # get reads overlapping each region
    all_reg_data1 = th.get_region_reads(
        plot_intervals, raw_read_coverage1, filter_no_cov=False, add_seq=False)
    all_reg_data2 = th.get_region_reads(
        plot_intervals, raw_read_coverage2, filter_no_cov=False, add_seq=False)

    # filter regions with no coverage in either read group
    merged_reg_data, all_reg_data1, all_reg_data2 = filter_and_merge_group_regs(
        all_reg_data1, all_reg_data2)
    if len(merged_reg_data) < len(plot_intervals):
        th._warning_message('Some selected regions contain no reads.')

    Titles, plot_types = get_plots_titles(
        all_reg_data1, all_reg_data2, overplot_type, overplot_thresh)

    merged_reg_data = th.add_reg_seq(merged_reg_data)
    BasesData = get_base_r_data(merged_reg_data)

    # get plotting data for either quantiles of raw signal
    SignalData1, QuantData1, BoxData1, EventData1 = get_plot_types_data(
        (all_reg_data1, plot_types, overplot_thresh, 'Group1'), 0.1)
    SignalData2, QuantData2, BoxData2, EventData2 = get_plot_types_data(
        (all_reg_data2, plot_types, overplot_thresh, 'Group2'), 0.5)

    if VERBOSE: th._status_message('Plotting.')
    r.r(resource_string(__name__, 'R_scripts/plotGroupComp.R').decode())
    r.r('pdf("' + pdf_fn + '", height=5, width=11)')
    r.globalenv[str('plotGroupComp')](
        r.DataFrame.rbind(SignalData1, SignalData2),
        r.DataFrame.rbind(QuantData1, QuantData2),
        r.DataFrame.rbind(BoxData1, BoxData2),
        r.DataFrame.rbind(EventData1, EventData2),
        BasesData, Titles, 0.4)
    r.r('dev.off()')

    if seqs_fn is not None:
        if VERBOSE: th._status_message('Outputting region seqeuences.')
        with io.open(seqs_fn, 'wt') as seqs_fp:
            for int_i in merged_reg_data:
                # get the interval from the base data struct
                reg_seq = int_i.seq if int_i.strand == '+' else th.rev_comp(
                    int_i.seq)
                seqs_fp.write('>{0}::{1:d}::{2} {3}\n{4}\n'.format(
                    int_i.chrm, int_i.start, int_i.strand, int_i.reg_text,
                    ''.join(reg_seq)))

    return

def get_reg_kmers(tb_model_fn, plot_intervals, raw_read_coverage,
                  min_reg_overlap=None, alt_model_fn=None):
    def get_reg_reads(reads, int_start, int_end):
        """ Filter reads obtained from expanded interval
        """
        return [r_data for r_data in reads
                if not (r_data.start >= int_end or r_data.end <= int_start)]
    std_ref = ts.TomboModel(tb_model_fn)
    # compute kmer values to make strand specific calculations easier
    dnstrm_bases = std_ref.kmer_width - std_ref.central_pos - 1
    expand_width = max(std_ref.central_pos, dnstrm_bases)
    filt_width = expand_width if min_reg_overlap is None else \
                 expand_width + min_reg_overlap
    if alt_model_fn is not None:
        alt_ref = ts.TomboModel(alt_model_fn)
        if (alt_ref.central_pos != std_ref.central_pos or
            alt_ref.kmer_width != alt_ref.kmer_width):
            th._error_message_and_exit(
                'Standard model not based on the same kmer position ' +
                'as alternative model.')

    # expand regions to get kmers at first and last positions
    expanded_intervals = [p_int._replace(start=p_int.start - expand_width,
                                         end=p_int.end + expand_width)
                          for p_int in plot_intervals]
    # get reads and region sequence
    expanded_intervals = th.get_region_reads(
        expanded_intervals, raw_read_coverage)
    expand_seqs = [int_i.seq for int_i in expanded_intervals]
    rev_expand_seqs = [th.rev_comp(int_i.seq) for int_i in expanded_intervals]
    # convert back to original plot_intervals with seq from exanded intervals
    all_reg_data = [
        int_i._replace(start=int_i.start + expand_width,
                       end=int_i.end - expand_width,
                       reads=get_reg_reads(int_i.reads, int_i.start + filt_width,
                                           int_i.end - filt_width),
                       seq=int_i.seq[expand_width:-expand_width])
        for int_i in expanded_intervals]

    all_reg_model_data, all_reg_alt_model_data = [], []
    for reg_data, reg_seq, rev_seq in zip(
            all_reg_data, expand_seqs, rev_expand_seqs):
        clipped_reg_seq = reg_seq
        clipped_rev_seq = rev_seq
        if std_ref.central_pos > dnstrm_bases:
            clipped_reg_seq = reg_seq[:dnstrm_bases-std_ref.central_pos]
            clipped_rev_seq = rev_seq[:dnstrm_bases-std_ref.central_pos]
        elif dnstrm_bases > std_ref.central_pos:
            clipped_reg_seq = reg_seq[dnstrm_bases-std_ref.central_pos:]
            clipped_rev_seq = rev_seq[dnstrm_bases-std_ref.central_pos:]
        fwd_kmers = [
            clipped_reg_seq[i:i + std_ref.kmer_width]
            for i in range(len(clipped_reg_seq) - std_ref.kmer_width + 1)]
        rev_kmers = [
            clipped_rev_seq[i:i + std_ref.kmer_width]
            for i in range(len(clipped_rev_seq) - std_ref.kmer_width + 1)]
        all_reg_model_data.append((
            reg_data.reg_id, reg_data.strand,
            [(reg_data.start + pos, std_ref.means[kmer],
              std_ref.sds[kmer])
             for pos, kmer in enumerate(fwd_kmers) if not th.invalid_seq(kmer)],
            [(reg_data.end - pos - 1, std_ref.means[kmer],
              std_ref.sds[kmer]) for pos, kmer in enumerate(rev_kmers)
             if not th.invalid_seq(kmer)]))
        # if alternative model is supplied add info
        if alt_model_fn is not None:
            all_reg_alt_model_data.append((
                reg_data.reg_id, reg_data.strand,
                [(reg_data.start + pos, alt_ref.means[kmer],
                  alt_ref.sds[kmer])
                 for pos, kmer in enumerate(fwd_kmers)
                 if not th.invalid_seq(kmer)],
                [(reg_data.end - pos - 1, alt_ref.means[kmer],
                  alt_ref.sds[kmer])
                 for pos, kmer in enumerate(rev_kmers)
                 if not th.invalid_seq(kmer)]))

    return all_reg_data, all_reg_model_data, all_reg_alt_model_data

def plot_motif_centered_with_stats(
        raw_read_coverage1, raw_read_coverage2, plot_intervals,
        stat_locs, overplot_thresh, pdf_fn, tb_model_fn, alt_model_fn=None):
    if VERBOSE: th._status_message('Preparing plot data.')

    ModelData = r.r('NULL')
    if raw_read_coverage2 is None:
        if tb_model_fn is None:
            merged_reg_data = th.get_region_reads(
                plot_intervals, raw_read_coverage1, filter_no_cov=False)
            plot_types = ['Downsample' for _ in merged_reg_data]
            SignalData, _, _, _ = get_plot_types_data(
                (merged_reg_data, plot_types, overplot_thresh, 'Group1'))
        else:
            (merged_reg_data, all_reg_model_data,
             all_reg_alt_model_data) = get_reg_kmers(
                 tb_model_fn, plot_intervals, raw_read_coverage1,
                 alt_model_fn=alt_model_fn)
            plot_types = ['Downsample' for _ in merged_reg_data]
            SignalData, _, _, _ = get_plot_types_data(
                (merged_reg_data, plot_types, overplot_thresh, 'Group1'))
            ModelData = get_model_r_data(all_reg_model_data)
            if alt_model_fn is not None:
                AltModelData = get_model_r_data(all_reg_alt_model_data)
    else:
        all_reg_data1 = th.get_region_reads(
            plot_intervals, raw_read_coverage1, filter_no_cov=False,
            add_seq=False)
        all_reg_data2 = th.get_region_reads(
            plot_intervals, raw_read_coverage2, filter_no_cov=False,
            add_seq=False)

        (merged_reg_data, all_reg_data1,
         all_reg_data2) = filter_and_merge_group_regs(
             all_reg_data1, all_reg_data2)
        plot_types = ['Downsample' for _ in merged_reg_data]
        merged_reg_data = th.add_reg_seq(merged_reg_data)

        # sigDat lists
        SignalData1, _, _, _ = get_plot_types_data(
            (all_reg_data1, plot_types, overplot_thresh, 'Group1'))
        SignalData2, _, _, _ = get_plot_types_data(
            (all_reg_data2, plot_types, overplot_thresh, 'Group2'))
        SignalData = r.DataFrame.rbind(SignalData1, SignalData2)

    BasesData = get_base_r_data(merged_reg_data)

    plot_poss, plot_stats = zip(*stat_locs)
    # stat lists
    StatsData = r.DataFrame({
        'Position':r.FloatVector(plot_poss),
        'Stat':r.FloatVector(plot_stats)})

    if VERBOSE: th._status_message('Plotting.')
    r.r(resource_string(__name__, 'R_scripts/plotMotifStats.R').decode())
    r.r('pdf("' + pdf_fn + '", height=5, width=8)')
    if alt_model_fn is None:
        r.globalenv[str('plotMotifStats')](
            SignalData, BasesData, StatsData, ModelData)
    else:
        r.globalenv[str('plotMotifStats')](
            SignalData, BasesData, StatsData, ModelData, AltModelData)
    r.r('dev.off()')

    return

def plot_model_single_sample(
        plot_intervals, raw_read_coverage, tb_model_fn,
        overplot_type, overplot_thresh, pdf_fn, alt_model_fn=None, seqs_fn=None):
    if VERBOSE: th._status_message('Preparing plot data.')
    # get reads overlapping each region along with all kmers
    all_reg_data, all_reg_model_data, all_reg_alt_model_data = get_reg_kmers(
        tb_model_fn, plot_intervals, raw_read_coverage,
        alt_model_fn=alt_model_fn)
    if len(all_reg_data) == 0:
        th._error_message_and_exit('No reads in any selected regions.')
    if len(all_reg_data) < len(plot_intervals):
        th._warning_message('Some selected regions contain no reads.')
    rna = th.is_rna(raw_read_coverage)

    Titles, plot_types = get_plots_titles(
        all_reg_data, None, overplot_type, overplot_thresh, True)
    ModelData = get_model_r_data(all_reg_model_data)
    if alt_model_fn is not None:
        AltModelData = get_model_r_data(all_reg_alt_model_data)
    BasesData = get_base_r_data(all_reg_data, is_rna=rna)
    SignalData, QuantData, BoxData, EventData = get_plot_types_data(
        (all_reg_data, plot_types, overplot_thresh, 'Group1'))

    if VERBOSE: th._status_message('Plotting.')
    r.r(resource_string(__name__, 'R_scripts/plotModelComp.R').decode())
    r.r('pdf("' + pdf_fn + '", height=5, width=11)')
    if alt_model_fn is None:
        r.globalenv[str('plotModelComp')](
            SignalData, QuantData, BoxData, EventData,
            BasesData, Titles, ModelData)
    else:
        r.globalenv[str('plotModelComp')](
            SignalData, QuantData, BoxData, EventData,
            BasesData, Titles, ModelData, AltModelData)
    r.r('dev.off()')

    if seqs_fn is not None:
        if VERBOSE: th._status_message('Outputting region seqeuences.')
        with io.open(seqs_fn, 'wt') as seqs_fp:
            for int_i in all_reg_data:
                reg_seq = int_i.seq if int_i.strand == '+' else th.rev_comp(
                    int_i.seq)
                seqs_fp.write('>{0}::{1:d}::{2} {3}\n{4}\n'.format(
                    int_i.chrm, int_i.start, int_i.strand, int_i.reg_text,
                    ''.join(reg_seq)))

    return

def plot_per_read_modification(
        all_reg_data, all_reg_stats, are_pvals, box_center, pdf_fn):
    if VERBOSE: th._status_message('Preparing plot data.')
    StatData, OrdData = get_reg_r_stats(all_reg_stats, are_pvals)
    BasesData = get_base_r_data(all_reg_data, zero_start=True)

    if VERBOSE: th._status_message('Plotting.')
    r.r(resource_string(__name__, 'R_scripts/plotPerReadStats.R').decode())
    r.r('pdf("' + pdf_fn + '", height=5, width=11)')
    r.globalenv[str('plotPerReadStats')](
        StatData, OrdData, BasesData, box_center, are_pvals)
    r.r('dev.off()')

    return


#################################
#### Plot processing methods ####
#################################

def get_valid_model_fns(
        tb_model_fn, plot_default_stnd, alt_model_fn,
        plot_default_alt, raw_read_coverage, f5_dirs2=None):
    # if no model was requested
    if (tb_model_fn is None and not plot_default_stnd and
        alt_model_fn is None and not plot_default_alt):
        return None, None

    if tb_model_fn is None:
        tb_model_fn, _ = ts.get_default_standard_ref(raw_read_coverage)
    if alt_model_fn is None and plot_default_alt is not None:
        alt_model_fn, _ = ts.get_default_alt_ref(
            plot_default_alt, raw_read_coverage)

    if f5_dirs2 is not None and tb_model_fn is not None:
        th._warning_message(
            'Both a second set of FAST5s and a tombo model were ' +
            'provided. Two samples with model plotting is not ' +
            'currently available. Models requested will be ignored.')

    return tb_model_fn, alt_model_fn

def plot_max_coverage(
        f5_dirs1, corrected_group, basecall_subgroups, pdf_fn,
        f5_dirs2, num_regions, num_bases, overplot_thresh, overplot_type,
        tb_model_fn, alt_model_fn, plot_default_stnd, plot_default_alt):
    raw_read_coverage = th.parse_fast5s(
        f5_dirs1, corrected_group, basecall_subgroups)

    tb_model_fn, alt_model_fn = get_valid_model_fns(
        tb_model_fn, plot_default_stnd, alt_model_fn, plot_default_alt,
        raw_read_coverage, f5_dirs2)
    if f5_dirs2 is None:
        coverage_regions = []
        for chrm, strand, cs_cov, cs_cov_starts in th.get_coverage_regions(
                raw_read_coverage):
            coverage_regions.extend(zip(
                cs_cov, cs_cov_starts, repeat(chrm), repeat(strand)))

        # max coverage plots both strands coverage
        plot_intervals = [
            th.intervalData('{:03d}'.format(rn), chrm, start, start + num_bases)
            for rn, (stat, start, chrm, strand) in
            enumerate(sorted(coverage_regions, reverse=True)[:num_regions])]

        if tb_model_fn is None:
            plot_single_sample(
                plot_intervals, raw_read_coverage, overplot_thresh,
                overplot_type, pdf_fn)
        else:
            plot_model_single_sample(
                plot_intervals, raw_read_coverage, tb_model_fn,
                overplot_type, overplot_thresh, pdf_fn, alt_model_fn)
    else:
        raw_read_coverage2 = th.parse_fast5s(
            f5_dirs2, corrected_group, basecall_subgroups)
        coverage_regions = []
        # only process chromosomes in both read groups
        for chrm, strand, cs_cov, cs_cov_starts in th.get_coverage_regions(
                raw_read_coverage, raw_read_coverage2):
            coverage_regions.extend(zip(
                cs_cov, cs_cov_starts, repeat(chrm), repeat(strand)))

        # max coverage plots both strands coverage
        plot_intervals = [
            th.intervalData('{:03d}'.format(rn), chrm, start, start + num_bases)
            for rn, (stat, start, chrm, strand) in
            enumerate(sorted(coverage_regions, reverse=True)[:num_regions])]

        plot_two_samples(
            plot_intervals, raw_read_coverage, raw_read_coverage2,
            overplot_thresh, overplot_type, pdf_fn)

    return

if _PROFILE_PLOT_MAX:
    _plot_max_wrapper = plot_max_coverage
    def plot_max_coverage(*args, **kwargs):
        import cProfile
        cProfile.runctx('_plot_max_wrapper(*args, **kwargs)', globals(), locals(),
                        filename='plot_max_cov.prof')
        return

def plot_genome_locations(
        f5_dirs1, corrected_group, basecall_subgroups, pdf_fn,
        f5_dirs2, num_bases, overplot_thresh, overplot_type,
        genome_locs, tb_model_fn, alt_model_fn, plot_default_stnd,
        plot_default_alt):
    if VERBOSE: th._status_message('Parsing genome locations.')
    # minus one here as all python internal coords are 0-based, but
    # genome is generally 1-based
    plot_intervals = []
    for i, (chrm, pos, strand) in enumerate(th.parse_genome_locations(
            genome_locs)):
        int_start = max(
            0, int(int(pos) - np.floor(num_bases / 2.0) - 1))
        plot_intervals.append(th.intervalData(
            '{:03d}'.format(i), chrm, int_start, int_start + num_bases, strand))

    raw_read_coverage = th.parse_fast5s(
        f5_dirs1, corrected_group, basecall_subgroups)
    tb_model_fn, alt_model_fn = get_valid_model_fns(
        tb_model_fn, plot_default_stnd, alt_model_fn, plot_default_alt,
        raw_read_coverage, f5_dirs2)

    if f5_dirs2 is None:
        if tb_model_fn is None:
            plot_single_sample(
                plot_intervals, raw_read_coverage, overplot_thresh,
                overplot_type, pdf_fn)
        else:
            plot_model_single_sample(
                plot_intervals, raw_read_coverage, tb_model_fn,
                overplot_type, overplot_thresh, pdf_fn, alt_model_fn)
    else:
        raw_read_coverage2 = th.parse_fast5s(
            f5_dirs2, corrected_group, basecall_subgroups)
        plot_two_samples(
            plot_intervals, raw_read_coverage, raw_read_coverage2,
            overplot_thresh, overplot_type, pdf_fn)

    return

def plot_per_read_mods_genome_location(
        f5_dirs, corrected_group, basecall_subgroups, pdf_fn,
        per_read_stats_fn, genome_locs, num_bases, num_reads, box_center,
        fasta_fn):
    if VERBOSE: th._status_message('Parsing genome locations.')
    plot_intervals = []
    for i, (chrm, pos, strand) in enumerate(th.parse_genome_locations(
            genome_locs, default_strand='+')):
        int_start = max(
            0, int(int(pos) - np.floor(num_bases / 2.0) - 1) + 1)
        plot_intervals.append(th.intervalData(
            '{:03d}'.format(i), chrm, int_start, int_start + num_bases, strand))

    # add sequence to each region if fast5s or fasta are provided
    if fasta_fn is not None:
        genome_index = th.Fasta(fasta_fn)
        plot_intervals_w_seq = []
        for int_data in plot_intervals:
            plot_intervals_w_seq.append(
                int_data._replace(seq=genome_index.get_seq(
                    int_data.chrm, int_data.start, int_data.end)))
        plot_intervals = plot_intervals_w_seq
    elif f5_dirs is not None:
        raw_read_coverage = th.parse_fast5s(
            f5_dirs, corrected_group, basecall_subgroups)
        plot_intervals = th.get_region_reads(plot_intervals, raw_read_coverage)
    else:
        th._warning_message(
            'No read FAST5 directory or genome FASTA file provided. ' +
            'Plotting without sequence.')

    if VERBOSE: th._status_message('Parsing per read statistics.')
    per_read_stats = ts.PerReadStats(per_read_stats_fn)
    interval_stats = []
    for int_data in plot_intervals:
        int_stats = per_read_stats.get_region_per_read_stats(int_data, num_reads)
        if int_stats is not None:
            # convert long form stats to matrix form (so they can be clustered)
            # regular sort doesn't seem to work for string (object) types
            # so using argsort
            int_stats = int_stats[np.argsort(int_stats['read_id'])]
            # use interval data instead of stats dimensions since regDat is
            # used to compute some window distances in R, so it must be full
            # matrix for the region with NAs
            int_len = int_data.end - int_data.start + 1
            all_read_stats = np.split(
                int_stats, np.where(int_stats['read_id'][:-1] !=
                                    int_stats['read_id'][1:])[0] + 1)
            read_stats_mat = np.empty((len(all_read_stats), int_len))
            read_stats_mat[:] = np.NAN
            for read_i, read_int_stats in enumerate(all_read_stats):
                #print(read_i, read_int_stats['read_id'][0])
                np.put(read_stats_mat[read_i,:],
                       read_int_stats['pos'] - int_data.start,
                       read_int_stats['stat'])
            interval_stats.append((int_data.reg_id, read_stats_mat))

    are_pvals = per_read_stats.are_pvals
    per_read_stats.close()

    plot_per_read_modification(
        plot_intervals, interval_stats, are_pvals, box_center, pdf_fn)

    return

def plot_motif_centered(
        f5_dirs1, corrected_group, basecall_subgroups, pdf_fn,
        f5_dirs2, num_regions, num_bases, overplot_thresh, overplot_type,
        motif, fasta_fn, deepest_coverage, tb_model_fn, alt_model_fn,
        plot_default_stnd, plot_default_alt):
    if VERBOSE: th._status_message('Identifying genomic k-mer locations.')
    genome_index = th.Fasta(fasta_fn)
    motif = th.TomboMotif(motif)

    def get_motif_locs(covered_chrms):
        motif_locs = []
        for chrm in genome_index.iter_chrms():
            if chrm not in covered_chrms: continue
            seq = genome_index.get_seq(chrm)
            for motif_loc in motif.motif_pat.finditer(seq):
                motif_locs.append((chrm, motif_loc.start(), '+'
                                   if not motif.is_palindrome else None))
            # search over negative strand as well if not palindromic
            if not motif.is_palindrome:
                for motif_loc in motif.rev_comp_pat.finditer(seq):
                    motif_locs.append((chrm, motif_loc.start(), '-'))

        if len(motif_locs) == 0:
            th._error_message_and_exit(
                'Motif (' + motif.raw_motif + ') not found in genome.')
        elif len(motif_locs) < num_regions:
            th._warning_message(
                'Motif (' + motif.raw_motif + ') only found ' +
                unicode(len(motif_locs)) + ' times in genome.')
            num_region = len(motif_locs)
        np.random.shuffle(motif_locs)

        return motif_locs

    def get_pos_cov(chrm, pos, strand, read_coverage, read_coverage2=None):
        def get_strand_cov(cov_strand):
            try:
                if read_coverage2 is None:
                    return read_coverage[(chrm, cov_strand)][pos]
                else:
                    return min(read_coverage[(chrm, cov_strand)][pos],
                               read_coverage2[(chrm, cov_strand)][pos])
            except (IndexError, KeyError):
                return 0

        # if strand is not specified get max coverage over both strands
        if strand is None:
            return max(get_strand_cov('+'), get_strand_cov('-'))
        # else get coverage for strand with motif
        return get_strand_cov(strand)


    raw_read_coverage = th.parse_fast5s(
        f5_dirs1, corrected_group, basecall_subgroups)
    tb_model_fn, alt_model_fn = get_valid_model_fns(
        tb_model_fn, plot_default_stnd, alt_model_fn, plot_default_alt,
        raw_read_coverage, f5_dirs2)

    if deepest_coverage:
        read_coverage = th.get_coverage(raw_read_coverage)
    if f5_dirs2 is None:
        covered_chrms = set(map(itemgetter(0), raw_read_coverage))
        # filter out motif_locs to chromosomes not covered
        motif_locs = get_motif_locs(covered_chrms)

        if deepest_coverage:
            if VERBOSE: th._status_message('Finding deepest coverage regions.')
            motif_locs_cov = sorted([
                (get_pos_cov(chrm, pos, strand, read_coverage),
                 chrm, pos, strand)
                for chrm, pos, strand in motif_locs], reverse=True)
            plot_intervals = []
            for i, (cov, chrm, pos, strand) in enumerate(motif_locs_cov):
                int_start = max(
                    0, pos - int((num_bases - motif.motif_len + 1) / 2.0))
                int_end = int_start + num_bases
                plot_intervals.append(th.intervalData(
                    '{:03d}'.format(i), chrm, int_start, int_end, strand))
                if len(plot_intervals) >= num_regions: break
        # plot random covered regions
        else:
            # iterate over regions and check if they have any coverage
            plot_intervals = []
            for i, (chrm, pos, strand) in enumerate(motif_locs):
                int_start = max(
                    0, pos - int((num_bases - motif.motif_len + 1) / 2.0))
                int_end = int_start + num_bases
                if strand is None and any(
                        ((chrm, s) in raw_read_coverage and
                         any(r_data.start < pos < r_data.end
                             for r_data in raw_read_coverage[(chrm, s)]))
                        for s in ('+', '-')):
                    plot_intervals.append(th.intervalData(
                        '{:03d}'.format(i), chrm, int_start, int_end, strand))
                elif ((chrm, strand) in raw_read_coverage and
                      any(r_data.start < pos < r_data.end
                          for r_data in raw_read_coverage[(chrm, strand)])):
                    plot_intervals.append(th.intervalData(
                        '{:03d}'.format(i), chrm, int_start, int_end, strand))
                if len(plot_intervals) >= num_regions: break

        if tb_model_fn is None:
            plot_single_sample(
                plot_intervals, raw_read_coverage, overplot_thresh,
                overplot_type, pdf_fn)
        else:
            plot_model_single_sample(
                plot_intervals, raw_read_coverage, tb_model_fn,
                overplot_type, overplot_thresh, pdf_fn, alt_model_fn)
    # two sample plot
    else:
        raw_read_coverage2 = th.parse_fast5s(
            f5_dirs2, corrected_group, basecall_subgroups)

        covered_chrms = set(map(itemgetter(0), raw_read_coverage)).intersection(
            map(itemgetter(0), raw_read_coverage2))
        # filter out motif_locs to chromosomes not covered
        motif_locs = get_motif_locs(covered_chrms)

        if deepest_coverage:
            read_coverage2 = th.get_coverage(raw_read_coverage2)
            if VERBOSE: th._status_message('Finding deepest coverage regions.')
            motif_locs_cov = sorted([
                (get_pos_cov(chrm, pos, strand, read_coverage, read_coverage2),
                 chrm, pos, strand)
                for chrm, pos, strand in motif_locs], reverse=True)
            if motif_locs_cov[0][0] == 0:
                th._error_message_and_exit(
                    'Motif not covered by both groups at any positions.')

            plot_intervals = []
            for i, (cov, chrm, pos, strand) in enumerate(motif_locs_cov):
                int_start = max(
                    0, pos - int((num_bases - motif.motif_len + 1) / 2.0))
                plot_intervals.append(th.intervalData(
                    '{:03d}'.format(i), chrm, int_start,
                    int_start + num_bases, strand))
                if len(plot_intervals) >= num_regions: break
        # plot random covered regions
        else:
            # iterate over regions and check if they have any coverage
            plot_intervals = []
            for i, (chrm, pos, strand) in enumerate(motif_locs):
                int_start = max(
                    0, pos - int((num_bases - motif.motif_len + 1) / 2.0))
                int_end = int_start + num_bases
                if strand is None and any((
                        (chrm, s) in raw_read_coverage and
                        (chrm, s) in raw_read_coverage2 and
                        any(r_data.start < pos < r_data.end
                            for r_data in raw_read_coverage[(chrm, s)]) and
                        any(r_data2.start < pos < r_data2.end
                            for r_data2 in raw_read_coverage2[(chrm, s)]))
                                          for s in ('+', '-')):
                    plot_intervals.append(th.intervalData(
                        '{:03d}'.format(i), chrm, int_start, int_end, strand))
                elif ((chrm, strand) in raw_read_coverage and
                      (chrm, strand) in raw_read_coverage2 and
                      any(r_data.start < pos < r_data.end
                          for r_data in raw_read_coverage[(chrm, strand)]) and
                      any(r_data2.start < pos < r_data2.end
                          for r_data2 in raw_read_coverage2[(chrm, strand)])):
                    plot_intervals.append(th.intervalData(
                        '{:03d}'.format(i), chrm, int_start, int_end, strand))

                if len(plot_intervals) >= num_regions: break

            if len(plot_intervals) == 0:
                th._error_message_and_exit(
                    'Motif not covered by both groups at any positions.')

        plot_two_samples(
            plot_intervals, raw_read_coverage, raw_read_coverage2,
            overplot_thresh, overplot_type, pdf_fn)

    return

def plot_max_diff(
        f5_dirs1, corrected_group, basecall_subgroups, pdf_fn,
        f5_dirs2, num_regions, num_bases, overplot_thresh, overplot_type,
        seqs_fn):
    raw_read_coverage1 = th.parse_fast5s(
        f5_dirs1, corrected_group, basecall_subgroups)
    raw_read_coverage2 = th.parse_fast5s(
        f5_dirs2, corrected_group, basecall_subgroups)

    if VERBOSE: th._status_message('Getting largest mean signal differences.')
    plot_intervals = [
        th.intervalData(
            '{:03d}'.format(rn), chrm, start, start + num_bases, strand,
            '(Mean diff: {:.2f})'.format(stat))
        for rn, (stat, start, chrm, strand) in
        enumerate(th.get_largest_signal_differences(
            raw_read_coverage1, raw_read_coverage2, num_regions, num_bases))]

    plot_two_samples(
        plot_intervals, raw_read_coverage1, raw_read_coverage2,
        overplot_thresh, overplot_type, pdf_fn, seqs_fn)

    return

def plot_most_signif(
        f5_dirs1, corrected_group, basecall_subgroups, pdf_fn,
        f5_dirs2, num_regions, overplot_thresh, seqs_fn, num_bases,
        overplot_type, stats_fn, tb_model_fn, alt_model_fn,
        plot_default_stnd, plot_default_alt, cov_damp_counts):
    if VERBOSE: th._status_message('Loading statistics from file.')
    plot_intervals = ts.TomboStats(stats_fn).get_most_signif_regions(
        num_bases, num_regions, cov_damp_counts=cov_damp_counts)

    raw_read_coverage = th.parse_fast5s(
        f5_dirs1, corrected_group, basecall_subgroups)
    tb_model_fn, alt_model_fn = get_valid_model_fns(
        tb_model_fn, plot_default_stnd, alt_model_fn, plot_default_alt,
        raw_read_coverage, f5_dirs2)

    if f5_dirs2 is None:
        if tb_model_fn is None:
            plot_single_sample(
                plot_intervals, raw_read_coverage, overplot_thresh,
                overplot_type, pdf_fn)
        else:
            plot_model_single_sample(
                plot_intervals, raw_read_coverage, tb_model_fn,
                overplot_type, overplot_thresh, pdf_fn, alt_model_fn)
    else:
        raw_read_coverage2 = th.parse_fast5s(
            f5_dirs2, corrected_group, basecall_subgroups)
        plot_two_samples(
            plot_intervals, raw_read_coverage, raw_read_coverage2,
            overplot_thresh, overplot_type, pdf_fn, seqs_fn)

    return

def get_unique_intervals(plot_intervals, covered_poss=None, num_regions=None):
    # unique genomic regions filter
    uniq_p_intervals = []
    used_intervals = defaultdict(set)
    for int_i in plot_intervals:
        # could have significant region immediately next to
        # beginning/end of reads
        interval_poss = list(range(int_i.start, int_i.end))
        if int_i.start not in used_intervals[(int_i.chrm, int_i.strand)] and (
                covered_poss is None or all(
                    pos in covered_poss[(int_i.chrm, int_i.strand)]
                    for pos in interval_poss)):
            uniq_p_intervals.append(int_i)
            used_intervals[(int_i.chrm, int_i.strand)].update(interval_poss)
        if num_regions is not None and len(uniq_p_intervals) >= num_regions:
            break

    return uniq_p_intervals

def plot_motif_centered_signif(
        f5_dirs1, corrected_group, basecall_subgroups, pdf_fn,
        f5_dirs2, num_regions, overplot_thresh, motif, stats_fn,
        context_width, num_stats, tb_model_fn, alt_model_fn,
        plot_default_stnd, plot_default_alt, fasta_fn, cov_damp_counts):
    try:
        importr(str('gridExtra'))
    except:
        th._error_message_and_exit(
            'Must have R packge `gridExtra` installed in order to ' +
            'create motif centered plots.')

    motif = th.TomboMotif(motif)
    genome_index = th.Fasta(fasta_fn)

    if VERBOSE: th._status_message('Loading statistics from file.')
    all_stats = ts.TomboStats(stats_fn)
    if VERBOSE: th._status_message('Sorting statistics.')
    all_stats.order_by_frac(cov_damp_counts)

    raw_read_coverage1 = th.parse_fast5s(
        f5_dirs1, corrected_group, basecall_subgroups)
    raw_read_coverage2 = th.parse_fast5s(
        f5_dirs2, corrected_group, basecall_subgroups) \
        if f5_dirs2 is not None else None

    tb_model_fn, alt_model_fn = get_valid_model_fns(
        tb_model_fn, plot_default_stnd, alt_model_fn, plot_default_alt,
        raw_read_coverage1, f5_dirs2)

    if VERBOSE: th._status_message('Finding most signficant regions with motif.')
    motif_regions_data = []
    search_width = ((context_width + motif.motif_len) * 2) - 1
    for reg_seq, chrm, strand, start, end in all_stats.iter_stat_seqs(
            genome_index, motif.motif_len + context_width - 1,
            motif.motif_len + context_width - 1, include_pos=True):
        reg_match = motif.motif_pat.search(reg_seq)
        if reg_match:
            offset = reg_match.start()
            if strand == '-':
                offset = search_width - offset - motif.motif_len
            reg_start = start + offset - context_width
            if reg_start < 0: continue
            # multiple signif stat locations might point to the same motif
            # location, so only take one of these
            if (reg_start, chrm, strand) not in motif_regions_data:
                motif_regions_data.append((reg_start, chrm, strand))
        if len(motif_regions_data) >= num_stats:
            break

    if len(motif_regions_data) == 0:
        th._error_message_and_exit(
            'No covered and tested sites contain motif of interest.')
    if len(motif_regions_data) < num_stats:
        th._warning_message(
            'Fewer covered and tested motif sites found than requested.')

    plot_width = motif.motif_len + (context_width * 2)
    def get_stat_pos(start, chrm, strand):
        # need to handle forward and reverse strand stats separately since
        # reverse strand stats are in reverse order wrt motif
        reg_pos_fracs = []
        for pos in range(start, start + plot_width):
            pos_frac = all_stats.get_pos_frac(chrm, strand, pos)
            if pos_frac is not None:
                reg_pos_fracs.append((pos - start if strand == '+' else
                                      -1 * (pos - start - plot_width + 1),
                                      pos_frac))

        return reg_pos_fracs

    if VERBOSE: th._status_message('Getting all regions statistics.')
    stat_locs = [
        loc_stat for motif_loc in motif_regions_data
        for loc_stat in get_stat_pos(*motif_loc)]

    # TODO: Fix so that negative strand reads are plotted too.
    # requires adding "don't reverse signal" option in getting plot data
    plot_intervals = []
    for i, (reg_start, chrm, strand) in enumerate(motif_regions_data):
        if strand == '-': continue
        plot_intervals.append(th.intervalData(
            '{:03d}'.format(i), chrm, reg_start, reg_start + plot_width,
            strand))
        if len(plot_intervals) >= num_regions:
            break

    plot_motif_centered_with_stats(
        raw_read_coverage1, raw_read_coverage2, plot_intervals,
        stat_locs, overplot_thresh, pdf_fn, tb_model_fn, alt_model_fn)

    return

def cluster_most_signif(
        f5_dirs1, corrected_group, basecall_subgroups, pdf_fn,
        f5_dirs2, num_regions, num_bases,
        r_struct_fn, num_processes, fasta_fn, stats_fn, slide_span):
    if VERBOSE: th._status_message('Loading statistics from file.')
    plot_intervals = ts.TomboStats(stats_fn).get_most_signif_regions(
        num_bases + (slide_span * 2), num_regions)

    raw_read_coverage1 = th.parse_fast5s(
        f5_dirs1, corrected_group, basecall_subgroups)
    raw_read_coverage2 = th.parse_fast5s(
        f5_dirs2, corrected_group, basecall_subgroups)

    # calculate positions covered by at least one read in both sets
    read_coverage1 = th.get_coverage(raw_read_coverage1)
    read_coverage2 = th.get_coverage(raw_read_coverage2)
    covered_poss = dict(
        (chrm_strand, set(
            np.where(read_coverage1[chrm_strand] > 0)[0]).intersection(
                np.where(read_coverage2[chrm_strand] > 0)[0]))
        for chrm_strand in set(read_coverage1).intersection(
                read_coverage2))

    # unique genomic regions filter
    plot_intervals = get_unique_intervals(plot_intervals, covered_poss)

    # get region data if outputting R data structure
    if r_struct_fn is not None:
        if VERBOSE: th._status_message('Getting sequences.')
        # expand regions for getting sequence by N in case motif is
        # the exact range found
        expand_pos = 2
        seq_intervals = [
            int_i._replace(
                start=int_i.start - expand_pos,
                end=int_i.start + expand_pos + num_bases + (slide_span * 2))
            for int_i in plot_intervals]
        if fasta_fn is None:
            # add region sequences to column names for saved dist matrix
            reg_seqs = [reg_data.seq for reg_data in th.get_region_sequences(
                seq_intervals, raw_read_coverage1, raw_read_coverage2)]
        else:
            genome_index = th.Fasta(fasta_fn)
            reg_seqs = [
                genome_index.get_seq(int_i.chrm, int_i.start, int_i.end)
                for int_i in seq_intervals]

    if VERBOSE: th._status_message('Getting region signal differences.')
    signal_diffs = th.get_signal_differences(
        raw_read_coverage1, raw_read_coverage2)
    slide_span_val = slide_span if slide_span else 0
    reg_sig_diffs = [signal_diffs[(int_i.chrm, int_i.strand)][
        int_i.start:int_i.start+num_bases+(slide_span_val*2)]
                     for int_i in plot_intervals]

    if VERBOSE: th._status_message('Getting distance between signals.')
    manager = mp.Manager()
    index_q = manager.Queue()
    dists_q = manager.Queue()

    for i in range(len(reg_sig_diffs)):
        index_q.put(i)

    args = (reg_sig_diffs, index_q, dists_q, slide_span)
    processes = []
    for p_id in range(num_processes):
        p = mp.Process(target=ts.get_pairwise_dists,
                       args=args)
        p.start()
        processes.append(p)

    reg_sig_diff_dists = []
    while any(p.is_alive() for p in processes):
        try:
            row_dists = dists_q.get(block=False)
            reg_sig_diff_dists.append(row_dists)
        except queue.Empty:
            sleep(1)
            continue
    # empty any entries left in queue after processes have finished
    while not dists_q.empty():
        row_dists = dists_q.get(block=False)
        reg_sig_diff_dists.append(row_dists)

    reg_sig_diff_dists = list(map(itemgetter(1), sorted(reg_sig_diff_dists)))

    reg_sig_diff_dists = r.r.matrix(
        r.FloatVector(np.concatenate(reg_sig_diff_dists)),
        ncol=len(reg_sig_diffs), byrow=True)

    if r_struct_fn is not None:
        reg_sig_diff_dists.colnames = r.StrVector(
            ['::'.join((seq, int_i.chrm, int_i.strand, unicode(int_i.start)))
             for seq, int_i in zip(reg_seqs, plot_intervals)])
        r_struct_fn = r.StrVector([r_struct_fn,])
    else:
        r_struct_fn = r.NA_Character

    if VERBOSE: th._status_message('Plotting (and saving data).')
    r.r(resource_string(__name__, 'R_scripts/plotSigMDS.R').decode())
    r.r('pdf("' + pdf_fn + '", height=7, width=7)')
    r.globalenv[str('plotSigMDS')](reg_sig_diff_dists, r_struct_fn)
    r.r('dev.off()')

    return


###########################
#### Test rpy2 install ####
###########################

def test_r_install():
    # first check for simple rpy2 install
    try:
        import rpy2
    except:
        th._error_message_and_exit(
            'Must have rpy2 installed in order to plot. Run ' +
            '`python -c "import rpy2"` to identify installation issues.')

    # try accessing previously loaded rpy2 modules (failed import indicates
    # that rpy2 and R failed to link during install)
    try:
        r
        importr
    except NameError:
        th._error_message_and_exit(
            'R and rpy2 must be linked during installation.\n\t\tRun ' +
            '`python -c "from rpy2 import robjects; from ' +
            'rpy2.robjects.packages import importr` to identify ' +
            'installation linking issues.\n\t\tIf using conda, check that ' +
            'R and rpy2 are installed from the same channel with ' +
            '`conda list | grep "r-base\|rpy2"` (last columns should match).')

    # now test for ggplot2 installation
    try:
        importr(str('ggplot2'))
    except:
        th._error_message_and_exit(
            'Must have R package ggplot2 installed in order to plot. ' +
            'Run `python -c "from rpy2.robjects.packages import importr; ' +
            'importr(str(\'ggplot2\'));"` to identify installation issues.')

    return


################################
#### Main plotting function ####
################################

def _plot_main(args):
    global VERBOSE
    VERBOSE = not args.quiet
    th.VERBOSE = VERBOSE
    ts.VERBOSE = VERBOSE

    test_r_install()

    # roc plotting doesn't use read dirs
    try:
        base_args = [args.fast5_basedirs, args.corrected_group,
                     args.basecall_subgroups, args.pdf_filename]
    except:
        pass
    try:
        genome_opts = [
            ('overplot_thresh', args.overplot_threshold),
            ('overplot_type', args.overplot_type)]
    except:
        pass
    nbase_opt = [('num_bases', args.num_bases
                  if 'num_bases' in args else None),]
    nreg_opt = [('num_regions', args.num_regions
                 if 'num_regions' in args else None),]
    nobs_opt = [('num_obs', args.num_obs
                 if 'num_obs' in args else None),]
    nread_opt = [('num_reads', args.num_reads
                  if 'num_reads' in args else None),]
    glocs_opt = [('genome_locs', args.genome_locations
                  if 'genome_locations' in args else None),]
    f5dirs2_opt = [('f5_dirs2', args.control_fast5_basedirs
                    if 'control_fast5_basedirs' in args else None),]
    rdata_opt = [('r_struct_fn', args.r_data_filename
                  if 'r_data_filename' in args else None),]
    tbmod_opt = [('tb_model_fn', args.tombo_model_filename \
                  if 'tombo_model_filename' in args else None),]
    dtbmod_opt = [('plot_default_stnd', args.plot_standard_model \
                  if 'plot_standard_model' in args else None),]
    atbmod_opt = [('alt_model_fn', args.alternate_model_filename \
                  if 'alternate_model_filename' in args else None),]
    datbmod_opt = [('plot_default_alt', args.plot_alternate_model \
                  if 'plot_alternate_model' in args else None),]
    fasta_opt = [('fasta_fn', args.genome_fasta
                  if 'genome_fasta' in args else None),]
    motif_opt = [('motif', args.motif if 'motif' in args else None),]
    seqfn_opt = [('seqs_fn', args.sequences_filename
                  if 'sequences_filename' in args else None),]
    statfn_opt = [('stats_fn', args.statistics_filename
                   if 'statistics_filename' in args else None),]
    covdamp_opt = [('cov_damp_counts', args.coverage_dampen_counts
                    if 'coverage_dampen_counts' in args else None),]

    if args.action_command == 'max_coverage':
        kwargs = dict(f5dirs2_opt + nreg_opt + nbase_opt + genome_opts +
                      tbmod_opt + atbmod_opt + dtbmod_opt + datbmod_opt)
        plot_max_coverage(*base_args, **kwargs)
    elif args.action_command == 'genome_locations':
        kwargs = dict(f5dirs2_opt + nbase_opt + genome_opts + glocs_opt +
                      tbmod_opt + atbmod_opt + dtbmod_opt + datbmod_opt)
        plot_genome_locations(*base_args, **kwargs)
    elif args.action_command == 'motif_centered':
        kwargs = dict(f5dirs2_opt + nreg_opt + nbase_opt + genome_opts +
                      fasta_opt + motif_opt + tbmod_opt + atbmod_opt +
                      dtbmod_opt + datbmod_opt +
                      [('deepest_coverage', args.deepest_coverage),])
        plot_motif_centered(*base_args, **kwargs)
    elif args.action_command == 'max_difference':
        kwargs = dict(f5dirs2_opt + nreg_opt + nbase_opt + genome_opts +
                      seqfn_opt)
        plot_max_diff(*base_args, **kwargs)
    elif args.action_command == 'most_significant':
        kwargs = dict(f5dirs2_opt + nreg_opt + nbase_opt + genome_opts +
                      seqfn_opt + statfn_opt + tbmod_opt + atbmod_opt +
                      dtbmod_opt + datbmod_opt + covdamp_opt)
        plot_most_signif(*base_args, **kwargs)
    elif args.action_command == 'motif_with_stats':
        kwargs = dict(f5dirs2_opt + nreg_opt + motif_opt + statfn_opt +
                      tbmod_opt + atbmod_opt + dtbmod_opt + datbmod_opt +
                      fasta_opt + covdamp_opt +
                      [('overplot_thresh', args.overplot_threshold),
                       ('context_width', args.num_context),
                       ('num_stats', args.num_statistics)])
        plot_motif_centered_signif(*base_args, **kwargs)
    elif args.action_command == 'cluster_most_significant':
        kwargs = dict(f5dirs2_opt + nreg_opt + nbase_opt +
                      fasta_opt + statfn_opt + rdata_opt +
                      [('num_processes', args.processes),
                       ('slide_span', args.slide_span)])
        cluster_most_signif(*base_args, **kwargs)
    elif args.action_command == 'per_read':
        kwargs = dict(glocs_opt + fasta_opt + nbase_opt +
                      [('per_read_stats_fn', args.per_read_statistics_filename),
                       ('num_reads', args.num_reads),
                       ('box_center', args.box_center)])
        plot_per_read_mods_genome_location(*base_args, **kwargs)
    elif args.action_command == 'kmer':
        kwargs = dict(nread_opt + rdata_opt +
                      [('read_mean', args.read_mean),
                       ('upstrm_bases', args.upstream_bases),
                       ('dnstrm_bases', args.downstream_bases),
                       ('kmer_thresh', args.num_kmer_threshold),
                       ('dont_plot', args.dont_plot)])
        plot_kmer_dist(*base_args, **kwargs)
    elif args.action_command == 'roc':
        kwargs = dict(fasta_opt + covdamp_opt +
                      [('pdf_fn', args.pdf_filename),
                       ('motif_descs', args.motif_descriptions),
                       ('stats_fns', args.statistics_filenames),
                       ('min_reads', args.minimum_test_reads)])
        plot_roc(**kwargs)
    elif args.action_command == 'per_read_roc':
        kwargs = dict(fasta_opt +
                      [('pdf_fn', args.pdf_filename),
                       ('motif_descs', args.motif_descriptions),
                       ('pr_stats_fns', args.per_read_statistics_filenames),
                       ('stats_per_block', args.statistics_per_block),
                       ('total_stats_limit', args.total_statistics_limit)])
        plot_per_read_roc(**kwargs)
    else:
        th._error_message_and_exit('Invalid tombo sub-command entered. ' +
                                   'Should have been caught by argparse.')

    return


if __name__ == '__main__':
    raise NotImplementedError(
        'This is a module. See commands with `tombo -h`')
