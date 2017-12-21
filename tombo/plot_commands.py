import os, sys

import re
import Queue

import numpy as np
import multiprocessing as mp

from time import sleep
from collections import defaultdict
from itertools import repeat, groupby
from pkg_resources import resource_string

# import tombo functions
import tombo_stats as ts
import tombo_helper as th

VERBOSE = False

# quantiles and especially density plots at leat 3 values
# to be meaningful
QUANT_MIN = 3

# plotting names for strands
FWD_STRAND = 'Forward Strand'
REV_STRAND = 'Reverse Strand'


###################################
#### ggplot via rpy2 functions ####
###################################

GG_LOAD_ERROR=(
    '*' * 60 + '\nERROR: Must have rpy2, R and ' +
    'R package ggplot2 installed in order to plot. If these ' +
    'packages are installed, run:\n\t\t`python -c "import rpy2.robjects; ' +
    'from rpy2.robjects.packages import importr; ' +
    'importr(\'ggplot2\');"`\n\t to see installation issues.\n' + \
    '*' * 60 + '\n\n')
try:
    import rpy2.robjects as r
    from rpy2.robjects.packages import importr
except:
    # pass here and raise error when main functions are actually called
    pass

############################################
#### Kmer signal distribution functions ####
############################################

def plot_kmer_dist(
        f5_dirs1, corrected_group, basecall_subgroups, pdf_fn,
        read_mean, upstrm_bases, dnstrm_bases, kmer_thresh, num_reads,
        r_struct_fn, dont_plot):
    if VERBOSE: sys.stderr.write('Parsing files and tabulating k-mers.\n')
    kmer_width = upstrm_bases + dnstrm_bases + 1
    reads_added = 0
    all_kmers = defaultdict(list)

    raw_read_coverage = th.parse_fast5s(
        f5_dirs1, corrected_group, basecall_subgroups)
    if VERBOSE: sys.stderr.write('Extracting read levels.\n')
    files = [r_data for cs_r_data in raw_read_coverage.values()
             for r_data in cs_r_data]
    np.random.shuffle(files)
    for r_data in files:
        means, seq = th.get_multiple_slots_read_centric(
            r_data, ['norm_mean', 'base'])
        if means is None: continue
        seq = ''.join(seq)

        read_kmers = defaultdict(list)
        for kmer, event_mean in zip(
                [seq[i:i+kmer_width]
                 for i in range(len(seq)-kmer_width+1)],
                means[upstrm_bases:]):
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
        sys.stderr.write(
            '****** ERROR ******\n\tOnly zero or one valid reads present. ' +
            'Check corrected group used in resquiggle as well as ' +
            '[--num-kmer-threshold] parameter especially if requested ' +
            'k-mer length is greater than 3 or 4. Consider setting ' +
            'to 0 for k-mer lengths > 4.\n')
        sys.exit()
    if reads_added < num_reads:
        sys.stderr.write(
            '****** WARNING ******\tFewer valid reads present than ' +
            'requested. Check corrected group used in ' +
            'resquiggle as well as [--num-kmer-threshold] ' +
            'parameter especially if requested k-mer length is ' +
            'greater than 3 or 4. Consider setting to 0 for k-mer ' +
            'legnths > 4.\n')

    if VERBOSE: sys.stderr.write('Preparing plot data.\n')
    kmer_levels = [kmer for means, kmer in sorted([
        (np.mean(zip(*means)[0]), kmer)
        for kmer, means in all_kmers.items()])]

    plot_data = [
        (kmer, kmer[upstrm_bases], sig_mean, read_i)
        for kmer in kmer_levels
        for sig_mean, read_i in all_kmers[kmer]]

    kmerDat = r.DataFrame({
        'Kmer':r.FactorVector(
            r.StrVector(zip(*plot_data)[0]),
            ordered=True, levels=r.StrVector(kmer_levels)),
        'Base':r.StrVector(zip(*plot_data)[1]),
        'Signal':r.FloatVector(zip(*plot_data)[2]),
        'Read':r.StrVector(zip(*plot_data)[3])})
    # df to plot kmers as tile of colors requires cowplot R package
    try:
        cowplot = importr("cowplot")
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
        sys.stderr.write(
            '********* WARNING: Install R package `cowplot` for ' +
            'visual kmer display. Using text kmer display. ********\n')
        baseDat = r.NA_Character

    if r_struct_fn is None:
        r_struct_fn = r.NA_Character
    else:
        r_struct_fn = r.StrVector([r_struct_fn,])
    dont_plot_r = r.BoolVector([dont_plot,])

    if VERBOSE: sys.stderr.write('Plotting.\n')
    r.r(resource_string(__name__, 'R_scripts/plotKmerDist.R'))
    if not dont_plot: r.r('pdf("' + pdf_fn + '", height=7, width=10)')
    if read_mean:
        r.globalenv['plotKmerDistWReadPath'](
            kmerDat, baseDat, r_struct_fn, dont_plot_r)
    else:
        r.globalenv['plotKmerDist'](
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
    (read_id, signal_data, raw_offset, shift, scale, lower_lim, upper_lim,
     old_segs, old_align_vals, new_align_vals, events_end,
     new_segs) = read_corr_data

    if reg_type == 'start':
        reg_start = 0
    elif reg_type == 'end':
        reg_start = events_end - num_obs
    elif reg_type == 'random':
        reg_start = np.random.randint(0, events_end - num_obs)
    else:
        # reg_type should be an integer which is the raw start position
        assert isinstance(reg_type, int)
        if r_strand == '+':
            reg_start = int(new_segs[reg_type - r_start] - (num_obs / 2) - 1)
        else:
            reg_start = int((new_segs[len(new_segs) - (reg_type - r_start) - 1]
                             - num_obs) + (num_obs / 2))

    norm_reg_signal, scale_values = th.normalize_raw_signal(
        signal_data, raw_offset + reg_start, num_obs,
        shift=shift, scale=scale, lower_lim=lower_lim,
        upper_lim=upper_lim)

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
                  for old_b, new_b in zip(
                          old_align_vals, new_align_vals)]
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

    old_bases, old_reg_segs = zip(*[
        (b, pos) for b, pos, in_reg in zip(*reg_align_vals)[0]
        if in_reg]) if (
                len(reg_align_vals) > 0 and
                sum(zip(*zip(*reg_align_vals)[0])[2]) > 0) else ([], [])
    new_bases, new_reg_segs = zip(*[
        (b, pos) for b, pos, in_reg in zip(*reg_align_vals)[1]
        if in_reg]) if (
                len(reg_align_vals) > 0 and
                sum(zip(*zip(*reg_align_vals)[1])[2]) > 0) else ([], [])

    # bring positions to zero start if aligning multiple sequences
    sig_range = range(reg_start, reg_start + num_obs)
    if start_at_zero:
        old_reg_segs = [
            old_seg_pos - reg_start for old_seg_pos in old_reg_segs]
        new_reg_segs = [
            new_seg_pos - reg_start for new_seg_pos in new_reg_segs]
        sig_range = range(0, num_obs)

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

def get_event_data(
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

def get_boxplot_data(
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

def get_quant_data(
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

def get_raw_signal_data(
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
                r_sig, overlap_seg_data, start_offset = th.get_raw_signal(
                    r_data, reg_data.start, reg_data.end)
            except:
                if not_warned:
                    not_warned = False
                    sys.stderr.write(
                        '********** WARNING *********\n\tGenome resolved ' +
                        'raw signal could not be retrieved for some reads. ' +
                        'Ensure that reads have been re-squiggled and that ' +
                        'all data slot corresponding accordingly.\n')
                continue
            for base_i, (b_start, b_end) in enumerate(zip(
                    overlap_seg_data[:-1], overlap_seg_data[1:])):
                Position.extend(
                    reg_data.start + base_i + start_offset +
                    np.linspace(0, 1, b_end - b_start, endpoint=False))
                Signal.extend(r_sig[b_start-overlap_seg_data[0]:
                                    b_end-overlap_seg_data[0]])
                Read.extend(list(repeat(
                    str(r_num) + '_' + group_num, b_end - b_start)))
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
    SignalData = get_raw_signal_data(*plot_args)
    QuantData = get_quant_data(*plot_args, pos_offest=quant_offset)
    BoxData = get_boxplot_data(*plot_args)
    EventData = get_event_data(*plot_args)

    return SignalData, QuantData, BoxData, EventData

def get_base_r_data(all_reg_data, zero_start=False, is_rna=False):
    BaseStart, Bases, BaseRegion, BaseStrand = [], [], [], []
    for reg_data in all_reg_data:
        if reg_data.strand == '+' or reg_data.strand is None:
            for i, base in enumerate(reg_data.seq):
                if is_rna and base == 'T':
                    base = 'U'
                if zero_start:
                    BaseStart.append(str(i))
                else:
                    BaseStart.append(str(i + reg_data.start))
                Bases.append(base)
                BaseRegion.append(reg_data.reg_id)
                BaseStrand.append(FWD_STRAND)
        if reg_data.strand == '-' or reg_data.strand is None:
            for i, base in enumerate(reg_data.seq):
                base = base.translate(th.COMP_BASES)
                if is_rna and base == 'T':
                    base = 'U'
                if zero_start:
                    BaseStart.append(str(i))
                else:
                    BaseStart.append(str(i + reg_data.start))
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
            for pos, (base_model_mean, base_model_sd) in fwd_model_data:
                Position.append(pos)
                Strand.append(FWD_STRAND)
                Mean.append(base_model_mean)
                SD.append(base_model_sd)
                Region.append(reg_id)
        if strand == '-' or strand is None:
            for pos, (base_model_mean, base_model_sd) in rev_model_data:
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



########################################
#### Base plotting linker functions ####
########################################

def plot_corrections(
        f5_dirs1, corrected_group, basecall_subgroups, pdf_fn,
        reg_type, num_obs, num_reads):
    if VERBOSE: sys.stderr.write('Preparing plot data.\n')
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
        sys.stderr.write(
            'ERROR: No reads were able to be processed. This command is ' +
            'only applicable to reads processed with event_resquiggle. ' +
            'Also check that --corrected-group and --basecall-subgroup ' +
            'match the event_resquiggle command.\n')
        sys.exit()
    if VERBOSE and len(OldSegDat) < num_reads:
        sys.stderr.write(
            'WARNING: Fewer reads than requested were able to ' +
            'be processed. Likely too few reads provided or ' +
            'those provided were not corrected.\n')
    OldSegDat = r.DataFrame.rbind(*OldSegDat)
    NewSegDat = r.DataFrame.rbind(*NewSegDat)
    SigDat = r.DataFrame.rbind(*SigDat)
    DiffDat = r.DataFrame.rbind(*DiffDat)

    if VERBOSE: sys.stderr.write('Plotting.\n')
    r.r(resource_string(__name__, 'R_scripts/plotReadCorr.R'))
    r.r('pdf("' + pdf_fn + '", height=7, width=11)')
    r.globalenv['plotReadCorr'](OldSegDat, NewSegDat, SigDat, DiffDat)
    r.r('dev.off()')

    return

def plot_multi_corrections(
        f5_dirs1, corrected_group, basecall_subgroups, pdf_fn,
        num_reads_per_plot, num_regions, num_obs, include_orig_bcs,
        genome_locations):
    num_regions = num_regions if num_regions % 2 == 0 else \
                  num_regions + 1
    raw_read_coverage = th.parse_fast5s(
        f5_dirs1, corrected_group, basecall_subgroups)
    read_coverage = th.get_coverage(raw_read_coverage)
    coverage_regions = []
    for chrom_strand, chrom_coverage in read_coverage.items():
        chrm_coverage_regions = [
            (x, len(list(y))) for x, y in groupby(chrom_coverage)]
        chrm_reg_starts = np.cumsum(np.insert(
            zip(*chrm_coverage_regions)[1], 0, 0))
        coverage_regions.extend(zip(
            zip(*chrm_coverage_regions)[0],
            [start + (reg_len / 2) for start, reg_len in
             zip(chrm_reg_starts, zip(*chrm_coverage_regions)[1])],
            repeat(chrom_strand[0]), repeat(chrom_strand[1])))

    if genome_locations is None:
        # randomly select regions with at least num_reads_to_plot regions
        coverage_regions = [
            (chrm, reg_center, strand)
            for stat, reg_center, chrm, strand in
            coverage_regions if stat >= num_reads_per_plot]
        np.random.shuffle(coverage_regions)
        plot_locs = zip(
            ['{:03d}'.format(rn) for rn in range(num_regions)],
            coverage_regions[:num_regions])
        if len(plot_locs) < num_regions:
            sys.stderr.write(
                '*' * 60 + '\nWarning: Fewer regions contain minimum ' +
                'number of reads than requested.\n' + '*' * 60 + '\n')
    else:
        if VERBOSE: sys.stderr.write('Parsing genome locations.\n')
        parsed_locations = []
        for chrm_pos_strand in genome_locations:
            split_vals = chrm_pos_strand.replace('"', '').replace(
                "'", "").split(':')[:3]
            # default to plus strand if not specified
            if len(split_vals) == 2:
                parsed_locations.append((
                    split_vals[0], split_vals[1], '+'))
            else:
                parsed_locations.append(split_vals)
        plot_locs = [
            ('{:03d}'.format(i), (chrm, int(pos) - 1, strand))
            for i, (chrm, pos, strand) in enumerate(parsed_locations)]
        # filter regions with no coverage
        plot_locs = [
            (reg_i, (chrm, start, strand))
            for (reg_i, (chrm, start, strand)) in plot_locs
            if (chrm, strand) in read_coverage and
            read_coverage[(chrm, strand)][start] > 0]
        if len(plot_locs) < len(parsed_locations):
            sys.stderr.write(
                '*' * 60 + '\nWarning: Some regions did not contain ' +
                'read coverage.\n' + '*' * 60 + '\n')

    if len(plot_locs) == 0:
        sys.stderr.write(
            '*' * 60 + '\nERROR: No regions contain minimum ' +
            'number of reads.\n' + '*' * 60 + '\n')
        sys.exit()

    if VERBOSE: sys.stderr.write('Preparing plot data.\n')
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

    if VERBOSE: sys.stderr.write('Plotting.\n')
    r.r(resource_string(__name__, 'R_scripts/plotMultiReadCorr.R'))
    r.r('pdf("' + pdf_fn + '", height=5, width=11)')
    if include_orig_bcs and OldSegDat is not None:
        r.globalenv['plotMultiReadCorr'](OldSegDat, NewSegDat, SigDat)
    else:
        r.globalenv['plotMultiReadCorrNoOrig'](NewSegDat, SigDat)
    r.r('dev.off()')

    return

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
                        " ::: Coverage: " + str(r_cov[0]) + r_ovp[0] + \
                        " + " + str(r_cov[1]) + r_ovp[1] + " -"
            else:
                cov_str = str(r_cov[0]) + r_ovp[0] if int_i.strand == '+' \
                          else str(r_cov[1]) + r_ovp[1]
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
                    str(r_cov[0]) + r_ovp[0] + " + " +
                    str(r_cov[1]) + r_ovp[1] + " -; Control (Black): " +
                    str(r_cov[2]) + r_ovp[2] + " + " +
                    str(r_cov[3]) + r_ovp[3] + " -")
            else:
                cov_str = (
                    'Sample (Red): ' + str(r_cov[0]) + r_ovp[0] +
                    '; Control (Black): ' + str(r_cov[2]) + r_ovp[2]
                ) if int_i.strand == '+' else (
                    'Sample (Red): ' + str(r_cov[1]) + r_ovp[1] +
                    '; Control (Black): ' + str(r_cov[3]) + r_ovp[3])
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
    if VERBOSE: sys.stderr.write('Preparing plot data.\n')
    all_reg_data = th.get_region_reads(plot_intervals, raw_read_coverage)
    if len(all_reg_data) == 0:
        sys.stderr.write(
            '*' * 60 + '\nERROR: No reads in any selected regions.\n'
            + '*' * 60 + '\n')
        sys.exit()
    rna = th.is_rna(raw_read_coverage)

    Titles, plot_types = get_plots_titles(
        all_reg_data, None, overplot_type, overplot_thresh)

    BasesData = get_base_r_data(all_reg_data, is_rna=rna)
    SignalData, QuantData, BoxData, EventData = get_plot_types_data(
        (all_reg_data, plot_types, overplot_thresh, 'Group1'))

    if VERBOSE: sys.stderr.write('Plotting.\n')
    r.r(resource_string(__name__, 'R_scripts/plotSingleRun.R'))
    r.r('pdf("' + pdf_fn + '", height=5, width=11)')
    r.globalenv['plotSingleRun'](SignalData, QuantData, BoxData,
                                 EventData, BasesData, Titles)
    r.r('dev.off()')

    return

def filter_and_merge_group_regs(g1_data, g2_data):
    filt_g1, filt_g2, merged_reg_data, both_no_cov = [], [], [], []
    for r1, r2 in zip(g1_data, g2_data):
        both_reads = r1.reads + r2.reads
        if len(both_reads) > 0:
            merged_reg_data.append(th.intervalData(
                r1.reg_id, r1.chrm, r1.start, r1.end,
                r1.strand, r1.reg_text, both_reads))
            filt_g1.append(r1)
            filt_g2.append(r2)
        else:
            both_no_cov.append(':'.join(map(str, (
                r1.chrm, str(r1.start) + '-' + str(r1.end), r1.strand))))

    if len(both_no_cov) > 0 and VERBOSE:
        sys.stderr.write(
            '*' * 60 + '\nWarning: Some regions include no reads: ' +
            '\t'.join(both_no_cov) + '\n' + '*' * 60 + '\n')

    if len(merged_reg_data) == 0:
        sys.stderr.write(
            '*' * 60 + '\nERROR: No reads in any selected regions.\n'
            + '*' * 60 + '\n')
        sys.exit()

    return merged_reg_data, filt_g1, filt_g2

def plot_two_samples(
        plot_intervals, raw_read_coverage1, raw_read_coverage2,
        overplot_thresh, overplot_type, pdf_fn, seqs_fn=None):
    if VERBOSE: sys.stderr.write('Preparing plot data.\n')
    # get reads overlapping each region
    all_reg_data1 = th.get_region_reads(
        plot_intervals, raw_read_coverage1, filter_no_cov=False, add_seq=False)
    all_reg_data2 = th.get_region_reads(
        plot_intervals, raw_read_coverage2, filter_no_cov=False, add_seq=False)

    # filter regions with no coverage in either read group
    merged_reg_data, all_reg_data1, all_reg_data2 = filter_and_merge_group_regs(
        all_reg_data1, all_reg_data2)

    Titles, plot_types = get_plots_titles(
        all_reg_data1, all_reg_data2, overplot_type, overplot_thresh)

    merged_reg_data = th.add_reg_seq(merged_reg_data)
    BasesData = get_base_r_data(merged_reg_data)

    # get plotting data for either quantiles of raw signal
    SignalData1, QuantData1, BoxData1, EventData1 = get_plot_types_data(
        (all_reg_data1, plot_types, overplot_thresh, 'Group1'), 0.1)
    SignalData2, QuantData2, BoxData2, EventData2 = get_plot_types_data(
        (all_reg_data2, plot_types, overplot_thresh, 'Group2'), 0.5)

    if VERBOSE: sys.stderr.write('Plotting.\n')
    r.r(resource_string(__name__, 'R_scripts/plotGroupComp.R'))
    r.r('pdf("' + pdf_fn + '", height=5, width=11)')
    r.globalenv['plotGroupComp'](
        r.DataFrame.rbind(SignalData1, SignalData2),
        r.DataFrame.rbind(QuantData1, QuantData2),
        r.DataFrame.rbind(BoxData1, BoxData2),
        r.DataFrame.rbind(EventData1, EventData2),
        BasesData, Titles, 0.4)
    r.r('dev.off()')

    if seqs_fn is not None:
        if VERBOSE: sys.stderr.write('Outputting region seqeuences.\n')
        with open(seqs_fn, 'w') as seqs_fp:
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
    def filter_reads(reads, int_start, int_end):
        """ Filter reads obtained from expanded interval
        """
        return [r_data for r_data in reads
                if not (r_data.start >= int_end or r_data.end <= int_start)]
    kmer_ref, upstrm_bases, _, _ = ts.parse_tombo_model(tb_model_fn)
    # compute kmer values to make strand specific calculations easier
    kmer_width = len(next(kmer_ref.iterkeys()))
    dnstrm_bases = kmer_width - upstrm_bases - 1
    expand_width = max(upstrm_bases, dnstrm_bases)
    filt_width = expand_width if min_reg_overlap is None else \
                 expand_width + min_reg_overlap
    if alt_model_fn is not None:
        alt_ref, alt_upstrm_bases, _, _ = ts.parse_tombo_model(
            alt_model_fn)
        assert (alt_upstrm_bases == upstrm_bases and
                kmer_width == len(next(alt_ref.iterkeys()))), (
                    '********* ERROR *********\n\tStandard model not based on ' +
                    'the same kmer position as alternative model.')

    # expand regions to get kmers at first and last positions
    expanded_intervals = [
        th.intervalData(p_int.reg_id, p_int.chrm, p_int.start - expand_width,
                        p_int.end + expand_width, p_int.strand, p_int.reg_text)
        for p_int in plot_intervals]
    # get reads and region sequence
    expanded_intervals = th.get_region_reads(
        expanded_intervals, raw_read_coverage)
    expand_seqs = [int_i.seq for int_i in expanded_intervals]
    rev_expand_seqs = [th.rev_comp(int_i.seq) for int_i in expanded_intervals]
    # convert back to original plot_intervals with seq from exanded intervals
    all_reg_data = [
        th.intervalData(
            int_i.reg_id, int_i.chrm, int_i.start + expand_width,
            int_i.end - expand_width, int_i.strand, int_i.reg_text,
            filter_reads(int_i.reads, int_i.start + filt_width,
                         int_i.end - filt_width),
            int_i.seq[expand_width:-expand_width])
        for int_i in expanded_intervals]

    all_reg_model_data, all_reg_alt_model_data = [], []
    for reg_data, reg_seq, rev_seq in zip(
            all_reg_data, expand_seqs, rev_expand_seqs):
        clipped_reg_seq = reg_seq
        clipped_rev_seq = rev_seq
        if upstrm_bases > dnstrm_bases:
            clipped_reg_seq = reg_seq[:dnstrm_bases-upstrm_bases]
            clipped_rev_seq = rev_seq[:dnstrm_bases-upstrm_bases]
        elif dnstrm_bases > upstrm_bases:
            clipped_reg_seq = reg_seq[dnstrm_bases-upstrm_bases:]
            clipped_rev_seq = rev_seq[dnstrm_bases-upstrm_bases:]
        all_reg_model_data.append((
            reg_data.reg_id, reg_data.strand,
            [(reg_data.start + pos, kmer_ref[''.join(bs)])
             for pos, bs in enumerate(zip(*[
                     clipped_reg_seq[i:] for i in range(kmer_width)]))
             if not any(b in ('-', 'N') for b in bs)],
            [(reg_data.end - pos - 1, kmer_ref[''.join(bs)])
             for pos, bs in enumerate(zip(*[
                     clipped_rev_seq[i:] for i in range(kmer_width)]))
             if not any(b in ('-', 'N') for b in bs)]))
        # if alternative model is supplied add info
        if alt_model_fn is not None:
            all_reg_alt_model_data.append((
                reg_data.reg_id, reg_data.strand,
                [(reg_data.start + pos, alt_ref[''.join(bs)])
                 for pos, bs in enumerate(zip(*[
                         clipped_reg_seq[i:] for i in range(kmer_width)]))
                 if not any(b in ('-', 'N') for b in bs)],
                [(reg_data.end - pos - 1, alt_ref[''.join(bs)])
                 for pos, bs in enumerate(zip(*[
                         clipped_rev_seq[i:] for i in range(kmer_width)]))
                 if not any(b in ('-', 'N') for b in bs)]))

    return all_reg_data, all_reg_model_data, all_reg_alt_model_data

def plot_motif_centered_with_stats(
        raw_read_coverage1, raw_read_coverage2, plot_intervals,
        stat_locs, overplot_thresh, pdf_fn, frac_order,
        tb_model_fn, alt_model_fn=None):
    if VERBOSE: sys.stderr.write('Preparing plot data.\n')

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

    # stat lists
    StatsData = r.DataFrame({
        'Position':r.FloatVector(zip(*stat_locs)[0]),
        'Stat':r.FloatVector(zip(*stat_locs)[1])})

    if VERBOSE: sys.stderr.write('Plotting.\n')
    r.r(resource_string(__name__, 'R_scripts/plotMotifStats.R'))
    r.r('pdf("' + pdf_fn + '", height=5, width=8)')
    r.globalenv['plotMotifStats'](
        SignalData, BasesData, StatsData, frac_order, ModelData)
    r.r('dev.off()')

    return

def plot_model_single_sample(
        plot_intervals, raw_read_coverage, tb_model_fn,
        overplot_type, overplot_thresh, pdf_fn, alt_model_fn=None, seqs_fn=None):
    if VERBOSE: sys.stderr.write('Preparing plot data.\n')
    # get reads overlapping each region along with all kmers
    all_reg_data, all_reg_model_data, all_reg_alt_model_data = get_reg_kmers(
        tb_model_fn, plot_intervals, raw_read_coverage,
        alt_model_fn=alt_model_fn)
    rna = th.is_rna(raw_read_coverage)

    Titles, plot_types = get_plots_titles(
        all_reg_data, None, overplot_type, overplot_thresh, True)
    ModelData = get_model_r_data(all_reg_model_data)
    if alt_model_fn is not None:
        AltModelData = get_model_r_data(all_reg_alt_model_data)
    BasesData = get_base_r_data(all_reg_data, is_rna=rna)
    SignalData, QuantData, BoxData, EventData = get_plot_types_data(
        (all_reg_data, plot_types, overplot_thresh, 'Group1'))

    if VERBOSE: sys.stderr.write('Plotting.\n')
    r.r(resource_string(__name__, 'R_scripts/plotModelComp.R'))
    r.r('pdf("' + pdf_fn + '", height=5, width=11)')
    if alt_model_fn is None:
        r.globalenv['plotModelComp'](SignalData, QuantData, BoxData,
                                     EventData, BasesData, Titles, ModelData)
    else:
        r.globalenv['plotModelComp'](SignalData, QuantData, BoxData, EventData,
                                     BasesData, Titles, ModelData, AltModelData)
    r.r('dev.off()')

    if seqs_fn is not None:
        if VERBOSE: sys.stderr.write('Outputting region seqeuences.\n')
        with open(seqs_fn, 'w') as seqs_fp:
            for int_i in all_reg_data:
                reg_seq = int_i.seq if int_i.strand == '+' else th.rev_comp(
                    int_i.seq)
                seqs_fp.write('>{0}::{1:d}::{2} {3}\n{4}\n'.format(
                    int_i.chrm, int_i.start, int_i.strand, int_i.reg_text,
                    ''.join(reg_seq)))

    return

def plot_per_read_modification(
        plot_intervals, raw_read_coverage, tb_model_fn, alt_model_fn, pdf_fn,
        fm_lag, num_reads, box_center):
    if alt_model_fn is not None:
        alt_ref, upstrm_bases, alt_base, _ = ts.parse_tombo_model(alt_model_fn)
        dnstrm_bases = len(next(alt_ref.iterkeys())) - upstrm_bases - 1

    def calc_read_lh(r_means, alt_base_locs, begin_lag, end_lag, reg_ref_means,
                     reg_ref_vars, reg_alt_means, reg_alt_vars):
        r_stats = np.empty(r_means.shape)
        r_stats[:] = np.NAN
        for alt_base_pos in alt_base_locs:
            alt_pos = alt_base_pos + end_lag
            r_stats[alt_pos] = ts.calc_llh_ratio(
                r_means[alt_pos-end_lag:alt_pos+begin_lag],
                reg_ref_means[alt_pos-end_lag:alt_pos+begin_lag],
                reg_ref_vars[alt_pos-end_lag:alt_pos+begin_lag],
                reg_alt_means[alt_pos-end_lag:alt_pos+begin_lag],
                reg_alt_vars[alt_pos-end_lag:alt_pos+begin_lag])

        return r_stats

    if VERBOSE: sys.stderr.write('Preparing plot data.\n')
    # get reads overlapping each region along with all kmers
    min_reg_overlap = (fm_lag * 2) + 1 if fm_lag > 0 else None
    all_reg_data, all_reg_model_data, all_reg_alt_model_data = get_reg_kmers(
        tb_model_fn, plot_intervals, raw_read_coverage, min_reg_overlap,
        alt_model_fn)
    all_reg_stats = []
    for reg_i, (reg_data, (reg_id, strand, fwd_model_data,
                           rev_model_data)) in enumerate(zip(
                               all_reg_data, all_reg_model_data)):
        model_data = (fwd_model_data if reg_data.strand == "+" else
                      rev_model_data[::-1])
        ref_means, ref_sds = zip(*zip(*model_data)[1])
        reg_events = get_reg_events(
            reg_data.reads, reg_data.start, reg_data.end, reg_data.strand,
            read_rows=True, num_reads=num_reads)
        if alt_model_fn is None:
            with np.errstate(invalid='ignore'):
                pvals = ts.z_score_to_p_value(
                    -np.abs(reg_events - ref_means) / ref_sds) * 2.0
                if fm_lag > 0:
                    all_reg_stats.append((
                        reg_data.reg_id,
                        ts.calc_window_fishers_method(pvals, fm_lag)))
                else:
                    all_reg_stats.append((reg_data.reg_id, pvals))
        else:
            fwd_alt_model, rev_alt_model = all_reg_alt_model_data[reg_i][2:]
            alt_model = (fwd_alt_model if reg_data.strand == "+" else
                         rev_alt_model[::-1])
            begin_lag = upstrm_bases if reg_data.strand =='+' else dnstrm_bases
            end_lag = dnstrm_bases if reg_data.strand =='+' else upstrm_bases
            alt_means, alt_sds = zip(*zip(*alt_model)[1])
            ref_vars = np.square(ref_sds)
            alt_vars = np.square(alt_sds)
            alt_base_locs = [alt_pos.start() for alt_pos in re.finditer(
                alt_base, reg_data.seq[end_lag:-begin_lag])]
            all_reg_stats.append((reg_data.reg_id, np.apply_along_axis(
                calc_read_lh, 1, reg_events,
                alt_base_locs, begin_lag, end_lag, ref_means,
                ref_vars, alt_means, alt_vars)))

    are_pvals = alt_model_fn is None
    StatData, OrdData = get_reg_r_stats(all_reg_stats, are_pvals)
    BasesData = get_base_r_data(all_reg_data, zero_start=True)

    if VERBOSE: sys.stderr.write('Plotting.\n')
    r.r(resource_string(__name__, 'R_scripts/plotPerReadStats.R'))
    r.r('pdf("' + pdf_fn + '", height=5, width=11)')
    r.globalenv['plotPerReadStats'](StatData, OrdData, BasesData, box_center,
                                    are_pvals)
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
        sys.stderr.write(
            '********* WARNING ******** Both a second set of FAST5s and a ' +
            'tombo model were provided. Two samples with model plotting is not ' +
            'currently available. Models requested will be ignored.\n')

    return tb_model_fn, alt_model_fn

def plot_max_coverage(
        f5_dirs1, corrected_group, basecall_subgroups, pdf_fn,
        f5_dirs2, num_regions, num_bases, overplot_thresh, overplot_type,
        tb_model_fn, alt_model_fn, plot_default_stnd, plot_default_alt):
    raw_read_coverage = th.parse_fast5s(
        f5_dirs1, corrected_group, basecall_subgroups)
    read_coverage = th.get_coverage(raw_read_coverage)

    tb_model_fn, alt_model_fn = get_valid_model_fns(
        tb_model_fn, plot_default_stnd, alt_model_fn, plot_default_alt,
        raw_read_coverage, f5_dirs2)
    if f5_dirs2 is None:
        coverage_regions = []
        for (chrom, strand), chrom_coverage in read_coverage.items():
            chrm_coverage_regions = [
                (x, len(list(y))) for x, y in groupby(chrom_coverage)]
            coverage_regions.extend(zip(
                zip(*chrm_coverage_regions)[0],
                np.cumsum(np.insert(
                    zip(*chrm_coverage_regions)[1], 0, 0)),
                repeat(chrom), repeat(strand)))

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
        read_coverage2 = th.get_coverage(raw_read_coverage2)
        coverage_regions = []
        # only process chromosomes in both read groups
        for (chrom, strand) in set(read_coverage).intersection(
                read_coverage2):
            chrom_coverage = read_coverage[(chrom, strand)]
            chrom_coverage2 = read_coverage2[(chrom, strand)]
            if chrom_coverage.shape[0] >= chrom_coverage2.shape[0]:
                merged_chrom_cov = np.pad(
                    chrom_coverage2, (0,chrom_coverage.shape[0] -
                                      chrom_coverage2.shape[0]),
                    'constant', constant_values=0) + chrom_coverage
            else:
                merged_chrom_cov = np.pad(
                    chrom_coverage, (0, chrom_coverage2.shape[0] -
                                     chrom_coverage.shape[0]),
                    'constant', constant_values=0) + chrom_coverage2

            chrm_coverage_regions = [
                (x, len(list(y))) for x, y in groupby(merged_chrom_cov)]
            coverage_regions.extend(zip(
                zip(*chrm_coverage_regions)[0],
                np.cumsum(np.insert(
                    zip(*chrm_coverage_regions)[1], 0, 0)),
                repeat(chrom), repeat(strand)))

        # max coverage plots both strands coverage
        plot_intervals = [
            th.intervalData('{:03d}'.format(rn), chrm, start, start + num_bases)
            for rn, (stat, start, chrm, strand) in
            enumerate(sorted(coverage_regions, reverse=True)[:num_regions])]

        plot_two_samples(
            plot_intervals, raw_read_coverage, raw_read_coverage2,
            overplot_thresh, overplot_type, pdf_fn)

    return

def plot_genome_locations(
        f5_dirs1, corrected_group, basecall_subgroups, pdf_fn,
        f5_dirs2, num_bases, overplot_thresh, overplot_type,
        genome_locations, tb_model_fn, alt_model_fn, plot_default_stnd,
        plot_default_alt):
    if VERBOSE: sys.stderr.write('Parsing genome locations.\n')
    # ignore strand for genome location plotting
    genome_locations = [
        chrm_pos.replace('"', '').replace("'", "").split(':')[:3]
        for chrm_pos in genome_locations]
    # minus one here as all python internal coords are 0-based, but
    # genome is generally 1-based
    plot_intervals = []
    for i, chrm_pos_strand in enumerate(genome_locations):
        if len(chrm_pos_strand) == 2:
            chrm, pos = chrm_pos_strand
            strand = None
        else:
            chrm, pos, strand = chrm_pos_strand
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
        f5_dirs1, corrected_group, basecall_subgroups, pdf_fn,
        num_bases, genome_locations, tb_model_fn, alt_model_fn,
        fm_lag, num_reads, box_center, plot_default_stnd,
        plot_default_alt):
    if (tb_model_fn is None and not plot_default_stnd and
        alt_model_fn is None and not plot_default_alt):
        sys.stderr.write('*********** WARNING ********\n\tNo model ' +
                         'indicated, so loading default standard model.\n')
        plot_default_stnd = True

    if VERBOSE: sys.stderr.write(
            'Parsing per read modifications at genome locations.\n')
    genome_locations = [
        chrm_pos.replace('"', '').replace("'", "").split(':')[:3]
        for chrm_pos in genome_locations]
    # minus one here as all python internal coords are 0-based, but
    # genome is generally 1-based
    plot_intervals = []
    for i, chrm_pos_strand in enumerate(genome_locations):
        if len(chrm_pos_strand) == 2:
            chrm, pos = chrm_pos_strand
            strand = '+'
        else:
            chrm, pos, strand = chrm_pos_strand
        int_start = max(
            0, int(int(pos) - np.floor(num_bases / 2.0) - 1) + 1)
        plot_intervals.append(th.intervalData(
            '{:03d}'.format(i), chrm, int_start, int_start + num_bases, strand))

    raw_read_coverage = th.parse_fast5s(
        f5_dirs1, corrected_group, basecall_subgroups)
    tb_model_fn, alt_model_fn = get_valid_model_fns(
        tb_model_fn, plot_default_stnd, alt_model_fn,
        plot_default_alt, raw_read_coverage)

    plot_per_read_modification(
        plot_intervals, raw_read_coverage, tb_model_fn, alt_model_fn, pdf_fn,
        fm_lag, num_reads, box_center)

    return

def plot_motif_centered(
        f5_dirs1, corrected_group, basecall_subgroups, pdf_fn,
        f5_dirs2, num_regions, num_bases, overplot_thresh, overplot_type,
        motif, fasta_fn, deepest_coverage, tb_model_fn, alt_model_fn,
        plot_default_stnd, plot_default_alt):
    if VERBOSE: sys.stderr.write('Identifying genomic k-mer locations.\n')
    fasta_records = th.parse_fasta(fasta_fn)
    motif_pat = th.parse_motif(motif)
    motif_len = len(motif)
    def get_motif_locs(covered_chrms):
        # TODO: search over negative strand as well
        motif_locs = []
        for chrm, seq in fasta_records.iteritems():
            if chrm not in covered_chrms: continue
            for motif_loc in motif_pat.finditer(seq):
                motif_locs.append((chrm, motif_loc.start()))

        if len(motif_locs) == 0:
            sys.stderr.write(
                'Motif (' + motif_pat.pattern +
                ') not found in genome.\n')
            sys.exit()
        elif len(motif_locs) < num_regions:
            sys.stderr.write(
                'WARNING: Motif (' + motif_pat.pattern +
                ') only found ' + str(len(motif_locs)) +
                'times in genome.\n')
            num_region = len(motif_locs)
        np.random.shuffle(motif_locs)

        return motif_locs

    raw_read_coverage = th.parse_fast5s(
        f5_dirs1, corrected_group, basecall_subgroups)
    tb_model_fn, alt_model_fn = get_valid_model_fns(
        tb_model_fn, plot_default_stnd, alt_model_fn, plot_default_alt,
        raw_read_coverage, f5_dirs2)

    if deepest_coverage:
        read_coverage = th.get_coverage(raw_read_coverage)
    if f5_dirs2 is not None:
        raw_read_coverage2 = th.parse_fast5s(
            f5_dirs2, corrected_group, basecall_subgroups)

        covered_chrms = set(zip(*raw_read_coverage)[0]).intersection(
            zip(*raw_read_coverage2)[0])
        # filter out motif_locs to chromosomes not covered
        motif_locs = get_motif_locs(covered_chrms)

        if deepest_coverage:
            read_coverage2 = th.get_coverage(raw_read_coverage2)
            if VERBOSE: sys.stderr.write(
                    'Finding deepest coverage regions.\n')
            def get_cov(chrm, pos):
                try:
                    plus_cov = min(read_coverage[(chrm, '+')][pos],
                                   read_coverage2[(chrm, '+')][pos])
                except (IndexError, KeyError):
                    plus_cov = 0
                try:
                    minus_cov = min(read_coverage[(chrm, '-')][pos],
                                    read_coverage2[(chrm, '-')][pos])
                except (IndexError, KeyError):
                    minus_cov = 0
                    return max(plus_cov, minus_cov)

            motif_locs_cov = sorted([
                (get_cov(chrm, pos), chrm, pos)
                for chrm, pos in motif_locs], reverse=True)
            if motif_locs_cov[0][0] == 0:
                sys.stderr.write(
                    '*' * 60 + '\nERROR: Motif not covered ' +
                    'by both groups at any positions.\n'
                    + '*' * 60 + '\n')
                sys.exit()

            plot_intervals = []
            for i, (cov, chrm, pos) in enumerate(motif_locs_cov):
                int_start = max(
                    0, pos - int((num_bases - motif_len + 1) / 2.0))
                plot_intervals.append(th.intervalData(
                    '{:03d}'.format(i), chrm, int_start,
                    int_start + num_bases, '+'))
                if len(plot_intervals) >= num_regions: break
        else:
            # iterate over regions and check if they have any coverage
            plot_intervals = []
            for i, (chrm, pos) in enumerate(motif_locs):
                int_start = max(
                    0, pos - int((num_bases - motif_len + 1) / 2.0))
                int_end = int_start + num_bases
                for strand in ('+', '-'):
                    if ((chrm, strand) in raw_read_coverage and
                        (chrm, strand) in raw_read_coverage2 and
                        any(r_data.start < pos < r_data.end
                            for r_data in raw_read_coverage[(chrm, strand)]) and
                        any(r_data2.start < pos < r_data2.end
                            for r_data2 in raw_read_coverage2[(chrm, strand)])):
                        plot_intervals.append(th.intervalData(
                            '{:03d}'.format(i), chrm, int_start,
                            int_end, strand))
                if len(plot_intervals) >= num_regions: break

        plot_two_samples(
            plot_intervals, raw_read_coverage, raw_read_coverage2,
            overplot_thresh, overplot_type, pdf_fn)
    else:
        covered_chrms = set(zip(*raw_read_coverage)[0])
        # filter out motif_locs to chromosomes not covered
        motif_locs = get_motif_locs(covered_chrms)

        if deepest_coverage:
            if VERBOSE: sys.stderr.write(
                    'Finding deepest coverage regions.\n')
            def get_cov(chrm, pos):
                try:
                    plus_cov = read_coverage[(chrm, '+')][pos]
                except (IndexError, KeyError):
                    plus_cov = 0
                try:
                    minus_cov = read_coverage[(chrm, '-')][pos]
                except (IndexError, KeyError):
                    minus_cov = 0
                return max(plus_cov, minus_cov)

            motif_locs_cov = sorted([
                (get_cov(chrm, pos), chrm, pos)
                for chrm, pos in motif_locs], reverse=True)

            plot_intervals = []
            for i, (cov, chrm, pos)in enumerate(motif_locs_cov):
                int_start = max(
                    0, pos - int((num_bases - motif_len + 1) / 2.0))
                int_end = int_start + num_bases
                plot_intervals.append(th.intervalData(
                    '{:03d}'.format(i), chrm, int_start, int_end, '+'))
                if len(plot_intervals) >= num_regions: break
        else:
            # iterate over regions and check if they have any coverage
            plot_intervals = []
            for i, (chrm, pos) in enumerate(motif_locs):
                int_start = max(
                    0, pos - int((num_bases - motif_len + 1) / 2.0))
                int_end = int_start + num_bases
                for strand in ('+', '-'):
                    if ((chrm, strand) in raw_read_coverage and
                        any(r_data.start < pos < r_data.end
                            for r_data in raw_read_coverage[(chrm, strand)])):
                        plot_intervals.append(th.intervalData(
                            '{:03d}'.format(i), chrm, int_start,
                            int_end, strand))
                if len(plot_intervals) >= num_regions: break

        if tb_model_fn is None:
            plot_single_sample(
                plot_intervals, raw_read_coverage, overplot_thresh,
                overplot_type, pdf_fn)
        else:
            plot_model_single_sample(
                plot_intervals, raw_read_coverage, tb_model_fn,
                overplot_type, overplot_thresh, pdf_fn, alt_model_fn)

    return

def plot_max_diff(
        f5_dirs1, corrected_group, basecall_subgroups, pdf_fn,
        f5_dirs2, num_regions, num_bases, overplot_thresh, overplot_type,
        seqs_fn):
    raw_read_coverage1 = th.parse_fast5s(
        f5_dirs1, corrected_group, basecall_subgroups)
    raw_read_coverage2 = th.parse_fast5s(
        f5_dirs2, corrected_group, basecall_subgroups)

    chrm_sizes = th.get_chrm_sizes(raw_read_coverage1, raw_read_coverage2)

    if VERBOSE: sys.stderr.write('Getting base signal.\n')
    base_means1 = th.get_all_mean_levels(raw_read_coverage1, chrm_sizes)
    base_means2 = th.get_all_mean_levels(raw_read_coverage2, chrm_sizes)

    if VERBOSE: sys.stderr.write(
            'Get differences between base signal.\n')
    # get num_region max diff regions from each chrm then find
    # global largest after
    largest_diff_indices = []
    for chrm, chrm_size in chrm_sizes.items():
        for strand in ('+', '-'):
            # calculate difference and set no coverage (nan) values
            # to zero
            chrm_diffs = np.nan_to_num(
                np.abs(base_means1[(chrm, strand)] -
                       base_means2[(chrm, strand)]))
            chrm_max_diff_regs = np.argsort(
                chrm_diffs)[::-1][:num_regions]
            largest_diff_indices.extend((
                chrm_diffs[pos], max(pos - int(num_bases / 2.0), 0),
                chrm, strand) for pos in chrm_max_diff_regs)

    plot_intervals = [
        th.intervalData(
            '{:03d}'.format(rn), chrm, start, start + num_bases, strand,
            '(Mean diff: {:.2f})'.format(stat))
        for rn, (stat, start, chrm, strand) in
        enumerate(sorted(largest_diff_indices, reverse=True)[:num_regions])]

    plot_two_samples(
        plot_intervals, raw_read_coverage1, raw_read_coverage2,
        overplot_thresh, overplot_type, pdf_fn, seqs_fn)

    return

def plot_most_signif(
        f5_dirs1, corrected_group, basecall_subgroups, pdf_fn,
        f5_dirs2, num_regions, overplot_thresh, seqs_fn, num_bases,
        overplot_type, qval_thresh, stats_fn, tb_model_fn, alt_model_fn,
        plot_default_stnd, plot_default_alt, stat_order):
    if VERBOSE: sys.stderr.write('Loading statistics from file.\n')
    all_stats, stat_type = ts.parse_stats(stats_fn)

    raw_read_coverage = th.parse_fast5s(
        f5_dirs1, corrected_group, basecall_subgroups)
    plot_intervals = ts.get_most_signif_regions(
        all_stats, num_bases, num_regions, qval_thresh,
        fraction_order=not stat_order)
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

def get_unique_intervals(plot_intervals, covered_poss, num_regions=None):
    # unique genomic regions filter
    uniq_p_intervals = []
    used_intervals = defaultdict(set)
    for int_i in plot_intervals:
        # could have significant region immediately next to
        # beginning/end of reads
        interval_poss = range(int_i.start, int_i.end)
        if int_i.start not in used_intervals[(int_i.chrm, int_i.strand)] and all(
                pos in covered_poss[(int_i.chrm, int_i.strand)]
                for pos in interval_poss):
            uniq_p_intervals.append(int_i)
            used_intervals[(int_i.chrm, int_i.strand)].update(interval_poss)
        if num_regions is not None and len(uniq_p_intervals) >= num_regions:
            break

    return uniq_p_intervals

def plot_motif_centered_signif(
        f5_dirs1, corrected_group, basecall_subgroups, pdf_fn,
        f5_dirs2, num_regions, overplot_thresh, motif, stats_fn,
        context_width, num_stats, stat_order, tb_model_fn):
    try:
        cowplot = importr("cowplot")
    except:
        sys.stderr.write(
            '*' * 60 + '\nERROR: Must have R packge `cowplot` ' +
            'installed in order to create motif centered plots ' +
            '(install via `install.packages(cowplot)` from ' +
            'an R prompt).\n' + '*' * 60 + '\n\n')
        sys.exit()

    motif_pat = th.parse_motif(motif)
    motif_len = len(motif)

    if VERBOSE: sys.stderr.write('Loading statistics from file.\n')
    all_stats, stat_type = ts.parse_stats(stats_fn)
    # check that fraction un-modified is included in stats data
    if not stat_order and all(np.isnan(all_stats['frac'][:50])):
        sys.stderr.write(
            '********** WARNING *********\tFraction requested for plotting, ' +
            'but does not appear to be included in statistics file, so ' +
            'plotting q-values.')
        stat_order = True
    if not stat_order:
        all_stats.sort(order='frac')

    raw_read_coverage1 = th.parse_fast5s(
        f5_dirs1, corrected_group, basecall_subgroups)
    raw_read_coverage2 = th.parse_fast5s(
        f5_dirs2, corrected_group, basecall_subgroups) \
        if f5_dirs2 is not None else None

    def log_max_stat(pval):
        return -np.log10(max(th.SMALLEST_PVAL, pval))
    def get_stats(stat):
        if not stat_order:
            return 1 - stat['frac']
        return log_max_stat(stat['mt_stat'])
    all_stats_dict = dict(
        ((str(stat['chrm']), str(stat['strand']), stat['pos']), get_stats(stat))
        for stat in all_stats)
    covered_poss = defaultdict(set)
    for stat in all_stats:
        covered_poss[(str(stat['chrm']), str(stat['strand']))].add(stat['pos'])

    if VERBOSE: sys.stderr.write(
            'Finding signficant regions with motif.\n')
    motif_regions_data = []
    for stat in all_stats:
        reg_data = th.get_region_sequences(
            [th.intervalData(
                '0', str(stat['chrm']), stat['pos'] - motif_len + 1,
                stat['pos'] + motif_len, str(stat['strand']))],
            raw_read_coverage1, raw_read_coverage2)[0]

        reg_match = motif_pat.search(reg_data.seq)
        if reg_match:
            motif_regions_data.append((
                stat['pos'], str(stat['chrm']), str(stat['strand']),
                reg_match.start()))
        if len(motif_regions_data) >= num_stats:
            break

    # get plot intervals for all stat regions then trim to
    # num_regions after getting all p-values for plotting
    plot_width = motif_len + (context_width * 2)
    plot_intervals = []
    for i, (pos, chrm, strand, offset) in enumerate(motif_regions_data):
        int_start = pos - motif_len + offset - context_width + 1
        plot_intervals.append(th.intervalData(
            '{:03d}'.format(i), chrm, int_start, int_start + plot_width,
            strand))
        if len(plot_intervals) >= num_stats: break
    # need to handle forward and reverse strand stats separately since
    # reverse strand stats are in reverse order wrt motif
    # note check for key in stats dict as significant position
    # may lie next to region with coverage below the threshold
    stat_locs = [
        (pos - int_i.start, all_stats_dict[(int_i.chrm, int_i.strand, pos)]
         if (chrm, strand, pos) in all_stats_dict else 0.0)
        for int_i in plot_intervals for pos in range(int_i.start, int_i.end)
        if int_i.strand == '+'] + [
                (-1 * (pos - int_i.end + 1),
                 all_stats_dict[(int_i.chrm, int_i.strand, pos)]
                 if (int_i.chrm, int_i.strand, pos) in all_stats_dict
                 else 0.0) for int_i in plot_intervals
                for pos in range(int_i.start, int_i.end) if int_i.strand == '-']
    # TODO: Fix so that negative strand reads are plotted too.
    # requires adding "don't reverse signal" option in getting plot
    # data
    plot_intervals = [
        int_i for int_i in plot_intervals if int_i.strand == '+']
    plot_intervals = get_unique_intervals(
        plot_intervals, covered_poss, num_regions)

    plot_motif_centered_with_stats(
        raw_read_coverage1, raw_read_coverage2, plot_intervals,
        stat_locs, overplot_thresh, pdf_fn, not stat_order, tb_model_fn)

    return

def cluster_most_signif(
        f5_dirs1, corrected_group, basecall_subgroups, pdf_fn,
        f5_dirs2, num_regions, qval_thresh, num_bases,
        r_struct_fn, num_processes, fasta_fn, stats_fn, slide_span):
    if VERBOSE: sys.stderr.write('Loading statistics from file.\n')
    all_stats, stat_type = ts.parse_stats(stats_fn)

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

    plot_intervals = ts.get_most_signif_regions(
        all_stats, num_bases + (slide_span * 2), num_regions, qval_thresh)

    # unique genomic regions filter
    plot_intervals = get_unique_intervals(plot_intervals, covered_poss)

    # get region data if outputting R data structure
    if r_struct_fn is not None:
        if VERBOSE: sys.stderr.write('Getting sequences.\n')
        # expand regions for getting sequence by N in case motif is
        # the exact range found
        expand_pos = 2
        seq_intervals = [
            th.intervalData(
                int_i.reg_id, int_i.chrm, int_i.start - expand_pos,
                int_i.start + expand_pos + num_bases + (slide_span * 2),
                int_i.strand, int_i.reg_text)
            for int_i in plot_intervals]
        if fasta_fn is None:
            # add region sequences to column names for saved dist matrix
            reg_seqs = [reg_data.seq for reg_data in th.get_region_sequences(
                seq_intervals, raw_read_coverage1, raw_read_coverage2)]
        else:
            fasta_records = th.parse_fasta(fasta_fn)
            reg_seqs = [
                fasta_records[int_i.chrm][int_i.start:int_i.end]
                for int_i in seq_intervals]

    if VERBOSE: sys.stderr.write('Getting base signal.\n')
    chrm_sizes = th.get_chrm_sizes(raw_read_coverage1, raw_read_coverage2)

    base_means1 = th.get_all_mean_levels(raw_read_coverage1, chrm_sizes)
    base_means2 = th.get_all_mean_levels(raw_read_coverage2, chrm_sizes)

    if VERBOSE: sys.stderr.write('Getting region signal difference.\n')
    slide_span_val = slide_span if slide_span else 0
    reg_sig_diffs = [
        np.nan_to_num(
            base_means1[(int_i.chrm, int_i.strand)][
                int_i.start:int_i.start+num_bases+(slide_span_val*2)] -
            base_means2[(int_i.chrm, int_i.strand)][
                int_i.start:int_i.start+num_bases+(slide_span_val*2)])
        for int_i in plot_intervals]

    # some code to output reg signal for discovery plotting
    """sys.stdout.write('\n'.join(
        '\t'.join(('\t'.join(map(str, reg_diff)), reg_seq,
                   int_i.chrm, str(int_i.start), int_i.strand))
        for reg_seq, int_i, reg_diff in
        zip(reg_seqs, plot_intervals reg_sig_diffs)) + '\n')
    sys.exit()"""

    if VERBOSE: sys.stderr.write('Getting distance between signals.\n')
    manager = mp.Manager()
    index_q = manager.Queue()
    dists_q = manager.Queue()

    for i in range(len(reg_sig_diffs)):
        index_q.put(i)

    args = (reg_sig_diffs, index_q, dists_q, slide_span)
    processes = []
    for p_id in xrange(num_processes):
        p = mp.Process(target=ts.get_pairwise_dists,
                       args=args)
        p.start()
        processes.append(p)

    reg_sig_diff_dists = []
    while any(p.is_alive() for p in processes):
        try:
            row_dists = dists_q.get(block=False)
            reg_sig_diff_dists.append(row_dists)
        except Queue.Empty:
            sleep(1)
            continue
    # empty any entries left in queue after processes have finished
    while not dists_q.empty():
        row_dists = dists_q.get(block=False)
        reg_sig_diff_dists.append(row_dists)

    reg_sig_diff_dists = zip(*sorted(reg_sig_diff_dists))[1]

    reg_sig_diff_dists = r.r.matrix(
        r.FloatVector(np.concatenate(reg_sig_diff_dists)),
        ncol=len(reg_sig_diffs), byrow=True)

    if r_struct_fn is not None:
        reg_sig_diff_dists.colnames = r.StrVector(
            ['::'.join((seq, int_i.chrm, int_i.strand, str(int_i.start)))
             for seq, int_i in zip(reg_seqs, plot_intervals)])
        r_struct_fn = r.StrVector([r_struct_fn,])
    else:
        r_struct_fn = r.NA_Character

    if VERBOSE: sys.stderr.write('Plotting (and saving data).\n')
    r.r(resource_string(__name__, 'R_scripts/plotSigMDS.R'))
    r.r('pdf("' + pdf_fn + '", height=7, width=7)')
    r.globalenv['plotSigMDS'](reg_sig_diff_dists, r_struct_fn)
    r.r('dev.off()')

    return




################################
#### Main plotting function ####
################################

def plot_main(args):
    global VERBOSE
    VERBOSE = not args.quiet
    th.VERBOSE = VERBOSE
    ts.VERBOSE = VERBOSE

    try:
        ggplot = importr("ggplot2")
    except:
        sys.stderr.write(GG_LOAD_ERROR)
        sys.exit()

    base_args = [args.fast5_basedirs, args.corrected_group,
                 args.basecall_subgroups, args.pdf_filename]
    try:
        genome_opts = [
            ('overplot_thresh', args.overplot_threshold),
            ('overplot_type', args.overplot_type)]
    except:
        pass
    nbase_opt = [('num_bases', args.num_bases if 'num_bases' in args else None),]
    nreg_opt = [('num_regions', args.num_regions
                 if 'num_regions' in args else None),]
    nobs_opt = [('num_obs', args.num_obs
                 if 'num_obs' in args else None),]
    nread_opt = [('num_reads', args.num_reads
                  if 'num_reads' in args else None),]
    glocs_opt = [('genome_locations', args.genome_locations
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
    qval_opt = [('qval_thresh', args.q_value_threshold
                 if 'q_value_threshold' in args else None),]
    statfn_opt = [('stats_fn', args.statistics_filename
                   if 'statistics_filename' in args else None),]
    if args.subcmd == 'plot_max_coverage':
        kwargs = dict(f5dirs2_opt + nreg_opt + nbase_opt + genome_opts +
                      tbmod_opt + atbmod_opt + dtbmod_opt + datbmod_opt)
        plot_max_coverage(*base_args, **kwargs)
    elif args.subcmd == 'plot_genome_location':
        kwargs = dict(f5dirs2_opt + nbase_opt + genome_opts + glocs_opt +
                      tbmod_opt + atbmod_opt + dtbmod_opt + datbmod_opt)
        plot_genome_locations(*base_args, **kwargs)
    elif args.subcmd == 'plot_motif_centered':
        kwargs = dict(f5dirs2_opt + nreg_opt + nbase_opt + genome_opts +
                      fasta_opt + motif_opt + tbmod_opt + atbmod_opt +
                      dtbmod_opt + datbmod_opt +
                      [('deepest_coverage', args.deepest_coverage),])
        plot_motif_centered(*base_args, **kwargs)
    elif args.subcmd == 'plot_max_difference':
        kwargs = dict(f5dirs2_opt + nreg_opt + nbase_opt + genome_opts +
                      seqfn_opt)
        plot_max_diff(*base_args, **kwargs)
    elif args.subcmd == 'plot_most_significant':
        kwargs = dict(f5dirs2_opt + nreg_opt + nbase_opt + genome_opts +
                      seqfn_opt + qval_opt + statfn_opt + tbmod_opt + atbmod_opt +
                      dtbmod_opt + datbmod_opt +
                      [('stat_order', args.statistic_order)])
        plot_most_signif(*base_args, **kwargs)
    elif args.subcmd == 'plot_motif_with_stats':
        kwargs = dict(f5dirs2_opt + nreg_opt + motif_opt + statfn_opt +
                      tbmod_opt +
                      [('overplot_thresh', args.overplot_threshold),
                       ('context_width', args.num_context),
                       ('stat_order', args.statistic_order),
                       ('num_stats', args.num_statistics)])
        plot_motif_centered_signif(*base_args, **kwargs)
    elif args.subcmd == 'plot_correction':
        kwargs = dict(nobs_opt + nread_opt + [('reg_type', args.region_type),])
        plot_corrections(*base_args, **kwargs)
    elif args.subcmd == 'plot_multi_correction':
        kwargs = dict(nreg_opt + nobs_opt + glocs_opt +
                      [('num_reads_per_plot', args.num_reads),
                       ('include_orig_bcs', args.include_original_basecalls)])
        plot_multi_corrections(*base_args, **kwargs)
    elif args.subcmd == 'cluster_most_significant':
        kwargs = dict(f5dirs2_opt + nreg_opt + qval_opt + nbase_opt +
                      fasta_opt + statfn_opt + rdata_opt +
                      [('num_processes', args.processes),
                       ('slide_span', args.slide_span)])
        cluster_most_signif(*base_args, **kwargs)
    elif args.subcmd == 'plot_kmer':
        kwargs = dict(nread_opt + rdata_opt +
                      [('read_mean', args.read_mean),
                       ('upstrm_bases', args.upstream_bases),
                       ('dnstrm_bases', args.downstream_bases),
                       ('kmer_thresh', args.num_kmer_threshold),
                       ('dont_plot', args.dont_plot)])
        plot_kmer_dist(*base_args, **kwargs)
    elif args.subcmd == 'plot_per_read':
        kwargs = dict(nbase_opt + glocs_opt + tbmod_opt + atbmod_opt +
                      dtbmod_opt + datbmod_opt +
                      [('fm_lag', args.fishers_method_context),
                       ('num_reads', args.num_reads),
                       ('box_center', args.box_center)])
        plot_per_read_mods_genome_location(*base_args, **kwargs)
    else:
        sys.stderr.write('ERROR: Invalid tombo sub-command entered. ' +
                         'Should have been caught by argparse.\n')

    return


if __name__ == '__main__':
    raise NotImplementedError, (
        'This is a module. See commands with `tombo -h`')
