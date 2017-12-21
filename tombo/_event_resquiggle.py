import os, sys

import re
import h5py
import Queue

import numpy as np
np.seterr(all='raise')
import multiprocessing as mp

from subprocess import call
from time import sleep, time
from itertools import repeat
from tempfile import NamedTemporaryFile
from distutils.version import LooseVersion
from collections import defaultdict, namedtuple, Counter

# import tombo functions
import tombo_helper as th

from c_helper import c_valid_cpts, c_valid_cpts_w_cap

VERBOSE = False

OPTIMIZE_RSQGL = False
OPTIMIZE_ALIGN = False

# allow this many times the alignment batch size into the queue of
# reads to be resquiggled
ALIGN_BATCH_MULTIPLIER = 5
PROGRESS_INTERVAL = 100

ALBACORE_TEXT = 'ONT Albacore Sequencing Software'

indelStats = namedtuple('indelStats', ('start', 'end', 'diff'))
indelGroupStats = namedtuple('indelGroupStats',
                             ('start', 'end', 'cpts', 'indels'))
readInfo = namedtuple(
    'readInfo',
    ('ID', 'Subgroup', 'ClipStart', 'ClipEnd',
     'Insertions', 'Deletions', 'Matches', 'Mismatches'))
mapperData = namedtuple('mapperData', ('exe', 'type', 'index'))
# set default index to None
mapperData.__new__.__defaults__ = (None,)

M5_FIELDS = (
    'qName', 'qLength', 'qStart', 'qEnd', 'qStrand',
    'tName', 'tLength', 'tStart', 'tEnd', 'tStrand',
    'score', 'numMatch', 'numMismatch', 'numIns', 'numDel',
    'mapQV', 'qAlignedSeq', 'matchPattern', 'tAlignedSeq')
SAM_FIELDS = (
    'qName', 'flag', 'rName', 'pos', 'mapq',
    'cigar', 'rNext', 'pNext', 'tLen', 'seq', 'qual')
CIGAR_PAT = re.compile('(\d+)([MIDNSHP=X])')
GAP_PAT = re.compile('-+')

#################################################
########## Raw Signal Re-squiggle Code ##########
#################################################

def get_valid_cpts(raw_signal, min_base_obs, num_cpts=None):
    if num_cpts is None:
        return c_valid_cpts(raw_signal, min_base_obs)
    return c_valid_cpts_w_cap(raw_signal, min_base_obs, num_cpts)

def get_indel_groups(
        alignVals, align_segs, raw_signal, min_base_obs, timeout,
        num_cpts_limit):
    def get_all_indels():
        # get genomic sequence for and between each indel
        read_align = ''.join(zip(*alignVals)[0])
        genome_align = ''.join(zip(*alignVals)[1])
        genome_gaps = [(m.start(), m.end()) for m in
                       GAP_PAT.finditer(genome_align)]
        read_gaps = [(m.start(), m.end())
                     for m in GAP_PAT.finditer(read_align)]
        all_indel_locs = sorted(
            genome_gaps + read_gaps +
            [(0,0), (len(read_align), len(read_align))])
        btwn_indel_seqs = [
            genome_align[m_start:m_end] for m_start, m_end in
            zip(zip(*all_indel_locs)[1][:-1],
                zip(*all_indel_locs)[0][1:])]
        # is each indel an ins(ertion) or deletion
        all_is_ins = [read_align[start:end].startswith('-')
                      for start, end in all_indel_locs[1:-1]]
        indel_seqs = [
            genome_align[start:end]
            if is_ins else read_align[start:end]
            for is_ins, (start, end) in
            zip(all_is_ins, all_indel_locs[1:-1])]

        # loop over indels along with sequence before and after in
        # order to check for ambiguous indels
        unambig_indels = []
        curr_read_len = len(btwn_indel_seqs[0])
        for indel_seq, before_seq, after_seq, is_ins in zip(
                indel_seqs, btwn_indel_seqs[:-1], btwn_indel_seqs[1:],
                all_is_ins):
            indel_len = len(indel_seq)
            # if this is an insertion then don't include indel in
            # length to end of indel
            # also extend indel end by 1 in order to check for new
            # breakpoints for neighboring segemnts
            indel_end = curr_read_len + 1 if is_ins else \
                        curr_read_len + indel_len + 1
            indel_diff = indel_len if is_ins else -1 * indel_len
            # indel without ambiguity correction
            # indelStats(curr_read_len - 1, indel_end, indel_diff)

            # extend ambiguous indels
            # only extend up to one position before beginning or end
            # as a one base pad is added outside of indel
            u, d = -1, 0
            while(d < len(after_seq) - 1 and
                  indel_seq[d%indel_len] == after_seq[d]):
                d += 1
            while(u * -1 <= len(before_seq) - 1 and
                  indel_seq[(u%indel_len)-indel_len] == before_seq[u]):
                u -= 1
            unambig_indels.append(indelStats(
                curr_read_len + u, indel_end + d, indel_diff))

            if not is_ins:
                curr_read_len += indel_len
            curr_read_len += len(after_seq)

        return unambig_indels

    def extend_group(indel_group):
        group_start = min(indel.start for indel in indel_group)
        group_stop = max(indel.end for indel in indel_group)
        num_cpts = sum(indel.diff for indel in indel_group
                       ) + group_stop - group_start - 1
        # check that there are enough points to split
        # add an extra set of values to ensure no zero changepoint
        while align_segs[group_stop] - align_segs[group_start] < (
                num_cpts + 2) * min_base_obs:
            num_cpts += int(group_start > 0) + int(
                group_stop < len(align_segs) - 1)
            group_start = max(0, group_start - 1)
            group_stop = min(len(align_segs) - 1, group_stop + 1)
        return group_start, group_stop, num_cpts
    def extend_and_join(indel_group):
        group_start, group_stop, num_cpts = extend_group(indel_group)
        # check if the extension hits the previous group
        while (len(indel_groups) > 0) and (
                group_start <= indel_groups[-1].end):
            indel_group = indel_groups[-1].indels + indel_group
            del indel_groups[-1]
            group_start, group_stop, num_cpts = extend_group(
                indel_group)
        return group_start, group_stop, num_cpts, indel_group
    def get_cpts(group_start, group_stop, num_cpts):
        """
        Get changepoints where the raw difference between min_base_obs
        obs to the left and min_base_obs obs to the right is largest
        while maintaining the min_base_obs between changepoints.
        Still need to test this function for off by one bugs etc.
        """
        if num_cpts_limit is not None and num_cpts > num_cpts_limit:
            raise RuntimeError, ('Reached maximum number of ' +
                                 'changepoints for a single indel')
        try:
            cpts = get_valid_cpts(
                raw_signal[align_segs[group_start]:align_segs[group_stop]],
                min_base_obs, num_cpts)
        # not implemented error returned when fewer cpts found than requested
        except NotImplementedError:
            return None
        cpts.sort()
        return cpts
    def extend_for_cpts(
            group_start, group_stop, num_cpts, indel_group):
        cpts = get_cpts(group_start, group_stop, num_cpts)
        # expand group until a valid set of changepoints can be identified
        while cpts is None:
            num_cpts += int(group_start > 0) + int(
                group_stop < len(align_segs) - 1)
            group_start = max(0, group_start - 1)
            group_stop = min(len(align_segs) - 1, group_stop + 1)
            while (len(indel_groups) > 0) and (
                    group_start <= indel_groups[-1].end):
                indel_group = indel_groups[-1].indels + indel_group
                del indel_groups[-1]
                group_start, group_stop, num_cpts = extend_group(
                    indel_group)
            cpts = get_cpts(group_start, group_stop, num_cpts)

        return (cpts + align_segs[group_start], group_start, group_stop,
                indel_group)

    if timeout is not None:
        timeout_start = time()

    # sort indels in order of start positions
    all_indels = get_all_indels()
    if len(all_indels) == 0:
        return []
    indel_groups = []
    curr_group = [all_indels[0],]
    for indel in all_indels[1:]:
        if timeout is not None and time() - timeout_start > timeout:
            raise RuntimeError, 'Read took too long to re-segment.'
        # check if indel hits current group
        if max(g_indel.end for g_indel in curr_group) >= indel.start:
            curr_group.append(indel)
        else:
            (curr_start, curr_stop, num_cpts,
             curr_group) = extend_and_join(curr_group)
            cpts, curr_start, curr_stop, curr_group = extend_for_cpts(
                curr_start, curr_stop, num_cpts, curr_group)
            # if the indel group still reaches the next indel
            if curr_stop >= indel.start:
                curr_group.append(indel)
            else:
                indel_groups.append(indelGroupStats(
                    curr_start, curr_stop, cpts, curr_group))
                curr_group = [indel,]

    # handle the last indel group if it is not yet included
    if len(indel_groups) == 0 or \
      indel_groups[-1].indels[-1] != all_indels[-1]:
        curr_start, curr_stop, num_cpts, curr_group = extend_and_join(
            curr_group)
        cpts, curr_start, curr_stop, curr_group = extend_for_cpts(
            curr_start, curr_stop, num_cpts, curr_group)
        indel_groups.append(indelGroupStats(
            curr_start, curr_stop, cpts, curr_group))

    return indel_groups

def find_read_start(
        norm_signal, starts_rel_to_read, min_base_obs,
        read_start_rel_to_raw, signal_length, num_obs=2000):
    # get only current starts up to requested num_obs
    begin_read_starts = starts_rel_to_read[
        :np.argmax(starts_rel_to_read >= num_obs)] \
        if starts_rel_to_read[-1] > num_obs else starts_rel_to_read
    if begin_read_starts.shape[0] <= 0:
        return norm_signal, starts_rel_to_read
    signal_cpts = get_valid_cpts(norm_signal[:num_obs], min_base_obs,
                                 begin_read_starts.shape[0])

    # identify the offset which aligns the most signal and read changepoints
    off_by_counts = Counter([
        signal_cpts[np.abs(np.subtract(
            signal_cpts, read_cpt)).argmin()] - read_cpt
        for read_cpt in begin_read_starts]).most_common()

    # if signal starts are ahead of read starts
    if off_by_counts[0][0] > 0:
        offset = off_by_counts[0][0]
        # don't let identified offset push past the end of the read signal
        if (offset + norm_signal.shape[0] +
            read_start_rel_to_raw) >= signal_length:
            offset = (signal_length - read_start_rel_to_raw -
                      norm_signal.shape[0])
        # add fake signal to the end of the read so the file does not have to
        # be queried again
        norm_signal = np.concatenate([norm_signal[offset:],
                                      [norm_signal[-1]] * offset])
        read_start_rel_to_raw += offset
    # if signal starts are behind of read starts
    elif off_by_counts[0][0] < 0:
        offset = off_by_counts[0][0] * -1
        # don't let identified start push start below 0
        if offset > read_start_rel_to_raw:
            offset = read_start_rel_to_raw
        # add fake signal to the end of the read so the file does not have to
        # be queried again
        norm_signal = np.concatenate([[norm_signal[0]] * offset,
                                      norm_signal[:-offset]])
        read_start_rel_to_raw -= offset

    return norm_signal, read_start_rel_to_raw

def resquiggle_read(
        fast5_fn, read_start_rel_to_raw, starts_rel_to_read,
        norm_type, outlier_thresh, alignVals, fix_read_start,
        timeout, num_cpts_limit, genome_loc, read_info,
        basecall_group, corrected_group, compute_sd, pore_model, obs_filter,
        min_base_obs=4, in_place=True, skip_index=False):
    # errors should not happen here since these slotes were checked
    # in alignment function, but old zombie processes might cause
    # problems here
    try:
        fast5_data = h5py.File(fast5_fn, 'r')
        channel_info = th.get_channel_info(fast5_data)

        # extract raw data for this read
        all_raw_signal = fast5_data[
            '/Raw/Reads/'].values()[0]['Signal'].value
        rna = th.is_read_rna(fast5_data)
        if rna:
            all_raw_signal = all_raw_signal[::-1]
        event_means, event_kmers = None, None
        if norm_type == 'pA':
            event_data = fast5_data[
                '/Analyses/' + basecall_group + '/' +
                read_info.Subgroup + '/Events'].value
            event_means = event_data['mean']
            event_kmers = event_data['model_state']
        fast5_data.close()
    except:
        raise NotImplementedError, (
            'Error opening file for re-squiggle. This should have ' +
            'been caught during the alignment phase. Check that there ' +
            'are no other tombo processes or processes accessing ' +
            'these HDF5 files running simultaneously.')

    # normalize signal
    # print read id for resquiggle shift and scale output
    #sys.stdout.write(read_info.ID + "\t")
    norm_signal, scale_values = th.normalize_raw_signal(
        all_raw_signal, read_start_rel_to_raw, starts_rel_to_read[-1],
        norm_type, channel_info, outlier_thresh, pore_model=pore_model,
        event_means=event_means, event_kmers=event_kmers)
    if fix_read_start:
        norm_signal, read_start_rel_to_raw = find_read_start(
            norm_signal, starts_rel_to_read, min_base_obs,
            read_start_rel_to_raw, all_raw_signal.shape[0])

    # group indels that are adjacent for re-segmentation
    indel_groups = get_indel_groups(
        alignVals, starts_rel_to_read, norm_signal, min_base_obs,
        timeout, num_cpts_limit)

    new_segs = []
    prev_stop = 0
    for group_start, group_stop, cpts, group_indels in indel_groups:
        ## add segments from last indel to this one and new segments
        new_segs.append(
            np.append(starts_rel_to_read[prev_stop:group_start+1],
                      cpts))
        prev_stop = group_stop
    # handle end of read
    new_segs.append(starts_rel_to_read[prev_stop:])
    new_segs = np.concatenate(new_segs).astype(np.int32)
    if np.diff(new_segs).min() < 1:
        raise NotImplementedError, (
            'New segments include zero length events.')
    if new_segs[0] < 0:
        raise NotImplementedError, (
            'New segments start with negative index.')
    if new_segs[-1] > norm_signal.shape[0]:
        raise NotImplementedError, (
            'New segments end past raw signal values.')

    # get just from alignVals
    align_seq = ''.join(zip(*alignVals)[1]).replace('-', '')
    if new_segs.shape[0] != len(align_seq) + 1:
        raise ValueError, ('Aligned sequence does not match number ' +
                           'of segments produced.')

    if in_place:
        # create new hdf5 file to hold new read signal
        th.write_new_fast5_group(
            fast5_fn, genome_loc, read_start_rel_to_raw, new_segs, align_seq,
            norm_signal, scale_values, corrected_group, read_info.Subgroup,
            norm_type, outlier_thresh, compute_sd, alignVals, read_info,
            starts_rel_to_read, rna)
    else:
        # create new hdf5 file to hold corrected read events
        pass

    if not skip_index:
        return th.prep_index_data(
            fast5_fn, genome_loc, read_start_rel_to_raw, new_segs,
            corrected_group, read_info.Subgroup, rna, obs_filter)

    return

def resquiggle_worker(
        basecalls_q, failed_reads_q, index_q, basecall_group, corrected_group,
        norm_type, outlier_thresh, timeout, num_cpts_limit, compute_sd,
        pore_model, obs_filter):
    num_processed = 0
    skip_index = index_q is None
    if not skip_index: proc_index_data = []
    while True:
        try:
            fast5_fn, sgs_align_data = basecalls_q.get(block=False)
            # None values placed in queue when all files have
            # been processed
            if fast5_fn is None: break
        except Queue.Empty:
            sleep(1)
            continue

        num_processed += 1
        if VERBOSE and num_processed % PROGRESS_INTERVAL == 0:
            sys.stderr.write('.')
            sys.stderr.flush()
        # process different read subgroups separately so that the same
        # file is never open simultaneously
        for align_data in sgs_align_data:
            (alignVals, genome_loc, starts_rel_to_read,
             read_start_rel_to_raw, read_info, fix_read_start) = align_data
            try:
                index_data = resquiggle_read(
                    fast5_fn, read_start_rel_to_raw, starts_rel_to_read,
                    norm_type, outlier_thresh, alignVals, fix_read_start,
                    timeout, num_cpts_limit, genome_loc, read_info,
                    basecall_group, corrected_group, compute_sd,
                    pore_model, obs_filter, skip_index=skip_index)
                if not skip_index:
                    proc_index_data.append(index_data)
            except Exception as e:
                # uncomment to identify mysterious errors
                #raise
                try:
                    th.write_error_status(
                        fast5_fn, corrected_group, read_info.Subgroup, str(e))
                except:
                    pass
                failed_reads_q.put((
                    str(e), read_info.Subgroup + ' :: ' + fast5_fn))

    if not skip_index: index_q.put(proc_index_data)

    return

if OPTIMIZE_RSQGL:
    resquiggle_wrapper = resquiggle_worker
    def resquiggle_worker(*args):
        import cProfile
        cProfile.runctx('resquiggle_wrapper(*args)', globals(), locals(),
                        filename='resquiggle_main.prof')
        return



############################################
########## Genomic Alignment Code ##########
############################################

def fix_raw_starts_for_clipped_bases(
        start_clipped_bases, end_clipped_bases, starts_rel_to_read,
        read_start_rel_to_raw):
    if start_clipped_bases > 0:
        start_clipped_obs = starts_rel_to_read[start_clipped_bases]
        starts_rel_to_read = starts_rel_to_read[
            start_clipped_bases:] - start_clipped_obs
        read_start_rel_to_raw += start_clipped_obs

    if end_clipped_bases > 0:
        starts_rel_to_read = starts_rel_to_read[
            :-1 * end_clipped_bases]

    return starts_rel_to_read, read_start_rel_to_raw

def fix_all_clipped_bases(batch_align_data, batch_reads_data):
    clip_fix_align_data = []
    for read_fn_sg, (
            alignVals, genome_loc, start_clipped_bases,
            end_clipped_bases) in batch_align_data.iteritems():
        (read_start_rel_to_raw, starts_rel_to_read, basecalls,
         channel_info, read_id, fix_read_start) = batch_reads_data[read_fn_sg]
        # fix raw start positions to match bases clipped in mapping
        starts_rel_to_read, read_start_rel_to_raw \
            = fix_raw_starts_for_clipped_bases(
                start_clipped_bases, end_clipped_bases,
                starts_rel_to_read, read_start_rel_to_raw)

        bc_subgroup, fast5_fn = read_fn_sg.split(th.FASTA_NAME_JOINER)
        num_ins, num_del, num_match, num_mismatch = 0, 0, 0, 0
        for rb, gb in alignVals:
            if rb == '-':
                num_del += 1
            elif gb == '-':
                num_ins += 1
            elif rb == gb:
                num_match += 1
            else:
                num_mismatch += 1
        read_info = readInfo(
            read_id, bc_subgroup, start_clipped_bases, end_clipped_bases,
            num_ins, num_del, num_match, num_mismatch)

        clip_fix_align_data.append((fast5_fn, (
            alignVals, genome_loc, starts_rel_to_read,
            read_start_rel_to_raw, read_info, fix_read_start)))

    return clip_fix_align_data

def clip_m5_alignment(alignVals, start, strand, chrm):
    # clip read to first matching bases
    start_clipped_read_bases = 0
    start_clipped_genome_bases = 0
    start_clipped_align_bases = 0
    r_base, g_base = alignVals[0]
    while r_base == '-' or g_base == '-':
        start_clipped_read_bases += int(r_base != '-')
        start_clipped_genome_bases += int(g_base != '-')
        start_clipped_align_bases += 1
        r_base, g_base = alignVals[start_clipped_align_bases]

    end_clipped_read_bases = 0
    end_clipped_genome_bases = 0
    end_clipped_align_bases = 0
    r_base, g_base = alignVals[-1]
    while r_base == '-' or g_base == '-':
        end_clipped_read_bases += int(r_base != '-')
        end_clipped_genome_bases += int(g_base != '-')
        end_clipped_align_bases += 1
        r_base, g_base = alignVals[-1 * (end_clipped_align_bases + 1)]

    alignVals = alignVals[start_clipped_align_bases:]
    if end_clipped_align_bases > 0:
        alignVals = alignVals[:-1*end_clipped_align_bases]

    if strand == '+' and start_clipped_genome_bases > 0:
        genome_loc = th.genomeLoc(
            start + start_clipped_genome_bases, '+', chrm)
    elif strand == '-' and end_clipped_genome_bases > 0:
        genome_loc = th.genomeLoc(
            start + end_clipped_genome_bases, '-', chrm)
    else:
        genome_loc = th.genomeLoc(start, strand, chrm)

    return alignVals, start_clipped_read_bases, \
        end_clipped_read_bases, genome_loc

def parse_m5_record(r_m5_record):
    if r_m5_record['tStrand'] != '+':
        raise NotImplementedError, (
            'Mapping indicates negative strand reference mapping.')

    if r_m5_record['qStrand'] == "+":
        alignVals = zip(r_m5_record['qAlignedSeq'],
                        r_m5_record['tAlignedSeq'])
    else:
        alignVals = zip(th.rev_comp(r_m5_record['qAlignedSeq']),
                        th.rev_comp(r_m5_record['tAlignedSeq']))

    alignVals, start_clipped_bases, end_clipped_bases, genome_loc \
        = clip_m5_alignment(
            alignVals, int(r_m5_record['tStart']),
            r_m5_record['qStrand'], r_m5_record['tName'])

    return (alignVals, genome_loc, start_clipped_bases,
            end_clipped_bases)

def parse_m5_output(align_output, batch_reads_data):
    alignments = dict(
        (read_fn_sg, None) for read_fn_sg in batch_reads_data.keys())
    for line in align_output:
        r_m5_record = dict(zip(M5_FIELDS, line.strip().split()))
        if len(r_m5_record) != len(M5_FIELDS):
            continue
        # store the alignment if none is stored for this read or
        # if this read has the lowest map quality thus far
        if alignments[r_m5_record['qName']] is None or \
           int(alignments[r_m5_record['qName']]['score']) < \
           int(r_m5_record['score']):
            alignments[r_m5_record['qName']] = r_m5_record

    batch_align_failed_reads = []
    batch_align_data = {}
    for read_fn_sg, r_m5_record in alignments.iteritems():
        if r_m5_record is None:
            batch_align_failed_reads.append(
                ('Alignment not produced.', read_fn_sg))
        else:
            try:
                batch_align_data[read_fn_sg] = parse_m5_record(
                    r_m5_record)
            except Exception as e:
                batch_align_failed_reads.append((str(e), read_fn_sg))

    return batch_align_failed_reads, batch_align_data

def parse_sam_record(r_sam_record, genome_index):
    def parse_cigar(strand):
        # parse cigar string
        cigar = [
            (int(reg_len), reg_type) for reg_len, reg_type in
            CIGAR_PAT.findall(r_sam_record['cigar'])]
        if len(cigar) < 1:
            raise RuntimeError, 'Invalid cigar string produced.'

        if strand == '-':
            cigar = cigar[::-1]

        return cigar

    def get_qseq(cigar, strand):
        # record clipped bases and remove from query seq as well as cigar
        qSeq = r_sam_record['seq'] if strand == '+' else th.rev_comp(
            r_sam_record['seq'])
        start_clipped_bases = 0
        end_clipped_bases = 0
        # handle clipping elements (H and S)
        if cigar[0][1] == 'H':
            start_clipped_bases += cigar[0][0]
            cigar = cigar[1:]
        if cigar[-1][1] == 'H':
            end_clipped_bases += cigar[-1][0]
            cigar = cigar[:-1]
        if cigar[0][1] == 'S':
            start_clipped_bases += cigar[0][0]
            qSeq = qSeq[cigar[0][0]:]
            cigar = cigar[1:]
        if cigar[-1][1] == 'S':
            end_clipped_bases += cigar[-1][0]
            qSeq = qSeq[:-cigar[-1][0]]
            cigar = cigar[:-1]

        return qSeq, start_clipped_bases, end_clipped_bases, cigar

    def get_tseq(qSeq, start_clipped_bases, end_clipped_bases, cigar, strand):
        tLen = sum([reg_len for reg_len, reg_type in cigar
                    if reg_type in 'MDN=X'])
        tSeq = genome_index[r_sam_record['rName']][
            int(r_sam_record['pos']) - 1:
            int(r_sam_record['pos']) + tLen - 1]
        if strand == '-': tSeq = th.rev_comp(tSeq)

        # check that cigar starts and ends with matched bases
        while cigar[0][1] not in 'M=X':
            if cigar[0][1] in 'ND':
                tSeq = tSeq[cigar[0][0]:]
            else:
                qSeq = qSeq[cigar[0][0]:]
                start_clipped_bases += cigar[0][0]
            cigar = cigar[1:]
        while cigar[-1][1] not in 'M=X':
            if cigar[-1][1] in 'ND':
                tSeq = tSeq[:-cigar[-1][0]]
            else:
                qSeq = qSeq[:-cigar[-1][0]]
                end_clipped_bases += cigar[-1][0]
            cigar = cigar[:-1]

        qLen = sum([reg_len for reg_len, reg_type in cigar
                    if reg_type in 'MI=X'])
        assert len(qSeq) == qLen, 'Read sequence from SAM and ' + \
            'cooresponding cigar string do not agree.'

        return tSeq, qSeq, start_clipped_bases, end_clipped_bases, cigar

    def get_align_vals(tSeq, qSeq, cigar, strand):
        alignVals = []
        tPos, qPos = 0, 0
        for reg_len, reg_type in cigar:
            if reg_type in 'M=X':
                alignVals.extend(zip(qSeq[qPos:qPos+reg_len],
                                     tSeq[tPos:tPos+reg_len]))
                tPos += reg_len
                qPos += reg_len
            elif reg_type in 'IP':
                alignVals.extend(zip(qSeq[qPos:qPos+reg_len], repeat('-')))
                qPos += reg_len
            else:
                alignVals.extend(zip(repeat('-'), tSeq[tPos:tPos+reg_len]))
                tPos += reg_len

        return alignVals

    strand = '-' if int(r_sam_record['flag']) & 0x10 else '+'
    cigar = parse_cigar(strand)
    qSeq, start_clipped_bases, end_clipped_bases, cigar = get_qseq(cigar, strand)
    tSeq, qSeq, start_clipped_bases, end_clipped_bases, cigar = get_tseq(
        qSeq, start_clipped_bases, end_clipped_bases, cigar, strand)
    alignVals = get_align_vals(tSeq, qSeq, cigar, strand)

    return (alignVals, th.genomeLoc(
        int(r_sam_record['pos']) - 1, strand, r_sam_record['rName']),
            start_clipped_bases, end_clipped_bases)

def parse_sam_output(align_output, batch_reads_data, genome_index):
    # create dictionary with empty slot to each read
    alignments = dict(
        (read_fn_sg, None) for read_fn_sg in batch_reads_data.keys())
    for line in align_output:
        if line.startswith('@'): continue
        r_sam_record = dict(zip(SAM_FIELDS, line.strip().split()))
        if len(r_sam_record) < len(SAM_FIELDS): continue
        if r_sam_record['rName'] == '*': continue
        # store the alignment if none is stored for this read or
        # if this read has the lowest map quality thus far
        qName = r_sam_record['qName'].replace(th.FN_SPACE_FILLER, ' ')
        if alignments[qName] is None or \
           int(alignments[qName]['mapq']) < \
           int(r_sam_record['mapq']):
            alignments[qName] = r_sam_record

    batch_align_failed_reads = []
    batch_align_data = {}
    for read_fn_sg, r_sam_record in alignments.iteritems():
        if r_sam_record is None:
            batch_align_failed_reads.append(
                ('Alignment not produced (if all reads failed ' +
                 'check for index files).', read_fn_sg))
        else:
            try:
                batch_align_data[read_fn_sg] = parse_sam_record(
                    r_sam_record, genome_index)
            except Exception as e:
                #raise
                batch_align_failed_reads.append((str(e), read_fn_sg))

    return batch_align_failed_reads, batch_align_data

def prep_graphmap_options(
        genome_fn, read_fn, out_fn, output_format, num_align_ps):
    return ['align', '-r', genome_fn, '-d', read_fn, '-o', out_fn,
            '-L', output_format, '-t', str(num_align_ps)]

def prep_bwa_mem_options(genome_fn, read_fn, num_align_ps):
    return ['mem', '-x', 'ont2d', '-v', '1', '-t', str(num_align_ps),
            genome_fn, read_fn]

def prep_minimap2_options(genome_fn, read_fn, num_align_ps, index_fn):
    mapper_genome = genome_fn if index_fn is None else index_fn
    return ['-ax', 'map-ont', '-t', str(num_align_ps), mapper_genome, read_fn]

def align_to_genome(batch_reads_data, genome_fn, mapper_data, genome_index,
                    num_align_ps, output_format='sam'):
    # prepare fasta text with batch reads
    batch_reads_fasta = ''
    for read_fn_sg, (_, _, basecalls, _, _, _) in \
        batch_reads_data.iteritems():
        # note spaces aren't allowed in read names so replace with
        # vertical bars and undo to retain file names
        batch_reads_fasta += ">" + read_fn_sg.replace(' ', th.FN_SPACE_FILLER) + \
                             '\n' + ''.join(basecalls) + '\n'

    read_fp = NamedTemporaryFile(suffix='.fasta')
    read_fp.write(batch_reads_fasta)
    read_fp.flush()
    out_fp = NamedTemporaryFile()

    # optionally suppress output from mapper with devnull sink
    with open(os.devnull, 'w') as FNULL:
        if mapper_data.type == 'graphmap':
            mapper_options = prep_graphmap_options(
                genome_fn, read_fp.name, out_fp.name,
                output_format, num_align_ps)
            stdout_sink = FNULL
        elif mapper_data.type == 'bwa_mem':
            mapper_options = prep_bwa_mem_options(
                genome_fn, read_fp.name, num_align_ps)
            stdout_sink = out_fp
        elif mapper_data.type == 'minimap2':
            mapper_options = prep_minimap2_options(
                genome_fn, read_fp.name, num_align_ps, mapper_data.index)
            stdout_sink = out_fp
        else:
            raise RuntimeError, 'Mapper not supported.'

        try:
            exitStatus = call([mapper_data.exe,] + mapper_options,
                              stdout=stdout_sink, stderr=FNULL)
            out_fp.seek(0)
            align_output = out_fp.readlines()
            # close files here so that they persist until
            # after basecalling is finished
            read_fp.close()
            out_fp.close()
        except:
            # whole mapping call failed so all reads failed
            return ([(
                'Problem running/parsing genome mapper. ' +
                'Ensure you have a compatible version installed.' +
                'Potentially failed to locate BWA index files.',
                read_fn_sg) for read_fn_sg
                     in batch_reads_data.keys()], [])

    if output_format == 'sam':
        batch_parse_failed_reads, batch_align_data = parse_sam_output(
            align_output, batch_reads_data, genome_index)
    elif output_format == 'm5':
        batch_parse_failed_reads, batch_align_data = parse_m5_output(
            align_output, batch_reads_data)
    else:
        raise RuntimeError, 'Mapper output type not supported.'

    clip_fix_align_data = fix_all_clipped_bases(
        batch_align_data, batch_reads_data)

    return batch_parse_failed_reads, clip_fix_align_data

def fix_stay_states(
        called_dat, starts_rel_to_read, basecalls,
        read_start_rel_to_raw, rna):
    move_states = called_dat['move'][1:] > 0
    if rna:
        move_states = move_states[::-1]
    start_clip = 0
    event_change_state = move_states[0]
    while not event_change_state:
        if start_clip >= len(move_states) - 2:
            raise RuntimeError, (
                'Read is composed entirely of stay model ' +
                'states and cannot be processed')
        start_clip += 1
        event_change_state = move_states[start_clip]
    end_clip = 0
    event_change_state = move_states[-1]
    while not event_change_state:
        end_clip += 1
        event_change_state = move_states[-(end_clip+1)]

    # clip all applicable data structures
    move_states = move_states[start_clip:]
    starts_rel_to_read = starts_rel_to_read[start_clip:]
    basecalls = basecalls[start_clip:]
    if end_clip > 0:
        move_states = move_states[:-end_clip]
        starts_rel_to_read = starts_rel_to_read[:-end_clip]
        basecalls = basecalls[:-end_clip]
    if start_clip > 0:
        start_clip_obs = starts_rel_to_read[0]
        starts_rel_to_read = starts_rel_to_read - start_clip_obs
        read_start_rel_to_raw += start_clip_obs

    # now actually remove internal stay states
    move_states = np.insert(
        move_states, (0, len(move_states) - 1), True)
    starts_rel_to_read = starts_rel_to_read[move_states]
    basecalls = basecalls[move_states[:-1]]

    return starts_rel_to_read, basecalls, read_start_rel_to_raw

def get_read_data(fast5_fn, basecall_group, basecall_subgroup):
    try:
        fast5_data = h5py.File(fast5_fn, 'r')
    except:
        raise NotImplementedError, (
            'Error opening file for alignment. This should have ' +
            'been caught during the HDF5 prep phase. Check that there ' +
            'are no other tombo processes or processes accessing ' +
            'these HDF5 files running simultaneously.')

    try:
        # get albacore version, or if not specified set to 0.0
        albacore_version = LooseVersion(fast5_data[
            '/Analyses/' + basecall_group].attrs['version']
            if 'version' in fast5_data['/Analyses/' +
                                       basecall_group].attrs else "0.0")
        called_dat = fast5_data[
            '/Analyses/' + basecall_group + '/' + basecall_subgroup +
            '/Events'].value
    except:
        raise RuntimeError, (
            'No events or corrupted events in file. Likely a ' +
            'segmentation error or mis-specified basecall-' +
            'subgroups (--2d?).')
    rna = th.is_read_rna(fast5_data)
    try:
        raw_attrs = dict(
            fast5_data['/Raw/Reads/'].values()[0].attrs.items())
        if rna:
            raw_len = fast5_data['/Raw/Reads/'].values()[0]['Signal'].shape[0]
    except:
        raise RuntimeError, (
            'Raw data is not stored in Raw/Reads/Read_[read#] so ' +
            'new segments cannot be identified.')

    try:
        channel_info = th.get_channel_info(fast5_data)
        fast5_data.close()
    except:
        raise RuntimeError, (
            'Error getting channel information and closing fast5 file.')

    read_id = raw_attrs['read_id']

    if albacore_version >= LooseVersion("2.0"):
        read_start_rel_to_raw = called_dat['start'][0].astype(np.int64)
    else:
        abs_event_start = np.round(
            called_dat['start'][0].astype(np.float64) *
            channel_info.sampling_rate).astype(np.uint64)
        if abs_event_start < raw_attrs['start_time']:
            # floating point errors can result in apparent read start before
            # the raw array, so set to zero
            read_start_rel_to_raw = 0
        else:
            read_start_rel_to_raw = int(
                abs_event_start - raw_attrs['start_time'])

    # check albacore version to determine which method to extract starts
    # relative to the read.
    # Before albacore version 1.0, events could be skipped so start times
    # should be used. Since then events are not removed so the length
    # slot is more reliable since it will have greater floating point
    # precision. Relevant discussion on community forum here:
    # https://community.nanoporetech.com/posts/albacore-zero-length-even
    if albacore_version < LooseVersion("1.0"):
        last_event = called_dat[-1]
        # convert starts to float64 to minimize floating point errors
        starts_rel_to_read = np.append(
            called_dat['start'], last_event['start'] +
            last_event['length']).astype(np.float64)
        # round to float64 to minimize floating point errors
        starts_rel_to_read = np.round(
            starts_rel_to_read *
            channel_info.sampling_rate).astype('int_') - abs_event_start
        kmer_reference_offset = 2
        fix_read_start = False
    elif albacore_version < LooseVersion("2.0"):
        # compute event starts from length slot as start slot is less
        # reliable due to storage as a float32
        starts_rel_to_read = np.cumsum(np.concatenate(
            [[0,], np.round(called_dat['length'] *
                            channel_info.sampling_rate).astype('int_')]))
        kmer_reference_offset = 2
        # Fix floating point errors in abs_event_start by comparing to
        # potential breakpoints using resquiggle criterion
        # don't actually fix here to avoid reading raw signal twice
        fix_read_start = True
    # assumed events format for other basecallers should be the most recent
    # albacore version.
    else:
        last_event = called_dat[-1]
        starts_rel_to_read = np.append(
            called_dat['start'],
            last_event['start'] + last_event['length']).astype('int_') - \
            read_start_rel_to_raw
        # raw basecalling caused the dominant kmer reference base to
        # move to the second position (from the third previously)
        # but raw was intorduced into rna basecalling one minor release later
        if rna and albacore_version < LooseVersion("2.1"):
            kmer_reference_offset = 2
        else:
            kmer_reference_offset = 1
        fix_read_start = False
    basecalls = np.array([event_state[kmer_reference_offset]
                          for event_state in called_dat['model_state']])

    if rna:
        starts_rel_to_read = -1 * (
            starts_rel_to_read[::-1] + read_start_rel_to_raw -
            raw_len)
        read_start_rel_to_raw = starts_rel_to_read[0]
        # due to floating point time value errors this value may extend beyond
        # the raw signal
        if read_start_rel_to_raw < 0:
            starts_rel_to_read -= read_start_rel_to_raw
            read_start_rel_to_raw = 0
        else:
            starts_rel_to_read = starts_rel_to_read - read_start_rel_to_raw
        basecalls = basecalls[::-1]

    if any(len(vals) <= 1 for vals in (
        starts_rel_to_read, basecalls,
        called_dat['model_state'])):
        raise NotImplementedError, (
            'One or no segments or signal present in read.')
    if min(np.diff(starts_rel_to_read)) < 1:
        raise NotImplementedError, (
            'Zero length event present in input data.')

    # remove stay states from the base caller
    (starts_rel_to_read, basecalls,
     read_start_rel_to_raw) = fix_stay_states(
         called_dat, starts_rel_to_read, basecalls,
         read_start_rel_to_raw, rna)

    return (read_start_rel_to_raw, starts_rel_to_read, basecalls,
            channel_info, read_id, fix_read_start)

def align_and_parse(
        fast5s_to_process, genome_fn, mapper_data, genome_index,
        basecall_group, basecall_subgroups, num_align_ps):
    batch_reads_data = {}
    batch_get_data_failed_reads = []
    for fast5_fn in fast5s_to_process:
        for bc_subgroup in basecall_subgroups:
            try:
                read_data = get_read_data(
                    fast5_fn, basecall_group, bc_subgroup)
                batch_reads_data[
                    bc_subgroup + th.FASTA_NAME_JOINER + fast5_fn] = read_data
            except Exception as e:
                # uncomment to identify mysterious errors
                #raise
                batch_get_data_failed_reads.append((
                    str(e), bc_subgroup + th.FASTA_NAME_JOINER + fast5_fn))

    batch_align_failed_reads, batch_align_data = align_to_genome(
        batch_reads_data, genome_fn, mapper_data,
        genome_index, num_align_ps)
    # regroup reads by filename (for 2D reads to be processed together
    # and avoid the same HDF5 file being opened simultaneuously)
    fn_batch_align_data = defaultdict(list)
    for fast5_fn, sg_align_data in batch_align_data:
        fn_batch_align_data[fast5_fn].append(sg_align_data)
    # uncomment to identify mysterious errors
    #print "Get data errors: " + str(batch_get_data_failed_reads)
    #print "Align read errors: " + str(batch_align_failed_reads)

    return (batch_get_data_failed_reads + batch_align_failed_reads,
            fn_batch_align_data)

def align_reads(
        fast5_batch, genome_fn, mapper_data, genome_index,
        basecall_group, basecall_subgroups, corrected_group,
        basecalls_q, overwrite, num_align_ps, in_place=True):
    batch_prep_failed_reads = []
    fast5s_to_process = []
    for fast5_fn in fast5_batch:
        prep_result = th.prep_fast5(
            fast5_fn, corrected_group, overwrite, in_place, basecall_group)
        if prep_result is None:
            fast5s_to_process.append(fast5_fn)
        else:
            batch_prep_failed_reads.append(prep_result)

    batch_align_failed_reads, batch_align_data = align_and_parse(
        fast5s_to_process, genome_fn, mapper_data, genome_index,
        basecall_group, basecall_subgroups, num_align_ps)
    for fast5_fn, sgs_align_data in batch_align_data.iteritems():
        basecalls_q.put((fast5_fn, sgs_align_data))
    # uncomment to identify mysterious errors
    #print "Prep reads fail: " + str(batch_prep_failed_reads)
    #print "Align reads fail: " + str(batch_align_failed_reads)

    return batch_prep_failed_reads + batch_align_failed_reads

def alignment_worker(
        fast5_q, basecalls_q, failed_reads_q, genome_fn,
        mapper_data, basecall_group, basecall_subgroups,
        corrected_group, overwrite, num_align_ps):
    # this is only needed for sam output format (not m5)
    genome_index = th.parse_fasta(genome_fn)
    while not fast5_q.empty():
        try:
            fast5_batch = fast5_q.get(block=False)
        except Queue.Empty:
            break

        batch_failed_reads = align_reads(
            fast5_batch, genome_fn, mapper_data,
            genome_index, basecall_group, basecall_subgroups,
            corrected_group, basecalls_q, overwrite, num_align_ps)
        for failed_read in batch_failed_reads:
            try:
                sg_fn = failed_read[1].split(th.FASTA_NAME_JOINER)
                if len(sg_fn) == 2:
                    subgroup, fast5_fn = sg_fn
                else:
                    subgroup, fast5_fn = None, sg_fn
                th.write_error_status(
                    fast5_fn, corrected_group, subgroup, failed_read[0])
            except:
                pass
            failed_reads_q.put(failed_read)

    return

if OPTIMIZE_ALIGN:
    alignment_wrapper = alignment_worker
    def alignment_worker(*args):
        import cProfile
        cProfile.runctx('alignment_wrapper(*args)', globals(), locals(),
                        filename='resquiggle_align.prof')
        return


def resquiggle_all_reads(
        fast5_fns, genome_fn, mapper_data,
        basecall_group, basecall_subgroups, corrected_group, norm_type,
        outlier_thresh, timeout, num_cpts_limit, overwrite,
        align_batch_size, num_align_ps, align_threads_per_proc,
        num_resquiggle_ps, compute_sd, pore_model, skip_index, obs_filter):
    manager = mp.Manager()
    fast5_q = manager.Queue()
    # set maximum number of parsed basecalls to sit in the middle queue
    basecalls_q = manager.Queue(
        align_batch_size * ALIGN_BATCH_MULTIPLIER)
    failed_reads_q = manager.Queue()
    index_q = manager.Queue() if not skip_index else None
    num_reads = 0
    fast5_batch = []
    for fast5_fn in fast5_fns:
        num_reads += 1
        fast5_batch.append(fast5_fn)
        # put batches of reads in queue
        if num_reads % align_batch_size == 0:
            fast5_q.put(fast5_batch)
            fast5_batch = []
    if len(fast5_batch) > 0:
        fast5_q.put(fast5_batch)

    align_args = (
        fast5_q, basecalls_q, failed_reads_q, genome_fn,
        mapper_data, basecall_group, basecall_subgroups,
        corrected_group, overwrite, align_threads_per_proc)
    align_ps = []
    for p_id in xrange(num_align_ps):
        p = mp.Process(target=alignment_worker, args=align_args)
        p.start()
        align_ps.append(p)

    rsqgl_args = (basecalls_q, failed_reads_q, index_q, basecall_group,
                  corrected_group, norm_type, outlier_thresh, timeout,
                  num_cpts_limit, compute_sd, pore_model, obs_filter)
    resquiggle_ps = []
    for p_id in xrange(num_resquiggle_ps):
        p = mp.Process(target=resquiggle_worker, args=rsqgl_args)
        p.start()
        resquiggle_ps.append(p)

    if VERBOSE: sys.stderr.write(
            'Correcting ' + str(num_reads) + ' files with ' +
            str(len(basecall_subgroups)) + ' subgroup(s)/read(s) ' +
            'each (Will print a dot for each ' + str(PROGRESS_INTERVAL) +
            ' reads completed).\n')
    failed_reads = defaultdict(list)
    all_index_data = []
    while any(p.is_alive() for p in align_ps):
        try:
            errorType, fn = failed_reads_q.get(block=False)
            failed_reads[errorType].append(fn)
        except Queue.Empty:
            sleep(1)
            continue

    # add None entried to basecalls_q to indicate that all reads have
    # been basecalled and processed
    for _ in xrange(num_resquiggle_ps):
        basecalls_q.put((None, None))

    while any(p.is_alive() for p in resquiggle_ps):
        try:
            errorType, fn = failed_reads_q.get(block=False)
            failed_reads[errorType].append(fn)
        except Queue.Empty:
            try:
                proc_index_data = index_q.get(block=False)
                all_index_data.extend(proc_index_data)
            except Queue.Empty:
                sleep(1)
                continue

    # empty any entries left in queue after processes have finished
    while not failed_reads_q.empty():
        errorType, fn = failed_reads_q.get(block=False)
        failed_reads[errorType].append(fn)
    while not index_q.empty():
        proc_index_data = index_q.get(block=False)
        all_index_data.extend(proc_index_data)

    # print newline after read progress dots
    if VERBOSE: sys.stderr.write('\n')

    return dict(failed_reads), all_index_data

def check_for_albacore(files, basecall_group, num_reads=50):
    has_albacore = False
    for fast5_fn in np.random.choice(files, num_reads):
        try:
            fast5_data = h5py.File(fast5_fn, 'r')
            if fast5_data['/Analyses/' + basecall_group].attrs['name'] == \
               ALBACORE_TEXT:
                has_albacore = True
                break
        except:
            continue

    if not has_albacore:
        sys.stderr.write(
            '******** WARNING ********* The provided FAST5 files do not ' +
            'appear to contain albacore basecalling events. ' +
            'tombo is only tested on albacore formatted results ' +
            'and other basecallers may not produce desired results.\n')

    return

def event_resquiggle_main(args):
    global VERBOSE
    VERBOSE = not args.quiet
    th.VERBOSE = VERBOSE

    # currently required, but adding new mappers shortly
    if all(map_exe is None for map_exe in (
            args.minimap2_executable, args.bwa_mem_executable,
            args.graphmap_executable)):
        sys.stderr.write(
            '*' * 60 + '\nERROR: Must provide either a ' + \
            'minimap2, graphmap or bwa-mem executable.\n' + '*' * 60 + '\n')
        sys.exit()
    if args.minimap2_executable is not None:
        mapper_data = mapperData(
            args.minimap2_executable, 'minimap2', args.minimap2_index)
    elif args.bwa_mem_executable is not None:
        mapper_data = mapperData(args.bwa_mem_executable, 'bwa_mem')
    else:
        mapper_data = mapperData(args.graphmap_executable, 'graphmap')

    if VERBOSE: sys.stderr.write('Getting file list.\n')
    try:
        if not os.path.isdir(args.fast5_basedir):
            sys.stderr.write(
                '*' * 60 + '\nERROR: Provided --fast5-basedir is ' +
                'not a directory.\n' + '*' * 60 + '\n')
            sys.exit()
        fast5_basedir = (
            args.fast5_basedir if args.fast5_basedir.endswith('/') else
            args.fast5_basedir + '/')
        files = th.get_files_list(fast5_basedir)
        if not args.skip_index:
            index_fn = th.get_index_fn(fast5_basedir, args.corrected_group)
            if os.path.exists(index_fn): os.remove(index_fn)
    except OSError:
        sys.stderr.write(
            '*' * 60 + '\nERROR: Reads base directory, a sub-directory ' +
            'or an old (hidden) index file does not appear to be ' +
            'accessible. Check directory permissions.\n' + '*' * 60 + '\n')
        sys.exit()
    if len(files) < 1:
        sys.stderr.write(
            '*' * 60 + '\nERROR: No files identified in the specified ' +
            'directory or within immediate subdirectories.\n' + '*' * 60 + '\n')
        sys.exit()

    check_for_albacore(files, args.basecall_group)

    outlier_thresh = args.outlier_threshold if (
        args.outlier_threshold > 0) else None

    # resolve processor and thread arguments
    num_proc = 2 if args.processes < 2 else args.processes
    align_threads_per_proc = \
      max(int(num_proc / (2 * args.align_processes)), 1) \
      if args.align_threads_per_process is None else \
      args.align_threads_per_process
    num_resquiggle_ps = int(num_proc / 2) \
      if args.resquiggle_processes is None \
      else args.resquiggle_processes

    # whether or not to skip SD calculation due to time
    compute_sd = args.include_event_stdev

    # parse pore model if k-mer conditional corrected pA values
    # are requested
    pore_model = None
    if args.normalization_type == 'pA':
        pore_model = th.parse_pore_model(args.pore_model_filename)

    obs_filter = th.parse_obs_filter(args.obs_per_base_filter) \
                 if 'obs_per_base_filter' in args else None

    failed_reads, all_index_data = resquiggle_all_reads(
        files, args.genome_fasta, mapper_data,
        args.basecall_group, args.basecall_subgroups,
        args.corrected_group, args.normalization_type, outlier_thresh,
        args.timeout, args.cpts_limit, args.overwrite,
        args.alignment_batch_size, args.align_processes,
        align_threads_per_proc, num_resquiggle_ps, compute_sd,
        pore_model, args.skip_index, obs_filter)
    if not args.skip_index:
        th.write_index_file(all_index_data, index_fn, fast5_basedir)
    fail_summary = [(err, len(fns)) for err, fns in failed_reads.items()]
    if len(fail_summary) > 0:
        total_num_failed = sum(zip(*fail_summary)[1])
        sys.stderr.write('Failed reads summary (' + str(total_num_failed) +
                         ' total failed):\n' + '\n'.join(
                             "\t" + err + " :\t" + str(n_fns)
                             for err, n_fns in sorted(fail_summary)) + '\n')
    else:
        sys.stderr.write('All reads successfully re-squiggled!\n')
    if args.failed_reads_filename is not None:
        with open(args.failed_reads_filename, 'w') as fp:
            fp.write('\n'.join((
                err + '\t' + ', '.join(fns)
                for err, fns in failed_reads.items())) + '\n')

    return

def args_and_main():
    import _option_parsers
    event_resquiggle_main(
        _option_parsers.get_resquiggle_parser().parse_args())
    return

if __name__ == '__main__':
    args_and_main()
