from __future__ import division, unicode_literals, absolute_import

from builtins import int, range, dict, map, zip

import io
import os
import re
import sys
import queue

# Future warning from cython in h5py
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
import h5py

from tqdm import tqdm
from time import sleep
from itertools import islice
from multiprocessing import Process, Queue, Pipe

if sys.version_info[0] > 2:
    unicode = str

from . import tombo_helper as th


VERBOSE = False

_MAX_QUEUE_SIZE = 1000
_ITER_QUEUE_LIMIT = 1000
_PROC_UPDATE_INTERVAL = 100

_MAX_FASTQ_QUEUE_SIZE = 10000
_SEQ_SUMMARY_FN_FIELD = 'filename'
_SEQ_SUMMARY_ID_FIELD = 'read_id'

# warning messages for annotate with fastqs over multiple processes,
# requiring passing warning codes to only print warning once.
_WARN_ID_VAL = 'ids'
_WARN_IO_VAL = 'io'
_WARN_MISMATCH_VAL = 'mismatch'
_WARN_OVRWRT_VAL = 'overwrite'
_WARN_UNIQ_VAL = 'uniq'
_WARN_CODES = (_WARN_ID_VAL, _WARN_IO_VAL, _WARN_MISMATCH_VAL, _WARN_OVRWRT_VAL)
_WARN_CODES_PREP = (_WARN_OVRWRT_VAL, _WARN_UNIQ_VAL)
_WARN_PREFIX = '****** WARNING ****** '


##########################
###### Annotate Raw ######
##########################

def _prep_fast5_for_fastq(fast5_data, bc_grp_name, bc_subgrp_name, overwrite):
    read_id = th.get_raw_read_slot(fast5_data).attrs.get('read_id')
    try:
        read_id = read_id.decode()
    except (AttributeError, TypeError):
        pass
    if read_id is None:
        return

    # if Analyses group doesn't exist yet, create it
    try:
        analyses_grp = fast5_data['/Analyses']
    except:
        analyses_grp = fast5_data.create_group('Analyses')

    # create Fastq slot, unless value exists and --overwrite is not set
    try:
        bc_grp = analyses_grp[bc_grp_name]
        bc_subgrp = analyses_grp[bc_subgrp_name]
    except:
        try:
            bc_grp = analyses_grp.create_group(bc_grp_name)
            bc_subgrp = bc_grp.create_group(bc_subgrp_name)
        except:
            if overwrite:
                del analyses_grp[bc_grp_name]
                bc_grp = analyses_grp.create_group(bc_grp_name)
                bc_subgrp = bc_grp.create_group(bc_subgrp_name)
            else:
                raise th.TomboError(
                    bc_grp_name + ' exists and --overwrite is not set.')

    return read_id

def _annotate_with_fastqs_worker(
        fastq_rec_q, fast5s_read_ids, fastq_slot, fq_slot_prepped,
        prog_q, warn_q, bc_grp_name, bc_subgrp_name, overwrite):
    been_warned = dict((warn_code, False) for warn_code in _WARN_CODES)
    num_recs_proc = 0
    while True:
        fastq_rec = fastq_rec_q.get()
        if fastq_rec is None:
            break

        # extract read_id from fastq (which should be the first text after
        # the "@" record delimiter up to the first white space or underscore
        read_id = fastq_rec[0].split()[0].split('_')[0][1:]
        if read_id not in fast5s_read_ids:
            if not been_warned[_WARN_ID_VAL]:
                been_warned[_WARN_ID_VAL] = True
                warn_q.put(_WARN_ID_VAL)
            continue

        try:
            with h5py.File(fast5s_read_ids[read_id], 'r+') as fast5_data:
                if not fq_slot_prepped:
                    try:
                        file_parsed_id = _prep_fast5_for_fastq(
                            fast5_data, bc_grp_name, bc_subgrp_name, overwrite)
                    except th.TomboError:
                        if not been_warned[_WARN_OVRWRT_VAL]:
                            been_warned[_WARN_OVRWRT_VAL] = True
                            warn_q.put(_WARN_OVRWRT_VAL)
                        continue
                    if read_id != file_parsed_id:
                        if not been_warned[_WARN_MISMATCH_VAL]:
                            been_warned[_WARN_MISMATCH_VAL] = True
                            warn_q.put(_WARN_MISMATCH_VAL)
                        continue
                bc_slot = fast5_data[fastq_slot]
                # add sequence to fastq slot
                bc_slot.create_dataset(
                    'Fastq', data=''.join(fastq_rec),
                    dtype=h5py.special_dtype(vlen=unicode))

                # progress q update
                num_recs_proc += 1
                if num_recs_proc % _PROC_UPDATE_INTERVAL == 0:
                    prog_q.put(_PROC_UPDATE_INTERVAL)
        except:
            if not been_warned[_WARN_IO_VAL]:
                been_warned[_WARN_IO_VAL] = True
                warn_q.put(_WARN_IO_VAL)
            continue

    # add last number of records reported from this process
    prog_q.put(num_recs_proc % _PROC_UPDATE_INTERVAL)

    return

def _feed_seq_records_worker(fastq_fns, fastq_rec_q, num_processes):
    for fastq_fn in fastq_fns:
        n_recs = 0
        with io.open(fastq_fn) as fastq_fp:
            while True:
                fastq_rec = list(islice(fastq_fp, 4))
                # if record contains fewer than 4 lines this indicates the
                # EOF, so move to next file
                if len(fastq_rec) != 4: break
                # if sequence identifier line does not start with "@" or quality
                # score line does not start with a "+" the file may be
                # corrupted, so don't process any more records
                if (re.match('@', fastq_rec[0]) is None or
                    re.match('\+', fastq_rec[2]) is None):
                    # TODO maybe send this as a warning code to avoid poorly
                    # formatted output
                    th.warning_message(
                        'Successfully parsed ' + unicode(n_recs) +
                        ' FASTQ records from ' + fastq_fn + ' before ' +
                        'encountering an invalid record. The rest of ' +
                        'this file will not be processed.')
                    break
                n_recs += 1
                fastq_rec_q.put(fastq_rec)

    # put none records to trigger annotation processes to exit
    for _ in range(num_processes):
        fastq_rec_q.put(None)

    return

def _get_ann_queues(prog_q, warn_q, num_read_ids, wp_conn):
    if VERBOSE: bar = tqdm(total=num_read_ids, smoothing=0)
    been_warned = dict((warn_code, False) for warn_code in _WARN_CODES)

    def update_warn(warn_val):
        if warn_val == _WARN_ID_VAL:
            if VERBOSE and not been_warned[_WARN_ID_VAL]:
                bar.write(
                    _WARN_PREFIX + 'Some FASTQ records contain read ' +
                    'identifiers not found in any FAST5 files or ' +
                    'sequencing summary files.',
                    file=sys.stderr)
            been_warned[_WARN_ID_VAL] = True
        elif warn_val == _WARN_IO_VAL:
            if VERBOSE and not been_warned[_WARN_IO_VAL]:
                bar.write(
                    _WARN_PREFIX + 'Some read files that could not be accessed.',
                    file=sys.stderr)
            been_warned[_WARN_IO_VAL] = True
        elif warn_val == _WARN_MISMATCH_VAL:
            if VERBOSE and not been_warned[_WARN_MISMATCH_VAL]:
                bar.write(
                    _WARN_PREFIX + 'Read ID(s) found in sequencing summary ' +
                    'and FAST5 file are discordant. Skipping such reads.',
                    file=sys.stderr)
            been_warned[_WARN_MISMATCH_VAL] = True
        elif warn_val == _WARN_OVRWRT_VAL:
            if VERBOSE and not been_warned[_WARN_OVRWRT_VAL]:
                bar.write(
                    _WARN_PREFIX + 'Basecalls exsit in specified slot for ' +
                    'some reads. Set --overwrite option to overwrite these ' +
                    'basecalls.', file=sys.stderr)
            been_warned[_WARN_OVRWRT_VAL] = True
        else:
            if VERBOSE: bar.write(
                    _WARN_PREFIX + 'Invalid warning code encountered.',
                    file=sys.stderr)

        return


    total_added_seqs = 0
    while True:
        try:
            iter_added = prog_q.get(block=False)
            total_added_seqs += iter_added
            if VERBOSE: bar.update(iter_added)
        except queue.Empty:
            try:
                warn_val = warn_q.get(block=False)
                update_warn(warn_val)
            except queue.Empty:
                sleep(0.1)
                # check if main thread has finished with all fastq records
                if wp_conn.poll():
                    break

    # collect all remaining warn and progress values
    while not prog_q.empty():
        iter_added = prog_q.get(block=False)
        total_added_seqs += iter_added
        if VERBOSE: bar.update(iter_added)
    while not warn_q.empty():
        warn_val = warn_q.get(block=False)
        update_warn(warn_val)

    if VERBOSE:
        bar.close()
        th.status_message('Added sequences to a total of ' +
                          str(total_added_seqs) + ' reads.')
        if total_added_seqs < num_read_ids:
            th.warning_message(
                'Not all read ids from FAST5s or sequencing summary files ' +
                'were found in FASTQs.\n\t\tThis can result from reads that ' +
                'failed basecalling or if full sets of FAST5s/sequence ' +
                'summaries are not processed with full sets of FASTQs.')

    return

def _annotate_with_fastqs(
        fastq_fns, fast5s_read_ids, fastq_slot, fq_slot_prepped, num_processes,
        bc_grp_name, bc_subgrp_name, overwrite):
    if VERBOSE: th.status_message('Annotating FAST5s with sequence from FASTQs.')
    fastq_rec_q = Queue(maxsize=_MAX_FASTQ_QUEUE_SIZE)
    prog_q = Queue()
    warn_q = Queue()

    # open a single process to read fastq files and feed the fastq record queue
    fq_feed_p = Process(target=_feed_seq_records_worker,
                        args=(fastq_fns, fastq_rec_q, num_processes))
    fq_feed_p.daemon = True
    fq_feed_p.start()

    # open fast5 annotation processes
    ann_args = (fastq_rec_q, fast5s_read_ids, fastq_slot, fq_slot_prepped,
                prog_q, warn_q, bc_grp_name, bc_subgrp_name, overwrite)
    ann_ps = []
    for p_id in range(num_processes):
        ann_p = Process(target=_annotate_with_fastqs_worker, args=ann_args)
        ann_p.daemon = True
        ann_p.start()
        ann_ps.append(ann_p)

    main_wp_conn, wp_conn = Pipe()
    warn_prog_p = Process(target=_get_ann_queues,
                          args=(prog_q, warn_q, len(fast5s_read_ids), wp_conn))
    warn_prog_p.daemon = True
    warn_prog_p.start()

    fq_feed_p.join()
    for ann_p in ann_ps:
        ann_p.join()
    # send signal to warn/progress queue that all other processes are complete
    main_wp_conn.send(True)
    warn_prog_p.join()

    return


##########################
#### Extract read_ids ####
##########################

def _get_prep_queue(read_ids_q, prog_q, warn_q, gp_conn, num_fast5s):
    """Process all records from all fast5 prep queues
    """
    ovrwrt_mess = (
        _WARN_PREFIX + 'Basecalls exsit in specified slot for some ' +
        'reads. Set --overwrite option to overwrite these basecalls.')
    fast5s_read_ids = {}
    # Warn non-unique read_ids in directory
    been_warned = dict((warn_code, False) for warn_code in _WARN_CODES_PREP)
    if VERBOSE: bar = tqdm(total=num_fast5s, smoothing=0)

    while True:
        try:
            read_id, fast5_fn = read_ids_q.get(block=False)
            if read_id in fast5s_read_ids:
                if VERBOSE and not been_warned[_WARN_UNIQ_VAL]:
                    bar.write(
                        _WARN_PREFIX + 'Multiple FAST5 files contain the ' +
                        'same read ID. Ensure that FAST5 files are from a ' +
                        'single run.', file=sys.stderr)
                    been_warned[_WARN_UNIQ_VAL] = True
                continue
            fast5s_read_ids[read_id] = fast5_fn
        except queue.Empty:
            try:
                warn_val = warn_q.get(block=False)
                if warn_val == _WARN_OVRWRT_VAL:
                    if VERBOSE and not been_warned[_WARN_OVRWRT_VAL]:
                        bar.write(ovrwrt_mess, file=sys.stderr)
                    been_warned[_WARN_OVRWRT_VAL] = True
                else:
                    bar.write(_WARN_PREFIX + 'Invalid warning code encountered.',
                              file=sys.stderr)
            except queue.Empty:
                try:
                    if VERBOSE: bar.update(prog_q.get(block=False))
                except queue.Empty:
                    sleep(0.1)
                    # check if main thread has finished with all FAST5s
                    if gp_conn.poll():
                        break

    while not read_ids_q.empty():
        read_id, fast5_fn = read_ids_q.get(block=False)
        fast5s_read_ids[read_id] = fast5_fn
    while not warn_q.empty():
        warn_val = warn_q.get(block=False)
        if warn_val == _WARN_OVRWRT_VAL:
            if VERBOSE and not been_warned[_WARN_OVRWRT_VAL]:
                bar.write(ovrwrt_mess, file=sys.stderr)
                been_warned[_WARN_OVRWRT_VAL] = True
            else:
                bar.write(_WARN_PREFIX + 'Invalid warning code encountered.',
                          file=sys.stderr)
    while not prog_q.empty():
        if VERBOSE: bar.update(prog_q.get(block=False))

    if VERBOSE: bar.close()
    gp_conn.send(fast5s_read_ids)

    return

def _prep_fastq_slot_worker(
        fast5_q, bc_grp, bc_subgrp, overwrite, read_ids_q, prog_q, warn_q):
    num_files_proc = 0
    been_warned_overwrite = False
    while True:
        try:
            fast5_fn = fast5_q.get(block=False)
        except queue.Empty:
            sleep(0.1)
            continue

        if fast5_fn is None:
            break

        num_files_proc += 1
        if num_files_proc % _PROC_UPDATE_INTERVAL == 0:
            prog_q.put(_PROC_UPDATE_INTERVAL)

        try:
            with h5py.File(fast5_fn) as fast5_data:
                try:
                    read_id = _prep_fast5_for_fastq(
                        fast5_data, bc_grp, bc_subgrp, overwrite)
                except th.TomboError:
                    # avoid the warn queue getting too large by sending overwite
                    # warnings for each read from each thread
                    if not been_warned_overwrite:
                        been_warned_overwrite = True
                        warn_q.put(_WARN_OVRWRT_VAL)
                    continue
        except:
            continue
        if read_id is None:
            continue

        read_ids_q.put((read_id, fast5_fn))

    prog_q.put(num_files_proc % _PROC_UPDATE_INTERVAL)

    return

def _fill_files_queue(fast5_q, fast5_fns, num_ps):
    for fast5_fn in fast5_fns:
        fast5_q.put(fast5_fn)
    for _ in range(num_ps):
        fast5_q.put(None)

    return

def _get_read_ids_and_prep_fastq_slot(
        fast5s_dir, bc_grp, bc_subgrp, overwrite, num_processes):
    """Extract read id from /Raw group and prep fastq slots for annotation with
    associated FASTQ files.
    """
    if VERBOSE: th.status_message(
            'Preparing reads and extracting read identifiers.')
    fast5_q = Queue(maxsize=_MAX_QUEUE_SIZE)
    read_ids_q = Queue()
    prog_q = Queue()
    warn_q = Queue()

    fast5_fns = th.get_files_list(fast5s_dir)
    files_p = Process(target=_fill_files_queue,
                      args=(fast5_q, fast5_fns, num_processes))
    files_p.daemon = True
    files_p.start()

    prep_args = (fast5_q, bc_grp, bc_subgrp, overwrite, read_ids_q,
                 prog_q, warn_q)
    prep_ps = []
    for p_id in range(num_processes):
        prep_p = Process(target=_prep_fastq_slot_worker, args=prep_args)
        prep_p.daemon = True
        prep_p.start()
        prep_ps.append(prep_p)

    main_gp_conn, gp_conn = Pipe()
    get_prep_p = Process(
        target=_get_prep_queue,
        args=(read_ids_q, prog_q, warn_q, gp_conn, len(fast5_fns)))
    get_prep_p.daemon = True
    get_prep_p.start()

    # join all processes into the main thread
    files_p.join()
    for prep_p in prep_ps:
        prep_p.join()
    # send signal to get_prep queue that all other processes are complete
    main_gp_conn.send(True)
    fast5s_read_ids = main_gp_conn.recv()

    return fast5s_read_ids

def _parse_sequencing_summary_files(fast5s_dir, seq_summary_fns):
    if VERBOSE: th.status_message('Getting read filenames.')
    full_fast5_fns = {}
    # walk through directory structure searching for fast5 files
    for root, _, fns in os.walk(fast5s_dir):
        for fn in fns:
            if not fn.endswith('.fast5'): continue
            full_fast5_fns[fn] = os.path.join(root, fn)

    if VERBOSE: th.status_message('Parsing sequencing summary files.')
    fast5s_read_ids = {}
    been_warned = False
    for seq_summary_fn in seq_summary_fns:
        with open(seq_summary_fn) as fp:
            try:
                header_fields = fp.readline().split()
                fn_field = next(i for i, h_field in enumerate(header_fields)
                                if re.match(_SEQ_SUMMARY_FN_FIELD, h_field))
                id_field = next(i for i, h_field in enumerate(header_fields)
                                if re.match(_SEQ_SUMMARY_ID_FIELD, h_field))
            except:
                th.warning_message(
                    'Could not extract header information for sequencing ' +
                    'summary file: ' + seq_summary_fn)
                continue
            try:
                for line in fp:
                    rec_fields = line.split()
                    rec_short_fn = rec_fields[fn_field]
                    try:
                        rec_full_fn = full_fast5_fns[rec_short_fn]
                    except KeyError:
                        if not been_warned:
                            th.warning_message(
                                'Some FASTQ records from sequencing summaries ' +
                                'do not appear to have a matching file.')
                        been_warned = True
                        continue
                    # convert filename to full filename and link to read id
                    fast5s_read_ids[rec_fields[id_field]] = rec_full_fn
            except:
                th.warning_message(
                    'Error parsing records for sequencing ' +
                    'summary file: ' + seq_summary_fn)

    return fast5s_read_ids


##################################
###### Annotate FAST5s Main ######
##################################

def annotate_reads_with_fastq_main(args):
    global VERBOSE
    VERBOSE = not args.quiet
    th.VERBOSE = VERBOSE

    fast5s_basedir = (
        args.fast5_basedir if args.fast5_basedir.endswith('/') else
        args.fast5_basedir + '/')
    if args.sequencing_summary_filenames:
        fast5s_read_ids = _parse_sequencing_summary_files(
            fast5s_basedir, args.sequencing_summary_filenames)
        fq_slot_prepped = False
    else:
        fast5s_read_ids = _get_read_ids_and_prep_fastq_slot(
            fast5s_basedir, args.basecall_group, args.basecall_subgroup,
            args.overwrite, args.processes)
        fq_slot_prepped = True
    fastq_slot = '/'.join(('/Analyses', args.basecall_group,
                           args.basecall_subgroup))
    _annotate_with_fastqs(
        args.fastq_filenames, fast5s_read_ids, fastq_slot, fq_slot_prepped,
        args.processes, args.basecall_group, args.basecall_subgroup,
        args.overwrite)

    return


if __name__ == '__main__':
    sys.stderr.write('This is a module. See commands with `tombo -h`')
    sys.exit(1)
