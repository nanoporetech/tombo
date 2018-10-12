"""
****************
Tombo Python API
****************

The Tombo python API is intended to give users access to some of the lower level processing functions used throughout Tombo analysis pipelines.

Primarily this includes inspection of per-read modified base detection results and observed (re-squiggle assigned) per-read signal levels.

.. note::

    Effort will be made to maintain this API interface introduced at Tombo version 1.4, but major structural changes to the Tombo framework may require changes to some API interface components. Such changes will be noted in github release notes where applicable.

-------------------
Python API Examples
-------------------

Import Tombo sub-modules::

    from tombo import tombo_helper, tombo_stats, resquiggle

Extract normalized raw signal assigned to each genomic base from each read in a region::

    # specify region of interest
    reg_data = tombo_helper.intervalData(
        chrm='chr20', start=10000, end=10100, strand='+')

    # parse Tombo index from previously re-squiggled set of reads
    reads_index = tombo_helper.TomboReads(['test_data/native_reads',])
    # extract reads that overlap this interval and then extract base signal
    # levels from 10 randomly selected reads
    reg_base_levels = reg_data.add_reads(
        reads_index).get_base_levels(num_reads=10)

Extracting per-read testing results::

    sample_per_read_stats = tombo_stats.PerReadStats(
        'test_stats.alt_model.5mC.tombo.per_read_stats')
    # reg_per_read_stats contains a numpy array containing per-read stats
    # over all reads covering the region of interest
    reg_per_read_stats = sample_per_read_stats.get_region_per_read_stats(
        reg_data)

Run standard resquiggle algorithm (may cause errors on some reads)::

    import h5py, mappy

    # set read values
    fast5_fn, reference_fn = 'path/to/read.fast5', 'genome.fasta'
    fast5_data = h5py.File(fast5_fn, 'r')
    seq_samp_type = tombo_helper.get_seq_sample_type(fast5_data)

    # prep aligner, signal model and parameters
    aligner = mappy.Aligner(reference_fn, preset=str('map-ont'), best_n=1)
    std_ref = tombo_stats.TomboModel(seq_samp_type=seq_samp_type)
    rsqgl_params = tombo_stats.load_resquiggle_parameters(seq_samp_type)

    # extract data from FAST5
    map_results = resquiggle.map_read(fast5_data, aligner, std_ref)
    all_raw_signal = tombo_helper.get_raw_read_slot(fast5_data)['Signal'][:]
    if seq_samp_type.rev_sig:
        all_raw_signal = all_raw_signal[::-1]
    map_results = map_results._replace(raw_signal=all_raw_signal)

    # run full re-squiggle
    rsqgl_results = resquiggle.resquiggle_read(
        map_results, std_ref, rsqgl_params, all_raw_signal=all_raw_signal)

    # or run individual steps
    num_events = tombo_stats.compute_num_events(
        all_raw_signal.shape[0], len(map_results.genome_seq),
        rsqgl_params.mean_obs_per_event)
    valid_cpts, norm_signal, scale_values = resquiggle.segment_signal(
        map_results, num_events, rsqgl_params)
    event_means = tombo_stats.compute_base_means(norm_signal, valid_cpts)
    dp_results = resquiggle.find_adaptive_base_assignment(
        valid_cpts, event_means, rsqgl_params, std_ref, map_results.genome_seq)
    norm_signal = norm_signal[
        dp_results.read_start_rel_to_raw:
        dp_results.read_start_rel_to_raw + dp_results.segs[-1]]
    segs = resquiggle.resolve_skipped_bases_with_raw(
        dp_results, norm_signal, rsqgl_params)
"""

from __future__ import unicode_literals, absolute_import
