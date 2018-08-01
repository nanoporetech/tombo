***********************
Read Filtering Commands
***********************

Read filtering commands can be useful to extract the most out out of a set of reads for modified base detection. Read filtering commands effect only the Tombo index file, and so filters can be cleared or applied iteratively without re-running the re-squiggle command. Five filters are currently made available (``genome_locations``, ``raw_signal_matching``, ``q_score``,  ``level_coverage`` and ``stuck``).

---------------------------------
``tombo filter genome_locations``
---------------------------------

The ``tombo filter genome_locations`` command filters out reads falling outside of a specified set of ``--include-regions``. These regions can either be whole chromosomes/sequence records or sub-regions within sequence records.

------------------------------------
``tombo filter raw_signal_matching``
------------------------------------

The ``tombo filter raw_signal_matching`` command filters out reads with poor matching between raw observed signal and expected signal levels from the canonical base model. Specify a new threshold to apply with the ``--signal-matching-score`` option. These scores are the mean half z-score (absolute value of z-score) taken over all bases of a read. A reasonable range for this threshold should be approxiamtely between 0.5 and 3. Reads with a larger fraction of modifications may require a larger value to process successfully.

------------------------
``tombo filter q_score``
------------------------

The ``tombo filter q_score`` command filters out reads with poor mean basecalling quality scores. This value can be indicative of low quality reads. Set this value with the ``--q-score`` option.

-------------------------------
``tombo filter level_coverage``
-------------------------------

The ``tombo filter level_coverage`` command aims to filter reads to achieve more even read depth across a genome/transcriptome. This may be useful in canonical and alternative model estimation. This filter may also help make test statistics more comparable across the genome.

This filter is applied by randomly selecting reads weighted by the approximate coverage at the mapped location of each read. The number of reads removed from downstream processing is defined by the ``--percent-to-filter`` option.

This filter is likely to be more useful for PCR'ed sample where duplicate locations are more likely to accumulate and cause large spikes in coverage.

----------------------
``tombo filter stuck``
----------------------

The ``tombo filter stuck`` command aims to remove reads where bases tend to get stuck in the pore for longer durations of time. These reads can be indicative of poor quality reads and thus negatively effect modified base detection.

This filter is based on the number of observations per genomic base along a read. The filter can be set on any number of percentiles of obervations per base. Reasonable values depend strongly on the sample type (DNA or RNA). A reasonable filter for DNA reads would be to filter reads with 99th percentile > 200 obs/base or a maximum base with > 5k obs/base. This filter would be set with the ``--obs-per-base-filter 99:200 100:5000`` option. Larger values should be used for RNA reads.

------------------------------
``tombo filter clear_filters``
------------------------------

The ``tombo filters clear_filters`` command removes any applied filters to this sample (including those applied during the ``resquiggle`` command; though reads that failed before signal to sequence assignment will not be included). New filters can then be applied to this set of reads.

All Tombo sub-commands will respect the filtered reads when parsed for processing.

.. hint::

   Save a set of filters for later use by copying the Tombo index file: ``cp path/to/native/rna/.fast5s.RawGenomeCorrected_000.tombo.index save.native.tombo.index``. To re-set to a set of saved filters after applying further filters simply replace the index file: ``cp save.native.tombo.index path/to/native/rna/.fast5s.RawGenomeCorrected_000.tombo.index``.
