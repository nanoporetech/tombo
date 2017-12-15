***********************
Read Filtering Commands
***********************

Read filtering commands can be useful to extract the most out out of a set of reads for modified base detection. Read filtering commands effect only the Tombo index file, and so filters can be cleared or applied iteratively without re-running any re-squiggle commands. Two filters are currently made available (``filter_stuck`` and ``filter_coverage``).

----------------
``filter_stuck``
----------------

The ``filter_stuck`` command aims to remove reads where bases tend to apparently get stuck in the pore for longer durations of time. These reads can be indicative of poor quality reads and thus negatively effect modified base detection.

This filter is based on the number of observations per genomic base along a read. The filter can be set on any number of percentiles of obervations per base. Reasonable values depend strongly on the sample type (DNA or RNA). A reasonable filter for DNA reads would be to filter reads with 99th percentile > 200 obs/base or a maximum base with > 5k obs/base. This filter would be set with the ``--obs-per-base-filter 99:200 100:5000`` option. Larger values should be used for RNA reads.

-------------------
``filter_coverage``
-------------------

The ``filter_coverage`` command aims to filter reads to achieve more even read depth across a genome. This may be useful particularly in canonical and particularly in alternative mode estimation. This filter may also help make some testing cases more comparable across the genome.

This filter is applied by randomly selecting reads weighted by the approximate coverage at the mapped location of each read. The number of reads removed from downstream processing is defined by the ``--percent-to-filter`` option.

This filter is likely to be more useful for PCR'ed sample where duplicate locations are more likely to accumulate and cause very large spikes in coverage.

-----------------
``clear_filters``
-----------------

The ``clear_filters`` simply removes any applied filters to this sample (failed reads from the re-squiggle command will still not be included). New filters can then be applied to this set of reads.

All Tombo commands will respect the filtered reads when they are parsed for procesing
