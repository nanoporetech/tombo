************
Text Outputs
************

Two main text outputs are available from Tombo:

1. Wiggle - genome browser statistics
2. Fasta - Sequence output surrounding most modified sites

``write_wiggles``
-----------------

The ``write_wiggles`` command takes in a set of reads (``--fast5-basedirs``) and potentially a control set of reads (``--control-fast5-basedirs``) or a pre-computed statistics file (``--statistics-filename``). Output wiggle files (`variableStep format <https://genome.ucsc.edu/goldenpath/help/wiggle.html>`_) will be produced for each requested statistic (both plus and minus strands).

Several statistics are available for output:

* ``coverage`` - The coverage level for mapped and validly re-squiggled reads
* ``fraction`` - The fraction of significantly modified reads
* ``signal`` - The mean signal level across all reads mapped to this location
* ``signal_sd`` - The mean signal standard deviation across all reads mapped to this location (not available unless ``--include-event-stdev`` was provided in resquiggle call)
* ``difference`` - The difference in normalized signal level between a sample and control set of reads

..

    Note that ``signal``, ``signal_sd`` and ``difference`` require each reads event level data to be queried and thus may be quite slow. ``coverage`` and ``fraction`` can be extracted simply from the tombo statistics file, which is much faster.

Potentially deprecated options:

* ``stat`` - Per base statistical testing -log10(p-values). Test values are quite dependent on read depth and thus this option may be deprecated at some point
* ``mt_stat`` - Multiple testing corrected statistical test values

Files will be output to individual wiggle files (two per statistic for plus and minus genomic strand) in the following format ``[wiggle-basename].[wiggle-type].[sample|control]?.[plus|minus].wig``

``write_most_significant_fasta``
--------------------------------

The ``write_most_significant_fasta`` command writes the genome sequence surrounding the most modified positions. This can be useful for several tasks related to modified base detection including motif discovery.

To run ``write_most_significant_fasta``, a ``--statistics-filename`` is required to extract the most significant locations and either a ``--fast5-basedirs`` or ``--genome-fasta`` is required to extract the genomic sequence. Several options are availble for selecting the sequence to be output:

* ``--num-regions`` - Defines the number of unique locations to be output
* ``--num-bases`` - Defines the number of bases to be output surrounding the significant locations

Potentially deprecated options:

* ``--statistic-order`` - Order regions by per-genomic base statistical testing values instead of fraction of reads with significant modified base results
* ``--q-value-threshold`` - Select the number of regions to output based on a q-value threhsold instead of a set number (This may produce very large files if not set carefully and so this option may be deprecated)

..

    Note that fraction of reads with a significant result at this location can produce non-optimal results with the alternative base comparison log likelihood ratio test. This may be replaced by an estimated fraction based on testing results instead of the current thresholding criterion.
