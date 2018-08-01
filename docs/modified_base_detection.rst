***********************
Modified Base Detection
***********************

Tombo enables three methods for detecting shifts in current signal level, indicative of non-canonical bases. These three methods allow researchers to investigate non-canonical bases given any sample type, while enabling more accurate detection of specific modifications when applicable.

----

.. figure::  _images/testing_method_comparison.png
   :align: center

   Tombo modified base testing methods.

----

All three methods are accessed by the ``tombo detect_modifications`` command group as described below.

**TL;DR**:

* To identify 5-methylcytosine (5mC; DNA or RNA) and N6-methyladenosine (6mA; DNA only), run ``tombo detect_modifications alternative_model`` with the ``--alternate-bases 5mC 6mA`` option
* For more experimental *de novo* modified base detection simply run ``tombo detect_modifications de_novo`` with just a set of reads
* For modified base detection via comparison to a control sample (e.g. PCR or IVT) run ``tombo detect_modifications sample_compare`` with a control set of reads (``--control-fast5-basedirs``)
* The ``tombo detect_modifications`` command will produce a binary file (not intended for use outside the Tombo framework)

  - To extract useful text files see the ``tombo text_output`` commands
  - To visualize raw signal around significant regions use the ``tombo plot most_significant`` command
  - To assess testing results around a known motif use the ``tombo plot motif_with_stats``, ``tombo plot roc``, and  ``tombo plot per_read_roc`` commands

.. hint::

   The ``tombo resquiggle`` command must be run on a set of reads before processing with ``tombo detect_modifications``.

Specific Alternate Base Detection (Recommended)
===============================================

In order to specifically detect 5mC and/or 6mA, use the ``tombo detect_modifications alternative_model`` command. This command computes a statistic similar to a log likelihood ratio (LLR) but dynamically scaled to be more robust to outlier signal levels. This statistic is computed for each "swap base" within each read provided (e.g. each cytosine for 5mC detection or each adenine for 6mA detection).

This statistic is computed by scaling the LLR by the normal likelihood function with the same variance and mean halfway between the canonical and alternative expected signal levels. Three additional scaling factors are added to this function in order to give greater weight to sequence contexts with larger differences between the canonical and alternative expected signal levels, which inherently provide more power to distinguish the canonical and alternative base signal levels. These parameters are also set so that values are on relatively the same scale as a log likelihood ratio for setting ``--single-read-threshold`` values. Default values for the scale factors below are :math:`S_f = 4`, :math:`S_{f2} = 3` and :math:`S_p = 0.3`, which produce the functions shown in the figure below. Users can experiment with the effect of these parameters with the provided ``scripts/outlier_robust_llr.R`` script.

.. math::

   \begin{align}
   ScaleDiff& = NormSignal - \frac{CanonicalMean + AltMean}{2}\\
   MeanDiffs& = |CanonicalMean - AltMean|\\
   OutlierRobustLlr& = \frac{e^{\frac{ScaleDiff^2}{S_f \cdot \sigma^2}} \cdot LLR}{\sigma^2 \cdot {MeanDiffs}^{S_p} \cdot S_{f2}}
   \end{align}

In order to compute a standard log likelihood ratio, use the ``--standard-log-likelihood-ratio`` option.

----

.. figure::  _images/outlier_robust_llr.gif
   :align: center

   Tombo outlier-robust versus standard likelihood ratio statistic over varied differences between canonical and alternative expected signal levels.

----

This statistic is computed and summed over all positions where the base of interest is included in the modeled k-mer. The default DNA model is a 6-mer, so the signal at the six surrounding genomic bases contribute to the resulting statistic at any one position. For example, for 5mC detection within in a TGGTA **C** GTCCG context, the signal will be tested against expected canonical and alternative 5mC levels at the following locations::

  TGGTA **C** GTCCG
  -----------------
  TGGTA **C**
   GGTA **C** G
    GTA **C** GT
     TA **C** GTC
      A **C** GTCC
        **C** GTCCG

New alternative base models will be added as they are trained and validated internally. This is the perferred method for modified base detection if a model is available for your biological sample of interest as the exact modification position is identified.

.. code-block:: bash

    tombo detect_modifications alternative_model --fast5-basedirs <fast5s-base-directory> \
        --alternate-bases 5mC 6mA --statistics-file-basename sample.alt_model

.. hint::

    Users may also train their own alternative base Tombo models and test against these with the advanced ``--alternate-model-filenames`` option. See more details in the :doc:`model_training` section.

De novo Non-canonical Base Method
=================================

In order to perform *de novo* non-canonical base detection, use the ``tombo detect_modifications de_novo`` command. This method is ideal for unknown modification motif detection when using in combination with the ``tombo text_output signif_sequence_context`` command and motif detection software (e.g. `MEME <http://meme-suite.org/tools/meme>`_).

For each read at each position, this method performs a hypothesis test against the canonical model based on the genomic sequence. Note that this method can be quite error prone and may result in a high false positive rate, especially on a per-read basis. This method also has the lowest barrier to entry, requiring only a set of reads and a reference sequence, allowing any nanopore researcher to start investigating potentially any type of modified base.

.. code-block:: bash

    tombo detect_modifications de_novo --fast5-basedirs <fast5s-base-directory> \
        --statistics-file-basename sample.de_novo

Canonical Sample Comparison Method
==================================

In order to perform *canonical sample comparison* modified base detection, use the ``tombo detect_modifications sample_compare`` command with a second set of reads from the same biological sample containing only canonical bases (e.g. PCR for DNA or IVT for RNA) via the ``--control-fast5-basedirs`` option.

For each sample read, this will perform a hypothesis test against a distribution estimated from the signal levels observed from the control sample reads at each genome position. As of version 1.4, this method uses the canonical base model as a prior for this estimated distribution improving results for low coverage regions (disable canonical prior with the ``--sample-only-estimates`` option or lower the prior impact on estimates by lowering the default ``--model-prior-weights`` values).

.. code-block:: bash

    tombo detect_modifications sample_compare --fast5-basedirs <fast5s-base-directory> \
        --control-fast5-basedirs <control-fast5s-base-directory> \
        --statistics-file-basename sample.compare_sample

.. note::

   Due to the nature of nanopore sequencing, the genomic context surrounding the read head effect that current at any position. Thus shifts in signal due to a modified base may occur at several positions to either side of the true modified location. In order to account for this the canonical sample and *de novo* modfied base detection methods accept the ``--fishers-method-context`` option which combines test values, using `Fisher's Method <https://en.wikipedia.org/wiki/Fisher%27s_method>`_, over a moving window across the genome. This  can help to center significant values on modified base positions. The default value for this parameter is 1, but reasonable results can be obtained for values between 0 and 3.

Aggregating Per-read Statistics
===============================

All of the above methods compute per-read, per-genome location test statistics. In order to facilitate research at the genomic location level, these per-read statistics are combined at each genomic location by applying a global threshold identifying each read as supporting a canonical or alternative base. This results in a fraction of reads indicating a modified base at each genomic location. This global threshold may consist of a single threshold value or a pair of values (where test statistics between the values do not contribute to the estimated fraction of modified reads).

All ``tombo detect_modifications`` methods enable output of per-read test statistics (``--per-read-statistics-basename``). Tombo also provides the ``tombo detect_modifications aggregate_per_read_stats`` command in order to apply different global threshold values to per-read statistics without re-computing these statistics. Note it is not possible to change other testing parameters from this command (e.g. ``--fishers-method-context``).

Dampened Fraction Estimates
===========================

At low coverage locations the fraction of modified reads estimates can be poor. Thus the ``--coverage-dampen-counts`` option is provided in order to dampen the estimated fraction of modified reads at low coverage locations. This allows easier use of the fraction statistic in downstream analysis.

  - The fraction estimate includes pseudo-counts added to the un-modified and modified read counts (as specified by the ``--coverage-dampen-counts`` option)
  - This is equivalent to using a beta prior when estimating the fraction of reads modified at each position
  - Test the effect of different dampen counts using the ``scripts/test_beta_priors.R`` (the default values are shown below)
  - The raw fraction is still included in the statistics file (access from python API)

----

.. figure::  _images/dampened_fraction.png
   :align: center

   Heatmap showing the resulting dampened farction of modified reads given the default ``--coverage-dampen-counts`` values over range of coverage and number of un-modified reads.

----

Multi-processing
================

Tombo statistical testing provides the option to perform testing spread across multiple processes. This also limits the memory requirement for modified base detection, as only signal levels within a multiprocess block are held in memory. For very high coverage samples, consider lowering the ``--multiprocess-region-size`` value to minimize computational memory usage.

Multi-processing is performed over batches delineated by regular intervals across chromosomes covered by at least one read. The interval size is determined by the ``--multiprocess-region-size`` option and processed by a number of processors indicated by the ``--processes`` option. The produced per-base (and per-reda) results are identical no matter the multi-processing options selected. These regions are also used as batches to store the pre-read statistics file.

Tombo Statistics File Format
============================

For all modified base detection methods, the result is a binary Tombo statistics file. This file contains statistics associated with each genomic base producing a valid result. This file is not intended for use outside of the Tombo framework. Several Tombo commands (e.g. ``tombo text_output browser_files``, ``tombo text_output signif_sequence_context`` and ``tombo plot most_significant``) take the binary statistics file as an input, accommodating many user pipelines downstream of modified base detection.

While the Tombo statistics file is meant to be a binary file not processed by outside tools its contents are described here for completeness. Access to this file is recommended through the ``tombo.tombo_helper.TomboStats`` object in the Tombo python API.

.. important::

   All other optional arguments to the ``tombo.tombo_stats.TomboStats`` constructor should be left as ``None``; setting these values will delete the file and construct a blank per-read statistics file.

The Tombo statistics file is in `HDF5 format <https://support.hdfgroup.org/HDF5/whatishdf5.html>`_. Attributes at the root level are 1) ``stat_type`` indicating which testing method was used (``model_compare``, ``de_novo`` or ``sample_compare``), 2) ``block_size`` indicating the number of genomic bases in each statistics block and 3) `Cov_Threshold`` containing the coverage threshold applied to this file.

Blocks of statistics are stored in the ``Statistic_Blocks`` group. Within this group, each block of statistics is found within a group named ``Group_NNN``. Each group contains attributes for the block ``start``, ``chrm`` and ``strand``. The ``block_stats`` data set contains the per-location statistics records. Each record contains the following attributes: ``frac``, ``pos``, ``chrm``, ``strand``, ``cov``, ``control_cov``, and ``valid_cov``.

``frac`` contains the fraction of valid (not including per-read statistics within the interval specified by ``--single_read_threshold``) reads at this genomic position identified as the standard base.

``cov``, ``control_cov``, and ``valid_cov`` contain the read coverage at the genomic position for the sample and control reads. ``control_cov`` is only applicable for the control sample comparison testing method. ``valid_cov`` contains the number of reads contributing to the ``frac`` of tested reads as defined by ``--single-read-threshold``.

Per-read Statistics File Format
===============================

Per-read statistics can be stored by setting the ``--per-read-statistics-basename`` option to any ``tombo detect_modifications`` command. This output file can then be used in downstream Tombo sub-commands (e.g. the ``tombo plot per_read`` and ``tombo detect_modifications aggregate_per_read_stats`` commands).

For advanced users, the Tombo per-read statsitics file can be accessed via the Tombo python API using the ``tombo.tombo_stats.PerReadStats`` class. This class provides initialization, simply taking the per-read statsitics filename. The ``PerReadStats`` class supports the ``get_region_stats`` function which takes a ``tombo.tombo_helper.intervalData`` object specifying an interval of interest. This will return a numpy array containing a record for each read (specified by the ``read_id`` field) and each tested genomic position (``pos`` field) along with the test statistic (``stat`` field) at that location.

.. important::

   All other optional arguments to the ``tombo.tombo_stats.PerReadStats`` constructor should be left as ``None``; setting these values will delete the file and construct a blank per-read statistics file.

The per-read statistics file is in the HDF5 format. All blocks are stored within the ``Statistic_Blocks`` slot. The size of the blocks is stored in the ``block_size`` attribute (defined by the ``--multiprocess-region-size`` option) and the type of statistical test applied is stored in the ``stat_type`` attribute.

Each genomic block is stored in a different ``Block_[##]`` slot. These slots do not have any particular order. Within each block the ``chrm``, ``strand`` and ``start`` of the block are stored. The block statistics are stored in the ``block_stats`` slot. Per-read statistics contain a record for each tested location within each read. Each record contains the genomic position (``pos``), the test statistic (``stat``; hypothesis test p-value or log likelihood ratio as indicated by the statistic type), and the ``read_id``. A single read spanning multiple blocks will contain statistics in more than one block. An individual read's statistics can be reconstructed using the ``read_id`` field.
