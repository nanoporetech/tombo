***********************
Modified Base Detection
***********************

Tombo enables three methods for detecting shifts in current signal level, indicative of non-canonical bases. These three methods allow researchers to investigate non-canonical bases given any sample type, while also enabling more accurate detection of specific modifications when applicable.

----

.. figure::  _images/testing_method_comparison.png
   :align: center
   :scale: 30%
   
   Tombo modified base testing methods.

----

All three methods are accessed by the ``detect_modifications`` Tombo command as described below.

**TL;DR**:

* To identify 5-methylcytosine (5mC) and N6-methyladenosine (6mA), run ``detect_modifications alternative_model`` with the ``--alternate-bases 5mC 6mA`` option
* For more experimental de novo modified base detection simply run ``detect_modifications de_novo`` with just a set of reads
* For modified base detection via comparison to a control sample (e.g. PCR) run ``detect_modifications sample_compare`` with a control set of reads (``--control-fast5-basedirs``)
* The ``detect_modifications`` command will produce a binary file (not intended for use outside the Tombo framework)
  
  - To extract useful text files see the ``text_output`` commands
  - To visualize raw signal around significant regions use the ``plot most_significant`` command
  - To assess testing results around a known motif use the ``plot motif_with_stats``, ``plot roc``, and  ``plot per_read_roc`` commands

.. hint::
   
   The ``resquiggle`` command must be run on a set of reads before processing with ``detect_modifications``.

-------------------
Statistical Testing
-------------------

For all statistical testing methods, the result is a binary Tombo statistics file. This file contains statistics associated with each genomic base producing a valid result. This file is not intended for use outside of the Tombo framework. Several Tombo commands (e.g. ``text_output browser_files``, ``text_output signif_sequence_context`` and ``plot most_significant``) take the statistics file as an input, accommodating many user pipelines downstream of modified base detection.

Of particular interest, the statistics file contains the fraction of reads at each genomic position passing a set threshold or falling outside of a set interval if 2 values are provided to the ``--single-read-threshold`` option. The default value for this parameter is set for each testing method and for DNA and RNA data types using the default settings. Note that changing testing parameters may require a new ``--single-read-threshold`` for optimal results. For example, changing the ``--fishers-method-context`` option value in either the ``de_novo`` or ``compare_sample`` methods is likely to require a new threshold value.

For ``--single-read-threshold`` values with an interval or for the ``alternative_model`` with values greater than 0, the number of reads falling outside of the threshold values is saved under the ``valid_cov`` column in the statistics file. These values can be output with the ``text_output browser_files --file-types valid_coverage`` command.

For the de novo and alternative model testing approaches a default canonical model is used (included with Tombo). Users may also train their own canonical Tombo model (possibly for an older chemistry version) and test against this model using the advanced ``--tombo-model-filename`` option. See more in the :doc:`model_training` section.

Another available output from the ``detect_modifications`` command is a per-read (and per-base) binary (HDF5) statistics file (via ``--per-read-statistics-basename`` option). This file is currently made available for research on per-read modified base detection including plotting via the ``plot per_read`` command and further computing via the ``detect_modifications aggregate_per_read_stats`` command. For advanced researchers, the per-read statistics data can be accessed (including random access to particular regions of the genome) using the ``tombo.tombo_stats.PerReadStats`` class from the Tombo python API.

Alternative Model Method
========================

In order to specifically detect 5mC and 6mA, use the ``detect_modifications alternative_model`` command. Users may also train their own alternative base Tombo models and test against these with the advanced ``--alternate-model-filenames`` option. See more details in the :doc:`model_training` section.

The ``detect_modifications alternative_model`` command will compute a statistic similar to a log likelihood ratio (LLR) but dynamically scaled to be more robust to outlier signal assignment. This statistic is computed for each "swap base" within each read provided (e.g. computed at each cytosine for 5mC detection and each adenine for 6mA detection).

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
   :scale: 30%
   
   Tombo outlier-robust versus standard likelihood ratio statistic over varied differences between canonical and alternative expected signal levels.

----

This statistic is computed and summed over all positions modeled. The default DNA model is a 6-mer, so the signal at the six surrounding genomic bases contribute to the resulting statistic at any one position. For example, for 5mC detection within in a TGGTA **C** GTCCG context, the signal will be tested against expected canonical and alternative 5mC levels at the following locations::

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

    # with user trained alternative base model
    tombo detect_modifications alternative_model --fast5-basedirs <fast5s-base-directory> \
        --alternate-model-filenames alternative_base.tombo.model \
        --statistics-file-basename sample.user_alt_model

De novo Non-canonical Base Method
=================================

In order to perform *de novo* non-canonical base detection, use the ``detect_modifications de_novo`` command.

For each read, this will perform a hypothesis test against the canonical model based on the genomic sequence at each position. Note that this method can be quite error prone and may result in a high false positive rate, but may be of use in a research and development setting. This method also has the lowest barrier to entry, requiring only a set of reads and a genome, allowing any nanopore researcher to start investigating potentially any type of modified base.

.. code-block:: bash

    tombo detect_modifications de_novo --fast5-basedirs <fast5s-base-directory> \
        --statistics-file-basename sample.de_novo

Canonical Sample Comparison Method
==================================

In order to perform *canonical sample comparison* modified base detection, use the ``detect_modifications sample_compare`` command with a second set of reads from the same biological sample containing only canonical bases (e.g. PCR for DNA or IVT for RNA) via the ``--control-fast5-basedirs``.

For each sample read, this will perform a hypothesis test against a normal distribution estimated from the signal level observed from the control sample reads at each genome position. This method does not always identify the exact modification position or the identity of the modified base as with the *de novo* method.

Note that no model is used in the application of this method. Instead the testing null distribution is estimated at each genomic location from the control set of reads.

For both this method, as well as the *de novo* method, the ``--fishers-method-context`` option will combine test values, using `Fisher's Method <https://en.wikipedia.org/wiki/Fisher%27s_method>`_, over a moving window extending a number of positions in either direction. Due to the nature of nanopore sequencing, the genomic context surrounding the read head effect that current at any position. Thus shifts in signal due to a modified base may occur at several positions to either side of the true modified location. Thus combining statistical test values across several genomic positions can help to center significant values on the truly modified position. The default value for this parameter is 1, but reasonable results can be obtained for values between 0 and 3.

.. code-block:: bash

    tombo detect_modifications sample_compare --fast5-basedirs <fast5s-base-directory> \
        --control-fast5-basedirs <control-fast5s-base-directory> \
        --statistics-file-basename sample.compare_sample

-----------------------------
Aggregate Per-read Statistics
-----------------------------

In order to facilitate research on the per-genomic base aggregation across reads, Tombo provides the ``detect_modifications aggregate_per_read_stats`` command. The primary utility for this command is to enable easier manipulation of the per-read threshold values. It is not possible to change other testing parameters from this command (e.g. ``--fishers-method-context`` or ``--tombo-model-filename``).

----------------
Multi-processing
----------------

Tombo statistical testing provides the option to perform testing spread across multiple processes. This also limits the memory requirement for modified base detection, as only signal levels within a multiprocess block are held in memory. For very high coverage samples, consider lowering the ``--multiprocess-region-size`` value to minimize computational memory usage.

Multi-processing is performed over batches delineated by regular intervals across chromosomes covered by at least one read. The interval size is determined by the ``--multiprocess-region-size`` option and processed by a number of processors indicated by the ``--processes`` option. The produced per-base (and per-reda) results are identical no matter the multi-processing options selected. These regions are also used as batches to store the pre-read statistics file.

----------------------------
Tombo Statistics File Format
----------------------------

While the Tombo statistics file is meant to be a binary file not processed by outside tools its contents are described here for completeness. The Tombo statistics file is `HDF5 format <https://support.hdfgroup.org/HDF5/whatishdf5.html>`_. There is one attribute at the root level, ``stat_type`` indicating which testing method was used (``model_compare``, ``de_novo`` or ``sample_compare``).

The per-base statistics are stored in a dataset, ``stats``, containing one record for each genomic base. Each record contains the following attributes: ``frac``, ``pos``, ``chrm``, ``strand``, ``cov``, ``control_cov``, and ``valid_cov``.

``pos``, ``chrm`` and ``strand`` define the zero-based genomic position for this record.

``frac`` contains the fraction of valid (not including per-read statistics within the interval specified by ``--single_read_threshold``) reads at this genomic position identified as the standard base.

``cov``, ``control_cov``, and ``valid_cov`` contain the read coverage at the genomic position for the sample and control reads. ``control_cov`` is only applicable for the control sample comparison testing method. ``valid_cov`` contains the number of reads contributing to the ``frac`` of tested reads as defined by ``--single-read-threshold``.

-------------------------------
Per-read Statistics File Format
-------------------------------

Per-read statistics can be stored by setting the ``--per-read-statistics-basename`` option to any ``detect_modifications`` command. This output file can then be used in downstream Tombo sub-commands (e.g. the ``plot per_read`` and ``detect_modifications aggregate_per_read_stats`` commands).

For advanced users, the Tombo per-read statsitics file can be accessed via the Tombo python API using the ``tombo.tombo_stats.PerReadStats`` class. This class provides initialization, simply taking the per-read statsitics filename. The ``PerReadStats`` class supports the ``get_region_stats`` function which takes a ``tombo.tombo_helper.intervalData`` object specifying an interval of interest. This will return a numpy array containing a record for each read (specified by the ``read_id`` field) and each tested genomic position (``pos`` field) along with the test statistic (``stat`` field) at that location.

.. important::
   
   All other optional arguments to the ``tombo.tombo_stats.PerReadStats`` constructor should be left as ``None``; setting these values will delete the file and construct a blank per-read statistics file.

The per-read statistics file is in the HDF5 format. All blocks are stored within the ``Statistic_Blocks`` slot. The size of the blocks is stored in the ``block_size`` attribute (defined by the ``--multiprocess-region-size`` option) and the type of statistical test applied is stored in the ``stat_type`` attribute.

Each genomic block is stored in a different ``Block_[##]`` slot. These slots do not have any particular order. Within each block the ``chrm``, ``strand`` and ``start`` of the block are stored. The block statistics are stored in the ``block_stats`` slot. Per-read statistics contain a record for each tested location within each read. Each record contains the genomic position (``pos``), the test statistic (``stat``; hypothesis test p-value or log likelihood ratio as indicated by the statistic type), and the ``read_id``. A single read spanning multiple blocks will contain statistics in more than one block. An individual read's statistics can be reconstructed using the ``read_id`` field.
