***********************
Modified Base Detection
***********************

Tombo enables three methods for detecting shifts in current signal level, indicative of non-canonical bases. These three methods allow researchers to investigate non-canonical bases given any sample type, while also enabling more accurate detection of specific known modifications (currently only 5-methyl cytosine, but more coming soon).

----

.. figure::  _images/testing_method_comparison.png
   :align: center
   :scale: 30%
   
   Tombo modified base testing methods.

----

All three methods are accessed by the ``test_significance`` tombo subcommand as described below.

TL;DR:

* To identify 5-methyl cytosine (5mC) run ``test_significance`` with the ``--alternate-bases 5mC`` option
* For more experimental de novo modified base detection simply run ``test_significance`` with a set of reads
* For modified base detection via comparison to a control sample (e.g. PCR) run ``test_significance`` with a control set of reads (``--control-fast5-basedirs``)
* The ``test_significance`` command will produce a binary file (not intended for direct use)
  - To extract useful text files use the ``write_wiggles`` command
  - To visualize raw signal around significant regions use the ``plot_most_significant`` command

..

   Note that the ``resquiggle`` command must be run on a set of reads before processing with ``test_significance``.

-------------------
Statistical Testing
-------------------

For all statistical testing methods, the result is a binary Tombo statistics file. This file contains statistics associated with each genomic base (per-read output is not currently supported, but may be in the future). This file is not intended for direct use, but several other Tombo commands (``write_wiggles``, ``write_most_significant_fasta``, ``plot_most_significant``, etc.) take the statistics file as an input, accomodating many user pipelines downstream of modified base detection.

Of particular interest, the statistics file contains the fraction of reads at each genomic position passing a set threshold (``--single-read-threshold``). This value is set to a default of 1% p-value for hypothesis testing methods (de novo and control sample comparison) and a log likelihood ratio of 2 for the alternative model testing method. Note that for likelihood ratio test fractions, some reads will fall between the +/- threshold values and so the sum of ``frac`` and ``alt_frac`` may not be 1.

For the de novo and alternative model testing approaches a default canonical model is used (included with Tombo code base). Users may also train their own canonical Tombo model (possibly for an older chemistry version) and test agaist this model using the ``--tombo-model-filename`` option. See more in the Model Training documentation section.

5mC (DNA) Detection
===================

In order to specifically detect 5-methyl cytosine (as well as other alternative bases soon), use the ``test_significance`` command with the ``--alternate-bases`` option. Users may also train their own alternative base Tombo models and test against these with the ``--alternate-model-filenames`` option. See more details in the Model Training documentation section.

This will perform a log likelihood ratio test using the default canonical and 5mC alternative models provided with Tombo. This likelihood ratio is computed over all positions modeled. The default model is a 6-mer, so the signal at the 6 surrounding genomic contexts, contribute to the log likelihood ratio test statistic at any one position.

For example with a **C** found in a TGGTA **C** GTCCG context, the signal will be tested against expected canonical and alternative distributions at the following locations::

  TGGTA **C** GTCCG
  -----------------
  TGGTA **C**
   GGTA **C** G
    GTA **C** GT
     TA **C** GTC
      A **C** GTCC
        **C** GTCCG

New alternative base models will be added as they are trained. This is the perferred method for modified base detection if a model is available for your biological sample of interest.

.. code-block:: bash

    tombo test_significance --fast5-basedirs <fast5s-base-directory> \
        --alternate-bases 5mC --statistics-file-basename sample_5mC_detection

    # if you have trained you own alternative base model
    tombo test_significance --fast5-basedirs <fast5s-base-directory> \
        --alternate-model-filename exotic_base.tombo.model \
        --statistics-file-basename sample_exotic_base_detection

Canonical Sample Comparison Detection
=====================================

In order to perform *canonical sample comparison* modified base detection, use the ``test_significance`` command with a second set of reads from the same biological sample containing only canonical bases (e.g. PCR) via the ``--control-fast5-basedirs``.

For each sample read, this will perform a hypothesis test against a normal distribution estimated from the the signal level observed from the control sample at each genomic position. This method provides the highest accuracy (as effects outside of the default modeled 6-mer are accounted for in the control sample), but does not always identify the exact modification position or the identity of the modified base.

Note that no global model is used in the application of this method. Instead the testing null distribution is estimated at each genomic location.

For both this method, as well as the canonical model method, the ``--fishers-method-context`` option will combine test values, using `Fisher's Method <https://en.wikipedia.org/wiki/Fisher%27s_method>`_, over a moving window extending a number of positions in either direction. Due to the nature of nanopore sequencing, the genomic context surrounding the read head effect that current at any position. Thus shifts in signal due to a modified base may occur at several positions to either side of the true modified location. Thus combining statistical test values across several genomic positions can help to center significant values on the truely modified position. The default value for this parameter is 1, but reasonable results can be obtained for values between 0 and 3.

.. code-block:: bash

    tombo test_significance --fast5-basedirs <fast5s-base-directory> \
        --control-fast5-basedirs <control-fast5s-base-directory> \
        --statistics-file-basename sample_canonical_compare

De novo Non-canonical Base Detection
====================================

In order to perform de novo non-canonical base detection, use the ``test_significance`` command with no other options (or a canonical Tombo model, ``--tombo-model-filename``, if not using the default canonical Tombo model).

For each read, this will perform a hypothesis test against the canonical model based on the genomic sequence at each position. Note that this method can be quite error prone and will likely result in a high number of false positives, but may be of use in a research and development setting. This method also has the least requirements of the three methods, requiring only a set of reads and a genome, allowing any nanopore researcher to start investigating modified bases.

.. code-block:: bash

    tombo test_significance --fast5-basedirs <fast5s-base-directory> \
        --statistics-file-basename sample_de_novo_detection

----------------
Multi-processing
----------------

Tombo statistical testing provides the option to perform testing spread acorss multiple processes. This also allows users to limit the memory requirement for processing statistical testing, as all events for all reads covering a region must be read into memory at once to perform testing. If the ``test_significance`` command seems to be using too much memory, consider setting the ``--multiprocess-region-size`` to a lower value.

Multi-processing is performed over batches delineated by regular intervals across chromosomes covered by at least one read. The interval size is determined by the ``--multiprocess-region-size`` option and processed by ``--processes`` individual processes independently. The produced per-base results are identical no matter the multi=processing options selected.

----------------------------
Tombo Statistics File Format
----------------------------

While the Tombo statistics file is meant to be a binary file not processed by outside tools its contents are described here. The Tombo statistics file is in the HDF5 format. There is one attribute at the root level, ``stat_type`` indicating which testing method was used (``model_compare``, ``de_novo`` or ``sample_compare``).

The per-base statistics are stored in a dataset, ``stats``, containing one record for each genomic base. Each record contains the following attributes: ``stat``, ``mt_stat``, ``frac``, ``alt_frac``, ``pos``, ``chrm``, ``strand``, ``cov`` and ``control_cov``.

``pos``, ``chrm`` and ``strand`` define the genomic position for this record. Note that position is 0-based indexing, and care should be taken if using this to compare to other genomic data sets.

``frac`` and ``alt_frac`` contain the fraction of reads at this genomic position producing a significant result as defined by ``--single-read-threshold``. Note that for log likelihood test results, via the alternative model comparison method, the sum of ``frac`` and ``alt_frac`` may not sum to 1.

``cov`` and ``control_cov`` contain the read coverage at the genomic position for the sample and control reads (control coverage is only applicable for the control sample comparison testing method).

``stat`` and ``mt_stat`` contain the "averaged" statistical testing values at each genomic location. ``mt_stat`` contains the `Benjamini-Hochberg <https://en.wikipedia.org/wiki/False_discovery_rate#Benjamini%E2%80%93Hochberg_procedure>`_ multiple testing corrected testing value. For the hypothesis testing based methods (de novo and control sample comparison) this "average" is taken as the `Fisher's method <https://en.wikipedia.org/wiki/Fisher%27s_method>`_ combination of p-values for each read at this genomic position. Note that this is highly dependent on read coverage and thus should be used with caution. For the alternative model comparison method, the "average" is taken as the average log likelihood test statistics over all reads covering this genomic position.
