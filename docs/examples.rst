**************
Tombo Commands
**************

Below are minimal use case examples. For more detail on each commands' options and further algorithm details, please see the corresponding documentation sections.

------------------------------------------
Re-squiggle (Raw Signal Genomic Alignment)
------------------------------------------

The re-squiggle algorithm aligns raw signal (electric current nanopore measurements) to genomic/transcriptomic sequence.

The ``resquiggle`` command will add infomation (the mapped genomic location and the raw signal to sequence assignment) to the read files provided (in FAST5 format), as well as producing an index file for more efficient file access in downstream commands.

.. important::
   
   The ``resquiggle`` command must be run before any further processing by Tombo commands.

**Note**: Tombo currenly includes default canonical models for both DNA or RNA data (including R9.4 and R9.5; 1D and 1D^2; R9.*.1 chemistries). Analysis of other nanopore data types is not supported at this time (e.g. R7 data). If DNA or RNA sample type is not explicitly specified (via ``--dna`` or ``--rna`` options) the sample type will be detected automatically from the raw read files.

For more details see the :doc:`re-squiggle documentation </resquiggle>`.

.. code-block:: bash

    # annotate raw FAST5s with FASTQ files produced from the same reads if the raw files do not contain FASTQ information
    tombo annotate_raw_with_fastqs --fast5-basedir <fast5s-base-directory> --fastq-filenames <reads.fastq>

    tombo resquiggle <fast5s-base-directory> <reference-fasta> --processes 4

-----------------------
Modified Base Detection
-----------------------

Tombo provides three methods for the investigation of modified bases (within the ``detect_modifications`` command group). Each method has different advantages and requirements.

All modified base detection methods poduce per-read, per-genomic position test statistics (which can be saved via the ``--per-read-statistics-basename`` option). A threshold is then applied to these statistics to produce a fraction of reads that appear modified at each genomic locaiton.

1. Specific alternative base detection

  - Run using ``tombo detect_modifications alternative_model`` command.
  - This method identifies signal that deviates from the canonical base expected signal level while matching a specific alternative base expected signal level.
  - This method produces a statistic similar to a log likelihood ratio, but scaled to be more robust to outlier signal assignments (similar to `Tukey's biweight function <http://mathworld.wolfram.com/TukeysBiweight.html>`_).
  - Alternative DNA models are currently available for 5-methylcytosine (5mC) and N6-methyladenosine (6mA) in all sequence contexts.
  - An alternative RNA model is available for 5mC.

2. *De novo* canonical model comparison

  - Run using ``tombo detect_modifications de_novo`` command.
  - This method compares re-squiggled signal to the default canonical model.
  - While this method may produce significant false positive and negative results per-read, it produces the best results for many statistical measures per-genomic location (fraction of modified bases across a set of reads).

3. Canonical (control) sample comparison

  - Run using ``tombo detect_modifications sample_compare`` command.
  - This method performs a hypothesis test against the distribution estimated from the control sample at each base.
  - This method requires the production of a second set of reads containing only the 4 canonical bases (e.g PCR for DNA of IVT for RNA).

..

    Both the control sample comparison and the *de novo* methods may not identify the exact modified base location (as the shifted signal does not always center on a modified base) and gives no information as to the identity of a modified base.

----

.. figure::  _images/testing_method_comparison.png
   :align: center
   :scale: 30%
   
   Tombo modified base testing methods.

----

The result of all ``detect_modifications`` calls will be a binary statistics file(s), which can be passed to other Tombo commands.

For more details see the :doc:`modified base detection documentation </modified_base_detection>`.

Specific Alternative Base Method
================================

In order to specifically detect 5mC and 6mA, use the ``detect_modifications alternative_model`` command.

This will compute a statistic similar to a log likelihood ratio using the default canonical model and the 5mC and 6mA alternative DNA models provided with Tombo.

This is the perferred method for modified base detection if a model is available for your biological modification of interest, as it identifies the exact location of the modified base and reduces false positives for spurious shifts in signal.

.. code-block:: bash

    tombo detect_modifications alternative_model --fast5-basedirs <fast5s-base-directory> \
        --alternate-bases 5mC 6mA --statistics-file-basename sample_alt_model

*De novo* Non-canonical Base Method
===================================

In order to perform *de novo* non-canonical base detection, use the ``detect_modifications de_novo`` command.

This will perform a hypothesis test against the default canonical base model provided with Tombo. Note that this method is quite error prone and may result in a high false positive rate on a per-read basis, but may be of use in a research and development setting. This method also has the lowest requirement, consisting of only a set of reads potentially containing modifications and a reference sequence.

.. code-block:: bash

    tombo detect_modifications de_novo --fast5-basedirs <fast5s-base-directory> \
        --statistics-file-basename sample_de_novo_detection

Canonical Sample Comparison Method
==================================

In order to execute the canonical sample comparison method, use the ``detect_modifications sample_compare`` command.

This will perform a hypothesis test against the signal level observed from the control sample (provided via ``--control-fast5-basedirs`` option) at each genomic position. This method currently performs the worst, but future updates to this method may increase the accuracy of this method. This method (like the ``de_novo`` method) does not always identify the exact modified base position.

.. code-block:: bash

    tombo detect_modifications sample_compare --fast5-basedirs <fast5s-base-directory> \
        --control-fast5-basedirs <control-fast5s-base-directory>  \
        --statistics-file-basename sample_canonical_compare

-----------
Text Output
-----------

Genome Browser File Output
==========================

In order to output the results of re-squiggling and statistical testing in a genome browser compatible format (either `wiggle format <https://genome.ucsc.edu/goldenpath/help/wiggle.html>`_ or `bedgraph format <https://genome.ucsc.edu/goldenpath/help/bedgraph.html>`_), the ``tombo text_output genome_browser`` command is provided.

.. code-block:: bash

    tombo text_output genome_browser --fast5-basedirs <fast5s-base-directory> \
        --statistics-filename sample_alt_model.5mC.tombo.stats \
         --browser-file-basename sample_alt_model --file-types dampened_fraction coverage

.. hint::

    Other ``--file-types`` available are ``fraction``, ``valid_coverage``, ``signal``, ``signal_sd``, ``dwell`` and ``difference``.

    The ``dampened_fraction`` option adds psuedo-counts to the detected number of un-modified and modified reads at each tested location (as specified by the ``--coverage-dampen-counts`` option), while the ``fraction`` option returns the raw fraction of modified reads at any genomic site from ``detect_modifications`` results. The ``dampen_fraction`` output is intended to allow the inclusion of low coverage regions in downstream analysis without causing potentially false site at the top of rank lists. Visualize different values of the ``--coverage-dampen-counts`` option with the included ``scripts/test_beta_priors.R`` script.

Genome Sequence Output
======================

For modified base analysis pipelines (e.g. motif detection), it may be useful to output the genomic sequence surrounding locations with the largest fraction of modified reads. The ``text_output signif_sequence_context`` command is provided for this purpose.

.. code-block:: bash

    tombo text_output signif_sequence_context --statistics-filename sample_alt_model.6mA.tombo.stats \
        --genome-fasta <reference-fasta> --sequences-filename sample_alt_model.6mA.most_signif.fasta

Example `meme <http://meme-suite.org/doc/meme.html>`_ command line modified base motif detection command.

.. code-block:: bash

   ./meme -oc motif_output.meme -dna -mod zoops sample_alt_model.6mA.most_signif.fasta
    
For more details see the :doc:`text output documentation </text_output>`.

-----------------
Plotting Commands
-----------------

Tombo provides many plotting functions for the visualization of modified bases and raw nanopore signal in general.

Most plotting commands are genome-anchored. That is the raw signal is plotted as the re-squiggle algorithm has assigned it to the genome. Thus each read may contain a different number of raw observations assigned to each genomic base. For regions with higher coverage, several over-plotting options are available. For those options producing a distribution, these are taken over each reads average signal assigned to a base. This requires extraction of these levels from all relevant FAST5 files and thus can be slow for very deep coverage regions.

Each genome anchored plotting command allows for the selection of genomic positions based on generally applicable criterion.

.. code-block:: bash

    tombo plot max_coverage --fast5-basedirs <fast5s-base-directory> --plot-standard-model
    
    tombo plot motif_centered --fast5-basedirs <fast5s-base-directory> --motif AWC \
        --genome-fasta genome.fasta --control-fast5-basedirs <control-fast5s-base-directory>
    
    tombo plot per_read --per-read-statistics-filename <per_read_statistics.tombo.stats> \
        --genome-locations chromosome:1000 chromosome:2000:- \
        --genome-fasta genome.fasta

For more details see the :doc:`plotting documentation </plotting>`.

.. tip::

    For additional command details, see the specific commands documentation section.
