**************
Tombo Examples
**************

Below are minimal use case examples. For more detail on each commands options and further algorithm details, please see the corresponding documentation sections.

------------------------------------------
Re-squiggle (Raw Signal Genomic Alignment)
------------------------------------------

The re-squiggle algorithm aligns raw signal (electric current nanopore measurements) to genomic sequence based on a genomic mapping.

This command will add infomation including the mapped genomic location and the raw signal to sequence assignment to the read files (in FAST5 format) provided, as well as producing an index file for more efficient file access in downstream commands.

The ``resquiggle`` command must be run before any further processing by Tombo commands.

**Important note**: Currently, only models for R9.4/5 (1D or 1D^2) DNA or RNA sequencing are included with Tombo. Analysis of other nanopore data types is not supported at this time. If DNA or RNA sample type is not explicitly specified (via ``--dna`` or ``--rna`` options) the sample type will be detected automatically for the set of reads.

For more details see the :doc:`re-squiggle documentation </resquiggle>`.

.. code-block:: bash

    # optionally annotate raw FAST5s with FASTQ files produced from the same reads
    tombo annotate_raw_with_fastqs --fast5-basedir <fast5s-base-directory> --fastq-filenames reads.fastq

    tombo resquiggle <fast5s-base-directory> <reference-fasta> --processes 4

-----------------------
Modified Base Detection
-----------------------

Tombo provides three methods for the investigation of modified bases. Each method has different advantages and requirements.

* The specific alternative base method is preferred. Alternative DNA models are currently available for 5-methylcytosine (5mA) and N6-methyladenosine (6mA) in all sequence contexts.
  
  - More modifications will continue to be added.
* The canonical (control) sample comparison method requires the production of a second set of reads containing only the 4 canonical bases (e.g PCR).
* The de novo method compares signal to the included canonical bases model.

  - This method is recommended only as a research tool and may produce high false positive rates.
* Both the control sample comparison and the de novo methods may not identify the exact modified base location and give no information as to the identity of a modified base.

----

.. figure::  _images/testing_method_comparison.png
   :align: center
   :scale: 30%
   
   Tombo modified base testing methods.

----

The result of all ``test_significance`` calls will be a binary statistics file(s), which can be passed to other Tombo sub-commands.

For more details see the :doc:`modified base detection documentation </modified_base_detection>`.

Specific Alternative Base Method
================================

In order to specifically detect 5mC and 6mA, use the ``test_significance`` command with the ``--alternate-bases 5mC 6mA`` option.

This will perform a log likelihood ratio test using the default canonical model and the 5mC and 6mA alternative models provided with Tombo.

New alternative base models will be added as they are trained. This is the perferred method for modified base detection if a model is available for your biological modification of interest as it identifies the exact location of the modified base and reduces false positives for spurious shifts in signal.

.. code-block:: bash

    tombo test_significance --fast5-basedirs <fast5s-base-directory> \
        --alternate-bases 5mC 6mA --statistics-file-basename sample_alt_model

Canonical Sample Comparison Method
==================================

In order to perform canonical-sample-comparison modified base detection, use the ``test_significance`` command with a second set of reads from the same biological sample containing only canonical bases (e.g. PCR) using the ``--control-fast5-basedirs`` option.

This will perform a hypothesis test against the signal level observed from the control sample at each genomic position. In some cases this method provides the highest accuracy, but does not always identify the exact modified base position.

.. code-block:: bash

    tombo test_significance --fast5-basedirs <fast5s-base-directory> \
        --control-fast5-basedirs <control-fast5s-base-directory>  \
        --statistics-file-basename sample_canonical_compare

De novo Non-canonical Base Method
=================================

In order to perform de novo non-canonical base detection, use the ``test_significance`` command without any other options (aside from the set of reads to test).

This will perform a hypothesis test against the default canonical base model provided with Tombo. Note that this method is quite error prone and may result in a high false positive rate, but may be of use in a research and development setting. This method also has the lowest requirement of only a set of reads and a genome.

.. code-block:: bash

    tombo test_significance --fast5-basedirs <fast5s-base-directory> \
        --statistics-file-basename sample_de_novo_detection

-----------
Text Output
-----------

Wiggle Format Output
====================

In order to output the results of re-squiggling and statistical testing in a genome browser compatible format (`wiggle format <https://genome.ucsc.edu/goldenpath/help/wiggle.html>`_), the ``write_wiggles`` sub-command is provided.

.. code-block:: bash

    tombo write_wiggles --fast5-basedirs <fast5s-base-directory> --wiggle-basename sample_alt_model \
        --statistics-filename sample_alt_model.5mC.tombo.stats --wiggle-types dampened_fraction coverage

.. hint::

    Other ``--wiggle-types`` available are ``fraction``, ``signal``, ``signal_sd``, ``dwell`` and ``difference``.

    The ``dampened_fraction`` option adds psuedo-counts to the detected number of un-modified and modified reads at each tested location (as specified by the ``--coverage-dampen-counts`` option), while the ``fraction`` option returns the raw fraction of modified reads at any genomic site. The ``dampen_fraction`` output is intended to allow the inclusion of low coverage regions in downstream analysis without causing potentially false site at the top of rank lists. Visualize different values of the ``--coverage-dampen-counts`` option with the included ``scripts/test_beta_priors.R`` script.

Genome Sequence Output
======================

For modified base analysis pipelines (e.g. motif detection), it may be useful to output the genomic sequence surrounding locations with the largest fraction of modified reads. The ``write_most_significant_fasta`` sub-command is provided for this purpose.

.. code-block:: bash

    tombo write_most_significant_fasta --statistics-filename sample_alt_model.6mA.tombo.stats \
        --genome-fasta <reference-fasta>

Example `meme <http://meme-suite.org/doc/meme.html>`_ command line modified base motif detection command.

.. code-block:: bash

   ./meme -oc motif_output.meme -dna -mod zoops tombo_results.significant_regions.fasta
    
For more details see the :doc:`text output documentation </text_output>`.

-----------------
Plotting Examples
-----------------

Tombo provides many plotting functions for the visualization of modified bases and raw nanopore signal in general.

Most plotting commands are genome-anchored. That is the raw signal is plotted as the re-squiggle algorithm has assigned it to the genome. Thus each read contains a different number of raw observations assigned to each genomic base. For summary distributions (overplotting optios not showing raw signal) the distributions are taken over each read's average signal level at the genomic position.

Each genome anchored plotting command allows for the selection of genomic positions based on generally applicable criterion.

.. code-block:: bash

    tombo plot_max_coverage --fast5-basedirs <fast5s-base-directory> --plot-standard-model
    
    tombo plot_motif_centered --fast5-basedirs <fast5s-base-directory> --motif AWC \
        --genome-fasta genome.fasta --control-fast5-basedirs <control-fast5s-base-directory>
    
    tombo plot_per_read --per-read-statistics-filename <per_read_statistics.tombo.stats> \
        --genome-locations chromosome:1000 chromosome:2000:- \
        --genome-fasta genome.fasta

For more details see the :doc:`plotting documentation </plotting>`.

.. tip::

    For additional command details, see the specific commands documentation section.
