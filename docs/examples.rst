**************
Tombo Commands
**************

This page contains breif descriptions of the most common Tombo commands. For more detail on each commands' options and further algorithm details, please see the corresponding documentation sections.

------------------------------------------
Re-squiggle (Raw Signal Genomic Alignment)
------------------------------------------

The re-squiggle algorithm aligns raw signal (electric current nanopore measurements) to reference genomic or transcriptomic sequence.

One of the major assumptions of the re-squiggle algorithm is that the provided reference sequence is correct. Thus for a poorly assembled reference or divergent sample, an assembly polishing step (possibly from the same data/sample) may improve results.

The ``resquiggle`` command will add infomation (the mapped reference location and the raw signal to sequence assignment) to the read files provided (in FAST5 format), as well as producing an index file for more efficient file access in downstream commands.

.. important::

   The ``resquiggle`` command must be run before any further processing by Tombo commands.

**Note**: Tombo currenly includes default canonical models for both DNA or RNA data (supporting R9.4 and R9.5; 1D and 1D^2; R9.*.1 chemistries). Analysis of other nanopore data types is not supported at this time (e.g. R7 data). If DNA or RNA sample type is not explicitly specified (via ``--dna`` or ``--rna`` options) the sample type will be detected automatically from the raw read files.

For more details see the :doc:`re-squiggle documentation </resquiggle>`.

.. code-block:: bash

    # annotate raw FAST5s with FASTQ files produced from the same reads
    # skip this step if raw read files already contain basecalls
    tombo preprocess annotate_raw_with_fastqs --fast5-basedir <fast5s-base-directory> --fastq-filenames <reads.fastq>

    tombo resquiggle <fast5s-base-directory> <reference-fasta> --processes 4 --num-most-common-errors 5

-----------------------
Modified Base Detection
-----------------------

Tombo provides four (including two types of sample comparison) methods for the investigation of modified bases (within the ``tombo detect_modifications`` command group). Each method has different advantages and requirements.

----

.. figure::  _images/testing_method_comparison.png
   :align: center

   Tombo modified base testing methods.

----

All modified base detection methods, except the ``level_sample_compare`` method, produce per-read, per-reference position test statistics (which can be saved via the ``--per-read-statistics-basename`` option). A threshold is then applied to these statistics to produce an estimate for the fraction of reads that appear modified at each reference location.

1. **Specific alternative base detection (recommended)**

  - Run using ``tombo detect_modifications alternative_model`` command.
  - This method identifies sites where signal matches a specific alternative base expected signal levels better than the canonical expected levels producing a statistic similar to a log likelihood ratio.
  - All-context alternative DNA models are currently available for 5-methylcytosine (5mC) and N6-methyladenosine (6mA; not currently available for RNA).
  - More accurate motif specific models are available for dam and dcm methylation (found in E. coli) and CpG methylation (found in human samples).

2. **De novo canonical model comparison**

  - Run using ``tombo detect_modifications de_novo`` command.
  - This method identifies sites where signal deviates from the expected canonical signal levels.
  - While this method has the highest error rates, it can be used effectively on any sample and is particularly useful for motif discovery for motif-specific modifications (e.g. bacterial samples).

3. **Sample comparison**

   - These two methods requires the production of a second set of reads for comparison.

   a. ``level_sample_compare``

      - Run using ``tombo detect_modifications level_sample_compare`` command.
      - This method performs either a ks-test (default), u-test or t-test to identify sites where signal levels deviate between two samples.
      - This method is ideal for high coverage samples (in order to accurately estimate the effect size measures) and comparison of two potentailly non-canonical samples (e.g. a KO experiment).

   b. ``model_sample_compare``

      - Run using ``tombo detect_modifications model_sample_compare`` command.
      - This uses a canonical control sample to adjust the expected signal levels due to un-modeled effects.
      - This adjusted model is then used to test for modifications as in the de novo method.
      - This was the ``sample_compare`` method prior to version 1.5.

..

    Both the sample comparison and the *de novo* methods may not identify the exact modified base location (as the shifted signal does not always center on a modified base) and gives no information as to the identity of a modified base.

The result of all ``tombo detect_modifications`` calls will be a binary statistics file(s), which can be passed to other Tombo commands.

For more details see the :doc:`modified base detection documentation </modified_base_detection>`.

-----------
Text Output
-----------

Genome Browser File Output
==========================

In order to output re-squiggle and/or modified base detection results in a genome browser compatible format (either `wiggle format <https://genome.ucsc.edu/goldenpath/help/wiggle.html>`_ or `bedgraph format <https://genome.ucsc.edu/goldenpath/help/bedgraph.html>`_), the ``tombo text_output genome_browser`` command is provided.

.. code-block:: bash

    tombo text_output browser_files --fast5-basedirs <fast5s-base-directory> \
        --statistics-filename sample_alt_model.CpG.tombo.stats \
        --browser-file-basename sample_alt_model --file-types dampened_fraction coverage

.. hint::

    All ``--file-types`` available are:

      - ``fraction``, ``dampened_fraction``, and ``valid_coverage`` derived from a (non-``level_sample_compare``) statistics file
      - ``statistic`` derived from a ``level_sample_compare`` statistics files
      - ``coverage`` derived from the Tombo index (fast)
      - ``signal``, ``signal_sd``, ``dwell``, and ``difference`` derived from read FAST5 files (slow)

Reference Sequence Output
=========================

For modified base analysis pipelines (e.g. motif detection), it may be useful to output the reference sequence surrounding the most likely modified sites. The ``text_output signif_sequence_context`` command is provided for this purpose.

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

Tombo provides many plotting functions for the visualization of modified bases within raw nanopore signal.

Most plotting commands are reference-anchored. That is the normalized raw signal is plotted as the re-squiggle algorithm has assigned it to the reference sequence.

Each reference anchored plotting command allows for the selection of reference positions based on generally applicable criterion.

.. code-block:: bash

    tombo plot max_coverage --fast5-basedirs <fast5s-base-directory> --plot-standard-model

    tombo plot motif_centered --fast5-basedirs <fast5s-base-directory> --motif AWC \
        --genome-fasta genome.fasta --control-fast5-basedirs <control-fast5s-base-directory>

    tombo plot per_read --per-read-statistics-filename <per_read_statistics.tombo.stats> \
        --genome-locations chromosome:1000 chromosome:2000:- \
        --genome-fasta genome.fasta

.. note::

   For regions with higher coverage, several over-plotting options are available. For those options producing a distribution, these are taken over each reads average signal assigned to a base. This requires extraction of these levels from all relevant FAST5 files and thus can be slow for very deep coverage regions.

For more details see the :doc:`plotting documentation </plotting>`.

--------------
Read Filtering
--------------

Read filtering commands can be useful to extract the most out out of a set of reads for modified base detection. Read filtering commands effect only the Tombo index file, and so filters can be cleared or applied iteratively without re-running the re-squiggle command. Five filters are currently made available (``genome_locations``, ``raw_signal_matching``, ``q_score``, ``level_coverage`` and ``stuck``).

.. code-block:: bash

    # filter reads to a specific genomic location
    tombo filter genome_locations --fast5-basedirs path/to/native/rna/fast5s/ \
        --include-regions chr1:0-10000000

    # apply a more strigent observed to expected signal score (default: 1.1 for DNA reads)
    tombo filter raw_signal_matching --fast5-basedirs path/to/native/rna/fast5s/ \
        --signal-matching-score 1.0

For more details see the :doc:`filter documentation </filtering>`.
