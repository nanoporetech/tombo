**************
Tombo Examples
**************

Below are minimal use case examples. For more detail on each commands options and further algorithm details, please see the corresponding documentation sections.

----------------------------------
Re-squiggle (Raw Signal Alignment)
----------------------------------

The re-squiggle algorithm aligns raw signal to genomic sequence based on a genomic mapping.

This command will add infomation including the mapped genomic location and the raw signal to sequence assignment into the FAST5 files provided as well as producing an index file for more efficient file access in downstream commands.

The ``resquiggle`` command must be run on any set of FAST5s before any further processing by Tombo.

**Important note**: Currently, only a models for R9.4/5 sequencing DNA or RNA are included with Tombo. Analysis of other nanopore data types is not supported at this time.

For more details see the :doc:`re-squiggle documentation </resquiggle>`.

.. code-block:: bash

    # optionally annotate raw FAST5s with FASTQ files produced from the same reads
    tombo annotate_raw_with_fastqs --fast5-basedir <fast5s-base-directory> --fastq-filenames reads.fastq

    tombo resquiggle [-h] <fast5s-base-directory> <reference-fasta> --minimap2-executable ./minimap2

-----------------------
Modified Base Detection
-----------------------

Tombo provides three methods for the investigation of modified bases. Each method has different advantages and requirements.

* The specific alternative base method is preferred, but is currently only available for 5mC in DNA (more modifications coming soon).
* The canonical (control) sample comparison method would be preferred next, but requires the producetion of a second set of reads containing only the 4 canonical bases.
* The de novo method is recommended only a research tool and may produce high false positive rates.
* Additionally, both the control sample and de novo methods may not identify the exact modified base location and give no information as to the identity of a modified base.

The result of all ``test_significance`` calls will be a binary statistics file, which can be passed to other Tombo sub-commands.

For more details see the :doc:`modified base detection documentation </modified_base_detection>`.

Specific Alternative Base Detection
===================================

In order to specifically detect 5-methyl cytosine (and other alternative bases in the future), use the ``test_significance`` command with the ``--alternate-bases 5mC`` option.

This will perform a log likelihood ratio test using the default canonical and 5mC alternative models provided with Tombo.

New alternative base models will be added as they are trained. This is the perferred method for modified base detection if a model is available for your biological sample of interest.

.. code-block:: bash

    tombo test_significance --fast5-basedirs <fast5s-base-directory> \
        --alternate-bases 5mC --statistics-file-basename sample_5mC_detection

Canonical Sample Comparison Detection
=====================================

In order to perform canonical-sample-comparison modified base detection, use the ``test_significance`` command with a second set of reads from the same biological sample containing only canonical bases (e.g. PCR) using the ``--control-fast5-basedirs`` option.

This will perform a hypothesis test against the signal level observed from the control sample at each genomic position. This method provides the highest accuracy, but does not always identify the exact modification position or the identity of the modified base.

.. code-block:: bash

    tombo test_significance --fast5-basedirs <fast5s-base-directory> \
        --control-fast5-basedirs <control-fast5s-base-directory>  \
        --statistics-file-basename sample_canonical_compare

De novo Non-canonical Base Detection
====================================

In order to perform de novo non-canonical base detection, use the ``test_significance`` command without any other options (aside from the set of reads to test).

This will perform a hypothesis test against the default canonical base model provided with Tombo. Note that this method is quite error prone and will likely result in a high number of false positives, but may be of use in a research and development setting. This method also has the lowest requirement of only a set of reads and a genome.

.. code-block:: bash

    tombo test_significance --fast5-basedirs <fast5s-base-directory> \
        --statistics-file-basename sample_de_novo_detection

-----------
Text Output
-----------

Wiggle Format Output
====================

In order to output the results of re-squiggling and statistical testing in a genome browser compatible format (`wiggle format <https://genome.ucsc.edu/goldenpath/help/wiggle.html>`_), the ``write_wiggles`` subcommand is provided.

.. code-block:: bash

    tombo write_wiggles --fast5-basedirs <fast5s-base-directory> --wiggle-basename sample_5mC_detection \
        --statistics-filename sample_5mC_detection.5mC.tombo.stats --wiggle-types fraction coverage

Genome Sequence Output
======================

For modified base analysis pipelines (e.g. motif detection), it may be of use to output the genomic sequence surrounding locations with the largest fraction of modified reads. The ``write_most_significant_fasta`` sub-command is provided for this purpose.

.. code-block:: bash

    tombo write_most_significant_fasta --statistics-filename sample_5mC_detection.5mC.tombo.stats \
        --genome-fasta <reference-fasta>

Example `meme <http://meme-suite.org/doc/meme.html>`_ command line modified base motif detection command.

.. code-block:: bash

   ./meme -oc motif_output.meme -dna -mod zoops tombo_results.significant_regions.fasta
    
For more details see the :doc:`text output documentation </text_output>`.

-----------------
Plotting Examples
-----------------

Tombo provides many plotting functions for the visualization of potentially modified bases and the raw nanopore signal in general.

Most plotting commands are genome-anchored. That is the raw signal is plotted as the re-squiggle algorithm has assigned it to the genome. Thus each read contain a different number of raw observations per genomic base. For summary distributions (not raw signal) the distributions are taken over each reads average signal level at the genomic position.

Each genome anchored plotting command allows for the selection of genomic positions based on generally applicable criterion.

.. code-block:: bash

    tombo plot_max_coverage --fast5-basedirs <fast5s-base-directory> --plot-standard-model
    
    tombo plot_motif_centered --fast5-basedirs <fast5s-base-directory> --motif AWC \
        --genome-fasta genome.fasta --control-fast5-basedirs <control-fast5s-base-directory>
        
    tombo plot_per_read --fast5-basedirs <fast5s-base-directory> \
        --genome-locations chromosome:1000 chromosome:2000:- --plot-alternate-model 5mC

For more details see the :doc:`plotting documentation </plotting>`.

..

    For additional command details, see the specific commands documentation section.
