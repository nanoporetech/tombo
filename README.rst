=============
Tombo Summary
=============

|travis_badge|

.. |travis_badge| image:: https://travis-ci.org/nanoporetech/tombo.svg?branch=master
    :target: https://travis-ci.org/nanoporetech/tombo

Tombo is a suite of tools primarily for the identification of modified nucleotides from nanopore sequencing data.

Tombo also provides tools for the analysis and visualization of raw nanopore signal.

============
Installation
============

|bioconda_badge| |pypi_badge|

.. |bioconda_badge| image:: https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat-square
    :target: http://bioconda.github.io/recipes/ont-tombo/README.html

.. |pypi_badge| image:: https://badge.fury.io/py/ont-tombo.svg
    :target: https://pypi.org/project/ont-tombo/

Basic tombo installation (python 2.7 and 3.4+ support)

::

    # install via bioconda environment
    conda install -c bioconda ont-tombo

    # or install pip package (numpy install required before tombo for cython optimization)
    pip install numpy
    pip install ont-tombo[full]

===========
Quick Start
===========

Call 5mC and 6mA sites from raw nanopore read files. Then output genome browser `wiggle format file <https://genome.ucsc.edu/goldenpath/help/wiggle.html>`_ for 5mA calls and plot raw signal around most significant 6mA sites.

::

   # skip this step if FAST5 files already contain basecalls
   tombo preprocess annotate_raw_with_fastqs --fast5-basedir path/to/fast5s/ \
       --fastq-filenames basecalls1.fastq basecalls2.fastq \
       --sequencing-summary-filenames seq_summary1.txt seq_summary2.txt \
       --processes 4
   
   tombo resquiggle path/to/fast5s/ genome.fasta --processes 4
   tombo detect_modifications alternative_model --fast5-basedirs path/to/fast5s/ \
       --statistics-file-basename sample.alt_modified_base_detection \
       --per-read-statistics-basename sample.alt_modified_base_detection \
       --alternate-bases 5mC 6mA --processes 4
   
   # produces "estimated fraction of modified reads" genome browser files
   # for 5mC testing
   tombo text_output browser_files --statistics-filename sample.alt_modified_base_detection.5mC.tombo.stats \
       --file-types dampened_fraction --browser-file-basename sample.alt_modified_base_detection.5mC
   # and 6mA testing (along with coverage bedgraphs)
   tombo text_output browser_files --statistics-filename sample.alt_modified_base_detection.6mA.tombo.stats \
       --fast5-basedirs path/to/fast5s/  --file-types dampened_fraction coverage\
       --browser-file-basename sample.alt_modified_base_detection.6mA
   
   # plot raw signal at most significant 6mA locations
   tombo plot most_significant --fast5-basedirs path/to/fast5s/ \
       --statistics-filename sample.alt_modified_base_detection.6mA.tombo.stats \
       --plot-standard-model --plot-alternate-model 6mA \
       --pdf-filename sample.most_significant_6mA_sites.pdf

Detect any deviations from expected signal levels for canonical bases to investigate any type of modification.

::

   tombo resquiggle path/to/fast5s/ genome.fasta --processes 4
   tombo detect_modifications de_novo --fast5-basedirs path/to/fast5s/ \
       --statistics-file-basename sample.de_novo_modified_base_detection \
       --per-read-statistics-basename sample.de_novo_modified_base_detection \
       --processes 4
   
   # produces "estimated fraction of modified reads" genome browser files from de novo testing
   tombo text_output browser_files --statistics-filename sample.de_novo_modified_base_detection.tombo.stats \
       --browser-file-basename sample.de_novo_modified_base_detection --file-types dampened_fraction

..
   
   All of these commands work for RNA data as well, but a transcriptome reference sequence must be provided for spliced transcripts.

=====================
Further Documentation
=====================

Run ``tombo -h`` to see all Tombo command groups and run ``tombo [command-group] -h`` to see all commands within each group.

Detailed documentation for all Tombo commands and algorithms can be found at https://nanoporetech.github.io/tombo/

==============
Tombo Commands
==============

Re-squiggle (Raw Data to Genome Alignment)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The ``resquiggle`` algorithm is the central point for the Tombo tookit. For each nanopore read, this command takes basecalled sequence and the raw nanopore signal values. The basecalled sequence is mapped to a genomic or transcriptomic reference. The raw nanopore signal is assigned to the mapped genomic or transcriptomic sequence based on expected signal levels from an included canonical base model. This anchors each raw signal observation from a read to a genomic position. This information is then leveraged to gain information about the potential location of modified nucleotides either within a single read or across a group of reads from a sample of interest.

::

    tombo resquiggle path/to/fast5s/ reference.fasta --processes 4

..

   - Only R9.4 and R9.5 data is supported at this time (including R9.*.1).
   - DNA or RNA sample type is automatically detected from FAST5s (set explicitly with ``--dna`` or ``--rna``).
   - FAST5 files need not contain ``Events`` data, but must contain ``Fastq`` slot containing basecalls. See ``preprocess annotate_raw_with_fastqs`` for pre-processing of raw FAST5s with basecalled reads.
   - The reference sequence file can be a genome/transcriptome FASTA file or a minimap2 index file.
   - The ``resquiggle`` command must be run before testing for modified bases.

Detect Modified Bases
^^^^^^^^^^^^^^^^^^^^^

There are three methods provided with Tombo to identify modified bases.

For more information on these methods see the `Tombo documentation here <https://nanoporetech.github.io/tombo/modified_base_detection.html>`_.

::

   # Identify deviations from the canoncial expected signal levels that specifically match the
   # expected levels from an alternative base e.g.5mC or 6mA (recommended method)
   tombo detect_modifications alternative_model --fast5-basedirs path/to/native/dna/fast5s/ \
       --alternate-bases 5mC 6mA --statistics-file-basename sample.alt_testing

   # Identify any deviations from the canonical base model
   tombo detect_modifications de_novo --fast5-basedirs path/to/native/dna/fast5s/ \
       --statistics-file-basename sample.de_novo_testing --processes 4

   # comparing to a control sample (e.g. PCR)
   tombo detect_modifications sample_compare --fast5-basedirs path/to/native/dna/fast5s/ \
       --control-fast5-basedirs path/to/amplified/dna/fast5s/ \
       --statistics-file-basename sample.compare_testing

..

    Must run ``resquiggle`` on reads before testing for modified bases.

    All ``detect_modifications`` commands produce a binary Tombo statistics file. For use in text output or plotting region selection see ``text_output browser_files`` or ``plot most_significant`` Tombo commands.

    Specify the ``--per-read-statistics-basename`` option to save per-read statistics for plotting or further processing (acces via the Tombo API).

Text Output
^^^^^^^^^^^

::

   # output estimated fraction  of reads modified at each genomic base and
   # valid coverage (after failed reads, filters and testing threshold are applied) in wiggle format
   tombo text_output browser_files --file-types dampened_fraction --statistics-filename sample.alt_testing.5mC.tombo.stats
   
   # output read coverage depth (after failed reads and filters are applied) in bedgraph format
   tombo text_output browser_files --file-types coverage --fast5-basedirs path/to/native/dna/fast5s/

..

    For more text output commands see the `Tombo text output documentation here <https://nanoporetech.github.io/tombo/text_output.html>`_.

Raw Signal Plotting
^^^^^^^^^^^^^^^^^^^

::

    # plot raw signal with standard model overlay at reions with maximal coverage
    tombo plot max_coverage --fast5-basedirs path/to/native/rna/fast5s/ --plot-standard-model
    
    # plot raw signal along with signal from a control (PCR) sample at locations with the AWC motif
    tombo plot motif_centered --fast5-basedirs path/to/native/rna/fast5s/ \
        --motif AWC --genome-fasta genome.fasta --control-fast5-basedirs path/to/amplified/dna/fast5s/
    
    # plot raw signal at genome locations with the most significantly/consistently modified bases
    tombo plot most_significant --fast5-basedirs path/to/native/rna/fast5s/ \
        --statistics-filename sample.alt_testing.5mC.tombo.stats --plot-alternate-model 5mC
    
    # plot per-read test statistics using the 6mA alternative model testing method
    tombo plot per_read --per-read-statistics-filename sample.alt_testing.6mA.tombo.per_read_stats \
        --genome-locations chromosome:1000 chromosome:2000:- --genome-fasta genome.fasta

..

    For more plotting commands see the `Tombo plotting documentation here <https://nanoporetech.github.io/tombo/plotting.html>`_.

Read Filtering
^^^^^^^^^^^^^^

::

    # filter reads to a specific genomic location
    tombo filter genome_locations --fast5-basedirs path/to/native/rna/fast5s/ \
        --include-regions chr1:0-10000000

    # apply a more strigent raw signal matching threshold
    tombo filter  --fast5-basedirs path/to/native/rna/fast5s/ \
        --signal-matching-score 1.0

..

    For more read filtering commands see the `Tombo filter documentation here <https://nanoporetech.github.io/tombo/filtering.html>`_.

    Hint: Save a set of filters for later use by copying the Tombo index file: ``cp path/to/native/rna/.fast5s.RawGenomeCorrected_000.tombo.index save.native.tombo.index``. To re-set to a set of saved filters after applying further filters simply replace the index file: ``cp save.native.tombo.index path/to/native/rna/.fast5s.RawGenomeCorrected_000.tombo.index``.

====================
Note on Tombo Models
====================

Tombo is currently provided with two canonical models (for DNA and RNA data) and three alternative models (DNA::5mC, DNA::6mA and RNA::5mC).

These models are used by default in the re-squiggle and modified base detection commands. The correct canonical model is automatically selected for DNA or RNA based on the contents of each FAST5 file and processed accordingly.

Additional models will be added in future releases.

=========================
Installation Requirements
=========================

python Requirements (handled by conda or pip):
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

-  numpy
-  scipy
-  h5py
-  cython
-  mappy>=2.10
-  tqdm

Optional packages (handled by conda, but not pip):
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

-  Plotting Packages (R and rpy2 must be linked during installation)
   
   +  R
   +  rpy2
   +  ggplot2
   +  gridExtra (required for ``plot_motif_with_stats`` and ``plot_kmer`` subcommands)

-  On-disk Random Fasta Access
   
   +  pyfaidx

Advanced Installation Instructions
----------------------------------

Minimal tombo installation without optional dependencies (enables re-squiggle, all modified base testing methods and text output)

::

    pip install ont-tombo

Install current github version of tombo

::

    pip install git+https://github.com/nanoporetech/tombo.git

Download and install github version of tombo

::

    git clone https://github.com/nanoporetech/tombo.git
    cd tombo
    pip install -e .

    # to update, run:
    git pull
    pip install -I --no-deps -e .

========
Citation
========

Stoiber, M.H. et al. De novo Identification of DNA Modifications Enabled by Genome-Guided Nanopore Signal Processing. bioRxiv (2016).

http://biorxiv.org/content/early/2017/04/10/094672

============
Known Issues
============

-  When running the ``detect_modifications`` commands on large genomes, the computational memory usage can become very high. It is currently recommended to processes smaller regions using the ``tombo filter genome_locations`` command (with saved Tombo index hint above). This problem is being addressed and will be resolved in a later release.

-  The Tombo conda environment (especially with python 2.7) may have installation issues.
   
   + Tombo works best in python 3.4+, so many problems can be solved by upgrading python.
   + If installed using conda:

      - Ensure the most recent version of conda is installed (``conda update -n root conda``).
      - It is recommended to set conda channels as described for `bioconda <https://bioconda.github.io>`_.
      - Run ``conda update --all``.
   + In python 2.7 there is an issue with the conda scipy.stats package. Down-grading to version 0.17 fixes this issue.
   + In python 2.7 there is an issue with the conda h5py package. Down-grading to version <=2.7.0 fixes this issue.
