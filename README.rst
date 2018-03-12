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

..

    Additional installation instructions options below

=============
Documentation
=============

Run ``tombo -h`` to see all Tombo sub-commands and run ``tombo [sub-command] -h`` to see the options for any Tombo sub-command.

Detailed documentation for all Tombo algorithms and commands can be found at https://nanoporetech.github.io/tombo/

==============
Tombo Examples
==============

Re-squiggle (Raw Data to Genome Alignment)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

::

    tombo resquiggle path/to/amplified/dna/fast5s/ genome.fasta --processes 4

..

    Only R9.4/5 data (including R9.[4/5].1) is supported at this time.

    DNA or RNA is automatically determined from FAST5s (set explicitly with ``--dna`` or ``--rna``).

    FAST5 files need not contain Events data, but must contain Fastq slot. See ``annotate_raw_with_fastqs`` for pre-processing of raw FAST5s.

Identify Modified Bases
^^^^^^^^^^^^^^^^^^^^^^^

::

    # comparing to an alternative 5mC and 6mA model (recommended method)
    tombo test_significance --fast5-basedirs path/to/native/dna/fast5s/ \
        --alternate-bases 5mC 6mA --statistics-file-basename sample

    # comparing to a control sample (e.g. PCR)
    tombo test_significance --fast5-basedirs path/to/native/dna/fast5s/ \
        --control-fast5-basedirs path/to/amplified/dna/fast5s/ --statistics-file-basename sample_compare

    # compare to the canonical base model
    tombo test_significance --fast5-basedirs path/to/native/dna/fast5s/ \
        --statistics-file-basename sample_de_novo --processes 4

..

    Must run ``resquiggle`` on reads before testing for modified bases.

    ``test_significance`` produces a binary file. See ``write_wiggles`` or ``plot_most_significant`` Tombo sub-commands for text output or genome region selection.

Text Output (Wiggle file format)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

::

    # extract fraction of reads modified at each genomic base in wiggle file format
    tombo write_wiggles --wiggle-types fraction --statistics-filename sample.5mC.tombo.stats

    # extract read depth from mapped and re-squiggled reads
    tombo write_wiggles --wiggle-types coverage --fast5-basedirs path/to/native/dna/fast5s/

Extract Sequences Surrounding Modified Positions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

::

    tombo write_most_significant_fasta --statistics-filename sample.6mA.tombo.stats \
        --genome-fasta genome.fasta

Plotting Examples
^^^^^^^^^^^^^^^^^

::

    # plot raw signal with standard model overlay at reions with maximal coverage
    tombo plot_max_coverage --fast5-basedirs path/to/native/rna/fast5s/ --plot-standard-model
    
    # plot raw signal along with signal from a control (PCR) sample at locations with the AWC motif
    tombo plot_motif_centered --fast5-basedirs path/to/native/rna/fast5s/ \
        --motif AWC --genome-fasta genome.fasta --control-fast5-basedirs path/to/amplified/dna/fast5s/
    
    # plot raw signal at genome locations with the most significantly/consistently modified bases
    tombo plot_most_significant --fast5-basedirs path/to/native/rna/fast5s/ \
        --statistics-filename sample.5mC.tombo.stats --plot-alternate-model 5mC
    
    # plot per-read test statistics using the 6mA alternative model testing method
    tombo plot_per_read --fast5-basedirs path/to/native/rna/fast5s/ \
        --genome-locations chromosome:1000 chromosome:2000:- --plot-alternate-model 6mA

===============
Common Commands
===============

Re-squiggle (Raw Data Alignment):
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

..

    Must be run before any other commands.

::

   resquiggle                    Re-annotate raw signal with genomic alignment from existing basecalls.

Modified Base Detection:
^^^^^^^^^^^^^^^^^^^^^^^^

::

   test_significance             Test for shifts in signal indicative of non-canonical bases.

Text Output Commands:
^^^^^^^^^^^^^^^^^^^^^

::

   write_wiggles                 Write text outputs for genome browser visualization and bioinformatic processing (wiggle file format).
   write_most_significant_fasta  Write sequence centered on most modified genomic locations.

Genome Anchored Plotting Commands:
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

::

   plot_max_coverage             Plot raw signal in regions with maximal coverage.
   plot_genome_location          Plot raw signal at defined genomic locations.
   plot_motif_centered           Plot raw signal at a specific motif.
   plot_max_difference           Plot raw signal where signal differs most between two read groups.
   plot_most_significant         Plot raw signal at most modified locations.
   plot_motif_with_stats         Plot example signal and statistic distributions around a motif of interst.
   plot_per_read                 Plot per read modified base probabilities.

Read Filtering:
^^^^^^^^^^^^^^^

::

   clear_filters                 Clear filters to process all successfully re-squiggled reads.
   filter_stuck                  Apply filter based on observations per base thresholds.
   filter_coverage               Apply filter to downsample for more even coverage.

..

    Get additional help for subcommands with ``tombo [command] -h``

====================
Note on Tombo Models
====================

Tombo is currently provided with two standard models (DNA and RNA) and two alternative models (DNA::5mC, DNA::6mA). These models are applicable only to R9.4/5 flowcells with 1D or 1D^2 kits (not 2D).

These models are used by default for the re-squiggle and testing commands. The correct model is automatically selected for DNA or RNA based on the contents of each FAST5 file and processed accordingly. Additional models will be added in future releases.

============
Requirements
============

python Requirements (handled by conda or pip):
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

-  numpy
-  scipy
-  h5py
-  cython
-  mappy

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

Install github version of tombo (versions on conda/pypi should be up-to-date)

::

    pip install git+https://github.com/nanoporetech/tombo.git

========
Citation
========

Stoiber, M.H. et al. De novo Identification of DNA Modifications Enabled by Genome-Guided Nanopore Signal Processing. bioRxiv (2016).

http://biorxiv.org/content/early/2017/04/10/094672

=======
Gotchas
=======

-  The Tombo conda environment (especially with python 2.7) may have installation issues.
   
   + Tombo works best in python 3.4+, so many problems can be solved by upgrading python.
   + If installed using conda:

      - Ensure the most recent version of conda is installed (``conda update conda``).
      - It is recommended to set conda channels as described for `bioconda <https://bioconda.github.io>`_.
      - Run ``conda update --all``.
   + In python 2.7 there is an issue with the conda scipy.stats package. Down-grading to version 0.17 fixes this issue.
