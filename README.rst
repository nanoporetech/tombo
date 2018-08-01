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

Re-squiggle raw nanopore read files and call 5mC and 6mA sites.

Then, for 5mA calls, output genome browser `wiggle format file <https://genome.ucsc.edu/goldenpath/help/wiggle.html>`_ and, for 6mA calls, plot raw signal around most significant locations.

::

   # skip this step if FAST5 files already contain basecalls
   tombo preprocess annotate_raw_with_fastqs --fast5-basedir path/to/fast5s/ \
       --fastq-filenames basecalls1.fastq basecalls2.fastq \
       --sequencing-summary-filenames seq_summary1.txt seq_summary2.txt \
       --processes 4

   tombo resquiggle path/to/fast5s/ genome.fasta --processes 4 --num-most-common-errors 5
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
       --fast5-basedirs path/to/fast5s/  --file-types dampened_fraction coverage \
       --browser-file-basename sample.alt_modified_base_detection.6mA

   # plot raw signal at most significant 6mA locations
   tombo plot most_significant --fast5-basedirs path/to/fast5s/ \
       --statistics-filename sample.alt_modified_base_detection.6mA.tombo.stats \
       --plot-standard-model --plot-alternate-model 6mA \
       --pdf-filename sample.most_significant_6mA_sites.pdf

Detect any deviations from expected signal levels for canonical bases to investigate any type of modification.

::

   tombo resquiggle path/to/fast5s/ genome.fasta --processes 4 --num-most-common-errors 5
   tombo detect_modifications de_novo --fast5-basedirs path/to/fast5s/ \
       --statistics-file-basename sample.de_novo_modified_base_detection \
       --per-read-statistics-basename sample.de_novo_modified_base_detection \
       --processes 4

   # produces "estimated fraction of modified reads" genome browser files from de novo testing
   tombo text_output browser_files --statistics-filename sample.de_novo_modified_base_detection.tombo.stats \
       --browser-file-basename sample.de_novo_modified_base_detection --file-types dampened_fraction

===
RNA
===

All Tombo commands work for RNA data as well, but a transcriptome reference sequence must be provided for spliced transcripts.

The reasons for this decision and other tips for processing RNA data within the Tombo framework can be found in the `RNA section <https://nanoporetech.github.io/tombo/rna.html>`_ of the detailed Tombo documentation.

=====================
Further Documentation
=====================

Run ``tombo -h`` to see all Tombo command groups and run ``tombo [command-group] -h`` to see all commands within each group.

Detailed documentation for all Tombo commands and algorithms can be found at https://nanoporetech.github.io/tombo/

========
Citation
========

Stoiber, M.H. et al. De novo Identification of DNA Modifications Enabled by Genome-Guided Nanopore Signal Processing. bioRxiv (2016).

http://biorxiv.org/content/early/2017/04/10/094672

============
Known Issues
============

-  The Tombo conda environment (especially with python 2.7) may have installation issues.

   + Tombo works best in python 3.4+, so many problems can be solved by upgrading python.
   + If installed using conda:

      - Ensure the most recent version of conda is installed (``conda update -n root conda``).
      - It is recommended to set conda channels as described for `bioconda <https://bioconda.github.io/#set-up-channels>`_.
      - Run ``conda update --all``.
   + In python 2.7 there is an issue with the conda scipy.stats package. Down-grading to version 0.17 fixes this issue.
   + In python 2.7 there is an issue with the conda h5py package. Down-grading to version <=2.7.0 fixes this issue.
