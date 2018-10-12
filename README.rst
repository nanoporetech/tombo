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

This quick start guides the steps to perform some common modified base detection analyses using the Tombo command line interface.

The first step in any Tombo analysis is to re-squiggle (raw signal to reference sequence alignment) raw nanopore reads. This creates an index and stores the information necessary to perform downstream analyses.

In this example, an E. coli sample is tested for dam and dcm methylation (CpG model also available for human analysis). Using these results, raw signal is plotted at the most significantly modified dcm positions and the dam modified base predictions are output to a `wiggle <https://genome.ucsc.edu/goldenpath/help/wiggle.html>`_ file for use in downstream processing or visualization in a genome browser.

::

   tombo resquiggle path/to/fast5s/ genome.fasta --processes 4 --num-most-common-errors 5
   tombo detect_modifications alternative_model --fast5-basedirs path/to/fast5s/ \
       --statistics-file-basename native.e_coli_sample \
       --alternate-bases dam dcm --processes 4

   # plot raw signal at most significant dcm locations
   tombo plot most_significant --fast5-basedirs path/to/fast5s/ \
       --statistics-filename native.e_coli_sample.dcm.tombo.stats \
       --plot-standard-model --plot-alternate-model dcm \
       --pdf-filename sample.most_significant_dcm_sites.pdf

   # produces wig file with estimated fraction of modified reads at each valid reference site
   tombo text_output browser_files --statistics-filename native.e_coli_sample.dam.tombo.stats \
       --file-types dampened_fraction --browser-file-basename native.e_coli_sample.dam
   # also produce successfully processed reads coverage file for reference
   tombo text_output browser_files --fast5-basedirs path/to/fast5s/ \
       --file-types coverage --browser-file-basename native.e_coli_sample

While motif models (``CpG``, ``dcm`` and ``dam``; most accurate) and all-context specific alternate base models (``5mC`` and ``6mA``; more accurate) are preferred, Tombo also allows users to investigate other or even unknown base modifications.

Here are two example commands running the ``de_novo`` method (detect deviations from expected cannonical signal levels) and the ``level_sample_compare`` method (detect deviation in signal levels between two samples of interest; works best with high coverage).

::

   tombo detect_modifications de_novo --fast5-basedirs path/to/fast5s/ \
       --statistics-file-basename sample.de_novo_detect --processes 4
   tombo text_output browser_files --statistics-filename sample.de_novo_detect.tombo.stats \
       --browser-file-basename sample.de_novo_detect --file-types dampened_fraction

   tombo detect_modifications level_sample_compare --fast5-basedirs path/to/fast5s/ \
       --control-fast5-basedirs path/to/control/fast5s/ --minimum-test-reads 50 \
       --processes 4 --statistics-file-basename sample.level_samp_comp_detect
   tombo text_output browser_files --statistics-filename sample.level_samp_comp_detect.tombo.stats \
       --browser-file-basename sample.level_samp_comp_detect --file-types statistic

..

   See more complete tutorials on the `documentation page <https://nanoporetech.github.io/tombo/tutorials.html>`_.

===
RNA
===

All Tombo commands work for direct RNA nanopore reads as well, but a transcriptome reference sequence must be provided for spliced transcripts.

Tips for processing direct RNA reads within the Tombo framework can be found in the `RNA section <https://nanoporetech.github.io/tombo/rna.html>`_ of the detailed Tombo documentation.

=====================
Further Documentation
=====================

Run ``tombo -h`` to see all Tombo command groups and run ``tombo [command-group] -h`` to see all commands within each group.

Detailed documentation for all Tombo commands and algorithms can be found on the `tombo documentation page <https://nanoporetech.github.io/tombo/>`_.

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
