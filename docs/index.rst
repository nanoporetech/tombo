*********************************
Welcome to Tombo's documentation!
*********************************

Tombo is a suite of tools primarily for the identification of modified nucleotides from nanopore sequencing data.

Tombo also provides tools for the analysis and visualization of raw nanopore signal.

------------
Installation
------------

|bioconda_badge| |pypi_badge|

.. |bioconda_badge| image:: https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat-square
    :target: http://bioconda.github.io/recipes/ont-tombo/README.html

.. |pypi_badge| image:: https://badge.fury.io/py/ont-tombo.svg
    :target: https://badge.fury.io/py/ont-tombo

Basic tombo installation (python 2.7 and 3.4+ support)

::

    # install via bioconda environment
    conda install -c bioconda ont-tombo

    # or install pip package (numpy install required before tombo for cython optimization)
    pip install numpy
    pip install ont-tombo[full]

See :doc:`tutorials` for common workflows.

===========
Quick Start
===========

This **quick start** guides the steps to perform some common modified base detection analyses using the Tombo command line interface.

The first step in any Tombo analysis is to re-squiggle (raw signal to reference sequence alignment) raw nanopore reads. This creates an index and stores the information necessary to perform downstream analyses.

In this example, an E. coli sample is tested for dam and dcm methylation (present in lab E. coli; CpG model also available for human analysis). Using these results, raw signal is plotted at the most significantly modified dcm positions and the dam results are output to a `wiggle <https://genome.ucsc.edu/goldenpath/help/wiggle.html>`_ file for use in downstream processing or visualization in a genome browser.

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

   # produces "estimated fraction of modified reads" genome browser files
   tombo text_output browser_files --statistics-filename native.e_coli_sample.dam.tombo.stats \
       --file-types dampened_fraction --browser-file-basename native.e_coli_sample.dam
   # also produce successfully processed reads coverage file for reference
   tombo text_output browser_files --fast5-basedirs path/to/fast5s/ \
       --file-types coverage --browser-file-basename native.e_coli_sample

While motif models (``CpG``, ``dcm`` and ``dam``; most accurate) and all-context specific alternate base models (``5mC`` and ``6mA``; more accurate) are preferred, Tombo also allows users to investigate other or even unknown base modifications.

Here are two example commands running the ``de_novo`` method (detect deviations from expected cannonical base signal levels) and the ``level_sample_compare`` method (detect deviation in signal levels between two samples of interest; works best with high >50X coverage).

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

.. note::

   All Tombo commands work for direct RNA nanopore reads as well, but a transcriptome reference sequence must be provided for spliced transcripts.

   Run ``tombo -h`` to see all Tombo command groups, run ``tombo [command-group] -h`` to see all commands within each group and run ``tombo [command-group] [comand] -h`` for help with arguments to each Tombo command.

   Detailed documentation for all Tombo algorithms and commands can be found through the links here.

------
Naming
------

Tombo Ahi is a Japanese name for albacore (the name of an Oxford Nanopore Technologies basecaller). So use albacore to identify canonical bases and then use Tombo to detect more exotic, non-canonical bases.

--------
Contents
--------

.. toctree::
   :maxdepth: 2

   tutorials
   examples
   resquiggle
   modified_base_detection
   text_output
   plotting
   filtering
   rna
   model_training

-------------------
Tombo API Reference
-------------------

.. toctree::
   :maxdepth: 2

   tombo

:ref:`Tombo Module Documentation Index <genindex>`
