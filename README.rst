.. image:: /ONT_logo.png
  :width: 800

******************

Tombo
"""""

Tombo is a suite of tools primarily for the identification of modified nucleotides from nanopore sequencing data. Tombo also provides tools for the analysis and visualization of raw nanopore signal.

Detailed documentation for all Tombo commands and algorithms can be found on the `tombo documentation page <https://nanoporetech.github.io/tombo/>`_.

Features
--------

- Modified Base Detection

  - Supports both DNA and direct RNA

    - `RNA processing details <https://nanoporetech.github.io/tombo/rna.html>`_
  - Three detection algorithms support broad range of applications

    - Alternative model (preferred)
    - Sample comparison
    - De novo
- Reference-anchored raw signal vizualization
- Raw signal analysis python API
- User-friendly model estimation methods with tutorial

*********************

Getting Started
"""""""""""""""

|bioconda_badge| |pypi_badge|

.. |bioconda_badge| image:: https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat-square
    :target: http://bioconda.github.io/recipes/ont-tombo/README.html

.. |pypi_badge| image:: https://badge.fury.io/py/ont-tombo.svg
    :target: https://pypi.org/project/ont-tombo/

Conda installation (preferred method)

::

    # install via bioconda environment (https://bioconda.github.io/#set-up-channels)
    conda install -c bioconda ont-tombo

The first step in any Tombo analysis is to re-squiggle (raw signal to reference sequence alignment) raw nanopore reads. This creates an index and stores the raw signal alignments necessary to perform downstream analyses.

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

Alternative Installation Methods
--------------------------------

Tombo is available for installation via pip, but requires an R installation as well as R package dependencies (ggplot2 and gridextra) for all visualization functions.

::

   # install pip package (numpy install required before tombo for cython optimization)
   pip install numpy
   pip install ont-tombo[full]

Tombo can also be installed directly from source (mostly for development) by running the following commands:

::

   git clone https://github.com/nanoporetech/tombo
   cd tombo
   pip install -e .

********

Known Issues
""""""""""""

Tombo does not support multi-read FAST5 format read data files. Please use the ``multi_to_single_fast5`` command from the `ont_fast5_api package <https://github.com/nanoporetech/ont_fast5_api>`_ in order to convert to single-read FAST5 format before processing with Tombo.

Help
""""

Licence and Copyright
---------------------

© 2017-18 Oxford Nanopore Technologies Ltd.

Tombo is distributed under the terms of the included MPL2 licence.

References and Supporting Information
-------------------------------------

Stoiber, M.H. et al. De novo Identification of DNA Modifications Enabled by Genome-Guided Nanopore Signal Processing. bioRxiv (2016).

http://biorxiv.org/content/early/2017/04/10/094672
