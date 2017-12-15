=======
Summary
=======

Tombo is a suite of tools primarily for the identification of modified nucleotides from nanopore sequencing data.

Tombo also provides tools for the analysis and visualization of raw nanopore signal.

============
Installation
============

Basic tombo installation (python2.7 support only)

::

    # install via bioconda environment
    conda install -c bioconda ont-tombo

    # or install pip package (numpy install required before tombo for cython optimization)
    pip install numpy
    pip install ont-tombo

..

    Additional installation instructions options below

==================
Full Documentation
==================

Detailed documentation can be found at https://nanoporetech.github.io/tombo/

==============
Tombo Examples
==============

Re-squiggle (Raw Data Alignment)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

::

    tombo resquiggle path/to/amplified/dna/fast5s/ genome.fasta --minimap2-executable ./minimap2 --processes 4

..

    FAST5 files need not contain Events data, but must contain Fastq slot.

    Only R9.4/5 data (DNA or RNA) is supported at this time. Processing of other samples may produce sub-optimal results.

Identify Modified Bases
^^^^^^^^^^^^^^^^^^^^^^^

::

    # comparing to an alternative 5mC model
    tombo test_significance --fast5-basedirs path/to/native/dna/fast5s/ \
        --alternate-bases 5mC --statistics-file-basename sample_compare

    # comparing to a control sample (e.g. PCR)
    tombo test_significance --fast5-basedirs path/to/native/dna/fast5s/ \
        --control-fast5-basedirs path/to/amplified/dna/fast5s/  --statistics-file-basename sample_compare

    # compare to the canonical base model
    tombo test_significance --fast5-basedirs path/to/native/dna/fast5s/ \
        --statistics-file-basename sample --processes 4

..

    Must run ``resquiggle`` on reads before testing for modified bases.
   
    ``test_significance`` produces a binary file. See ``write_wiggles`` for several text outputs or ``plot_most_significant`` to use for genome region selection.

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

    tombo write_most_significant_fasta --statistics-filename sample_compare.5mC.tombo.stats \
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
        --statistics-filename sample_compare.5mC.tombo.stats --plot-alternate-model 5mC
    
    # plot per-read test statistics using the 5mC alternative model testing method
    tombo plot_per_read --fast5-basedirs path/to/native/rna/fast5s/ \
        --genome-locations chromosome:1000 chromosome:2000:- --plot-alternate-model 5mC

===============
Common Commands
===============

::

   # get tombo help
   tombo -h
   # run tombo sub-commands
   tombo [command] [options]

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

Tombo is currently provided with two standard models (DNA and RNA) and one alternative model (DNA::5mC). These models are applicable only to R9.4/5 flowcells with 1D or 1D^2 kits (not 2D).

These models are used by default for the re-squiggle and testing commands. The correct model is automatically selected for DNA or RNA based on the contents of each FAST5 file and processed accordingly. Additional models will be added in future releases.

============
Requirements
============

At least one supported mapper:
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

-  minimap2 (https://github.com/lh3/minimap2)
-  BWA-MEM (http://bio-bwa.sourceforge.net/)
-  graphmap (https://github.com/isovic/graphmap)

-  HDF5 (http://micro.stanford.edu/wiki/Install_HDF5#Install)

python Requirements (handled by pip):
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

-  numpy (must be installed before installing tombo)
-  scipy
-  h5py
-  cython

Optional packages for plotting (install R packages with ``install.packages([package_name])`` from an R prompt):
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

-  rpy2 (along with an R installation)
-  ggplot2 (required for any plotting subcommands)
-  cowplot (required for plot_motif_with_stats subcommand)

Optional packages for alternative model estimation:
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

-  sklearn

Advanced Installation Instructions
----------------------------------

Install tombo with all optional dependencies (for plotting and model estimation)

::

    pip install ont-tombo[full]

Install tombo with plotting dependencies (requires separate installation
of R packages ggplot2 and cowplot)

::

    pip install ont-tombo[plot]

Install tombo with alternative model estimation dependencies

::

    pip install ont-tombo[alt_est]

Install github version of tombo (most versions on pypi should be up-to-date)

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

-  If plotting commands fail referencing rpy2 images, shared object files, etc., this may be an issue with the version of libraries installed by conda. In order to resolve this issue, remove the conda-forge channel and re-install ont-tombo.
