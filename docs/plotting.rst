*****************
Plotting Commands
*****************

In order to enhance modified base detection and give users a better understanding of raw nanopore data, Tombo provides a number of plotting commands.

------------------------
Genome Anchored Plotting
------------------------

Plot Region Selection
^^^^^^^^^^^^^^^^^^^^^

Most Tombo plotting functions are genome anchored. These commands create plots analogous to a genome browser, but with all raw signal within a region. The available commands each differ in their mode of genome region selection. This allows users to plot regions of interest for many research contexts.

* ``plot_max_coverage`` - Select regions with maximal coverage
* ``plot_genome_location`` - Select specified genomic locations
* ``plot_motif_centered`` - Select regions with a specific motif (follows `NEB single letter codes <https://www.neb.com/tools-and-resources/usage-guidelines/single-letter-codes>`_)
* ``plot_max_difference`` - Select regions where two samples' average signal differs most
* ``plot_most_significant`` - Select most consistently/significantly mofidied locations

These plotting commands produce raw signal level plots such at the example below. Options are available for each of these plots to logically select genomic regions based on the given criterion.

----

.. figure::  _images/single_samp.png
   :align: center
   :scale: 30%
   
   Single sample raw signal plot

----

Model Plotting
^^^^^^^^^^^^^^

Plots are also enabled to visualize the different testing frameworks available in Tombo. Thses plots would additionally include a control sample, the standard model or any non-standard base model, visualizing the control sample comparison, de novo testing and log likelihood ratio tests respectively.

Control these plots with these options: ``--control-fast5-basedirs``, ``--tombo-model-filename``, ``--alternate-model-filename``, ``--plot-standard-model``, and ``--plot-alternate-model``

----

.. figure::  _images/sample_comp.png
   :align: center
   :scale: 30%
   
   Control sample comparison plot

.. figure::  _images/model_comp.png
   :align: center
   :scale: 30%
   
   Standard model plot

.. figure::  _images/alt_model_comp.png
   :align: center
   :scale: 30%
   
   Alternate model plot

----

Over-Plotting
^^^^^^^^^^^^^

When high coverage regions are plotted the raw signal plots can become less interpretable. By default, when read coverage exceeds 50X reads are randomly downsampled (change this option with ``--overplot-threshold``). Three additional over-plotting options (boxplot, quantile and density) are available as shown below (chose which over-plotting type to use with the ``--overplot-type`` option).

----

.. figure::  _images/boxplot.png
   :align: center
   :scale: 30%
   
   Boxplot over-plotting option

.. figure::  _images/quantile.png
   :align: center
   :scale: 30%
   
   Quantile over-plotting option

.. figure::  _images/density.png
   :align: center
   :scale: 30%
   
   Density over-plotting option

----

Per-read Plotting
^^^^^^^^^^^^^^^^^

As testing is applied on a per-read setting, per-read statistics plots are also available. As per-read statistics are not stored in the current Tombo framework, test values are re-computed for this plotting command (and the control sample comparison method is not currently enabled). Create these plots with the ``plot_per_read`` command.

----

.. figure::  _images/pre_read_5mC.png
   :align: center
   :scale: 30%
   
   Alternative 5mC model testing

.. figure::  _images/per_read_do_novo.png
   :align: center
   :scale: 30%
   
   De novo, standard model, per-read testing

----

Motif-centered Statistic Plotting
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In several biological contexts base modifications occur at specific motifs. In order to visualize the distribution of Tombo statistical test results centered on a motif of biolgical interest (or a discovered motif) the ``plot_motif_with_stats`` command is provided.

This command identifies a number (defined by ``--num-statistics``) of genomic regions centered on this motif with the highest significance testing values. Importantly, the identified highest testing values need not be found within the actual motif, but simply within a region containing the motif defined by ``--num-context``. In this way, non-interesting motifs (motifs which don't direct modifications) will not contain more significant statistics centered on a specific position within the provided motif. A number (defined by ``--num-regions``) of example regions with the highest test statistics centered on the motif of interest are plotted as well.

----

.. figure::  _images/stat_dist.png
   :align: center
   :scale: 30%
   
   Example statistics distribution around `biologically relevant CCWGG motif in E. coli <https://www.neb.com/tools-and-resources/usage-guidelines/dam-and-dcm-methylases-of-e-coli>`_

----

-----------------------
Other Plotting Commands
-----------------------

K-mer Level Distributions
^^^^^^^^^^^^^^^^^^^^^^^^^

In order to investigate the k-mer signal current levels of a particular set of reads, the ``plot_kmer`` command is provided.

----

.. figure::  _images/kmer_levels.png
   :align: center
   :scale: 30%
   
   Example k-mer current level distribution plot

----

Correction Plotting
^^^^^^^^^^^^^^^^^^^

Plotting commands, ``plot_correction`` and ``plot_multi_correction``, are provided to visualize the old event-based re-squiggle process. These commands are thus only applocable on reads that have been processed with ``event_reqsuiggle``.
