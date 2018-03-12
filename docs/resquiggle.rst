*********************
Re-squiggle Algorithm
*********************

The signal level data produced from a nanopore read is referred to as a squiggle. Base calling this squiggle information generally contains some errors compared to a refernce genome. The re-squiggle algorithm defines a new squiggle to genomic sequence assignment, hence a re-squiggle.

The re-squiggle algorithm is the basis for the Tombo framework. The re-squiggle algorithm takes as input a read file (in FAST5 format) containing raw signal and associated base calls. The base calls are mapped to a genome reference and then the raw signal is assigned to the genomic sequence based on an expected current level model.

**TL;DR**:

*  Re-squiggle must be run before any other Tombo command (aside from the ``annotate_raw_with_fastqs`` pre-processing sub-command).
*  Minimally the command takes a directory containing FAST5 files and a genome/transcriptome reference.
   
   - Genome reference may be previously known or discovered from this sample.

*  FAST5 files must contain basecalls (as produced by albacore in fast5 mode or added with ``annotate_raw_with_fastqs``), but need not contain the "Events" table.
*  Tombo currently only supports R9.4 and R9.5 data (via included default models). R9.4.1 and R9.5.1 are supported. Other data may produce sub-optimal results.
*  DNA and RNA reads will be detected automatically and processed accordingly (set explicitly with ``--dna`` or ``--rna``).
   
   -  Tombo does not perform spliced mapping. Thus a transcriptime reference must be passed to the re-squiggle command for RNA samples. For futher details on Tombo RNA processing see the :doc:`rna` section.

*  Run ``resquiggle`` over multiple cores with the ``--processes`` option.

   - The ``--threads-per-process`` is also provided, but it is generally recommended that this option remains set to the default of 1, though it may improve results on some computing environments.

-----------------
Algorithm Details
-----------------

The re-squiggle algorithm occurs in four main steps described below.

* Genome Mapping
* Event Detection
* Sequence to Signal Assignment
* Resolve Skipped Bases

Genome Mapping
--------------

The genome mapping is performed via the python API to ``minimap2`` (`mappy python package <https://pypi.python.org/pypi/mappy>`_).

Read base called sequence location within the FAST5 file is defined by the ``--basecall-group`` and ``--basecall-subgroups`` command line options. The default values of these parameters point to the default location for base calls from albacore or ``annotate_raw_with_fastqs``.

The genomic sequence for successfully mapped reads are then passed on to the sequence to signal assignment stage.

.. tip::
   
   Unless the optional dependency ``pyfaidx`` is installed (included in default conda installation), each process reads the whole reference genome into memory in order to extract genomic seqeunce. Take care when running Tombo on larger genomes to avoid overflowing a systems memory. This is true even if the optional ``--minimap2-index`` parameter is provided. The minimap2 index parameter only effects the mapping call itself.

Event Detection
---------------

The Tombo algorithm does not require the "Events" table (raw signal assignment to base calls). Instead, Tombo discovers events from the raw signal. This segmented signal makes downstream processing steps more efficient and stable. This event detection algorithm is different from the event detection performed in previous versions of albacore, but produces similar results.

Events are determined by identifying large shifts in current level, by taking the running difference between neighboring windows of raw signal (explicitly set this parameter with the ``--segmentation-parameters`` option). The largest jumps (or most significant via a t-test for RNA) are chosen as the breakpoints between events. The mean of normalized raw signal is then computed for each event.

Raw signal normalization estimates a median shift parameter and a median absolute deviation (MAD) scale parameter. By default, a global scale value is taken as the mean of MAD computed from a random sample of reads and used to scale all reads. This behaviour can be overriden with ``--fit-scale-per-read`` option or the ``--fixed-scale`` option to manually set the global scaling value (advanced users only). Raw signal is also windsorized, ``--outlier-threshold`` parameter. These scaling parameters are stored in the Tombo FAST5 slot for access in later commands. Note that only median signal normalization is available within the Tombo framework.

The ``--segmentation-parameters`` values have been optimized for DNA and RNA data types, so DNA and RNA read types should not be mixed in processing.

Sequence to Signal Assignment
-----------------------------

Given the mapped genomic sequence and segmented signal, the sequence to signal assignment algorithm finds the most likely matching between these two.

The algorithm first uses a large bandwidth (5000 events over the first 500 genomic bps) to identify the start of the genomic sequence within the events (see figure below). This is necessary as some portion at the beginning of a read is not base called and some additional sequence may have been trimmed from the alignment. The matching is determined by applying a dynamic programming/dynamic time warpping algorithm to find the most likely matching between the event signal levels and the expected signal levels given the genomic sequence.

----

.. figure::  _images/begin_half_z_scores.png
   :align: center
   :scale: 110%
   
   Read start shifted half-normal scores

.. figure::  _images/begin_forward_pass.png
   :align: center
   :scale: 110%
   
   Read start forward pass scores and traceback path

----

A static banded matrix is constructed by computing the z-score for event level (x-axis) against genomic positions (y-axis). The negative absolute value z-score is shifted to an expected value of zero to fill the banded matrix (see figure **a** above). A forward pass computes the maximal cummulative score up to each matched event to genome position (see figure **b** above).

At each iteration the maximal score is taken over three possibilities 1) staying in the same genomic position, and accumulating the shifted z-score 2) matching an event with a genomic position (with score bonus) 3) skipping this genomic position (with a score penalty). The score match and skip penalties are definied by the ``--signal-align-parameters``. The default values have been optimized for DNA and RNA data types. From this forward pass, the maximal score along the last genomic position is taken and traced back to obtain the starting position of matching sequence and signal.

If a read is short enough (less than 5500 events or less than 500 bps of called sequence), then the whole sequence to signal matching will be performed with a single run with an appropriate static bandwidth.

For longer reads, the above computed start matching position is taken and then the same dynamic programming solution is applied except a smaller adaptive band is now used (see figure below). The bandwidth is definied by the ``--signal-align-parameters`` option and again has been optimized for DNA and RNA data types. At each genomic position, the band position is defined to center on the maximal score of the forward pass from the previous base. This aims to ensure that the traceback path will remain within the adaptive window. There are edge cases where the valid matching leaves the adaptive band. These reads are filtered out and included in the failed read group ``Read event to sequence alignment extends beyond --bandwidth``.

----

.. figure::  _images/adaptive_half_z_scores.png
   :align: center
   :scale: 80%
   
   Full read adaptive banded shifted half-normal scores

.. figure::  _images/adaptive_forward_pass.png
   :align: center
   :scale: 80%
   
   Full read adaptive banded forward pass scores

----

Resolve Skipped Bases
---------------------

After the dynamic programming step, skipped bases must be resolved using the raw signal to obtain a matching of each genomic base to a bit of raw signal. A window around each skipped genomic base is identified. If a window does not contain enough raw signal to perform a raw signal search the window is expanded until enough signal is found. Overlapping windows are collapsed into a single window.

After deletion windows are identified, a dynamic programming algorithm very similar to the last step is performed. Importantly, the raw signal is used instead of events and the skip move is no longer allowed. Additionally, each genomic base is forced to contain a minimal number of raw observations to produce more robust assignments (explicitly set this value with the ``--segmentation-parameters`` option). This completes the re-squiggle procedure producing a matching of a read's raw signal to the mapped genomic sequence.

-------------------------------
Common Failed Read Descriptions
-------------------------------

``Fastq slot not present in --basecall-group``
``Raw data is not found in Raw/Reads/Read_[read#]``

*  These error indicates that a necessary peice of information for Tombo to run was not found in the FAST5 file.

``Alignment not produced``

*  This error indicates that minimap2 (via mappy API) did not produce a valid mapping.

``Could not close FAST5 file. Possibly corrupted file``

*  This error indicates that an unexpected error occurred trying to open or close a read file. This can happen if the reads are being actively accessed by another program or if a file has been corrupted.

``Read contains too many potential genomic deletions``
``Not enough raw signal around potential genomic deletion(s)``

*  These errors indicate that the sequence to signal matching algorithm was unable to identify a valid path.

``Poor raw to expected signal matching``
``Poor raw to expected signal matching at read start``

*  These errors indicate that the dynamic programming algorithm produce a poorly scored matching of genomic sequence to raw signal. Some potential sources for these errors include incorrect primary genomic mapping, incorrect genome sequence (compared to the biological sample), poor quality raw signal or an incompatible flowcell/library with included canonical models (only R9.5/4 flowcells currently supported; 2D reads are not supported; DNA and RNA are supported).

------------------
Tombo FAST5 Format
------------------

The result of the re-squiggle algorithm writes the sequence to signal assignment back into the read FAST5 files (found in the ``--corrected-group`` slot; the default value is the default for all other Tombo commands to read in this data). When running the re-squiggle algorithm a second time on a set of reads, the ``--overwrite`` option is required in order to write over the previous Tombo results.

The ``--corrected-group`` slot contains attributes for the signal normalization (shift, scale, upper_limit, lower_limit and outlier_threshold) as well as a boolean flag indicating whether the read is DNA or RNA. Within the ``Alignment`` group, the gemomic mapped start, end, strand and chromosome as well as mapping statistics (number clipped start and end bases, matching, mismatching, inserted and deleted bases) are stored.

The ``Events`` slot contains a matrix with the matching of raw signal to genomic sequence. This slot contains a single attribute (``read_start_rel_to_raw``) giving the zero-based offset indicating the beginning of the read genomic sequence within the raw signal. Each entry in the ``Events`` table indicates the normalized mean signal level (``norm_mean``), optionally (triggered by the ``--include-event-stdev`` option) the normalized signal standard deviation (``norm_stdev``), the start position of this base (``start``), the length of this event in raw signal values (``length``) and the genomic base (``base``).

This information is accessed as needed for down-stream Tombo processing commands.

This data generally adds ~75% to the memory footprint of a minimal FAST5 file (containing raw and sequence data; not including a basecalling Events table). This may vary across files and sample types.

**Important RNA note**: Tombo performs only un-spliced mapping. As such, for potentially spliced transcripts a transcriptome file must be provided. While this makes Tombo RNA processing annotation dependent the transcriptome is the more appropriate setting for modified base detection and thus this path has been chosen for Tomob RNA processing. More details about RNA processing can be found in the :doc:`rna` section.

**Minor RNA note**: RNA reads pass through the pore in the 3' to 5' direction during sequencing. As such, the raw signal and albacore events are stored in the reverse direction from the genome. Tombo events for RNA data are stored in the opposite direction (corresponding to the genome sequence direction, not sequencing time direction) for several practical reasons. Thus if events are to be compared to the raw signal, the raw signal must be reversed. Tombo RNA models are stored in the same direction and thus may be considered inverted as compared to some other RNA HMM signal level models.

----------------
Tombo Index File
----------------

By default, Tombo will create a hidden file containing the essential genome mapping location for each validly processed read. This file will be located alongside the base fast5s directory. For example, if resquiggle is called on this directory ``/path/to/fast5/reads/`` then the following file will be created ``/path/to/fast5/.reads.RawGenomeCorrected_000.tombo.index``. This file should generally be about 1,000 times smaller than the corresponding FAST5s and so should not incur significantly more disk usage. If desired, the index can be skipped with the ``--skip-index`` option. Note that this will make most all downstream commands much slower and filtering cannot be applited without this index file.

-------------------------------
Additional Command Line Options
-------------------------------

``--failed-reads-filename``

*  This option outputs the filenames for each read that failed via each failure mode. This can be useful for tracking down bothersome errors.

``--obs-per-base-filter``

*  This option applies a filter to "stuck" reads (too many observations per genomic base). This filter is applied only to the Tombo index file and can be cleared later. See the :doc:`filtering` section for more details.

``--ignore-read-locks``

*  Multiple independent ``resquiggle`` commands on the same set of reads should NOT be run simultaneously. This can cause hard to track errors and read file corruption. To protect against this, Tombo adds a lock (only acknowledged by Tombo) to each directory being processed. If a previous ``resquiggle`` command fails in a very unexpected fashion these locks can be left in place. In this case the ``--ignore-read-locks`` option is provided. This is the only intended use for this option.

---------------------
Pre-process Raw Reads
---------------------

Nanopore raw signal-space data consumes more memory than sequence-space data. As such, many users will produce only FASTQ basecalls initially from a set of raw reads in FAST5 format. The Tombo framework requires the linking of these basecalls with the raw signal-space data. The ``annotate_raw_with_fastqs`` sub-command is provided to assist with this workflow.

Given a directory (or nested directories) of FAST5 raw read files and a set of FASTQ format basecalls from these reads, the ``annotate_raw_with_fastqs`` adds the sequence information from the FASTQ files to the appropriate FAST5 files. This command generally adds 15-25% to the disk usage for the raw reads.

This functionality requires that the FASTQ seqeunce header lines begin with the read identifier from the FAST5 file. This has been tested with the Oxford Nanopore Technologies supported basecaller, albacore. Third-party basecallers may be not be processed correctly.

.. code-block:: bash

    tombo annotate_raw_with_fastqs --fast5-basedir <fast5s-base-directory> --fastq-filenames reads.fastq

---------------------------
Non-standard Data Locations
---------------------------

In the Tombo framework, it is possible to access and store basecalls and genome-anchored re-squiggle results in custom locations within FAST5 files.

For example, basecalls can be found in the ``Basecall_1D_001`` slot in a set of reads that have been basecalled more than one time. In this case the basecalls can be accessed in Tombo by specifying the ``--basecall-group`` option to the ``resquiggle`` command.

It can also be adventageous to store re-squiggle results in a non-standard locations. For example, if one would like to test multiple sets of re-squiggle parameters or reference versions without having to overwrite old results and re-run the ``resquiggle`` command, the ``--corrected-group`` option can be specified. This will store the re-squiggle results in a new slot within the FAST5 file as well as creating a new Tombo index file.

.. important::
   
   If the ``--corrected-group`` is specified in the ``resquiggle`` command, this same value must be passed to all other Tombo sub-commands in order to access these results. This inlcudes all filtering, plotting, significance testing, and text output commands.
