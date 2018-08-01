**************
RNA Processing
**************

.. important::

   Tombo cannot currently process spliced alignments. Thus processing RNA data requires that a transcriptome (NOT genome) reference provided for organisms with spliced transcription products.

Processing RNA data within the Tombo framework requires some extra care. The major item to consider when performing RNA processing is that a transcriptome reference must be supplied as spliced mapping is not supported. The lack of spliced mapping support within the Tombo framework is a conscious decision for identification of modified RNA bases. This is because the transcriptome is the natural setting for the detection of modified RNA bases. When modified RNA bases are projected onto the genome reference any potential transcript isoform-specfic modification information is lost or the signal diluted. Leaving open the potential for isoform-specific modified base detection is one reason for the choice to force mapping modified bases to a transcriptome. Regions at the edge of alternative exons also have divergent expected signal levels and thus genome statistics computed at these positions would be very difficult to process into a logical output. Processing would also be very sensetive to shifts in the mapped splice boundaries which can be variable with nanopore reads.

Tools to investigate isoform-specific modified bases is a future goal within the Tombo framework. This does pose some informatic challenges for downstream processing of Tombo RNA data. A recommended Tombo RNA processing pipeline will be posted here soon to help make integrative modified RNA processing more streamlined with other genome bioinformatic tools.

A second minor note is that since RNA is currently sequenced in the 3' to 5' direction; thus special care must be taken when accessing Tombo re-squiggled raw signal data. The raw signal (from MinKNOW) and albacore basecalled events are stored in the reverse direction from the genome (3' to 5' for reads mapping to the plus genome strand). Tombo events for RNA data are stored in the opposite direction (corresponding to the genome strand sequence direction, not sequencing time direction) for several practical reasons. Thus if Tombo events are to be compared to the raw signal, the raw signal must be reversed. Tombo RNA models are stored in this direction as well and thus may be considered inverted as compared to some other RNA HMM signal level models processing data in the sequencing time direction.

-----------------------
RNA Processing Workflow
-----------------------

As Tombo RNA processing presents unique informatic challenges, a recommended processing pipeline will be posted here soon.

This pipeline is for users looking to process a sample from a genome seqeuence reference and a gene annotation file (GTF or GFF). For users successfully processing data from a transcriptome reference this processing workflow will not be applicable.

This pipeline aims to address the majority of use cases for RNA modified base detection, namely porting Tombo results to a genome browser compatible format. Please check back soon for the recommended Tombo RNA processing pipeline!

This pipeline will likely be built on the `R/bioconductor ensembldb package <https://bioconductor.org/packages/release/bioc/html/ensembldb.html>`_. This package allows the creation of custom data bases and `mapping between genome and transcriptome coordinates <https://bioconductor.org/packages/release/bioc/vignettes/ensembldb/inst/doc/coordinate-mapping.html>`_. The functions from this software are recommended in order to project Tombo RNA results into a genome coordinate space. A full tutorial/example script for a full Tombo pipeline based on this package will be provided soon.
