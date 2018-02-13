**************
RNA Processing
**************

Processing RNA data within the Tombo framework requires some extra care. The major item to consider when performing RNA processing is that a transcriptome reference must be supplied as spliced mapping is not supported within the Tombo framework. The lack of spliced mapping support within the Tombo framework is a conscious decision for identification of modified RNA bases. This is because the transcriptome is the natural setting for the detection of modified RNA bases. When modified RNA bases are projected onto the genome reference any potential transcript isoform-specfic modification information is lost. Leaving open the potential for isoform-specific modified base detection is a main reason for the choice to force mapping modified bases to a transcriptome. Tools to investigate isoform-specific modified bases is a future goal within the Tombo framework. This does pose some informatic challenges for downstream processing of Tombo RNA data. A recommended Tombo RNA processing pipeline will be posted here soon.

A second minor note is that since RNA is currently sequenced in the 3' to 5' direction; thus special care must be taken when accessing Tombo re-squiggled binary data. The raw signal (from MinKNOW) and albacore basecalled events are stored in the reverse direction from the genome (3' to 5' for reads mapping to the plus genome strand). Tombo events for RNA data are stored in the opposite direction (corresponding to the genome sequence direction, not sequencing time direction) for several practical reasons. Thus if Tombo events are to be compared to the raw signal, the raw signal must be reversed. Tombo RNA models are stored in this direction as well and thus may be considered inverted as compared to some other RNA HMM signal level models processing data in the sequencing time direction.

-----------------------
RNA Processing Workflow
-----------------------

As Tombo RNA processing presents unique informatic challenges a recommended processing pipeline will be posted here soon. This pipeline aims to address the majority of use cases for RNA modified base detection including porting Tombo results to a genome browser compatible format. Please check back soon for the recommended Tombo RNA processing pipeline!
