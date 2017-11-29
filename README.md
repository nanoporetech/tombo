## Summary
Tombo is a suite of tools primarily for the identification of modified nucleotides from nanopore sequencing data. Tombo also provides tools for the analysis and visualization of raw nanopore signal.


## Installation
Basic tombo installation (includes resquiggle, statistical testing and text file outputs)
```
# numpy install required before tombo for cython optimization
pip install numpy
pip install git+https://github.com/nanoporetech/tombo.git
```

> Full installation instructions below


## Use Cases

#### Perform re-squiggle algorithm (once per root directory of FAST5 files)
###### * Need not contain Events data, but must contain Fastq slot
```
tombo resquiggle path/to/amplified/dna/fast5s/ genome.fasta --minimap2-executable ./minimap2 --processes 4
```

#### Identify modified nucleotides comparing to an alternative 5mC model
```
tombo test_significance --fast5-basedirs path/to/native/dna/fast5s/ --alternate-bases 5mC --statistics-file-basename sample_compare
tombo write_wiggles --fast5-basedirs path/to/native/dna/fast5s/ --wiggle-basename sample_compare.5mC --statistics-filename sample_compare.5mC.tombo.stats
```

#### Identify modified nucleotides from a single sample
```
tombo test_significance --fast5-basedirs path/to/native/dna/fast5s/ --statistics-file-basename sample --processes 4
tombo write_wiggles --fast5-basedirs path/to/native/dna/fast5s/ --wiggle-basename sample --statistics-filename sample.tombo.stats
```

#### Identify modified nucleotides comparing two samples
```
tombo test_significance --fast5-basedirs path/to/native/dna/fast5s/ --control-fast5-basedirs path/to/amplified/dna/fast5s/  --statistics-file-basename sample_compare
tombo write_wiggles --fast5-basedirs path/to/native/dna/fast5s/ --control-fast5-basedirs path/to/amplified/dna/fast5s/ --wiggle-basename sample_compare --statistics-filename sample_compare.tombo.stats
```

#### Identify modified RNA nucleotides from a single sample
```
tombo resquiggle path/to/native/rna/fast5s/ transcriptome.fasta --minimap2-executable ./minimap2 --processes 4
tombo test_significance --fast5-basedirs path/to/native/rna/fast5s/ --statistics-file-basename rna_sample --processes 4
tombo write_wiggles --fast5-basedirs path/to/native/rna/fast5s/ --wiggle-basename rna_sample --statistics-filename rna_sample.tombo.stats
```

#### Some plotting examples
```
tombo plot_max_coverage --fast5-basedirs path/to/native/rna/fast5s/ --tombo-model-filename /path/to/tombo/root/tombo/tombo_models/tombo.DNA.model
tombo plot_motif_centered --fast5-basedirs path/to/native/rna/fast5s/ --motif AWC --genome-fasta genome.fasta --control-fast5-basedirs path/to/amplified/dna/fast5s/ --deepest-coverage
tombo plot_per_read --fast5-basedirs path/to/native/rna/fast5s/ --genome-locations chromosome:1000 chromosome:2000:- --alternate-model-filename /path/to/tombo/root/tombo/tombo_models/tombo.DNA.5mC.model
```

#### Extract sequences surrounding modified positions
```
tombo write_most_significant_fasta --statistics-filename sample_compare.5mC.tombo.stats --genome-fasta genome.fasta
```


## Usage

```
tombo -h
tombo [command] [options]
```

#### Resquiggle (Must be run before any other commands):
     resquiggle               	   Re-annotate raw signal with genomic aignement of existing basecalls.

#### Statistical Testing Command:
     test_significance             Test for shifts in signal against a reference or against another set of reads.

#### Text Output Commands:
     write_wiggles                 Write wiggle files for nanopore signal values, coverage, and statistics.
     write_most_significant_fasta  Write sequence where signal differs the most significantly between two groups.

#### Genome Anchored Plotting Commands:
     plot_max_coverage             Plot signal in regions with maximal coverage.
     plot_genome_location          Plot signal at defined genomic locations.
     plot_motif_centered           Plot locations centered on a specific motif.
     plot_max_difference           Plot locations where signal differs the most between two groups.
     plot_most_significant         Plot locations where signal differs most significantly between two groups.
     plot_model_most_significant   Plot locations where signal differs most significantly from the kmer model.
     plot_motif_with_stats         Plot signal from several regions and test statistics centered on a k-mer of interst.
     plot_per_read                 Plot per read modified base predictions.

#### Sequencing Time Anchored Plotting Commands:
     plot_correction               Plot segmentation before and after correction.
     plot_multi_correction         Plot multiple raw signals anchored by genomic location.

#### Other Plotting Commands:
     cluster_most_significant      Clustering traces at bases with significant differences.
     plot_kmer                     Plot signal quantiles acorss kmers.

#### Read Filtering (Only effects tombo index file):
     clear_filters                 Clear filters to process all successfully re-squiggled reads.
     filter_stuck                  Apply filter based on observations per base thresholds.

####  Event-based Re-squiggle (Primarily for producing new models):
     event_resquiggle              Re-annotate raw signal with genomic alignment from existing basecalls.
     model_resquiggle              Re-annotate raw signal with genomic bases by shifting the signal to more closely match a k-mer model.
     estimate_kmer_reference       Estimate reference k-mer model derived from the provided reads.
     estimate_alt_reference        Estimate alternative tombo model from a sample containing standard bases spiked with a single non-standard base at random positions.

> Get additional help for subcommands with `tombo [command] -h`


## Requirements
At least one supported mapper:
- minimap2 (<https://github.com/lh3/minimap2>)
- BWA-MEM (<http://bio-bwa.sourceforge.net/>)
- graphmap (<https://github.com/isovic/graphmap>)

- HDF5 (<http://micro.stanford.edu/wiki/Install_HDF5#Install>)

#### python Requirements (handled by pip):
- numpy (must be installed before installing tombo)
- scipy
- h5py
- cython

#### Optional packages for plotting (install R packages with `install.packages([package_name])` from an R prompt):
- rpy2 (along with an R installation)
- ggplot2 (required for any plotting subcommands)
- cowplot (required for plot_motif_with_stats subcommand)

#### Optional packages for alternative model estimation:
- sklearn


## Installation along with additional dependencies

Install tombo with all optional dependencies (for plotting and model estimation)
```
pip install git+https://github.com/nanoporetech/tombo.git[full]
```

Install tombo with plotting dependencies (requires separate installation of R packages ggplot2 and cowplot)
```
pip install git+https://github.com/nanoporetech/tombo.git[plot]
```

Install tombo with alternative model estimation dependencies
```
pip install git+https://github.com/nanoporetech/tombo.git[alt_est]
```


## Citation
Stoiber, M.H. et al. De novo Identification of DNA Modifications Enabled by Genome-Guided Nanopore Signal Processing. bioRxiv (2016).

http://biorxiv.org/content/early/2017/04/10/094672
