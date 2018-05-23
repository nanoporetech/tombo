natDir='test_data/native_reads/'
ampDir='test_data/amplified_reads/'
rcsvDir='test_data/recursive_test/'
natFsq="test_data/fastqs.native.fasta"
natFqDir='test_data/native_reads.for_fastq_ann/'

nrModFn='../tombo/tombo/tombo_models/tombo.DNA.model'
altModFn='../tombo/tombo/tombo_models/tombo.DNA.5mC.model'
poreModel="r9_250bps.nucleotide.5mer.template.model"

genomeFn="e_coli.K12.NEB5alpha.fasta"
mmiFn="e_coli.K12.NEB5alpha.mmi"

genomeLocs='"CP017100.1:1505285" "CP017100.1:1504705"'
strandGenomeLocs='"CP017100.1:1505285:+" "CP017100.1:1504705:+"'

runHelps=false
runResquiggle=true

printf "********* Testing help commands **********\n"
tombo -h || { echo 'Main tombo help failed. Likely a syntax error.' ; exit 1; }

if [ $runHelps == true ]
then
tombo resquiggle -h

tombo preprocess annotate_raw_with_fastqs -h

tombo filter clear_filters -h
tombo filter stuck -h
tombo filter level_coverage -h
tombo filter q_score -h
tombo filter raw_signal_matching -h
tombo filter genome_locations -h

tombo detect_modifications de_novo -h
tombo detect_modifications alternative_model -h
tombo detect_modifications sample_compare -h
tombo detect_modifications aggregate_per_read_stats -h

tombo text_output browser_files -h
tombo text_output signif_sequence_context -h

tombo plot max_coverage -h
tombo plot genome_locations -h
tombo plot motif_centered -h
tombo plot max_difference -h
tombo plot most_significant -h
tombo plot motif_with_stats -h
tombo plot per_read -h

tombo plot roc -h
tombo plot per_read_roc -h
tombo plot kmer -h
tombo plot cluster_most_significant -h

tombo build_model estimate_scale -h
tombo build_model event_resquiggle -h
tombo build_model estimate_reference -h
tombo build_model estimate_alt_reference -h
fi

if [ $runResquiggle == true ]
then
printf "\n\n********* Testing re-squiggle command **********\n"
tombo resquiggle \
      $natDir $genomeFn \
      --failed-reads-filename testing.native.failed_read.txt \
      --processes 4 --overwrite --include-event-stdev
tombo resquiggle \
      $ampDir $genomeFn \
      --failed-reads-filename testing.amplified.failed_read.txt \
      --processes 4 --overwrite --include-event-stdev

printf "\n\n********* Testing FASTQ annotation and re-squiggle **********\n"
tombo preprocess annotate_raw_with_fastqs --fast5-basedir $natFqDir \
      --fastq-filenames $natFsq --overwrite
tombo resquiggle \
      $natFqDir $genomeFn \
      --corrected-group FastqAnnotation \
      --failed-reads-filename testing.native.fastq_ann.failed_read.txt \
      --processes 4 --overwrite

printf "\n\n********* Testing re-squiggle command with filename **********\n"
tombo resquiggle \
      $natDir $genomeFn --tombo-model-filename $nrModFn \
      --corrected-group RawWFilenameCorrected \
      --processes 4 --overwrite \
      --failed-reads-filename testing.native.fn_model.failed_read.txt
tombo resquiggle \
      $ampDir $genomeFn --tombo-model-filename $nrModFn \
      --corrected-group RawWFilenameCorrected \
      --processes 4 --overwrite \
      --failed-reads-filename testing.amplified.fn_model.failed_read.txt

printf "\n\n********* Testing event-based resquiggle **********\n"
tombo build_model event_resquiggle \
      $natDir $genomeFn --minimap2-executable ./minimap2 \
      --corrected-group RawEventCorrected --processes 4 --overwrite \
      --failed-reads-filename testing.native.failed_read.event.txt

printf "\n\n********* Testing minimap2 index **********\n"
tombo resquiggle \
      $natDir $mmiFn \
      --corrected-group RawMinimapIndexCorrected \
      --processes 4 --overwrite \
      --failed-reads-filename testing.native.failed_read.txt

printf "\n\n********* Testing pA normalization **********\n"
tombo build_model event_resquiggle --minimap2-executable ./minimap2 \
      $natDir $genomeFn \
      --normalization-type pA_raw --processes 4 \
      --corrected-group RawGenomeCorrected_pA_raw_000 --overwrite \
      --failed-reads-filename testing.native.pA_raw.failed_read.txt
tombo build_model event_resquiggle \
        $natDir $genomeFn --minimap2-executable ./minimap2 \
        --normalization-type pA --pore-model-filename $poreModel \
        --corrected-group RawGenomeCorrected_pA_000 --overwrite \
        --failed-reads-filename testing.native.pA.failed_read.txt \
        --processes 4

printf "\n\n********* Testing recursive resquiggle **********\n"
tombo resquiggle \
        $rcsvDir $genomeFn \
        --failed-reads-filename testing.recursive.failed_read.txt \
        --processes 4 --overwrite
fi

printf "\n\n********* Testing filter functions **********\n"
tombo filter clear_filters --fast5-basedirs $natDir
tombo filter stuck --fast5-basedirs $natDir \
      --obs-per-base-filter 99:200 100:5000
tombo filter level_coverage --fast5-basedirs $natDir \
      --percent-to-filter 10
tombo filter q_score --fast5-basedirs $natDir --q-score 21
tombo filter raw_signal_matching --fast5-basedirs $natDir \
      --signal-matching-score 0.75
tombo filter clear_filters --fast5-basedirs $natDir
tombo filter genome_locations --fast5-basedirs $natDir \
      --include-regions CP017100.1:1,458,474-1,558,736
tombo filter clear_filters --fast5-basedirs $natDir


printf "\n\n********* Testing estimate global scale function **********\n"
tombo build_model estimate_scale $natDir

printf "\n\n********* Testing single sample genome-anchored plotting functions **********\n"
tombo plot max_coverage --fast5-basedirs $natDir \
        --num-bases 21 --overplot-threshold 1000 \
        --pdf-filename testing.max_cov.1_samp.pdf
tombo plot genome_locations --fast5-basedirs $ampDir \
        --genome-locations $genomeLocs \
        --num-bases 21 --overplot-threshold 1000 \
        --pdf-filename testing.genome_loc.1_samp.pdf
tombo plot motif_centered --fast5-basedirs $natDir --motif ATC \
        --genome-fasta $genomeFn \
        --num-bases 21 --overplot-threshold 1000 \
        --pdf-filename testing.motif_centered.1_samp.pdf
tombo plot motif_centered --fast5-basedirs $natDir --motif TWA \
        --genome-fasta $genomeFn \
        --num-bases 21 --overplot-threshold 1000 \
        --pdf-filename testing.motif_centered.palindrome.1_samp.pdf
tombo plot motif_centered --fast5-basedirs $natDir --motif ATC \
        --genome-fasta $genomeFn \
        --num-bases 21 --overplot-threshold 1000 --deepest-coverage \
        --pdf-filename testing.motif_centered.deepest.1_samp.pdf
tombo plot max_coverage --fast5-basedirs $rcsvDir \
        --num-bases 21 --overplot-threshold 1000 \
        --pdf-filename testing.max_cov.1_samp.recursive.pdf


printf "\n\n********* Testing multiple sample genome-anchored plotting functions **********\n"
tombo plot max_coverage --fast5-basedirs $natDir \
        --control-fast5-basedirs $ampDir \
        --num-bases 21 --overplot-threshold 1000 \
        --pdf-filename testing.max_cov.2_samp.pdf
tombo plot genome_locations --fast5-basedirs $natDir \
        --control-fast5-basedirs $ampDir \
        --genome-locations $genomeLocs \
        --num-bases 21 --overplot-threshold 1000 \
        --pdf-filename testing.genome_loc.2_samp.pdf
tombo plot motif_centered --fast5-basedirs $natDir --motif ATC \
        --genome-fasta $genomeFn \
        --control-fast5-basedirs $ampDir \
        --num-bases 21 --overplot-threshold 1000 \
        --pdf-filename testing.motif_centered.2_samp.pdf
tombo plot motif_centered --fast5-basedirs $natDir --motif TWA \
        --genome-fasta $genomeFn \
        --control-fast5-basedirs $ampDir \
        --num-bases 21 --overplot-threshold 1000 \
        --pdf-filename testing.motif_centered.palindrome.2_samp.pdf
tombo plot motif_centered --fast5-basedirs $natDir --motif ATC \
        --genome-fasta $genomeFn \
        --control-fast5-basedirs $ampDir \
        --num-bases 21 --overplot-threshold 1000 --deepest-coverage \
        --pdf-filename testing.motif_centered.deepest.2_samp.pdf

printf "\n\n********* Testing statistical testing. **********\n"
rm test_stats.de_novo.tombo.stats test_stats.2samp.tombo.stats \
   test_stats.alt_model.5mC.tombo.stats \
   test_stats.alt_model.6mA.tombo.stats \
   test_stats.alt_default_model.5mC.tombo.stats \
   test_stats.alt_default_model.6mA.tombo.stats \
   test_stats.de_novo.tombo.per_read_stats test_stats.2samp.tombo.per_read_stats \
   test_stats.alt_model.5mC.tombo.per_read_stats \
   test_stats.alt_model.6mA.tombo.per_read_stats \
   test_stats.alt_default_model.5mC.tombo.per_read_stats \
   test_stats.alt_default_model.6mA.tombo.per_read_stats \
   test_standard.model test_stats.de_novo.new_thresh.tombo.stats \
   test_alt_est.alt_C.tombo_model
tombo detect_modifications de_novo --fast5-basedirs $natDir \
      --minimum-test-reads 5 \
      --statistics-file-basename test_stats.de_novo \
      --per-read-statistics-basename test_stats.de_novo
tombo detect_modifications de_novo --fast5-basedirs $natDir \
      --minimum-test-reads 5 --single-read-threshold 0.1 0.75 \
      --statistics-file-basename test_stats.de_novo.two_way_thresh \
      --per-read-statistics-basename test_stats.de_novo.two_way_thresh
tombo detect_modifications sample_compare --fast5-basedirs $natDir \
      --control-fast5-basedirs $ampDir \
      --minimum-test-reads 5 \
      --statistics-file-basename test_stats.2samp \
      --per-read-statistics-basename test_stats.2samp
tombo detect_modifications alternative_model --fast5-basedirs $natDir \
      --alternate-bases 5mC 6mA \
      --statistics-file-basename test_stats.alt_default_model \
      --per-read-statistics-basename test_stats.alt_default_model
tombo detect_modifications alternative_model --fast5-basedirs $natDir \
      --tombo-model-filename $nrModFn \
      --alternate-model-filenames $altModFn \
      --statistics-file-basename test_stats.alt_model \
       --per-read-statistics-basename test_stats.alt_model
tombo build_model estimate_reference --fast5-basedirs $natDir \
      --tombo-model-filename test_standard.model \
      --upstream-bases 1 --downstream-bases 1 --minimum-kmer-observations 1
tombo build_model estimate_alt_reference --fast5-basedirs $natDir \
      --control-fast5-basedirs $ampDir \
      --tombo-model-filename test_standard.model \
      --alternate-model-filename test_alt.model \
      --alternate-model-name 5mC --alternate-model-base C \
      --minimum-kmer-observations 1 --save-density-basename test_save_dens
tombo build_model estimate_alt_reference \
      --alternate-density-filename test_save_dens.alternate_density.txt \
      --control-density-filename test_save_dens.control_density.txt \
      --tombo-model-filename test_standard.model \
      --alternate-model-filename test_alt.use_densities.model \
      --alternate-model-name 5mC --alternate-model-base C \
      --minimum-kmer-observations 1

printf "\n\n********* Testing aggregate per-read stats **********\n"
tombo detect_modifications aggregate_per_read_stats --minimum-test-reads 5 \
      --single-read-threshold 0.4 \
      --statistics-file-basename test_stats.de_novo.new_thresh \
      --per-read-statistics-filename test_stats.de_novo.tombo.per_read_stats

printf "\n\n********* Testing ROC and Precision-Recall plotting **********\n"
tombo plot roc --genome-fasta e_coli.K12.NEB5alpha.fasta \
      --statistics-filenames test_stats.2samp.tombo.stats \
      test_stats.alt_default_model.5mC.tombo.stats \
      test_stats.alt_default_model.6mA.tombo.stats \
      test_stats.de_novo.tombo.stats test_stats.de_novo.new_thresh.tombo.stats \
      --motif-descriptions \
      CCWGG:2:"dcm 5mC Samp Comp"::GATC:2:"dam 6mA Samp Comp" \
      CCWGG:2:"dcm 5mC Alt Test" GATC:2:"dam 6mA Alt Test" \
      CCWGG:2:"dcm 5mC De Novo"::GATC:2:"dam 6mA De Novo" \
      CCWGG:2:"dcm 5mC De Novo New Thresh"::GATC:2:"dam 6mA De Novo New Thresh"

tombo plot per_read_roc --genome-fasta e_coli.K12.NEB5alpha.fasta \
      --per-read-statistics-filenames test_stats.2samp.tombo.per_read_stats \
      test_stats.alt_default_model.5mC.tombo.per_read_stats \
      test_stats.alt_default_model.6mA.tombo.per_read_stats \
      test_stats.de_novo.tombo.per_read_stats --motif-descriptions \
      CCWGG:2:"dcm 5mC Samp Comp"::GATC:2:"dam 6mA Samp Comp" \
      CCWGG:2:"dcm 5mC Alt Test" GATC:2:"dam 6mA Alt Test" \
      CCWGG:2:"dcm 5mC De Novo"::GATC:2:"dam 6mA De Novo"

printf "\n\n********* Testing mutliple sample statistical testing genome-anchored plotting functions **********\n"
tombo plot max_difference --fast5-basedirs $natDir \
      --control-fast5-basedirs $ampDir \
      --num-bases 21 --overplot-threshold 1000 \
      --pdf-filename testing.max_diff.pdf
tombo plot most_significant --fast5-basedirs $natDir \
      --control-fast5-basedirs $ampDir \
      --num-bases 21 --overplot-threshold 1000 \
      --statistics-filename test_stats.2samp.tombo.stats \
      --pdf-filename testing.most_signif.2samp.pdf
tombo plot most_significant --fast5-basedirs $natDir \
      --control-fast5-basedirs $ampDir \
      --num-bases 21 --overplot-threshold 1000 \
      --statistics-filename test_stats.de_novo.tombo.stats \
      --pdf-filename testing.most_signif.model.pdf
tombo plot most_significant --fast5-basedirs $natDir \
      --control-fast5-basedirs $ampDir \
      --num-bases 21 --overplot-threshold 1000 \
      --plot-standard-model \
      --statistics-filename test_stats.alt_model.5mC.tombo.stats \
      --pdf-filename testing.most_signif.alt_model_5mC.pdf
tombo plot motif_with_stats --fast5-basedirs $natDir \
      --motif CAW --genome-fasta $genomeFn --overplot-threshold 1000  \
      --plot-standard-model --statistics-filename test_stats.de_novo.tombo.stats \
      --pdf-filename testing.motif_w_stats.pdf
tombo plot motif_with_stats --fast5-basedirs $natDir \
      --tombo-model-filename $nrModFn --motif CCWGG \
      --genome-fasta $genomeFn --overplot-threshold 1000 \
      --statistics-filename test_stats.de_novo.tombo.stats \
      --pdf-filename testing.motif_w_stats.model.pdf
tombo plot motif_with_stats --fast5-basedirs $natDir \
      --control-fast5-basedirs $ampDir --motif CCWGG \
      --genome-fasta $genomeFn --overplot-threshold 1000 \
      --statistics-filename test_stats.2samp.tombo.stats \
      --pdf-filename testing.motif_w_stats.2samp.pdf
tombo plot motif_with_stats --fast5-basedirs $natDir \
      --plot-alternate-model 5mC --motif CCWGG --genome-fasta $genomeFn \
      --statistics-filename test_stats.alt_model.5mC.tombo.stats \
      --pdf-filename testing.motif_w_stats.alt_model_5mC.pdf
tombo plot motif_with_stats --fast5-basedirs $natDir \
      --plot-alternate-model 6mA --motif GATC --genome-fasta $genomeFn \
      --statistics-filename test_stats.alt_default_model.6mA.tombo.stats \
      --pdf-filename testing.motif_w_stats.alt_model_6mA.alt_dist.pdf

printf "\n\n********* Testing overplotting options **********\n"
tombo plot max_coverage --fast5-basedirs $natDir \
        --num-bases 21 --overplot-threshold 1 --overplot-type Downsample \
        --pdf-filename testing.max_coverage.Downsample.pdf
tombo plot max_coverage --fast5-basedirs $natDir \
        --num-bases 21 --overplot-threshold 1 --overplot-type Boxplot \
        --pdf-filename testing.max_coverage.Boxplot.pdf
tombo plot max_coverage --fast5-basedirs $natDir \
        --num-bases 21 --overplot-threshold 1 --overplot-type Quantile \
        --pdf-filename testing.max_coverage.Quantile.pdf
tombo plot max_coverage --fast5-basedirs $natDir \
        --num-bases 21 --overplot-threshold 1 --overplot-type Density \
        --pdf-filename testing.max_coverage.Density.pdf
tombo plot max_coverage --fast5-basedirs $natDir \
        --control-fast5-basedirs $ampDir \
        --num-bases 21 --overplot-threshold 1 --overplot-type Downsample \
        --pdf-filename testing.max_coverage.2samp.Downsample.pdf
tombo plot max_coverage --fast5-basedirs $natDir \
        --control-fast5-basedirs $ampDir \
        --num-bases 21 --overplot-threshold 1 --overplot-type Boxplot \
        --pdf-filename testing.max_coverage.2samp.Boxplot.pdf
tombo plot max_coverage --fast5-basedirs $natDir \
        --control-fast5-basedirs $ampDir \
        --num-bases 21 --overplot-threshold 1 --overplot-type Quantile \
        --pdf-filename testing.max_coverage.2samp.Quantile.pdf
tombo plot max_coverage --fast5-basedirs $natDir \
        --control-fast5-basedirs $ampDir \
        --num-bases 21 --overplot-threshold 1 --overplot-type Density \
        --pdf-filename testing.max_coverage.2samp.Density.pdf

printf "\n\n********* Testing model-based plotting **********\n"
tombo plot max_coverage --fast5-basedirs $natDir \
        --tombo-model-filename $nrModFn \
        --pdf-filename testing.max_cov.1_samp.model.pdf
tombo plot most_significant --fast5-basedirs $natDir \
        --tombo-model-filename $nrModFn \
        --statistics-filename test_stats.de_novo.tombo.stats \
        --pdf-filename testing.model_plotting.pdf
tombo plot most_significant --fast5-basedirs $natDir \
        --tombo-model-filename $nrModFn \
        --statistics-filename test_stats.de_novo.tombo.stats \
        --overplot-threshold 1 --overplot-type Downsample \
        --pdf-filename testing.model_plotting.downsample.pdf
tombo plot most_significant --fast5-basedirs $natDir \
        --tombo-model-filename $nrModFn \
        --statistics-filename test_stats.de_novo.tombo.stats \
        --overplot-threshold 1 --overplot-type Boxplot \
        --pdf-filename testing.model_plotting.boxplot.pdf
tombo plot most_significant --fast5-basedirs $natDir \
        --tombo-model-filename $nrModFn \
        --statistics-filename test_stats.de_novo.tombo.stats \
        --overplot-threshold 1 --overplot-type Quantile \
        --pdf-filename testing.model_plotting.quant.pdf
tombo plot most_significant --fast5-basedirs $natDir \
        --tombo-model-filename $nrModFn \
        --statistics-filename test_stats.de_novo.tombo.stats \
        --overplot-threshold 1 --overplot-type Density \
        --pdf-filename testing.model_plotting.density.pdf
tombo plot genome_locations --fast5-basedirs $ampDir \
        --tombo-model-filename $nrModFn \
        --alternate-model-filename $altModFn \
        --genome-locations $genomeLocs \
        --num-bases 21 --overplot-threshold 1000 \
        --pdf-filename testing.genome_loc.two_model_comp.pdf

printf "\n\n********* Testing event-resquiggled plotting **********\n"
tombo plot max_coverage --fast5-basedirs $natDir \
        --tombo-model-filename $nrModFn \
        --pdf-filename testing.max_cov.1_samp.model.model_resq.pdf \
        --corrected-group RawEventCorrected
tombo plot most_significant --fast5-basedirs $natDir \
        --tombo-model-filename $nrModFn \
        --corrected-group RawEventCorrected \
        --statistics-filename test_stats.de_novo.tombo.stats \
        --pdf-filename testing.model_plotting.resq_most_signif.pdf
tombo plot max_coverage --fast5-basedirs $natDir \
        --tombo-model-filename $nrModFn \
        --pdf-filename testing.max_cov.1_samp.model.model_resq.pdf \
        --corrected-group RawEventCorrected

printf "\n\n********* Testing per-read testing plot **********\n"
tombo plot per_read --genome-locations $genomeLocs --num-bases 101 \
      --per-read-statistics-filename test_stats.2samp.tombo.per_read_stats \
      --genome-fasta $genomeFn --pdf-filename testing.per_read.pdf
tombo plot per_read --genome-locations $genomeLocs --num-bases 101 \
      --per-read-statistics-filename test_stats.de_novo.tombo.per_read_stats \
      --genome-fasta $genomeFn --pdf-filename testing.de_novo.per_read.pdf
tombo plot per_read --fast5-basedirs $natDir \
      --genome-locations $genomeLocs --num-bases 101 \
      --per-read-statistics-filename test_stats.alt_model.5mC.tombo.per_read_stats \
      --pdf-filename testing.per_read.w_alt.pdf
tombo plot per_read --genome-locations $genomeLocs --num-bases 101 \
      --per-read-statistics-filename test_stats.alt_model.5mC.tombo.per_read_stats \
      --pdf-filename testing.per_read.wo_seq.pdf

printf "\n\n********* Testing auxilliary commands **********\n"
tombo text_output signif_sequence_context --fast5-basedirs $natDir $ampDir \
        --statistics-filename test_stats.de_novo.tombo.stats \
        --sequences-filename testing_signif_regions.from_fast5s.fasta
tombo text_output signif_sequence_context \
        --statistics-filename test_stats.de_novo.tombo.stats \
        --sequences-filename testing_signif_regions.from_fasta.fasta \
        --genome-fasta $genomeFn
tombo text_output browser_files --fast5-basedirs $natDir \
        --control-fast5-basedirs $ampDir \
        --file-types coverage signal signal_sd dwell difference
tombo text_output browser_files --fast5-basedirs $natDir \
        --control-fast5-basedirs $ampDir \
        --file-types coverage fraction signal signal_sd dwell difference \
        --statistics-filename test_stats.2samp.tombo.stats
tombo text_output browser_files --file-types fraction dampened_fraction \
        valid_coverage --statistics-filename \
        test_stats.de_novo.two_way_thresh.tombo.stats

printf "\n\n********* Testing other plotting commands **********\n"
tombo plot cluster_most_significant --fast5-basedirs $natDir \
        --control-fast5-basedirs $ampDir \
        --genome-fasta $genomeFn --num-regions 100 \
        --statistics-filename test_stats.2samp.tombo.stats
tombo plot cluster_most_significant --fast5-basedirs $natDir \
        --control-fast5-basedirs $ampDir \
        --genome-fasta $genomeFn --num-regions 100 \
        --r-data-filename testing.cluster_data.RData \
        --statistics-filename test_stats.2samp.tombo.stats
tombo plot kmer --fast5-basedirs $natDir \
        --num-kmer-threshold 0 \
        --pdf-filename testing.kmer_dist.median.all_events.pdf
tombo plot kmer --fast5-basedirs $natDir --read-mean \
        --num-kmer-threshold 1 \
        --pdf-filename testing.kmer_dist.median.pdf
tombo plot kmer --fast5-basedirs $natDir --read-mean \
        --corrected-group RawGenomeCorrected_pA_raw_000 \
        --num-kmer-threshold 1 \
        --pdf-filename testing.kmer_dist.pA_raw.pdf
tombo plot kmer --fast5-basedirs $natDir --read-mean \
        --corrected-group RawGenomeCorrected_pA_000 \
        --num-kmer-threshold 1 \
        --pdf-filename testing.kmer_dist.pA.pdf
