suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggridges))

## set _DEBUG_PARAMS = True in resquiggle.py
## example run for min_obs_per_base testing:
##for i in {0..2}; do
##    testParam=`echo $i | awk '{print ($1 * 1) + 1}'`
##    tombo resquiggle param_test_reads/ genome.fasta \
##        --segmentation-parameters 5 3 $testParam 5 \
##        --signal-align-parameters 4.2 4.2 300 1500 20.0 40 750 2500 250 \
##        --processes 4
##done > param_values.txt


dat <- read.table('param_values.txt', header=TRUE)

## update with output from tombo.resquiggle._write_params_debug
colnames(dat) <- c(
    'running_window', 'min_obs_per_base', 'raw_min_obs_per_base',
    'mean_obs_per_event', 'match_evalue', 'skip_pen', 'bandwidth',
    'read_name', 'mean_score')

## filter out save bandwidth reads if bandwidth was not tested here
if(length(unique(dat$bandwidth))  == 2){
    dat <- dat %>% filter(bandwidth == min(as.numeric(dat$bandwidth)))
}

#3 convert params to factors and get stat that was investigated in this run
param_names <- setdiff(colnames(dat), c('mean_score', 'read_name'))
for(param_name in param_names){
    dat[,param_name] <- factor(dat[,param_name])
}
stat <- param_names[which.max(sapply(param_names, function(param_name)
    length(unique(dat[,param_name]))))]

## take min score over same read with same params (over re-scalings)
dat <- dat %>% group_by_at(c(param_names, 'read_name')) %>%
    summarize(mean_score=min(mean_score))

## filter for reads included in all parameter groups
rdat <- dat %>% group_by(read_name) %>% summarize(nreads=n())
maxNReads <- rdat$read_name[which(rdat$nreads == max(rdat$nreads))]
fdat <- dat %>% filter(read_name %in% maxNReads)

## compute and print mean and median stats for plotting
sumDat <- dat %>% group_by_at(stat) %>%
    summarize(med=median(mean_score), mean=mean(mean_score))
sumFDat <- fdat %>% group_by_at(stat) %>%
    summarize(med=median(mean_score), mean=mean(mean_score))

sumDat %>% print.data.frame(digits=6)
sumFDat %>% print.data.frame(digits=6)

pdf(paste0('param_values.', stat, '.pdf'), width=10)
ggplot(dat, aes_string(y=stat, x='mean_score', fill=stat)) +
    geom_vline(aes(xintercept=min(sumDat$med))) +
    geom_density_ridges(alpha=0.3, cex=0.5) +
    geom_point(aes_string(x='med', y=stat),
               color='red', size=5, data=sumDat) +
    geom_point(aes_string(x='mean', y=stat),
               color='orange', size=5, data=sumDat) +
    theme_minimal() + theme(legend.position="none") +
    coord_cartesian(xlim=quantile(dat$mean_score, c(0.01, 0.98)))
ggplot(fdat, aes_string(y=stat, x='mean_score', fill=stat)) +
    geom_vline(aes(xintercept=min(sumFDat$med))) +
    geom_density_ridges(alpha=0.3, cex=0.5) +
    geom_point(aes_string(x='med', y=stat),
               color='red', size=5, data=sumFDat) +
    geom_point(aes_string(x='mean', y=stat),
               color='orange', size=5, data=sumFDat) +
    theme_minimal() + theme(legend.position="none") +
    coord_cartesian(xlim=quantile(fdat$mean_score, c(0.01, 0.98)))
foo <- dev.off()
