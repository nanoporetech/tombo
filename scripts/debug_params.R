library(dplyr)
library(ggplot2)
library(ggbeeswarm)

## set _DEBUG_PARAMS = True in resquiggle.py
## example run for min_obs_per_base testing:
##for i in {0..6}; do
##    testParam=`echo $i | awk '{print ($1 * 1) + 2}'`
##    tombo resquiggle param_test_reads/ genome.fasta --segmentation-parameters 5 $testParam 5 --signal-align-parameters 4.2 4.2 1200 1.75 5.0 --processes 4
##done > param_values.txt

stat <- 'min_obs_per_base'

dat <- read.table('param_values.txt')
colnames(dat) <- c('running_window', 'min_obs_per_base', 'mean_obs_per_event',
                   'match_evalue', 'skip_pen', 'bandwidth',
                   'read_name', 'mean_score')
dat$mean_obs_per_event <- factor(dat$mean_obs_per_event)
dat$running_window <- factor(dat$running_window)
dat$min_obs_per_base <- factor(dat$min_obs_per_base)
dat$match_evalue <- factor(dat$match_evalue)
dat$skip_pen <- factor(dat$skip_pen)
dat$bandwidth <- factor(dat$bandwidth)

dat <- dat %>% group_by(mean_obs_per_event, min_obs_per_base, running_window,
                        match_evalue, skip_pen, bandwidth, read_name) %>% summarize(mean_score=min(mean_score))

rdat <- dat %>% group_by(read_name) %>% summarize(nreads=n())
maxNReads <- rdat$read_name[which(rdat$nreads == max(rdat$nreads))]
fdat <- dat %>% filter(read_name %in% maxNReads)

minMed <- dat %>% group_by_at(stat) %>% summarize(med=median(mean_score)) %>% summarize(min(med))
minMedF <- fdat %>% group_by_at(stat) %>% summarize(med=median(mean_score)) %>% summarize(min(med))

pdf(paste0('param_values.', stat, '.pdf'), width=10)
ggplot(dat, aes_string(x=stat, y='mean_score', color=stat)) +
    geom_hline(aes(yintercept=minMed)) +
    geom_beeswarm(alpha=0.3, cex=0.5) +
    stat_summary(fun.y=median, color='red', geom='point', size=2) +
    stat_summary(fun.y=mean, color='orange', geom='point', size=2) +
    theme_bw() + theme(axis.text.x=element_text(angle=60, hjust=1))
ggplot(fdat, aes_string(x=stat, y='mean_score', color=stat)) +
    geom_hline(aes(yintercept=minMedF)) +
    geom_beeswarm(alpha=0.3, cex=0.5) +
    stat_summary(fun.y=median, color='red', geom='point', size=2) +
    stat_summary(fun.y=mean, color='orange', geom='point', size=2) +
    theme_bw() + theme(axis.text.x=element_text(angle=60, hjust=1))
foo <- dev.off()
