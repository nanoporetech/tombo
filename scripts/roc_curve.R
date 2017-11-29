library(dplyr)

dat <- read.table('roc_output.txt', header=TRUE, sep="\t", stringsAsFactors=FALSE)
dat %>% group <- by(Comparison) %>% summarize(auc=sum(TP * (FP - lag(FP)), na.rm=TRUE))

dat$Comparison[dat$Comparison == "Alternative 5mC Model"] <- "Alternative\n5mC Model"
dat$Comparison[dat$Comparison == 'PCR or Standard Model Comparison'] <- 'PCR or Standard\nModel Comparison'

library(ggplot2)
pdf('comp_model_roc.pdf', height=4.5, width=6)
ggplot(dat) + geom <- abline(slope=1, intercept=0) +
        geom <- path(aes(x=FP, y=TP, color=Comparison)) + theme <- bw() +
                xlab('False Positive Rate') + ylab('True Positive Rate')
foo <- dev.off()
