library(ggplot2)
library(stringr)
library(ggridges)


## density basename and alternative base
densBase <- 'debug_est_alt'
altBase <- 'C'


## parse density data
densDat <- read.table(paste0(densBase, '.alternate_density.txt'), header=TRUE)
standardDensDat <- read.table(paste0(densBase, '.control_density.txt'), header=TRUE)
densDat$Sample <- "Alternative"
standardDensDat$Sample <- "Standard"

allDensDat <- rbind.data.frame(densDat, standardDensDat)
sAllDat <- split(allDensDat, allDensDat$Kmer)
sDiffs <- sort(unlist(lapply(sAllDat, function(x)
    weighted.mean(x[x$Sample == 'Alternative','Signal'], x[x$Sample == 'Alternative','Density']) -
    weighted.mean(x[x$Sample == 'Standard','Signal'], x[x$Sample == 'Standard','Density']))))


## plot k-mers with the largest shifts in average signal level
upDat <- do.call(rbind.data.frame, lapply(names(head(sDiffs, 20)), function(kmer) sAllDat[[kmer]]))
dnDat <- do.call(rbind.data.frame, lapply(names(tail(sDiffs, 20)), function(kmer) sAllDat[[kmer]]))

pdf(paste0(densBase, '.density.pdf'), width=10)
ggplot(upDat, aes(x=Signal, y=Kmer, height=Density, fill=Sample)) +
    geom_ridgeline(alpha=0.4, size=0, color='white') +
    scale_fill_discrete(name='Contains\nAlternative\nBase') +
    theme_ridges() + theme(axis.text.y=element_text(family="mono"))
ggplot(dnDat, aes(x=Signal, y=Kmer, height=Density, fill=Sample)) +
    geom_ridgeline(alpha=0.4, size=0, color='white') +
    scale_fill_discrete(name='Contains\nAlternative\nBase') +
    theme_ridges() + theme(axis.text.y=element_text(family="mono"))
foo <- dev.off()


## plot alternative densities in k-mers without the alternative base of interest
noAltBaseDiffs <- sDiffs[!str_detect(names(sDiffs), altBase)]

upDat <- do.call(rbind.data.frame, lapply(names(head(noAltBaseDiffs, 20)), function(kmer) sAllDat[[kmer]]))
dnDat <- do.call(rbind.data.frame, lapply(names(tail(noAltBaseDiffs, 20)), function(kmer) sAllDat[[kmer]]))

pdf(paste0(densBase, '.noAlt.density.pdf'), width=10)
ggplot(upDat, aes(x=Signal, y=Kmer, height=Density, fill=Sample)) +
    geom_ridgeline(alpha=0.4, size=0, color='white') +
    scale_fill_discrete(name='Contains\nAlternative\nBase') +
    theme_ridges() + theme(axis.text.y=element_text(family="mono"))
ggplot(dnDat, aes(x=Signal, y=Kmer, height=Density, fill=Sample)) +
    geom_ridgeline(alpha=0.4, size=0, color='white') +
    scale_fill_discrete(name='Contains\nAlternative\nBase') +
    theme_ridges() + theme(axis.text.y=element_text(family="mono"))
foo <- dev.off()


singleAltBaseDiffs <- sDiffs[str_count(names(sDiffs), altBase) == 1]

upDat <- do.call(rbind.data.frame, lapply(names(head(singleAltBaseDiffs, 20)), function(kmer) sAllDat[[kmer]]))
dnDat <- do.call(rbind.data.frame, lapply(names(tail(singleAltBaseDiffs, 20)), function(kmer) sAllDat[[kmer]]))

pdf(paste0(densBase, '.singleAlt.density.pdf'), width=10)
ggplot(upDat, aes(x=Signal, y=Kmer, height=Density, fill=Sample)) +
    geom_ridgeline(alpha=0.4, size=0, color='white') +
    scale_fill_discrete(name='Contains\nAlternative\nBase') +
    theme_ridges() + theme(axis.text.y=element_text(family="mono"))
ggplot(dnDat, aes(x=Signal, y=Kmer, height=Density, fill=Sample)) +
    geom_ridgeline(alpha=0.4, size=0, color='white') +
    scale_fill_discrete(name='Contains\nAlternative\nBase') +
    theme_ridges() + theme(axis.text.y=element_text(family="mono"))
foo <- dev.off()


onePlusAltBaseDiffs <- sDiffs[str_count(names(sDiffs), altBase) > 0]

upDat <- do.call(rbind.data.frame, lapply(names(head(onePlusAltBaseDiffs, 20)), function(kmer) sAllDat[[kmer]]))
dnDat <- do.call(rbind.data.frame, lapply(names(tail(onePlusAltBaseDiffs, 20)), function(kmer) sAllDat[[kmer]]))

pdf(paste0(densBase, '.onePlusAlt.density.pdf'), width=10)
ggplot(upDat, aes(x=Signal, y=Kmer, height=Density, fill=Sample)) +
    geom_ridgeline(alpha=0.4, size=0, color='white') +
    scale_fill_discrete(name='Contains\nAlternative\nBase') +
    theme_ridges() + theme(axis.text.y=element_text(family="mono"))
ggplot(dnDat, aes(x=Signal, y=Kmer, height=Density, fill=Sample)) +
    geom_ridgeline(alpha=0.4, size=0, color='white') +
    scale_fill_discrete(name='Contains\nAlternative\nBase') +
    theme_ridges() + theme(axis.text.y=element_text(family="mono"))
foo <- dev.off()


## plot estimated shift correction due to high mod content samples
getMedMean <- function(x, sampType){
    sampDat <- x[x$Sample == sampType,]
    densCum <- cumsum(sampDat$Density)
    return(c(sampDat$Signal[which.min(abs(densCum - (tail(densCum, 1) / 2)))],
             weighted.mean(sampDat$Signal, sampDat$Density)))
}

standDiffs <- do.call(rbind.data.frame, lapply(sAllDat, function(x){
    stdMedMean <- getMedMean(x, 'Standard')
    altMedMean <- getMedMean(x, 'Alternative')
    data.frame(stdMed=stdMedMean[1],
               altMed=altMedMean[1],
               stdMean=stdMedMean[2],
               altMean=altMedMean[2],
               kmer=x$Kmer[1],
               hasAltBase=str_detect(x$Kmer[1], altBase))}))
standDiffs$MedDiff <- standDiffs$stdMed - standDiffs$altMed
standDiffs$MeanDiff <- standDiffs$stdMean - standDiffs$altMean

pdf(paste0('signal_shifts.', densBase, '.pdf'))
ggplot(standDiffs[!standDiffs$hasAltBase,]) + geom_point(aes(x=stdMed, y=MedDiff), alpha=0.3) +
    geom_smooth(aes(x=stdMed, y=MedDiff), color='red', method = "lm", formula = y ~ x + I(x^2)) + theme_bw()
ggplot(standDiffs[!standDiffs$hasAltBase,]) + geom_point(aes(x=stdMean, y=MeanDiff), alpha=0.3) +
    geom_smooth(aes(x=stdMean, y=MeanDiff), color='red', method = "lm", formula = y ~ x + I(x^2)) + theme_bw()
ggplot(standDiffs) + geom_point(aes(x=stdMed, y=MedDiff, color=hasAltBase), alpha=0.3) +
    geom_smooth(aes(x=stdMed, y=MedDiff, color=hasAltBase), method = "lm", formula = y ~ x + I(x^2)) +  theme_bw()
ggplot(standDiffs) + geom_point(aes(x=stdMean, y=MeanDiff, color=hasAltBase), alpha=0.3) +
    geom_smooth(aes(x=stdMean, y=MeanDiff, color=hasAltBase), method = "lm", formula = y ~ x + I(x^2)) +  theme_bw()
foo <- dev.off()
