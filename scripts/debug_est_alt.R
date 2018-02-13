library(ggplot2)
library(ggridges)

densBase <- 'debug_est_alt'
altBase <- 'C'

densDat <- read.table(paste0(densBase, '.alternate_density.txt'), header=TRUE)
standardDensDat <- read.table(paste0(densBase, '.control_density.txt'), header=TRUE)
densDat$Sample <- "Alternative"
standardDensDat$Sample <- "Standard"

allDensDat <- rbind.data.frame(densDat, standardDensDat)
sAllDat <- split(allDensDat, allDensDat$Kmer)
sDiffs <- sort(unlist(lapply(sAllDat, function(x)
    weighted.mean(x[x$Sample == 'Alternative','Signal'], x[x$Sample == 'Alternative','Density']) -
    weighted.mean(x[x$Sample == 'Standard','Signal'], x[x$Sample == 'Standard','Density']))))

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
