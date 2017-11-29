library(ggplot2)
library(ggridges)

densDat <- read.table('debug_est_alt.C.density.txt', header=TRUE)
standardDensDat <- read.table('debug_est_standard_ref.density.txt', header=TRUE)
densDat$Sample <- "Alternative"
standardDensDat$Sample <- "Standard"

allDensDat <- rbind.data.frame(densDat, standardDensDat)
sAllDat <- split(allDensDat, allDensDat$Kmer)
sDiffs <- sort(unlist(lapply(sAllDat, function(x)
    weighted.mean(x[x$Sample == 'Alternative','Signal'], x[x$Sample == 'Alternative','Density']) -
    weighted.mean(x[x$Sample == 'Standard','Signal'], x[x$Sample == 'Standard','Density']))))

upDat <- do.call(rbind.data.frame, lapply(names(head(sDiffs, 1)), function(kmer) sAllDat[[kmer]]))
dnDat <- do.call(rbind.data.frame, lapply(names(tail(sDiffs, 1)), function(kmer) sAllDat[[kmer]]))

pdf('alternate_model_estimation.density.C.one_kmer.pdf', width=15)
ggplot(upDat, aes(x=Signal, y=Kmer, height=Density, fill=Sample)) + geom_ridgeline(alpha=0.4, size=0, color='white') +
    scale_fill_discrete(name='Contains\nAlternative\nBase') + theme_ridges()
ggplot(dnDat, aes(x=Signal, y=Kmer, height=Density, fill=Sample)) + geom_ridgeline(alpha=0.4, size=0, color='white') +
    scale_fill_discrete(name='Contains\nAlternative\nBase') + theme_ridges()
foo <- dev.off()



dat <- read.table('debug_est_alt.C.txt', header=TRUE,
                  stringsAsFactors=FALSE, row.names=NULL, sep='\t')

sd_width <- dat[1,'sd_width']
min_frac <- dat[1,'min_frac']

pdf('alternate_model_estimation.C.pdf', width=10)
ggplot(dat) +
    geom_point(aes(x=peaks_diff, color=contains_alt, y=minor_frac), alpha=0.4) +
    geom_rect(aes(xmin=xs, xmax=xe, ymin=ys, ymax=ye),
              data.frame(xs=c(-3,-sd_width),
                         xe=c(3,sd_width),
                         ys=c(0,min_frac), ye=c(min_frac,0.5)), alpha=0.2) +
    xlab('Difference Between Estimated Means') +
    ylab('Alternative Disribution Proportion') +
    scale_color_discrete(name='Contains\nAlternative\nBase') + theme_bw()
foo <- dev.off()
