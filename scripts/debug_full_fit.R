library(ggplot2)
library(cowplot)

bandwidth <- 500

dat <- read.table('debug_event_align.full_fit.txt', header=TRUE)
failedDat <- read.table('debug_event_align.full_failed.txt', header=TRUE, sep="\t")

pdf('debug_event_align.full.pdf', width=15, height=5)
for(reg in unique(dat$Region)){
    regDat <- dat[dat$Region == reg,]
    p1 <- ggplot(regDat) + geom_line(aes(x=EventPos, y=EventScore), size=0.1) +
        theme_minimal() + ylim(c(-15,2)) +
        geom_hline(aes(yintercept=yint), data=data.frame(yint=c(0,0.5)),
                   color='red')
    p2 <- ggplot(regDat) +
        geom_hline(aes(yintercept=bandwidth/2), color='blue') +
        geom_line(aes(x=EventPos, y=BandPos), size=0.1) +
        theme_minimal() + ylim(c(0,bandwidth))
    title <- ggdraw() + draw_label(as.character(
                            failedDat[failedDat$Region == reg,'DidFail']))
    print(plot_grid(title, p1, p2, align='v', ncol=1, rel_heights=c(1,6,6)))
}
foo <- dev.off()
