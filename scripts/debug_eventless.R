library(ggplot2)

dat <- read.table('debug_event_align.txt', header=TRUE)
tbdat <- read.table('debug_event_align.traceback.txt', header=TRUE)

pdf('debug_event_align.pdf', height=4.5, width=6)
for(reg in unique(dat$Region)){
    regDat <- dat[dat$Region == reg,]
    mp <- ifelse(grepl('fwd_end', reg),
                 mean(regDat$Score[grepl('fwd_end', regDat$Region)]), 0)
    if (grepl('fwd', reg)){
        tbReg <- tbdat[tbdat$Region == reg,]
        print(ggplot(regDat) + geom_tile(aes(x=EventPos, y=SeqPos, fill=Score)) +
              scale_fill_gradient2(high='#67001f', mid='#ffffbf', low='#1a1a1a',
                                   midpoint=mp) +
              geom_line(aes(x=EventPos, y=SeqPos),
                        data=tbReg, color='steelblue') +
              theme_minimal() + xlab('Segmented Signal') +
              ylab('Genomic Sequence'))
    } else {
        print(ggplot(regDat) + geom_tile(aes(x=EventPos, y=SeqPos, fill=Score)) +
              scale_fill_gradient2(high='#67001f', mid='#ffffbf', low='#1a1a1a',
                                   midpoint=mp) +
              theme_minimal() + xlab('Segmented Signal') +
              ylab('Genomic Sequence'))
    }
}
foo <- dev.off()
