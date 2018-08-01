minDPScore <- -250

plotDP <- function(dpDat, tbDat){
    for(reg in unique(dpDat$Region)){
        regDat <- dpDat[dpDat$Region == reg,]
        regDat$Score[regDat$Score < minDPScore] <- minDPScore
        mp <- mean(range(regDat$Score[grepl(reg, regDat$Region)]))
        if (grepl('fwd', reg)){
            tbReg <- tbDat[tbDat$Region == reg,]
            print(ggplot(regDat) +
                  geom_tile(aes(x=EventPos, y=SeqPos, fill=Score)) +
                  scale_fill_gradient2(high='#67001f', mid='#ffffbf',
                                       low='#1a1a1a', midpoint=mp) +
                  geom_line(aes(x=EventPos, y=SeqPos),
                            data=tbReg, color='steelblue') +
                  theme_minimal() + xlab('Segmented Signal') +
                  ylab('Genomic Sequence') + ggtitle(tbDat$Region[1]))
        } else {
            print(ggplot(regDat) +
                  geom_tile(aes(x=EventPos, y=SeqPos, fill=Score)) +
                  scale_fill_gradient2(high='#67001f', mid='#ffffbf',
                                       low='#1a1a1a', midpoint=mp) +
                  theme_minimal() + xlab('Segmented Signal') +
                  ylab('Genomic Sequence') + ggtitle(regDat$Region[1]))
        }
    }
}
