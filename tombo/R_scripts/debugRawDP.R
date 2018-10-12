minDPScore <- -250

plotRawDP <- function(zDat, fwdDat, tbDat, sigDat){
    for(reg in unique(zDat$Region)){
        tbReg <- tbDat[tbDat$Region == reg,]

        regZDat <- zDat[zDat$Region == reg,]
        regZDat$Score[regZDat$Score < minDPScore] <- minDPScore
        zP <- ggplot(regZDat) +
            geom_tile(aes(x=EventPos, y=SeqPos, fill=Score)) +
            scale_fill_gradient2(
                high='#67001f', mid='#ffffbf', low='#1a1a1a',
                midpoint=mean(range(regZDat$Score))) +
            geom_line(aes(x=EventPos, y=SeqPos),
                  data=tbReg, color='steelblue') +
            theme_minimal() + ylab('Genomic Sequence') +
            theme(axis.title.x=element_blank(), axis.text.x=element_blank()) +
            ggtitle(tbDat$Region[1]) +
            xlim(min(sigDat$Pos) - 1, max(sigDat$Pos) + 1)

        regFwdDat <- fwdDat[fwdDat$Region == reg,]
        regFwdDat$Score[regFwdDat$Score < minDPScore] <- minDPScore
        fP <- ggplot(regFwdDat) +
            geom_tile(aes(x=EventPos, y=SeqPos, fill=Score)) +
            scale_fill_gradient2(
                high='#67001f', mid='#ffffbf', low='#1a1a1a',
                midpoint=mean(range(regFwdDat$Score))) +
            geom_line(aes(x=EventPos, y=SeqPos),
                      data=tbReg, color='steelblue') +
            theme_minimal() + ylab('Genomic Sequence') +
            theme(axis.title.x=element_blank(), axis.text.x=element_blank()) +
            ggtitle(tbDat$Region[1]) +
            xlim(min(sigDat$Pos) - 1, max(sigDat$Pos) + 1)

        sP <- ggplot(sigDat) + geom_line(aes(x=Pos, y=Signal, color=0)) +
            theme_minimal() + xlab('Position') +
            xlim(min(sigDat$Pos) - 1, max(sigDat$Pos) + 1)

        print(plot_grid(zP, sP, align='v', ncol=1, rel_heights=c(5,1)))
        print(plot_grid(fP, sP, align='v', ncol=1, rel_heights=c(5,1)))
    }
}
