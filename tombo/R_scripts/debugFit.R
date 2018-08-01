plotFit <- function(fitDat, readBw){
    for(reg in unique(fitDat$Region)){
        regDat <- fitDat[fitDat$Region == reg,]
        p1 <- ggplot(regDat) +
            geom_line(aes(x=EventPos, y=EventScore), size=0.1) +
            ggtitle(regDat$Region) + theme_minimal() +
            #geom_hline(aes(yintercept=yint), data=data.frame(yint=c(-4)),
            #           color='red') +
            #ylim(c(-15,5)) +
            theme(axis.title.x=element_blank())
        p2 <- ggplot(regDat) +
            geom_hline(aes(yintercept=readBw / 2), color='blue') +
            geom_line(aes(x=EventPos, y=BandPos), size=0.5) +
            theme_minimal() + ylim(c(0, readBw)) +
            theme(axis.title.x=element_blank())
        p3 <- ggplot(regDat) +
            geom_line(aes(x=EventPos, y=ModelMean), size=0.1, color='blue') +
            geom_line(aes(x=EventPos, y=EventMean), size=0.1) +
            theme_minimal()
        print(plot_grid(p1, p2, p3, align='v', ncol=1, rel_heights=c(7,6,7)))
    }
}
