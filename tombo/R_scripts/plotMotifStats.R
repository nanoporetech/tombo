numModelVals <- 20
pseudoQuants <- seq(1/numModelVals,1-(1/numModelVals),
                    1/numModelVals)

plotMotifStats <- function(PlotDat, BaseDat, StatsDat,
                           ModelDat, AltModelDat=NULL){
    ylim <- 4
    regions <- unique(PlotDat$Region)
    midReg <- regions[(length(regions) + 1) / 2]
    ## "ggplot_gtable(ggplot_build(" writes a page so sink it to dev null
    pdf('/dev/null')
    ps <- lapply(regions, function(region){
        rBaseDat <- BaseDat[BaseDat$Region==region,]
        rPlotDat <- PlotDat[PlotDat$Region==region,]
        ## randomize so all of one group isn't on top
        sRPlotDat <- split(rPlotDat, rPlotDat$Read)
        rPlotDat <- do.call(rbind.data.frame, sample(sRPlotDat))
        # plot signal and distribution
        p <- ggplot(rPlotDat)
        if(! is.null(ModelDat)){
            rModelDat <- ModelDat[ModelDat$Region==region,]
            modDensDat <- lapply(split(
                rModelDat, paste0(rModelDat$Position,
                                  rModelDat$Strand)),
                function(psDat){
                    psDens <- density(qnorm(
                        pseudoQuants, mean=psDat$Mean[1], sd=psDat$SD[1]))
                    nDens <- length(psDens$x)
                    data.frame(Position=psDens$y,
                               Signal=psDens$x,
                               Strand=rep(psDat$Strand[1], nDens),
                               gPos=rep(psDat$Position[1], nDens),
                               Group=rep(psDat$Region[1], nDens))
                })
            maxDens <- max(unlist(lapply(modDensDat, function(x) x$Position)))
            if(! is.null(AltModelDat)) {
                rAltModelDat <- AltModelDat[AltModelDat$Region==region,]
                altModDensDat <- lapply(split(
                    rAltModelDat, paste0(rAltModelDat$Position,
                                         rAltModelDat$Strand)),
                    function(psDat){
                        psDens <- density(qnorm(
                            pseudoQuants, mean=psDat$Mean[1], sd=psDat$SD[1]))
                        nDens <- length(psDens$x)
                        data.frame(Position=psDens$y,
                                   Signal=psDens$x,
                                   Strand=rep(psDat$Strand[1], nDens),
                                   gPos=rep(psDat$Position[1], nDens),
                                   Group=rep(psDat$Region[1], nDens))
                    })
                altMaxDens <- max(unlist(lapply(
                    altModDensDat, function(x) x$Position)))
                maxDens <- max(altMaxDens, maxDens)
            }
            normDensDat <- do.call(
                rbind.data.frame,
                lapply(modDensDat, function(posDens){
                    data.frame(Position=(posDens$Position / maxDens) +
                                   posDens$gPos[1],
                               Signal=posDens$Signal,
                               Strand=posDens$Strand,
                               gPos=posDens$gPos,
                               Group=posDens$Group)
                }))
            p <- p + geom_polygon(aes(x=Position, y=Signal, group=gPos),
                                  data=normDensDat, fill='black', alpha=0.5,
                                  size=0, show.legend=FALSE)
            if(! is.null(AltModelDat)) {
                normAltDensDat <- do.call(
                    rbind.data.frame,
                    lapply(altModDensDat, function(posDens){
                        data.frame(Position=(posDens$Position / maxDens) +
                                       posDens$gPos[1],
                                   Signal=posDens$Signal,
                                   Strand=posDens$Strand,
                                   gPos=posDens$gPos,
                                   Group=posDens$Group)
                    }))
                p <- p + geom_polygon(aes(x=Position, y=Signal, group=gPos),
                                      data=normAltDensDat, fill='red', alpha=0.5,
                                      size=0, show.legend=FALSE)
            }
        }
        p <- p + geom_path(aes(x=Position, y=Signal, color=Group, group=Read),
                           alpha=0.5, size=0.1, show.legend=FALSE)
        p <- p + geom_text(aes(x=Position+0.5, y=-ylim,
                               label=Base, color=Base),
                           data=rBaseDat,
                           hjust=0.5, vjust=0, size=3, show.legend=FALSE) +
            scale_color_manual(
                values=c(
                    'A'='#00CC00', 'C'='#0000CC', 'G'='#FFB300',
                    'T'='#CC0000', '-'='black', 'N'='black',
                    'Group1'='red', 'Group2'='black')) +
            geom_vline(
                xintercept=min(rBaseDat$Position):
                (max(rBaseDat$Position) + 1), size=0.01) +
           scale_x_continuous(expand=c(0,0)) +
           coord_cartesian(ylim=c(-ylim, ylim)) +
           theme_bw() +
           theme(axis.text.x=element_blank(),
                 axis.text.y=element_blank(),
                 axis.title.x=element_blank(),
                 axis.ticks.x=element_blank(),
                 axis.ticks.y=element_blank(),
                 plot.margin=margin(0,0,0,0,'lines'),
                 panel.grid.major.x=element_blank(),
                 panel.grid.minor.x=element_blank(),
                 panel.grid.major.y=element_blank(),
                 panel.grid.minor.y=element_blank())
        if(region != midReg){
            p <- p + theme(axis.title.y=element_blank())
        }
        return(ggplot_gtable(ggplot_build(p)))
    })

    maxStat <- max(StatsDat$Stat)
    if(maxStat <= 1){ tickVals <- c(0,0.2,0.4,0.6,0.8,1)
    } else if(maxStat < 10){ tickVals <- seq(0,10,by=2)
    } else { tickVals <- seq(0,100,by=5) }
    ps[[length(ps) + 1]] <- ggplot_gtable(ggplot_build(
        ggplot(StatsDat) +
        geom_violin(aes(
            x=Position+0.5, y=Stat,
            group=cut_width(Position, 0.9999)), size=0.1, fill='black') +
        scale_x_continuous(expand=c(0,0)) +
        scale_y_continuous(breaks=tickVals) +
        theme_bw() +
        theme(axis.text.x=element_blank(),
              axis.ticks.x=element_blank(),
              axis.title.x=element_blank(),
              panel.grid.minor.y=element_blank()) +
        ylab('Fraction Modified')))
    maxWidth <- do.call(grid::unit.pmax,
                        sapply(ps, function(x) x$widths[2:3]))
    ps <- lapply(ps, function(p){
        p$widths[2:3] <- maxWidth
        return(p)})
    # close dev null sink
    foo <- dev.off()
    do.call(
        grid.arrange,
        c(ps, list(ncol=1, heights=c(rep(1, length(regions)), 3))))
    ##print(do.call(
    ##    plot_grid,
    ##    c(ps, list(ncol=1, align='v',
    ##               rel_heights=c(rep(1, length(regions)), 3)))))
}
