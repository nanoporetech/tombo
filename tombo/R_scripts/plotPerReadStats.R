# set thresholds for plotting tile
pointMaxReads <- 40
pointMaxBases <- 200
textLim <- 150

plotPerReadStats <- function(StatData, OrdData, baseDat, boxCenter, arePvals){
    all_reg_ids <- unique(StatData$Region)
    last_reg_id <- tail(all_reg_ids, 1)
    for(reg_i in all_reg_ids){
        regDat <- StatData[StatData$Region == reg_i,]
        regOrd <- OrdData[OrdData$Region == reg_i,'Read']
        regDat$Read <- factor(regDat$Read, ordered=TRUE, levels=regOrd)
        boxDat <- data.frame(xS=mean(range(regDat$Position))-1,
                             xE=mean(range(regDat$Position))+1,
                             yS=0.5, yE=length(unique(regDat$Read))+0.5)
        regDat <- regDat[!is.na(regDat$Stats),]
        reg_base_dat <- baseDat[baseDat$Region==reg_i,]
        p <- ggplot(regDat)
        ## add stat values
        if(arePvals){
            if(length(unique(regDat$Read)) > pointMaxReads ||
               length(unique(regDat$Position)) > pointMaxBases){
                p <- p + geom_tile(aes(x=Position, y=Read, fill=Stats))
            } else {
                p <- p + geom_point(
                             aes(x=Position, y=Read, fill=Stats),
                             stroke=0, color='#969696', size=5, shape=21)
            }
            p <- p +
                scale_fill_gradient2(
                    low="#ffffff", mid="#ffffff", high='#cb181d',
                    midpoint=0.1, name='-Log10\nP-Value')
        } else {
            if(length(unique(regDat$Read)) > pointMaxReads ||
               length(unique(regDat$Position)) > pointMaxBases){
                p <- p + geom_tile(aes(x=Position, y=Read, fill=Stats))
            } else {
                p <- p + geom_point(aes(x=Position, y=Read, fill=Stats),
                                    stroke=0, color='#969696', size=5, shape=21,
                                    alpha=0.5)
            }
            lhRatioMax <- max(abs(regDat$Stats))
            breaks <- seq(-lhRatioMax, lhRatioMax, length.out=5)
            p <- p +
                scale_fill_gradient2(
                    low="#b2182b", mid='#ffffff', high='#4d4d4d', midpoint=0,
                    name='Log\nLikelihood\nRatio\n', breaks=breaks,
                    labels=c('Alternative\nBase', breaks[2], '0',
                             breaks[4], 'Standard\nBase'))
        }
        ## add either text or tile-like base data
        if(nrow(reg_base_dat) > 0){
            if(nrow(reg_base_dat) < textLim){
                p <- p + geom_text(
                             aes(x=Position, y=0.5, label=Base, color=Base),
                             data=reg_base_dat, hjust=0.5, size=3,
                             show.legend=FALSE, vjust=1.2, angle=0)
            } else {
                p <- p + geom_point(
                             aes(x=Position, y=0, color=Base),
                             data=reg_base_dat, show.legend=FALSE, shape=15)
            }
        }
        if(boxCenter){
            p <- p + geom_rect(aes(xmin=xS, xmax=xE, ymin=yS, ymax=yE),
                               data=boxDat, fill=NA, color='black', size=0.2)
        }
        p <- p +
            scale_color_manual(
                values=c('A'='#00CC00', 'C'='#0000CC', 'G'='#FFB300',
                         'T'='#CC0000', '-'='black', 'N'='black')) +
            scale_x_continuous(expand=c(0, 0)) + ylab('Reads') +
            theme_minimal() +
            theme(panel.grid=element_blank(), axis.text.y=element_blank(),
                  axis.text.x=element_blank(), legend.text.align=0.5,
                  legend.title.align=0.5)
        ## need to set clip to off so bases aren't cut off
        gt <- ggplot_gtable(ggplot_build(p))
        gt$layout$clip[gt$layout$name == "panel"] <- "off"
        grid::grid.draw(gt)
        if(reg_i != last_reg_id){ grid::grid.newpage() }
    }
}
