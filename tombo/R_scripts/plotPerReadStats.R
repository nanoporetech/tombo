ngpValMax <- 20
lhRatioMax <- 6

plotPerReadStats <- function(StatData, OrdData, baseDat, boxCenter, arePvals){
    all_reg_ids <- unique(StatData$Region)
    last_reg_id <- tail(all_reg_ids, 1)
    for(reg_i in all_reg_ids){
        regDat <- StatData[StatData$Region == reg_i,]
        regOrd <- OrdData[OrdData$Region == reg_i,'Read']
        if(arePvals){
            regDat$Stats[regDat$Stats > ngpValMax] <- ngpValMax
        } else {
            regDat$Stats[regDat$Stats > lhRatioMax] <- lhRatioMax
            regDat$Stats[regDat$Stats < -lhRatioMax] <- -lhRatioMax
        }
        regDat$Read <- factor(regDat$Read, ordered=TRUE, levels=regOrd)
        boxDat <- data.frame(xS=mean(range(regDat$Position))-1.5,
                             xE=mean(range(regDat$Position))+0.5,
                             yS=0.5, yE=length(unique(regDat$Read))+0.5)
        regDat <- regDat[!is.na(regDat$Stats),]
        reg_base_dat <- baseDat[baseDat$Region==reg_i,]
        p <- ggplot(regDat)
        if(arePvals){
            p <- p + geom_point(aes(x=Position, y=Read, fill=Stats),
                                stroke=0, color='#969696', size=5, shape=21) +
                scale_fill_gradient(low="#fff7ec", high='#7f0000',
                                    name='-Log10\nP-Value')
        } else {
            p <- p + geom_point(aes(x=Position, y=Read, fill=Stats),
                                stroke=0, color='#969696', size=5, shape=21) +
                scale_fill_gradient2(
                    low="#b2182b", mid='#ffffff', high='#4d4d4d', midpoint=0,
                    name='Log\nLikelihood\nRatio\n', breaks=c(-6,-3,0,3,6),
                    labels=c('Alternative\nBase', '-3','0','3', 'Standard\nBase'))
        }
        if(nrow(reg_base_dat) > 0){
            p <- p + geom_text(
                         aes(x=Position, y=0.5, label=Base, color=Base),
                         data=reg_base_dat, hjust=0.5, size=3, show.legend=FALSE,
                         vjust=1.2, angle=0)
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
