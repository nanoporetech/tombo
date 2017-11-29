plotKmerDist <- function(dat, baseDat, saveDatFn, dontPlot){
    if (!is.na(saveDatFn)){
        save(dat, file=saveDatFn)
    }
    if (dontPlot){ return() }
    mainP <- ggplot(dat) +
        ##geom_boxplot(aes(x=Kmer, y=Signal, color=Base)) +
        ##scale_color_manual(
        ##    values=c('#00CC00', '#0000CC', '#FFB300', '#CC0000')) +
        geom_violin(aes(x=Kmer, y=Signal, fill=Base), size=0) +
        scale_fill_manual(
            values=c('#00CC00', '#0000CC', '#FFB300', '#CC0000')) +
        theme_bw() +
        theme(
            axis.text.x=element_text(angle=60, hjust=1, size=8),
            legend.position='bottom',
            panel.grid.major.x=element_line(size=0.01))
    if (is.na(baseDat)){
        print(mainP)
    } else {
        mainL <- get_legend(mainP)
        mainP <- mainP + theme(legend.position='none')
        baseP <- ggplot(baseDat) +
            geom_tile(aes(x=Kmer, y=Position, fill=Base)) +
            scale_fill_manual(
                values=c('A'='#00CC00', 'C'='#0000CC',
                         'G'='#FFB300', 'T'='#CC0000')) +
            theme_bw() +
            theme(
                axis.text.x=element_text(angle=60, hjust=1, size=8),
                legend.position='none')
        mainP <- mainP + theme(axis.text.x=element_blank(),
                               axis.title.x=element_blank())
        if (nchar(as.character(dat$Kmer[1])) > 3){
            mainP <- mainP + theme(axis.text.x=element_blank())
            baseP <- baseP + theme(axis.text.x=element_blank())
        }
        print(plot_grid(plot_grid(mainP, baseP, ncol=1,
                                  rel_heights=c(5,1), align='v'),
                        mainL, ncol=1, rel_heights=c(10,1)))
    }}

plotKmerDistWReadPath <- function(dat, baseDat, saveDatFn, dontPlot){
    if (!is.na(saveDatFn)){
        save(dat, file=saveDatFn)
    }
    if (dontPlot){ return() }
    maxVal <- max(as.numeric(dat$Kmer))
    mainP <- ggplot(dat) +
        #geom_boxplot(aes(x=Kmer, y=Signal, color=Base)) +
        #scale_color_manual(
        #    values=c('#00CC00', '#0000CC', '#FFB300', '#CC0000')) +
        geom_violin(aes(x=Kmer, y=Signal, fill=Base), size=0) +
        scale_fill_manual(
            values=c('#00CC00', '#0000CC', '#FFB300', '#CC0000')) +
        theme_bw() +
        theme(axis.text.x=element_text(angle=60, hjust=1, size=8),
              legend.position='bottom',
              panel.grid.major.x=element_line(size=0.01))
    readP <- ggplot(dat) +
        geom_path(aes(x=Kmer, y=Signal, group=Read),
                  alpha=0.05, size=0.5) +
        theme_bw() +
        theme(axis.text.x=element_text(angle=60, hjust=1, size=8),
              panel.grid.major.x=element_line(size=0.01)) +
        scale_color_manual(
            values=c('#00CC00', '#0000CC', '#FFB300', '#CC0000'))
    if (is.na(baseDat)){
        print(mainP)
        print(readP)
    } else {
    mainL <- get_legend(mainP)
    mainP <- mainP + theme(legend.position='none')
    baseP <- ggplot(baseDat) +
        geom_tile(aes(x=Kmer, y=Position, fill=Base)) +
        scale_fill_manual(
            values=c('A'='#00CC00', 'C'='#0000CC',
                     'G'='#FFB300', 'T'='#CC0000')) +
        theme_bw() + theme(
            axis.text.x=element_text(angle=60, hjust=1, size=8),
            legend.position='none')
    mainP <- mainP + theme(axis.text.x=element_blank(),
                           axis.title.x=element_blank())
    readP <- readP + theme(axis.text.x=element_blank(),
                           axis.title.x=element_blank())
    if (nchar(as.character(dat$Kmer[1])) > 3){
        mainP <- mainP + theme(axis.text.x=element_blank())
        readP <- readP + theme(axis.text.x=element_blank())
        baseP <- baseP + theme(axis.text.x=element_blank())
    }
    print(plot_grid(plot_grid(mainP, baseP, ncol=1,
                              rel_heights=c(5,1), align='v'),
          mainL, ncol=1, rel_heights=c(10,1)))
    print(plot_grid(plot_grid(readP, baseP, ncol=1,
                              rel_heights=c(5,1), align='v'),
          mainL, ncol=1, rel_heights=c(10,1)))
}}
