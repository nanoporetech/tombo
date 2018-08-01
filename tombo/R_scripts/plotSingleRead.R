plotSingleRead <- function(sigDat, vlineDat, hDat, hrDat, stDat){
    sigDat$Position <- sigDat$Position / 1000
    vlineDat$Position <- vlineDat$Position / 1000
    p <- ggplot(sigDat) +
        geom_vline(aes(xintercept=Position),
                   data=vlineDat, color='red', size=5) +
        geom_path(aes(x=Position, y=Signal), size=0.3) +
        theme_bw() +
        theme(axis.text=element_text(size=24),
              axis.title=element_text(size=28)) +
        xlab('Position (1000 raw obs.)') + ylim(-5,5)
    if(! is.null(hDat)){
        hDat$Position <- hDat$Position / 1000
        p <- p + geom_point(aes(x=Position, y=Signal),
                            data=hDat, color='red', size=2)
    }
    if(! is.null(hrDat)){
        hrDat$Position <- hrDat$Position / 1000
        hrDat$PositionEnd <- hrDat$PositionEnd / 1000
        p <- p +
            geom_rect(aes(xmin=Position, xmax=PositionEnd, ymin=-Inf, ymax=Inf),
                      data=hrDat, fill='red', alpha=0.3, color=NA)
    }
    if(is.null(stDat)){
        print(p)
    } else {
        stDat$Position <- stDat$Position / 1000
        p2 <- ggplot(stDat) + geom_line(aes(x=Position, y=Value)) + theme_bw()
        if(! is.null(hrDat)){
            p2 <- p2 +
                geom_rect(aes(xmin=Position, xmax=PositionEnd,
                              ymin=-Inf, ymax=Inf),
                          data=hrDat, fill='red', alpha=0.3, color=NA)
        }
        print(plot_grid(p, p2, ncol=1, align='v'))
    }
}
