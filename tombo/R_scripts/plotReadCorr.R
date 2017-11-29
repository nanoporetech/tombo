plotReadCorr <- function(OldSegDat, NewSegDat, SigDat, DiffDat){
    OldSegDat <- cbind.data.frame(OldSegDat, Type='Signal')
    NewSegDat <- cbind.data.frame(NewSegDat, Type='Signal')

    for(readId in unique(SigDat$Read)){
    rOldSegDat <- OldSegDat[OldSegDat$Read == readId,]
    rNewSegDat <- NewSegDat[NewSegDat$Read == readId,]
    rSigDat <- SigDat[SigDat$Read == readId,]
    rDiffDat <- DiffDat[DiffDat$Read == readId,]
    rSigDiffDat <- rbind.data.frame(
        cbind.data.frame(rSigDat, Type='Signal'),
        cbind.data.frame(rDiffDat, Type='Running Difference'))

    sig_max <- max(rSigDat$Signal)
    sig_min <- min(rSigDat$Signal)
    sig_range <- sig_max - sig_min
    p <- ggplot(rSigDiffDat) +
        geom_line(aes(x=Position, y=Signal), size=0.2) +
        facet_grid(Type ~ ., scales='free') +
        scale_color_manual(values=c('FALSE'='black', 'TRUE'='red')) +
        theme_bw() + theme(legend.position='none',
                           axis.title.y=element_blank())
    if(nrow(rOldSegDat) > 0){
        p <- p +
        geom_segment(
            data=rOldSegDat,
            aes(x=Position, xend=Position, y=sig_max,
                yend=sig_max - (sig_range * 0.3), color=IsDel)) +
        geom_text(
            data=rOldSegDat,
            aes(x=Position, y=sig_max, label=Base, color=IsMismatch),
            hjust=0, vjust=1, size=5)
    }
    if(nrow(rNewSegDat) > 0){
        p <- p +
        geom_segment(
            data=rNewSegDat,
            aes(x=Position, xend=Position, y=sig_min,
                yend=sig_min + (sig_range * 0.3), color=IsIns)) +
        geom_text(
            data=rNewSegDat,
            aes(x=Position, y=sig_min, label=Base),
            hjust=0, vjust=0, size=5)
    }
    print(p)
}}
