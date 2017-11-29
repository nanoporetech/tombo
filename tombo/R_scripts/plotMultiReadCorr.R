plotMultiReadCorr <- function(OldSegDat, NewSegDat, SigDat){
    for(regId in unique(OldSegDat$Region)){
        rOldSegDat <- OldSegDat[OldSegDat$Region == regId,]
        rNewSegDat <- NewSegDat[NewSegDat$Region == regId,]
        rSigDat <- SigDat[SigDat$Region == regId,]

        regCenter <- median(SigDat$Position) + 1
        sig_max <- max(rSigDat$Signal)
        sig_min <- min(rSigDat$Signal)
        sig_range <- sig_max - sig_min
        print(ggplot(rSigDat) +
              geom_line(aes(x=Position, y=Signal), size=0.2) +
              geom_vline(aes(xintercept=regCenter), color='red') +
              geom_segment(
                  data=rOldSegDat,
                  aes(x=Position, xend=Position, y=sig_max,
                      yend=sig_max - (sig_range * 0.3), color=IsDel)) +
              geom_text(
                  data=rOldSegDat,
                  aes(x=Position, y=sig_max, label=Base,
                      color=IsMismatch), hjust=0, vjust=1, size=3) +
              geom_segment(
                  data=rNewSegDat,
                  aes(x=Position, xend=Position, y=sig_min,
                      yend=sig_min + (sig_range * 0.3), color=IsIns)) +
              geom_text(
                  data=rNewSegDat,
                  aes(x=Position, y=sig_min, label=Base),
                  hjust=0, vjust=0, size=3) +
              facet_grid(Read ~ .) +
              scale_color_manual(
                  values=c('FALSE'='black', 'TRUE'='red')) +
              theme_bw() + theme(legend.position='none',
                                 panel.grid.minor.y=element_blank()))
    }}

plotMultiReadCorrNoOrig <- function(NewSegDat, SigDat){
    for(regId in unique(NewSegDat$Region)){
        rNewSegDat <- NewSegDat[NewSegDat$Region == regId,]
        rSigDat <- SigDat[SigDat$Region == regId,]

        regCenter <- median(SigDat$Position) + 1
        sig_max <- max(rSigDat$Signal)
        sig_min <- min(rSigDat$Signal)
        sig_range <- sig_max - sig_min
        print(ggplot(rSigDat) +
              geom_line(aes(x=Position, y=Signal), size=0.2) +
              geom_vline(aes(xintercept=regCenter), color='red') +
              geom_segment(
                  data=rNewSegDat,
                  aes(x=Position, xend=Position, y=sig_min,
                      yend=sig_min + (sig_range * 0.3), color=IsIns)) +
              geom_text(
                  data=rNewSegDat,
                  aes(x=Position, y=sig_min, label=Base, color=Base),
                  hjust=0, vjust=0, size=3) +
              facet_grid(Read ~ .) +
              scale_color_manual(
                  values=c(
                      'A'='#00CC00', 'C'='#0000CC', 'G'='#FFB300',
                      'T'='#CC0000', '-'='black', 'N'='black',
                      'FALSE'='black', 'TRUE'='red')) +
              theme_bw() + theme(legend.position='none',
                                 panel.grid.minor.y=element_blank()))
    }}
