plotGroupComp <- function(sigDat, quantDat, boxDat, eventDat,
                          baseDat, TitleDat, QuantWidth){
    ## fix 0 baased coordinates passed in
    sigDat$Position <- sigDat$Position + 1
    quantDat$Position <- quantDat$Position + 1
    boxDat$Position <- boxDat$Position + 1
    eventDat$Position <- eventDat$Position + 1
    baseDat$Position <- baseDat$Position + 1
    regions <- sort(c(unique(as.character(sigDat$Region)),
                      unique(as.character(quantDat$Region)),
                      unique(as.character(boxDat$Region)),
                      unique(as.character(eventDat$Region))))
    for(reg_i in regions){
        reg_base_dat <- baseDat[baseDat$Region==reg_i,]
        title <- TitleDat[TitleDat$Region==reg_i,'Title']
        if(reg_i %in% sigDat$Region){
            reg_sig_dat <- sigDat[sigDat$Region == reg_i,]
            ## randomize so all of one group isn't on top
            s_reg_sig_dat <- split(reg_sig_dat, reg_sig_dat$Read)
            reg_sig_dat <- do.call(rbind.data.frame, sample(s_reg_sig_dat))
            base_pos <- min(reg_sig_dat$Signal)
            p <- ggplot(reg_sig_dat) +
                geom_path(
                    aes(x=Position, y=Signal, color=Group, group=Read),
                    alpha=0.3, size=0.5, show.legend=FALSE)
        } else if(reg_i %in% quantDat$Region) {
            reg_quant_dat <- quantDat[quantDat$Region == reg_i,]
            base_pos <- min(reg_quant_dat$Lower)
            p <- ggplot(reg_quant_dat) +
                geom_rect(aes(xmin=Position, xmax=Position + QuantWidth,
                              ymin=Lower, ymax=Upper, fill=Group),
                          alpha=0.1, show.legend=FALSE) +
                ylab('Signal') + xlab('Position')
        } else if (reg_i %in% boxDat$Region) {
            reg_box_dat <- boxDat[boxDat$Region == reg_i,]
            base_pos <- min(reg_box_dat$SigMin)
            p <- ggplot(reg_box_dat) +
                geom_boxplot(
                    aes(Position + 0.5, ymin=SigMin, lower=Sig25,
                        middle=SigMed, upper=Sig75, ymax=SigMax,
                        fill=Group), size=0.2, alpha=0.3,
                    stat="identity", show.legend=FALSE) +
                ylab('Signal') + xlab('Position')
        } else {
            reg_event_dat <- eventDat[eventDat$Region == reg_i,]
            base_pos <- min(reg_event_dat$Signal)
            ## convert event data to be back to back densities instead
            ## of side by side violin
            g1 <- reg_event_dat$Group[1]
            densDat <- lapply(split(
                reg_event_dat, paste0(reg_event_dat$Group,
                                      reg_event_dat$Position,
                                      reg_event_dat$Strand)),
                              function(pgDat){
                                  pgDens <- density(pgDat$Signal)
                                  nDens <- length(pgDens$x)
                                  data.frame(Position=pgDens$y,
                                             Signal=pgDens$x,
                                             Strand=rep(pgDat$Strand[1], nDens),
                                             gPos=rep(pgDat$Position[1], nDens),
                                             Group=rep(pgDat$Group[1], nDens))
                              })
            maxDens <- max(unlist(lapply(densDat, function(x) x$Position)))
            normDensDat <- do.call(
                rbind.data.frame,
                lapply(densDat, function(posDens){
                    pos <- ifelse(posDens$Group == g1, posDens$Position * -1,
                                  posDens$Position) / (maxDens * 2)
                    data.frame(Position=pos + posDens$gPos[1],
                               Signal=posDens$Signal,
                               Strand=posDens$Strand,
                               gPos=posDens$gPos,
                               Group=posDens$Group)
                }))
            p <- ggplot(normDensDat) +
                geom_polygon(aes(x=Position + 0.5, y=Signal, fill=Group,
                                 group=paste0(Group, gPos)),
                             size=0, show.legend=FALSE) +
                ylab('Signal') + xlab('Position')
        }
        print(p + facet_grid(Strand ~ .) +
              geom_text(
                  aes(x=Position+0.5, y=base_pos, label=Base, color=Base),
                  data=reg_base_dat, hjust=0.5, size=3, show.legend=FALSE,
                  vjust=ifelse(reg_base_dat$Strand=='Forward Strand',0,1),
                  angle=ifelse(reg_base_dat$Strand=='Forward Strand',0,180)) +
              scale_color_manual(
                  values=c(
                      'A'='#00CC00', 'C'='#0000CC', 'G'='#FFB300',
                      'T'='#CC0000', '-'='black', 'N'='black',
                      'Group1'='red', 'Group2'='black')) +
              scale_fill_manual(
                  values=c('Group1'='red', 'Group2'='black')) +
              geom_vline(
                  xintercept=
                      min(reg_base_dat$Position):
                  (max(reg_base_dat$Position) + 1),
                  size=0.01) +
              ggtitle(title) +
              theme_bw() + theme(axis.text.x=element_text(hjust=0),
                                 panel.grid.major.x=element_blank(),
                                 panel.grid.minor.x=element_blank(),
                                 panel.grid.minor.y=element_blank()))
    }}
