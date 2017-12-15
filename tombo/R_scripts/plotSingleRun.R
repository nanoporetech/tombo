plotSingleRun <- function(sigDat, quantDat, boxDat, eventDat,
                          baseDat, TitleDat){
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
            base_pos <- min(reg_sig_dat$Signal)
            p <- ggplot(reg_sig_dat) +
                geom_path(aes(x=Position, y=Signal, group=Read),
                          alpha=0.3, size=0.2, show.legend=FALSE)
        } else if(reg_i %in% quantDat$Region) {
            reg_quant_dat <- quantDat[quantDat$Region == reg_i,]
            base_pos <- min(reg_quant_dat$Lower)
            p <- ggplot(reg_quant_dat) +
                geom_rect(aes(xmin=Position, xmax=Position+1,
                              ymin=Lower, ymax=Upper),
                          alpha=0.1, show.legend=FALSE) +
                ylab('Signal')
        } else if(reg_i %in% boxDat$Region) {
            reg_box_dat <- boxDat[boxDat$Region == reg_i,]
            base_pos <- min(reg_box_dat$SigMin)
            p <- ggplot(reg_box_dat) +
                geom_boxplot(
                    aes(Position + 0.5, ymin=SigMin, lower=Sig25,
                        middle=SigMed, upper=Sig75, ymax=SigMax),
                    size=0.2, alpha=0.5, stat="identity",
                    show.legend=FALSE) +
                ylab('Signal') + xlab('Position')
        } else {
            reg_event_dat <- eventDat[eventDat$Region == reg_i,]
            base_pos <- min(reg_event_dat$Signal)
            p <- ggplot(reg_event_dat) +
                geom_violin(aes(
                    x=Position + 0.5, y=Signal, group=Position),
                    size=0, show.legend=FALSE, fill='black') +
                ylab('Signal') + xlab('Position')
        }
        print(p + facet_grid(Strand ~ .) +
              geom_text(
                  aes(x=Position+0.5, y=base_pos, label=Base, color=Base),
                  data=reg_base_dat, hjust=0.5, size=3, show.legend=FALSE,
                  vjust=ifelse(reg_base_dat$Strand=='Forward Strand',0,1),
                  angle=ifelse(reg_base_dat$Strand=='Forward Strand',0,180)) +
              scale_color_manual(
                  values=c('A'='#00CC00', 'C'='#0000CC', 'G'='#FFB300',
                           'T'='#CC0000', 'U'='#CC0000', '-'='black', 'N'='black')) +
              geom_vline(
                  xintercept=min(reg_base_dat$Position):
                  (max(reg_base_dat$Position) + 1),
                  size=0.01) + ggtitle(title) +
              theme_bw() + theme(axis.text.x=element_text(hjust=0),
                                 panel.grid.major.x=element_blank(),
                                 panel.grid.minor.x=element_blank(),
                                 panel.grid.minor.y=element_blank()))
    }}
