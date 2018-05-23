boxQuants <- c(0.01,0.25,0.5,0.75,0.99)
quantQuants <- c(0.01,0.1,0.2,0.3,0.4)

plotModelComp <- function(sigDat, quantDat, boxDat, eventDat,
                          baseDat, TitleDat, modelDat, altModelDat=NULL,
                          numModelVals=20){
    pseudoQuants <- seq(1/numModelVals,1-(1/numModelVals),
                        1/numModelVals)
    ## fix 0 baased coordinates passed in
    sigDat$Position <- sigDat$Position + 1
    quantDat$Position <- quantDat$Position + 1
    boxDat$Position <- boxDat$Position + 1
    eventDat$Position <- eventDat$Position + 1
    baseDat$Position <- baseDat$Position + 1
    modelDat$Position <- modelDat$Position + 1
    if(!is.null(altModelDat)){
        altModelDat$Position <- altModelDat$Position + 1
    }
    regions <- sort(c(unique(as.character(sigDat$Region)),
                      unique(as.character(quantDat$Region)),
                      unique(as.character(boxDat$Region)),
                      unique(as.character(eventDat$Region))))
    for(reg_i in regions){
        reg_model_dat <- modelDat[modelDat$Region==reg_i,]
        if(!is.null(altModelDat)){
            reg_alt_model_dat <- altModelDat[modelDat$Region==reg_i,]
        }
        reg_base_dat <- baseDat[baseDat$Region==reg_i,]
        title <- TitleDat[TitleDat$Region==reg_i,'Title']
        if(reg_i %in% sigDat$Region){
            modDensDat <- lapply(split(
                reg_model_dat, paste0(reg_model_dat$Position,
                                      reg_model_dat$Strand)),
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
            if(!is.null(altModelDat)){
                altModDensDat <- lapply(split(
                    reg_alt_model_dat, paste0(reg_alt_model_dat$Position,
                                              reg_alt_model_dat$Strand)),
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
                maxDens <- max(unlist(c(
                    lapply(modDensDat, function(x) x$Position),
                    lapply(altModDensDat, function(x) x$Position))))
            } else {
                maxDens <- max(unlist(lapply(
                    modDensDat, function(x) x$Position)))
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
            reg_sig_dat <- sigDat[sigDat$Region == reg_i,]
            base_pos <- min(reg_sig_dat$Signal)
            p <- ggplot(reg_sig_dat) +
                geom_polygon(aes(x=Position, y=Signal, group=gPos),
                             data=normDensDat, fill='black', alpha=0.4,
                             size=0, show.legend=FALSE)
            if(!is.null(altModelDat)){
                altNormDensDat <- do.call(
                    rbind.data.frame,
                    lapply(altModDensDat, function(posDens){
                        data.frame(Position=(posDens$Position / maxDens) +
                                       posDens$gPos[1],
                                   Signal=posDens$Signal,
                                   Strand=posDens$Strand,
                                   gPos=posDens$gPos,
                                   Group=posDens$Group)
                    }))
                p <- p +
                    geom_polygon(aes(x=Position, y=Signal, group=gPos),
                                 data=altNormDensDat, fill='red', alpha=0.4,
                                 size=0, show.legend=FALSE)
            }
            p <- p + geom_path(
                         aes(x=Position, y=Signal, color=Group, group=Read),
                         alpha=0.3, size=0.2, show.legend=FALSE)
        } else if(reg_i %in% quantDat$Region) {
            reg_model_dat <- do.call(
                rbind.data.frame,
                lapply(split(
                    reg_model_dat, paste0(reg_model_dat$Position,
                                          reg_model_dat$Strand)),
                    function(psDat){
                        do.call(rbind.data.frame,
                                lapply(quantQuants, function(x){
                                    psQuantVals <- qnorm(
                                        c(x,1-x), mean=psDat$Mean[1],
                                        sd=psDat$SD[1])
                                    data.frame(
                                        Position=psDat$Position[1] + 0.4,
                                        Strand=psDat$Strand[1],
                                        Group='Model',
                                        Region=psDat$Region[1],
                                        Lower=psQuantVals[1],
                                        Upper=psQuantVals[2])}))
                    }))
            reg_quant_dat <- rbind.data.frame(
                quantDat[quantDat$Region == reg_i,],
                reg_model_dat)
            base_pos <- min(reg_quant_dat$Lower)
            p <- ggplot(reg_quant_dat) +
                geom_rect(aes(xmin=Position+0.1, xmax=Position + 0.5,
                              ymin=Lower, ymax=Upper, fill=Group),
                          alpha=0.1, show.legend=FALSE) +
                ylab('Signal') + xlab('Position')
        } else if (reg_i %in% boxDat$Region) {
            reg_model_dat <- do.call(
                rbind.data.frame,
                lapply(
                    split(
                        reg_model_dat, paste0(reg_model_dat$Position,
                                              reg_model_dat$Strand)),
                    function(psDat){
                        psBoxVals <- qnorm(
                            boxQuants, mean=psDat$Mean[1], sd=psDat$SD[1])
                        data.frame(Position=psDat$Position[1],
                                   Strand=psDat$Strand[1],
                                   Group='Model',
                                   Region=psDat$Region[1],
                                   SigMin=psBoxVals[1],
                                   Sig25=psBoxVals[2],
                                   SigMed=psBoxVals[3],
                                   Sig75=psBoxVals[4],
                                   SigMax=psBoxVals[5])
                    }))
            reg_box_dat <- rbind.data.frame(
                boxDat[boxDat$Region == reg_i,],
                reg_model_dat)
            base_pos <- min(reg_box_dat$SigMin)
            p <- ggplot(reg_box_dat) +
                geom_boxplot(
                    aes(Position + 0.5, ymin=SigMin, lower=Sig25,
                        middle=SigMed, upper=Sig75, ymax=SigMax,
                        fill=Group), size=0.2, alpha=0.3,
                    stat="identity", show.legend=FALSE) +
                ylab('Signal') + xlab('Position')
        } else {
            reg_model_dat <- do.call(
                rbind.data.frame,
                lapply(
                    split(
                        reg_model_dat, paste0(reg_model_dat$Position,
                                              reg_model_dat$Strand)),
                    function(psDat){
                        psPsuedoEvents <- qnorm(
                            pseudoQuants, mean=psDat$Mean[1], sd=psDat$SD[1])
                        nObs <- length(psPsuedoEvents)
                        data.frame(Position=rep(psDat$Position[1], nObs),
                                   Signal=psPsuedoEvents,
                                   Strand=rep(psDat$Strand[1], nObs),
                                   Group=rep('Model', nObs),
                                   Region=rep(psDat$Region[1], nObs))
                    }))
            reg_event_dat <- rbind.data.frame(
                eventDat[eventDat$Region == reg_i,],
                reg_model_dat)
            base_pos <- min(reg_event_dat$Signal)
            ## convert event data to be back to back densities instead
            ## of side by side violins
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
                      'T'='#CC0000', 'U'='#CC0000', '-'='black', 'N'='black',
                      'Group1'='red', 'Model'='black')) +
              scale_fill_manual(
                  values=c('Group1'='red', 'Model'='black')) +
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
