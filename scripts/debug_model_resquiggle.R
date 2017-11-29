library(dplyr)
library(ggplot2)
library(cowplot)

dat <- read.table('debug_signal_space.window_z_scores.txt', header=TRUE)
tbDat <- read.table('debug_signal_space.window_traceback.txt', header=TRUE)
maxPathDat <- read.table('debug_signal_space.window_max_path.txt', header=TRUE)
sigMaxPathDat <- read.table(
    'debug_signal_space.window_signal_max_path.txt', header=TRUE)
origPathDat <- read.table('debug_signal_space.window_orig_path.txt', header=TRUE)
switchDat <- read.table('debug_signal_space.window_switch_points.txt', header=TRUE)
sig <- read.table('debug_signal_space.signal.txt', header=TRUE)

diagDat <- read.table('debug_signal_space.window_last_diag.txt', header=TRUE)
diagDat$LastDiagCount[diagDat$LastDiagCount > 3] <- 3
diagDat$LastDiagCount <- factor(diagDat$LastDiagCount)

pdf('debug_signal_space.window.pdf', width=11)
for(reg_i in unique(dat$Region)){
    for(iter_i in 0:2){
        print(iter_i)
        if(sum(dat$Region == reg_i && dat$Iteration == iter_i) == 0){ next }
        print('.')
    regDat <- dat %>% filter(Region==reg_i, Iteration==iter_i)
    regPath <- tbDat %>% filter(Region==reg_i, Iteration==iter_i)
    regDiag <- diagDat %>% filter(Region==reg_i, Iteration==iter_i)
    zMean <- mean(regDat$ZScore)
    regPath$pathVal <- regPath$pathVal - (zMean * (regPath$SignalPos + 1))
    regPath$pathVal <- regPath$pathVal * max(regPath$SignalPos + 1) /
        ((regPath$SignalPos + 1) * max(regPath$pathVal))

    pCols <- c('Original Path'='#ffffbf',
               'Signal-based\nTraceback Path'="#f16913",
               'Max Probability\nTraceback Path'="#cb181d")
    zscrP <- ggplot(regDat) +
        geom_tile(aes(x=SignalPos, y=BasePos, fill=ZScore)) +
        scale_fill_gradient2(
            low="#9970ab", mid='#f7f7f7', high="#00441b",
            midpoint=mean(range(regDat$ZScore, na.rm=TRUE)),
            name='Lower Tail\nZ-Score') + ylab('Base') +
        geom_line(aes(x=SignalPos, y=BasePos), color=pCols[1], size=0.3,
                  data=origPathDat %>% filter(Region==reg_i, Iteration==iter_i)) +
        geom_line(aes(x=SignalPos, y=BasePos), color=pCols[2], size=0.5,
                  data=sigMaxPathDat %>% filter(Region==reg_i, Iteration==iter_i)) +
        geom_line(aes(x=SignalPos, y=BasePos), color=pCols[3], size=0.5,
                  data=maxPathDat %>% filter(Region==reg_i, Iteration==iter_i)) +
        geom_point(aes(x=SignalPos, y=BasePos + 1), size=0.1,
                   data=switchDat %>% filter(Region==reg_i, Iteration==iter_i)) +
        theme_minimal() + theme(axis.title.x=element_blank(),
                                axis.text.x=element_blank())

    diagP <- ggplot(regDiag) +
        geom_tile(aes(x=SignalPos, y=BasePos, fill=LastDiagCount)) +
        geom_line(aes(x=SignalPos, y=BasePos), color=pCols[1], size=0.3,
                  data=origPathDat %>% filter(Region==reg_i, Iteration==iter_i)) +
        geom_line(aes(x=SignalPos, y=BasePos), color=pCols[2], size=0.5,
                  data=sigMaxPathDat %>% filter(Region==reg_i, Iteration==iter_i)) +
        geom_line(aes(x=SignalPos, y=BasePos), color=pCols[3], size=0.5,
                  data=maxPathDat %>% filter(Region==reg_i, Iteration==iter_i)) + ylab('Base') +
        geom_point(aes(x=SignalPos, y=BasePos + 1), size=0.1,
                   data=switchDat %>% filter(Region==reg_i, Iteration==iter_i)) +
        theme_minimal() + theme(axis.title.x=element_blank(),
                                axis.text.x=element_blank())

    pathP <- ggplot(regPath) +
        geom_tile(aes(x=SignalPos, y=BasePos, fill=pathVal)) +
        geom_line(aes(x=SignalPos, y=BasePos, color=names(pCols)[1]), size=0.3,
                  data=origPathDat %>% filter(Region==reg_i, Iteration==iter_i)) +
        geom_line(aes(x=SignalPos, y=BasePos, color=names(pCols)[2]), size=0.5,
                  data=sigMaxPathDat %>% filter(Region==reg_i, Iteration==iter_i)) +
        geom_line(aes(x=SignalPos, y=BasePos, color=names(pCols)[3]), size=0.5,
                  data=maxPathDat %>% filter(Region==reg_i, Iteration==iter_i)) +
        scale_fill_gradient2(
            low="#9970ab", mid='#f7f7f7', high="#00441b",
            midpoint=mean(range(regPath$pathVal, na.rm=TRUE)),
            name='Normalized\nCumulative\nZ-Score') +
        scale_color_manual(values=pCols, name='') +
        ylab('Base') + theme_minimal() + theme(axis.title.x=element_blank(),
                                               axis.text.x=element_blank())
    sigP <- ggplot(sig %>% filter(Region==reg_i, Iteration==iter_i)) +
        geom_path(aes(x=SignalPos + 0.5, y=Signal, color='')) +
        scale_color_manual(values=c("#cb181d"), name='') +
        xlab('Position') + theme_minimal()
    print(plot_grid(zscrP, diagP, pathP, sigP,
                    labels=c(paste0(reg_i, " ", iter_i),"","",""),
                    rel_heights=c(3,3,3,1.5), ncol=1, align='v'))
}}
foo <- dev.off()
