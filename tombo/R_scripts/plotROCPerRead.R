## should check if the stat is p-values, but this won't effect that
## so not fixing it now
lhRatioMax <- 25

plotROCPerRead <- function(rocDat, denStats){
    print(ggplot(rocDat) + geom_abline(slope=1, intercept=0) +
          geom_path(aes(x=FP, y=TP, color=Comparison)) + theme_bw() +
          xlab('False Positive Rate') + ylab('True Positive Rate'))
    print(ggplot(rocDat) +
          geom_path(aes(x=TP, y=Precision, color=Comparison)) + theme_bw() +
          xlab('Recall') + ylab('Precision'))
    for(modName in names(denStats)){
        denStats[[modName]]$stat[denStats[[modName]]$stat >
                                 lhRatioMax] <- lhRatioMax
        denStats[[modName]]$stat[denStats[[modName]]$stat <
                                 -lhRatioMax] <- -lhRatioMax
        print(ggplot(denStats[[modName]]) +
              geom_density(aes(x=stat, fill=motif_match),
                           alpha=0.5, color='white', size=0.01) +
              theme_bw() + ggtitle(modName) +
              scale_fill_discrete(name="Ground Truth\nModified"))
    }
}
