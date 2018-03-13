plotROCPerRead <- function(rocDat, denStats){
    print(ggplot(rocDat) + geom_abline(slope=1, intercept=0) +
          geom_path(aes(x=FP, y=TP, color=Comparison)) + theme_bw() +
          xlab('False Positive Rate') + ylab('True Positive Rate'))
    print(ggplot(rocDat) +
          geom_path(aes(x=Precision, y=TP, color=Comparison)) + theme_bw() +
          xlab('Precision') + ylab('Recall'))
    for(modName in names(denStats)){
        print(ggplot(denStats[[modName]]) +
              geom_density(aes(x=stat, fill=motif_match),
                           alpha=0.5, color='white', size=0.01) +
              theme_bw() + ggtitle(modName))
    }
}