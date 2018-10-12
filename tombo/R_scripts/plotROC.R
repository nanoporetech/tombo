plotROC <- function(rocDat){
    print(ggplot(rocDat) + geom_abline(slope=1, intercept=0) +
          geom_path(aes(x=FP, y=TP, color=Comparison)) + theme_bw() +
          xlab('False Positive Rate') + ylab('True Positive Rate'))
    print(ggplot(rocDat) +
          geom_path(aes(x=TP, y=Precision, color=Comparison)) + theme_bw() +
          xlab('Recall') + ylab('Precision'))
}
