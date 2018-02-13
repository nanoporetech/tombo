plotROC <- function(rocDat){
    print(ggplot(rocDat) + geom_abline(slope=1, intercept=0) +
          geom_path(aes(x=FP, y=TP, color=Comparison)) + theme_bw() +
          xlab('False Positive Rate') + ylab('True Positive Rate'))
    print(ggplot(rocDat) +
          geom_path(aes(x=Precision, y=TP, color=Comparison)) + theme_bw() +
          xlab('Precision') + ylab('Recall'))
}
