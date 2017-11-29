plotSigMDS <- function(sigDiffDists, matFn){
    if(!is.na(matFn)){ save(sigDiffDists, file=matFn) }
    fit <- cmdscale(as.dist(sigDiffDists), eig=TRUE, k=2)
    gdat <- as.data.frame(fit$points)
    colnames(gdat) <- c('Coordinate1', 'Coordinate2')
    print(ggplot(gdat, aes(x=Coordinate1, y=Coordinate2)) +
          geom_point(alpha=0.3, size=0.5) + theme_bw())
}
