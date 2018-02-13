library(ggplot2)

a <- 0.5
b <- 2

mDat <- do.call(rbind.data.frame, lapply(1:20, function(i)
    data.frame(postProb=sapply(0:i, function(j) (a+i-j) / (a+i+b)), coverage=i, notMod=0:i)))

pdf('test_priors.pdf', width=8)
ggplot(mDat) + geom_tile(aes(x=coverage, y=notMod, fill=postProb)) +
    theme_bw() + scale_fill_gradient2(low='#67001f', mid='#f7f7f7', high='#2166ac', midpoint=0.5)
dev.off()
