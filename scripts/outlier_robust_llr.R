library(ggplot2)

scaleFactor <- 4
scaleFactor2 <- 3
scalePower <- 0.2

firstMean <- 0
secondMeans <- seq(0.05, 0.8, 0.05)
constVar <- 0.12

xlims <- c(-2,2)

pdf('llr_test.pdf', width=10)
for (secondMean in secondMeans){
    secondMean <- firstMean + secondMean
    test_vals <- seq(xlims[1], xlims[2], 0.05)

    ref_diffs <- test_vals - firstMean
    alt_diffs <- test_vals - secondMean
    scale_diffs <- test_vals - ((firstMean + secondMean) / 2)
    space_btwn <- (firstMean - secondMean)^2

    dat <- rbind.data.frame(
        data.frame(value=log(dnorm(test_vals, firstMean) / dnorm(test_vals, secondMean)), x=test_vals,
                   type='Log Likelihood Ratio'),
        data.frame(value=exp(-(scale_diffs^2) / (scaleFactor * constVar)) * ((alt_diffs)^2 - (ref_diffs)^2) /
                       (constVar * space_btwn^scalePower * scaleFactor2), x=test_vals,
                   type='Outlier-Robust LLR'),
        data.frame(value=dnorm(test_vals, firstMean, sqrt(constVar)), x=test_vals, type='Canonical Expected\nSignal Level'),
        data.frame(value=dnorm(test_vals, secondMean, sqrt(constVar)), x=test_vals, type='Alternative Expected\nSignal Level')
        )

    print(ggplot(dat) + geom_density(aes(x=x, y=value, fill=type), stat='identity', color='white', size=0, alpha=0.3) +
          theme_bw() + ggtitle(paste('Expected Signal Level Difference:', secondMean)) + ylim(-2,2) + xlim(xlims[1],xlims[2]))
}
foo <- dev.off()
