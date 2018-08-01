library(dplyr)
library(ggplot2)

## set _DEBUG_BANDWIDTH = True (or _DEBUG_START_BANDWIDTH) in resquiggle.py
##   tombo resquiggle bandwidth_test_reads/ genome.fasta --processes 8 > band_boundary_tuning.txt


dat <- read.table('band_boundary_tuning.txt', header=TRUE)

std_bw <- names(sort(table(dat$bandwidth), decreasing=TRUE))[1]
print(std_bw)

pdf('band_boundary_tuning.pdf', width=11)
ggplot(dat %>% filter(bandwidth == std_bw)) +
    geom_density(aes(x=min_bw_edge_buffer), fill='black') + theme_bw()
foo <- dev.off()


## print bandwidth values for each percentile of reads included in that bandwidth
print(((as.numeric(std_bw) / 2) -
 quantile(dat$min_bw_edge_buffer,
          c(0.01, 0.05, 0.1, 0.2, 0.25, 0.5))) * 2)

## inverse for selected bw
sel_bw <- 200
1 - ecdf(dat$min_bw_edge_buffer)((as.numeric(std_bw) / 2) - (sel_bw / 2))




## for start benadwidth debugging
## print bandwidth values for each percentile of reads included in that bandwidth
print(
    quantile(dat$min_bw_edge_buffer,
             1 - c(0.01, 0.05, 0.1, 0.2, 0.25, 0.5)))

## inverse for selected bw
sel_bw <- 750
ecdf(dat$min_bw_edge_buffer)(sel_bw)
