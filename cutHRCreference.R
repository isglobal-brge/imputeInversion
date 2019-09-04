#'#################################################################################
#'#################################################################################
#' Create short version of references with SNPs close to inversion regions to speed up computation
#'#################################################################################
#'#################################################################################

## URL to download: https://www.well.ox.ac.uk/~wrayner/tools/1000GP_Phase3_combined.legend.gz

library(data.table)
library(GenomicRanges)
library(scoreInvHap)

dt <- fread("1000GP_Phase3_combined.legend")

SNPGR <- makeGRangesFromDataFrame(dt, start.field = "position", end.field = "position")

seqlevelsStyle(inversionGR) <- "NCBI"
seqlevels(inversionGR)[seqlevels(inversionGR) == "X"] <- "23"

ov <- findOverlaps(SNPGR, inversionGR + 5e5)
selSNPs <- from(ov)

dtFilt <- dt[selSNPs, ]
write.table(dtFilt, file = "1000GP_Phase3_combined.legend.imputeInversion", 
            col.names = TRUE, quote = FALSE, row.names = FALSE)