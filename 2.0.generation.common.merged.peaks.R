# how many peaks are in common between different datasets
#https://satijalab.org/signac/articles/merging.html
library(Signac)
library(Seurat)
library(GenomicRanges)
library(future)

#generate a set of peaks among 4 datasets across 1x and 2x
peaks.6798.1x <- read.table(
  file = "path.../peaks.bed",
  col.names = c("chr", "start", "end")
)

gr.6798.1x <- makeGRangesFromDataFrame(peaks.6798.1x)

peaks.6798.2x <- read.table(
  file = "path.../peaks.bed",
  col.names = c("chr", "start", "end")
)

gr.6798.2x <- makeGRangesFromDataFrame(peaks.6798.2x)

##merge 6798.1x and 6798.2x--> merged.6798, merge 6800.1x and 6800.2x--> merged.6800, then keep common peaks between merged.6798 and  merged.6800
gr.6798.merged <- reduce(x = c(gr.6798.1x, gr.6798.2x))
# GRanges object with 114843 ranges and 0 metadata columns:
  
peaks.6800.1x <- read.table(
  file = "path.../peaks.bed",
  col.names = c("chr", "start", "end")
)
# dim(peaks.6800.1x)
# [1] 102154      3

gr.6800.1x <- makeGRangesFromDataFrame(peaks.6800.1x)

# length(gr.6800.1x)
# [1] 102154

peaks.6800.2x <- read.table(
  file = "path.../peaks.bed",
  col.names = c("chr", "start", "end")
)
# > dim(peaks.6800.2x)
# [1] 104991      3

gr.6800.2x <- makeGRangesFromDataFrame(peaks.6800.2x)
# > length(gr.6800.2x)
# [1] 104991

gr.6800.merged <- reduce(x = c(gr.6800.1x, gr.6800.2x))
# GRanges object with 114051 ranges and 0 metadata columns:
  
gr.merged.peaks.6798.overlap <- subsetByOverlaps(gr.6798.merged, gr.6800.merged, ignore.strand=TRUE)#return the peaks in the gr.6798.merged that has overlapping with the peaks in gr.6800.merged
#GRanges object with 101674 ranges and 0 metadata columns:

common.merged.peaks.6798.6800 <- reduce(x = c(gr.merged.peaks.6798.overlap, gr.6800.merged))
#GRanges object with 111636 ranges and 0 metadata columns:

saveRDS(common.merged.peaks.6798.6800, file = "common.merged.peaks.6798.6800.rds")

peakwidths.common.merged.peaks.6798.6800 <- width(common.merged.peaks.6798.6800)
hist(peakwidths.common.merged.peaks.6798.6800)

common.merged.peaks.6798.6800 <- common.merged.peaks.6798.6800[peakwidths.common.merged.peaks.6798.6800  < 10000 & peakwidths.common.merged.peaks.6798.6800 > 20]
common.merged.peaks.6798.6800

saveRDS(common.merged.peaks.6798.6800, file = "common.merged.peaks.6798.6800.filtered.rds")#110444

peakwidths.common.merged.peaks.6798.6800.new <- width(common.merged.peaks.6798.6800)


common.merged.peaks.6798.6800 <- readRDS("/work/abg/pyang19/opt/pig.6798.6800.PBMC.Satija.pipeline.result/pbmc.1x.2x.cellrangeratac1.2.0.wd/common.merged.peaks.6798.6800.filtered.rds")
save(common.merged.peaks.6798.6800, file = "common.merged.peaks.6798.6800.RData")
load("common.merged.peaks.6798.6800.RData")







