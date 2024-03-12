library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(ggplot2)
library(patchwork)
# library(AnnotationHub)
# library(harmony)
#library(chromVAR)
library(dplyr)#top_n
library(rtracklayer)
set.seed(1234)

integrated <- readRDS("integrated.four.with.cluster.2.4.2.30.rds")#resolution 2.4,2:30, Idents(integrated) is cluster


DE.pig.resol2.4.2.30.findallmarker.onlypos <- FindAllMarkers(integrated,
                         logfc.threshold = 0.25,
                         min.pct = 0.2,
                         only.pos = TRUE)

DE.pig.resol2.4.2.30.findallmarker.onlypos.p0.05 <- subset(DE.pig.resol2.4.2.30.findallmarker.onlypos, p_val_adj < 0.05)

saveRDS(DE.pig.resol2.4.2.30.findallmarker.onlypos.p0.05, file = "DE.pig.resol2.4.2.30.findallmarker.onlypos.p0.05.rds")
DE.pig.resol2.4.2.30.findallmarker.onlypos.p0.05 <- readRDS("DE.pig.resol2.4.2.30.findallmarker.onlypos.p0.05.rds")
write.table(DE.pig.resol2.4.2.30.findallmarker.onlypos.p0.05, file= "DE.pig.resol2.4.2.30.findallmarker.onlypos.p0.05.txt", sep="\t", quote = FALSE)
scp  pyang19@novadtn.its.iastate.edu:/work/abg/pyang19/opt/pig.6798.6800.PBMC.Satija.pipeline.result/pbmc.1x.2x.cellrangeratac1.2.0.wd/strict.filter.new.annotation/strict.filter.new.annotation.common.peaks/commonpeaks.110444.new.directory/DE.pig.resol2.4.2.30.findallmarker.onlypos.p0.05.txt /Users/Pengxin/Desktop/scATAC.manuscript/manuscript.writing/supplementary.table

# write.table(DE.pig.resol2.findallmarker.onlypos, file= "integrated2.findallmarker.only.pos.txt", sep="\t", quote = FALSE)


DE.pig.resol2.4.2.30.findallmarker.pos.neg <- FindAllMarkers(integrated,
                                                      logfc.threshold = 0.25,
                                                      min.pct = 0.2,
                                                      only.pos = FALSE)
#assay = "peaks")
DE.pig.resol2.4.2.30.findallmarker.pos.neg.p0.05 <- subset(DE.pig.resol2.4.2.30.findallmarker.pos.neg, p_val_adj < 0.05)
DE.pig.resol2.4.2.30.findallmarker.pos.neg.p0.05 <- DE.pig.resol2.4.2.30.findallmarker.pos.neg.p0.05[order(DE.pig.resol2.4.2.30.findallmarker.pos.neg.p0.05$p_val_adj),]

saveRDS(DE.pig.resol2.4.2.30.findallmarker.pos.neg.p0.05, file = "DE.pig.resol2.4.2.30.findallmarker.pos.neg.p0.05.rds")
DE.pig.resol2.4.2.30.findallmarker.pos.neg.p0.05 <- readRDS("DE.pig.resol2.4.2.30.findallmarker.pos.neg.p0.05.rds")

# write.table(DE.pig.resol2.findallmarker.pos.neg, file= "DE.pig.resol2.findallmarker.pos.neg.txt", sep="\t", quote = FALSE)
