#4.cell type annotation new wd on by using scATAC data, in dependent of scRNA dataset

setwd("/work/abg/pyang19/opt/pig.6798.6800.PBMC.Satija.pipeline.result/pbmc.1x.2x.cellrangeratac1.2.0.wd/strict.filter.new.annotation/strict.filter.new.annotation.common.peaks/commonpeaks.110444.new.directory")
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
library(GenomicRanges)
set.seed(1234)


integrated <- readRDS("integrated.four.with.cluster.2.4.2.30.rds")#resolution 2.4,2:30, 3.integration.last.strict.non.batch.common.peaks.new.directory.2.30.R
# An object of class Seurat 
# 110444 features across 17207 samples within 1 assay 
# Active assay: peaks (110444 features, 110272 variable features)
# 3 dimensional reductions calculated: integrated_lsi, umap, lsi
Idents(integrated)#35 Levels: 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 ... 34

gene.activities <- GeneActivity(integrated, extend.upstream = 2000, extend.downstream = 2000) #by default, extend.upstream = 2000, extend.downstream = 0,

dim(gene.activities)#[1] 20035 17207
head(gene.activities)
# 6 x 17207 sparse Matrix of class "dgCMatrix"
# [[ suppressing 34 column names ‘ACGTTAGCAGTCAGCC-1_1’, ‘GCACCTTTCGGATGTT-1_1’, ‘TGGCAATTCTGGGCGT-1_1’ ... ]]
# 
# NA.1    1 1 . 2 2 1 1 1 2 1 2 2 1 1 3 2 1 3 2 1 1 . 2 3 2 2 . 3 1 4 1 1 4 .
# PSMB1   1 1 1 1 2 1 2 2 2 1 2 2 1 1 3 2 1 3 2 1 1 2 2 3 2 2 . 3 1 4 1 1 4 .
# FAM120B . 4 5 1 . 2 . 4 2 1 2 2 . 3 1 1 1 1 . . 1 . 1 1 3 2 2 1 3 1 1 1 4 .
# DLL1    . . . . 2 3 . 1 2 . . 4 . . 2 1 1 1 . 2 . . 1 2 2 . . . . . . . 1 3
# NA.11   1 . . . . . . . . . . . . . . . 1 . . . . . . . . . . . . . . . . .
# ERMARD  3 . . 2 2 8 4 6 8 1 1 5 2 1 3 1 5 4 2 3 . 1 4 7 4 2 3 3 2 3 4 1 4 3

length(unique(rownames(gene.activities)))#20028
# dim(gene.activities.upstream)#

saveRDS(gene.activities, file = "gene.activities.resol2.4.2.30.4000.rds")#resolution 2.4, dim 2:30,20035 17207
gene.activities <- readRDS("gene.activities.resol2.4.2.30.4000.rds")#resolution 2.4,2:30,20035 17207


integrated[['RNA']] <- CreateAssayObject(counts = gene.activities)

DefaultAssay(integrated) <- 'RNA'

integrated <- NormalizeData(
  object = integrated,
  assay = 'RNA',
  normalization.method = 'LogNormalize',#default scale.factor = 10000 which is the same as the one used for RNA assay in scRNA dataset
  scale.factor = median(integrated$nCount_RNA)
)
integrated <- ScaleData(integrated, assay = "RNA")

saveRDS(integrated, file = "integratedscATAC.withRNA.assay2.4.2.30.4000.rds")#resoltion2.4,2:30, +-2000, Active assay: RNA
integrated <- readRDS("integratedscATAC.withRNA.assay2.4.2.30.4000.rds")#resoltion2.4,2:30, +-2000, Active assay: RNA

marker.genes1 <- c("CD19", "CD79A", "PAX5", "IRF4","PRDM1")#B,ASC,
marker.genes2 <- c("CD14", "CD163", "NLRP3","TLR4","CSF1R", "CD86")#mono  
marker.genes3 <- c("FLT3","FCER1A","SLA-DRB1","SLA-DRA")#DC
marker.genes4 <- c("TCF4","IRF8","CD4","CD8B","XBP1","CLEC12A","CD93")#pDC
marker.genes5 <- c("CD3E","TRDC","CD4","CD8A","CD8B","CD2","PRF1")#T
marker.genes6 <- c("KLRB1","KLRK1","CD5","CD6","TYROBP","HCST")#T

p <- FeaturePlot(
  object = integrated,
  features = marker.genes1,
  pt.size = 0.03,
  label = FALSE,
  label.size = 1, ncol = 2, combine = TRUE, coord.fixed = FALSE,
  cols = c("lightgrey", "red"))
p
ggsave("FeaturePlot.integrated.2.4.2.30.4000.all.fragments.allpeaks.marker.genes1.pdf",path = "./plot.manuscript.final")

DefaultAssay(integrated)#RNA
#CD4, CD3E, CD19, CSF1R, CD14
list <-list("CD4","CD3E", "CD19", "CSF1R", "CD14", "CD8A", "CD8B")
for (i in 1:length(list))
{
  p <- FeaturePlot(
    object = integrated,
    features = list[[i]],
    pt.size = 0.03,
    label = FALSE,
    label.size = 1, ncol = 2, combine = TRUE, coord.fixed = FALSE,
    cols = c("lightgrey", "red"))
  p
  ggsave(paste0("FeaturePlot.gene.activity.",list[[i]],".pdf"), path = "./plot.manuscript.final/findallmarkers.celltype.annotation")
}

DefaultAssay(integrated) <- 'peaks'
p <- CoveragePlot(
  object = integrated,
  region = "1-237781345-237790424",#0,3,4,8,11,18,20 findallmarkers only pos
  features = "PAX5",
  annotation = TRUE,
  peaks = TRUE,
  tile = TRUE
) + ggtitle('integrated.2.4.2.30.annotated.withoutscRNA.1-237781345-237790424.PAX5.tile')+ NoLegend() + FontSize(x.title = 10, y.title = 10)
p
ggsave("integrated.2.4.2.30.annotated.withoutscRNA.PAX5.1.237781345.237790424.tile.pdf", path = "./plot.manuscript.final/findallmarkers.celltype.annotation")

Idents(integrated)# 35 clusters
DefaultAssay(integrated) <- 'peaks'
p <- FeaturePlot(
  object = integrated,
  features = "3-57984874-57986859",
  pt.size = 0.03,
  label = FALSE,
  label.size = 1, ncol = 2, combine = TRUE, coord.fixed = FALSE,
  cols = c("lightgrey", "red")) + ggtitle("CD8B:3-57984874-57986859")
p
ggsave(paste0("cluster21.CD8B.3-57984874-57986859.pdf"),path = "./plot.manuscript.final/findallmarkers.celltype.annotation" )

DefaultAssay(integrated) <- 'peaks'
p <- FeaturePlot(
  object = integrated,
  features = "3-57990004-57991024",
  pt.size = 0.03,
  label = FALSE,
  label.size = 1, ncol = 2, combine = TRUE, coord.fixed = FALSE,
  cols = c("lightgrey", "red")) + ggtitle("CD8B:3-57990004-57991024")
p
ggsave(paste0("cluster23.CD8B.3-57990004-57991024.pdf"),path = "./plot.manuscript.final/findallmarkers.celltype.annotation" )


#annotate our scATAC-seq-derived clusters
DefaultAssay(integrated) <- 'RNA'
Idents(integrated)

integrated <- RenameIdents(
  object = integrated,
  '0' = 'B',
  '1' = 'CD2negGD',
  '2' = 'CD4posab',
  '3' = 'B',
  '4' = 'B',
  '5' = 'Monocytes',
  '6' = 'CD4posab',
  '7' = 'Monocytes',
  '8' = 'B',
  '9' = 'Monocytes',
  '10' = 'NK',
  '11' = 'B',
  '12' = 'CD4posab',
  '13' = 'CD8abPOSab',
  '14' = 'CD2negGD',
  '15' = 'unknown',
  '16' = 'CD2negGD',
  '17' = 'CD4posab',
  '18' = 'B',
  '19' = 'CD8abPOSab',
  '20' = 'B',
  '21' = 'CD4posCD8pos',# shall this cluster 21 be annotated as CD4CD8ACD8
  '22' = 'CD2posGD',#diff
  '23' = 'CD8abPOSab',#diff, originally its annotated as CD4posab
  '24' = 'CD8abPOSab',#diff
  '25' = 'CD8aPOSabT_NK',#diff
  '26' = 'cDCs',
  '27' = 'ASC',
  '28' = 'CD2negGD',
  '29' = 'pDCs',
  '30' = 'CD8abPOSab',
  '31' = 'unknown',
  '32' = 'CD8abPOSab',
  '33' = 'CD8aPOSabT_NK',#diff
  '34' = 'CD4posab'
)

integrated <- subset(integrated, idents = "unknown", invert = TRUE)
# 130479 features across 16604 samples within 2 assays 
# An object of class Seurat 
# 130479 features across 16604 samples within 2 assays 
# Active assay: RNA (20035 features, 0 variable features)
# 1 other assay present: peaks
# 3 dimensional reductions calculated: integrated_lsi, umap, lsi

#new annotation below
# [1] cDCs          CD8abPOSab    Monocytes     CD4posab      B            
# [6] CD2negGD      ASC           CD2posGD      pDCs          CD8aPOSabT_NK
# [11] NK            CD4posCD8pos 
# 12 Levels: B CD2negGD CD4posab Monocytes NK CD8abPOSab ... pDCs

integrated$celltype_only_scATAC <- Idents(integrated)

Idents(integrated) <- integrated$seurat_clusters

saveRDS(integrated, file = "integratedscATAC.withRNA.assay2.4.2.30.4000.annotated.without.scRNAintegration.newest.rds")#resoltion2.4,2:30, +-2000, newest annotation, 15 unknown, able to distinguish CD2- and CD2+, albe to identify CD8aPOSabT_NK etc, 24 is CD8
integrated <- readRDS("integratedscATAC.withRNA.assay2.4.2.30.4000.annotated.without.scRNAintegration.newest.rds")#resoltion2.4,2:30, +-2000, newest annotation, 15 unknown, able to distinguish CD2- and CD2+, albe to identify CD8aPOSabT_NK etc, 24 is CD8


Idents(integrated) <- integrated$celltype_only_scATAC

# > table(Idents(integrated)), the order of the cell types here match with the sequence of the colors
# B      CD2negGD      CD4posab     Monocytes            NK 
# 5011          2219          2865          2341           682 
# CD8abPOSab  CD4posCD8pos      CD2posGD CD8aPOSabT_NK          cDCs 
# 1843           376           334           344           215 
# ASC          pDCs 
# 215           159 

levels(integrated)# this order of cell type matches with the color order
# [1] "B"             "CD2negGD"      "CD4posab"      "Monocytes"    
# [5] "NK"            "CD8abPOSab"    "CD4posCD8pos"  "CD2posGD"     
# [9] "CD8aPOSabT_NK" "cDCs"          "ASC"           "pDCs"   

cols.scATAC <- c('lightpink', 'red',   #B    CD2negGD    
                 'plum',#CD4posab
                 'orange', 'gold', 'darkgreen', #Monocytes            NK  CD8abPOSab
                 'darkviolet', #CD4posCD8pos
                 'hot pink', 'brown', 'black',#CD2posGD     CD8aPOSabT_NK       cDC 
                 'mediumseagreen', 'skyblue2')#ASC          pDCs      without integration


p <- DimPlot(
  object = integrated, cols = cols.scATAC)#even if do not mention group.by, it would by default use identity 
p
ggsave("celltype.annotation.WO.scRNA.no.15.31.new.annotation.newest.pdf", path = "./plot.manuscript.final")
# ggsave("celltype.annotation.WO.scRNA.no.15.31.new.annotation.cluster23.pdf", path = "./plot.manuscript.final")

