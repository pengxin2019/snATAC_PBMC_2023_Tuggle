
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
set.seed(1234)
 

pbmc_sc_rna <- readRDS("/work/abg/pyang19/opt/scRNA/scRNA.rds")

Idents(pbmc_sc_rna) <- pbmc_sc_rna$celltypes

# Idents(pbmc_sc_rna)
# # 13 Levels: Monocytes cDCs pDCs B cells ASC ... Erythrocytes

DefaultAssay(pbmc_sc_rna) <- 'RNA'
# An object of class Seurat 
# 30796 features across 28810 samples within 3 assays 
# Active assay: RNA (16140 features, 0 variable features)
# 2 other assays present: SCT, integrated
# 3 dimensional reductions calculated: pca, umap, tsne
pbmc_sc_rna <- FindVariableFeatures(pbmc_sc_rna)

Idents(pbmc_sc_rna)
# CD2- GD T cells            NK cells 
# 13 Levels: Monocytes cDCs pDCs B cells ASC ... Erythrocytes

pbmc_sc_rna <- RenameIdents(
  object = pbmc_sc_rna,
  'CD4+ ab T cells' = 'CD4posab',
  'B cells' = 'B',
  'CD2- GD T cells' = 'CD2negGD',
  'pDCs' = 'pDCs',
  'CD8a+ ab T/NK cells' = 'CD8aPOSabT_NK',
  'Monocytes' = 'Monocytes',
  'CD2+ GD T cells' = 'CD2posGD',
  'NK cells' = 'NK',
  'CD8ab+ ab T cells' = 'CD8abPOSab',
  'ASC' = 'ASC',
  'cDCs' = 'cDCs',
  'Unknown' = 'Unknown',
  'Erythrocytes' = 'Erythrocytes')

saveRDS(pbmc_sc_rna, file = "pbmc_sc_rna.with.shortened.cell.type.with.variable.rds")
pbmc_sc_rna <- readRDS("pbmc_sc_rna.with.shortened.cell.type.with.variable.rds")

# An object of class Seurat 
# 30796 features across 28810 samples within 3 assays 
# Active assay: RNA (16140 features, 2000 variable features)
# 2 other assays present: SCT, integrated
# 3 dimensional reductions calculated: pca, umap, tsne

Idents(pbmc_sc_rna)
# 13 Levels: CD4posab B CD2negGD pDCs CD8aPOSabT_NK Monocytes CD2posGD ... Erythrocytes

DE.pig.scRNA <- FindAllMarkers(pbmc_sc_rna,logfc.threshold = 0.25, min.pct = 0.2,only.pos = TRUE)
DE.pig.scRNA <- subset(DE.pig.scRNA, p_val_adj < 0.05)
DE.pig.scRNA <- DE.pig.scRNA[order(DE.pig.scRNA$p_val_adj),]#

saveRDS(DE.pig.scRNA, file = "DE.pig.PBMC.scRNA.by.shortened.cell.type.using.integrated.assay.shortened.rds")#4423    7
DE.pig.scRNA <- readRDS("DE.pig.PBMC.scRNA.by.shortened.cell.type.using.integrated.assay.shortened.rds")# 4423    7

all.cell.types.list <- c("ASC", "B","CD2negGD","CD2posGD","CD4posab","CD8abPOSab","CD8aPOSabT_NK","cDCs","Monocytes","NK","pDCs") # it is a factor


for (i in 1:11)
{
  DE.sub <- subset(DE.pig.scRNA, cluster == all.cell.types.list[i], ) 
  DE.gene.df <- as.data.frame(DE.sub$gene)
  names(DE.gene.df) <- "gene"
  
  write.table(DE.gene.df,file = paste0(all.cell.types.list[i], ".scRNA.DE.gene.df.", dim(DE.gene.df)[1], ".genes.txt"), row.names = FALSE, sep="\t", quote = FALSE)#coaccess > 0.05
}


