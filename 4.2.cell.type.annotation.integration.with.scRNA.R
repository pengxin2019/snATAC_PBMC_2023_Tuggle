#4.cell type annotation new wd on by using scATAC data, in dependent of scRNA dataset

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

integrated <- readRDS("integratedscATAC.withRNA.assay2.4.2.30.4000.rds")#resoltion2.4,2:30, +-2000
# An object of class Seurat 
# 130479 features across 17207 samples within 2 assays 
# Active assay: RNA (20035 features, 0 variable features)
# 1 other assay present: peaks
# 3 dimensional reductions calculated: integrated_lsi, umap, lsi

Idents(integrated)
# 35 Levels: 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 ... 34
DefaultAssay(integrated)#RNA
# Load the pre-processed scRNA-seq data for PBMCs
pbmc_sc_rna <- readRDS("/work/abg/pyang19/opt/scRNA/scRNA.rds")
# pbmc_sc_rna
# An object of class Seurat 
# 30796 features across 28810 samples within 3 assays 
# Active assay: RNA (16140 features, 0 variable features)
# 2 other assays present: SCT, integrated
# 3 dimensional reductions calculated: pca, umap, tsne

Idents(pbmc_sc_rna) <- pbmc_sc_rna$celltypes

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
pbmc_sc_rna <- readRDS("pbmc_sc_rna.with.shortened.cell.type.with.variable.rds")#from 5.2 script

DefaultAssay(pbmc_sc_rna)#RNA

# pbmc_sc_rna
# An object of class Seurat 
# 30796 features across 28810 samples within 3 assays 
# Active assay: RNA (16140 features, 0 variable features)
# 2 other assays present: SCT, integrated
# 3 dimensional reductions calculated: pca, umap, tsne

pbmc_sc_rna <- FindVariableFeatures(pbmc_sc_rna)
# pbmc_sc_rna
# An object of class Seurat 
# 30796 features across 28810 samples within 3 assays 
# Active assay: RNA (16140 features, 2000 variable features)
# 2 other assays present: SCT, integrated
# 3 dimensional reductions calculated: pca, umap, tsne
pbmc_sc_rna$shorted_celltype <- Idents(pbmc_sc_rna)

#FindTransferAnchors: features, Features to use for dimensional reduction. If not specified, set as variable features of the reference object which are also present in the query.
# if use RNA assay for both dataset, Error: No features to use in finding transfer anchors. To troubleshoot, try explicitly providing features to the features parameter and ensure that they are present in both reference and query assays.

transfer.anchors.integration.RNA <- FindTransferAnchors(
  reference = pbmc_sc_rna,
  query = integrated,
  reduction = 'cca'#3Run a CCA on the reference and query
)

saveRDS(transfer.anchors.integration.RNA, file = "transfer.anchors.integration.RNA.2.4.2.30.4000.RNA2RNA.rds")#resol 2.4, +-2000, 2:30, use RNA assay of scRNA
transfer.anchors.integration.RNA <- readRDS("transfer.anchors.integration.RNA.2.4.2.30.4000.RNA2RNA.rds")#resol 2.4, +-2000, 2:30, use RNA assay of scRNA


predicted.labels <- TransferData(
  anchorset = transfer.anchors.integration.RNA,
  refdata = pbmc_sc_rna$shorted_celltype,
  weight.reduction = integrated[['integrated_lsi']],
  dims = 2:30
)

integrated <- AddMetaData(object = integrated, metadata = predicted.labels)

saveRDS(integrated, file = "integrated.withRNA.assay.annotated.integration.scRNA.celltype2.4.2.30.RNA2RNA.rds")#resol 2.4, +-2000, 2.30, 
integrated <- readRDS("integrated.withRNA.assay.annotated.integration.scRNA.celltype2.4.2.30.RNA2RNA.rds")#resol 2.4, +-2000, 2.30, 

# An object of class Seurat 
# 130479 features across 17207 samples within 2 assays 
# Active assay: RNA (20035 features, 0 variable features)
# 1 other assay present: peaks
# 3 dimensional reductions calculated: integrated_lsi, umap, lsi


unique(integrated$predicted.id)
# [1] "CD8abPOSab"    "cDCs"          "B"             "CD2posGD"     
# [5] "CD2negGD"      "Monocytes"     "CD4posab"      "ASC"          
# [9] "CD8aPOSabT_NK" "NK"            "pDCs"          "Erythrocytes" 

table(integrated$predicted.id) # 

# ASC             B      CD2negGD      CD2posGD      CD4posab 
# 176          5202          2738           339          3023 
# CD8abPOSab CD8aPOSabT_NK          cDCs  Erythrocytes     Monocytes 
# 1875           427           166             1          2335 
# NK          pDCs 
# 777           148 

all.predicted.cell.type <- list("prediction.score.CD4posab",
                                "prediction.score.B",
                                "prediction.score.CD2negGD",
                                "prediction.score.pDCs",
                                "prediction.score.CD8aPOSabT_NK",
                                "prediction.score.Monocytes",
                                "prediction.score.CD2posGD",
                                "prediction.score.NK",
                                "prediction.score.CD8abPOSab",
                                "prediction.score.ASC",
                                "prediction.score.cDCs",
                                "prediction.score.Unknown",
                                "prediction.score.Erythrocytes",
                                "prediction.score.max")


for (i in 1:length(all.predicted.cell.type))#1:14
{
  # sub.cell.type <- subset(x = integrated, idents = all.predicted.cell.type[[i]])
  
  p <- FeaturePlot(
    object = integrated,
    features = all.predicted.cell.type[[i]],
    pt.size = 0.03,
    label = TRUE,
    label.size = 1, ncol = 2, combine = TRUE, coord.fixed = FALSE,
    cols = c("lightgrey", "red")) + ggtitle( all.predicted.cell.type[[i]])
  p
  ggsave(paste0("prediction.score.RNA2RNA.", all.predicted.cell.type[[i]], ".pdf"), path = "./plot.manuscript.final")
  # ggsave(paste0("prediction.score.old.doublet.RNA2RNA.", all.predicted.cell.type[[i]], ".pdf"), path = "./plots.new.doublet.removal.manuscript")
}


unique(integrated$predicted.id)
# [1] "CD8abPOSab"    "cDCs"          "B"             "CD2posGD"     
# [5] "CD2negGD"      "Monocytes"     "CD4posab"      "ASC"          
# [9] "CD8aPOSabT_NK" "NK"            "pDCs"          "Erythrocytes" 
cols.scRNA <- c('lightpink', 'red', 'orange', 'gold', 'darkgreen', 'mediumseagreen', 'skyblue2', 'steelblue', 'navy', 'plum3', 'darkmagenta', 'black', 'grey')#13
cols.scATAC <- c('lightpink', 'red', 'orange', 'gold', 'darkgreen', 'mediumseagreen', 'skyblue2', 'steelblue', 'navy', 'plum3','darkmagenta','black')#12
cols.scATAC <- c('lightpink', 'red', 'orange', 'gold', 'darkgreen', 'mediumseagreen', 'skyblue2', 'steelblue', 'navy', 'plum3','darkmagenta')#11

p <- DimPlot(
  object = pbmc_sc_rna,
  group.by = 'shorted_celltype',
  cols = cols.scRNA)  + ggtitle('scRNA-seq')
p
ggsave("scRNA.cell.type.pdf", path = "./plot.manuscript.final")

p <- DimPlot(
  object = integrated,
  group.by = 'predicted.id',
  cols = cols.scATAC) + ggtitle('scATAC.integration.with.scRNA.before.filter')
p
ggsave("scATAC.celltypescATAC.integration.RNA2RNA.pdf", path = "./plot.manuscript.final")



hist(integrated$prediction.score.max)
abline(v = 0.5, col = "red")# the line on the x axis is 0.5

integrated <- subset(integrated, subset = prediction.score.max > 0.5)
# An object of class Seurat 
# 130479 features across 16043 samples within 2 assays 
# Active assay: RNA (20035 features, 0 variable features)
# 1 other assay present: peaks
# 3 dimensional reductions calculated: integrated_lsi, umap, lsi

# tst <- subset(integrated, subset = predicted.id == "Erythrocytes")
# Error: No cells found

saveRDS(integrated, file = "integrated.withRNA.assay.annotated.integration.scRNA.celltype2.4.2.30.after.filter.RNA2RNA.rds")#resol 2.4, +-2000, 2.30
integrated <- readRDS("integrated.withRNA.assay.annotated.integration.scRNA.celltype2.4.2.30.after.filter.RNA2RNA.rds")#resol 2.4, +-2000, 2.30

cols.scATAC <- c('lightpink', 'red', 'orange', 'gold', 'darkgreen', 'mediumseagreen', 'skyblue2', 'steelblue', 'navy', 'plum3','darkmagenta')#11
plot2 <- DimPlot(
  object = integrated,
  group.by = 'predicted.id',
  cols = cols.scATAC) + ggtitle('scATAC.integration.with.scRNA.after.filter')
plot2
ggsave("scATAC.celltypescATAC.seq.2.4.2.30.RNA2RNA.after.filter.pdf", path = "./plot.manuscript.final")

Idents(integrated) <- integrated$seurat_clusters

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
  '15' = 'unknown',# 
  '16' = 'CD2negGD',# 
  '17' = 'CD4posab',#CD2negGD: 181; CD4posab: 191
  '18' = 'B',
  '19' = 'CD8abPOSab',
  '20' = 'B',
  '21' = 'CD8abPOSab',
  '22' = 'CD2posGD',
  '23' = 'CD8abPOSab',
  '24' = 'CD8abPOSab',
  '25' = 'CD8aPOSabT_NK',
  '26' = 'cDCs',
  '27' = 'ASC',
  '28' = 'CD2negGD',
  '29' = 'pDCs',
  '30' = 'CD8abPOSab',
  '31' = 'unknown',
  '32' = 'CD8abPOSab',
  '33' = 'CD8aPOSabT_NK',
  '34' = 'CD4posab'
)

integrated$cell_type_integration_scRNA <- Idents(integrated)

saveRDS(integrated, file = "integrated.withRNA.assay.annotated.integration.scRNA.celltype2.4.2.30.after.filter.new.idents.RNA2NRA.rds")#resol 2.4, +-2000, 2.30
integrated <- readRDS("integrated.withRNA.assay.annotated.integration.scRNA.celltype2.4.2.30.after.filter.new.idents.RNA2NRA.rds")#resol 2.4, +-2000, 2.30

integrated <- subset(integrated, idents = "unknown", invert = TRUE)#
# An object of class Seurat 
# 130479 features across 15626 samples within 2 assays 
# Active assay: RNA (20035 features, 0 variable features)
# 1 other assay present: peaks
# 3 dimensional reductions calculated: integrated_lsi, umap, lsi
saveRDS(integrated, file = "integrated.withRNA.assay.annotated.integration.scRNA.celltype2.4.2.30.after.filter.new.idents.RNA2NRA.wo.unknown.rds")#resol 2.4, +-2000, 2.30
integrated <- readRDS("integrated.withRNA.assay.annotated.integration.scRNA.celltype2.4.2.30.after.filter.new.idents.RNA2NRA.wo.unknown.rds")#resol 2.4, +-2000, 2.30

cols.scATAC <- c('lightpink', 'red',   #B    CD2negGD    
                 'plum',#CD4posab
                 'orange', 'gold', 'darkgreen', #Monocytes            NK  CD8abPOSab
                 'hot pink', #  CD2posGD
                 'brown', #CD8aPOSabT_NK
                 'black', #cDCs
                 'mediumseagreen', 'skyblue2')#ASC          pDCs      with integration

p <- DimPlot(
  object = integrated,
  cols = cols.scATAC)# + ggtitle('scATAC-seq.2.4.2.30.annotated.integration.with.scRNA.after.filter.prediction.score.0.5.no.31.new.ident.RNA2NRA')
p
ggsave("celltype.annotation.WITH.scRNA.prediction.score.0.5.no.15.31.RNA2RNA.pdf", path = "./plot.manuscript.final")



## Fig 3 make a dotplot for chromatin accessibility to visualize DAPs
DefaultAssay(integrated)#RNA
DefaultAssay(integrated) <- 'peaks'
Idents(integrated) <- integrated$cell_type_integration_scRNA

levels(integrated) <- c('Monocytes',
                        'cDCs',
                        'pDCs',
                        'B',
                        'ASC',
                        'CD4posab',
                        'CD8abPOSab',
                        'CD8aPOSabT_NK',
                        'NK',
                        'CD2posGD',
                        'CD2negGD') # # this order will be the order from UP to BOTTOM of plot Reorder the clusters b

p <- DotPlot(integrated, features = c('2-151108982-151113383',# Monocytes, CSF1R
                                '2-142344198-142344875', # Monocytes, CD14
                                '11-5414316-5415903', #cDCs FLT3 p_val_adj= 1.000466e-33
                                # '11-5414316-5415903', #pDCs FLT3 p_val_adj= 1.015993e-87
                                '1-237781345-237790424', # B PAX5 tss
                                '3-18591541-18594649', # B CD19
                                '9-45618596-45622334', # T, CD3E
                                '5-63915302-63918975', # CD4posab, CD4
                                '3-57998628-58000477', # CD8abPOSab, CD8A, TSS
                                '3-57967427-57970939', # CD8abPOSab, CD8B
                                '14-73519353-73521831', # NK PRF1
                                '5-61640161-61642007', # NK KLRK1
                                #'3-58027313-58029236', # NK CD8A, lowest p adj value, rna dotplot will have 2 CD8As
                                '7-76544172-76549080' #CD2posGD  NO CD2, TRDC region is same as that of CD2negGD, beow is 7-76544172-76549080  CD2negGD TRDC 
                                ), 
        cols = c('yellow', 'red')) + RotatedAxis()+ theme(axis.text.x = element_text(size = 5)) 
p
ggsave(paste0("dotplot.CA.celltype.11DAP.scATAC.by.celltype.pdf"), path = "./plot.manuscript.final")

### scRNA dot plot for gene markers across different cell types
cell.marker.list <- list('CSF1R','CD14', 'FLT3','PAX5', 'CD19','CD3E','CD4','CD8A','CD8B','PRF1','KLRK1','CD8A','TRDC')

levels(pbmc_sc_rna) <- c('Unknown',
                         'Erythrocytes',
                         'Monocytes',
                         'cDCs',
                         'pDCs',
                         'B',
                         'ASC',
                         'CD4posab',
                         'CD8abPOSab',
                         'CD8aPOSabT_NK',
                         'NK',
                         'CD2posGD',
                         'CD2negGD'
                         ) # # this order will be the order from BORROM to UP of plot Reorder the clusters b

pbmc_sc_rna <- subset(x = pbmc_sc_rna, idents = c("Unknown", "Erythrocytes"), invert = TRUE)#28570 cells
# An object of class Seurat 
# 30796 features across 28570 samples within 3 assays 
# Active assay: RNA (16140 features, 2000 variable features)
# 2 other assays present: SCT, integrated
# 3 dimensional reductions calculated: pca, umap, tsne

levels(pbmc_sc_rna) <- c('Monocytes',
                         'cDCs',
                         'pDCs',
                         'B',
                         'ASC',
                         'CD4posab',
                         'CD8abPOSab',
                         'CD8aPOSabT_NK',
                         'NK',
                         'CD2posGD',
                         'CD2negGD'
) # # this order will be the order from BORROM to UP of plot Reorder the clusters b

p <- DotPlot(pbmc_sc_rna, features = c('CSF1R','CD14', 'FLT3','PAX5', 'CD19','CD3E','CD4','CD8A','CD8B','PRF1','KLRK1','TRDC'), cols = c('yellow', 'red')) + RotatedAxis()
p
ggsave("dotplot.scRNA.11celltypes..pdf", path = "./plot.manuscript.final")

### scRNA Ridge plot for gene markers across different cell types
#CD4 CSF1R PAX5
# Ridge plots - from ggridges. Visualize single cell expression distributions in each cluster

for (i in 1:length(cell.marker.list))
{
  p <- RidgePlot(pbmc_sc_rna, features = cell.marker.list[[i]],  slot = "counts")
  p
  ggsave(paste0("ridgeplot", cell.marker.list[[i]], ".scRNA.pdf"), path = "./plot.manuscript.final")
}
