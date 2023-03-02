
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
library(VennDiagram)
set.seed(1234)

DE.pig.scRNA <- readRDS("DE.pig.PBMC.scRNA.by.shortened.cell.type.using.integrated.assay.shortened.rds")# 4423    7, script 5.2


closest_genes_result <- readRDS("DAP.pig.2.4.2.30.findallmarkers.pos.integration.with.scRNA.closest.genes.result.RNA2RNA.avg_log2FC.0.25.no.restriction.distance.rds")#11, script 5.0

####### only identify overlaped genes who has a gene name #################################################

DE.pig.scRNA.remove <- DE.pig.scRNA[grep("ENSSSCG000000", DE.pig.scRNA$gene), ]#751   7
gene.use <- setdiff(DE.pig.scRNA$gene, DE.pig.scRNA.remove$gene)#2100
DE.pig.scRNA.use <-  DE.pig.scRNA[DE.pig.scRNA$gene %in% gene.use,]# 3672    7



all.cell.types.list <- list("cDCs", "Monocytes", "CD4posab","CD2negGD","B","ASC","CD2posGD","CD8abPOSab","NK","pDCs","CD8aPOSabT_NK")
overlapped_genes_result <- list()

for (i in 1:length(all.cell.types.list))
{
  DAP.genes.sub <- na.omit(closest_genes_result[[i]]$gene_name)  
  
  DE.pig.scRNA.use.sub <- subset(DE.pig.scRNA.use,cluster == all.cell.types.list[[i]])
  DE.pig.scRNA.use.sub.genes <- DE.pig.scRNA.use.sub$gene
  
  overlapped.genes <- intersect(DAP.genes.sub,DE.pig.scRNA.use.sub.genes)
  overlapped_genes_result[[i]] <- overlapped.genes
  
  write.table(overlapped.genes, file = paste0(all.cell.types.list[[i]], ".overlapped.genes.DE.DAP.txt"), sep = "\t",row.names = FALSE,col.names=FALSE,quote = FALSE )
  
  fraction <- length(overlapped.genes)/length(DE.pig.scRNA.use.sub.genes)*100
  
  grid.newpage()
  draw.pairwise.venn(length(DAP.genes.sub), length(DE.pig.scRNA.use.sub.genes), length(overlapped.genes), category = c(paste0(all.cell.types.list[[i]],".DAP.genes"), paste0(all.cell.types.list[[i]], ".DE.genes")), lty = rep("blank", 
                                                                                          2), fill = c("light blue", "pink"), alpha = rep(0.5, 2), cat.pos = c(0, 
                                                                                                                                                               0), cat.dist = rep(0.025, 2))

  # ggsave(paste0(all.cell.types.list[[i]],".venn.overlapped.genes.pdf"), path = "./plots.strict.filter.new.annotation.common.peaks.new.directory/manuscript_plot")
  
  print(paste0("The number of overlapped genes in  ", all.cell.types.list[[i]],  " : ", length(overlapped.genes), " number of DAP.genes.sub is: ",length(DAP.genes.sub) ,"; number of DE.genes is ", length(DE.pig.scRNA.use.sub.genes), ";   overlapped.genes#/genesDE# is : ", fraction))
}

names(overlapped_genes_result) <- list("cDCs", "Monocytes", "CD4posab","CD2negGD","B","ASC","CD2posGD","CD8abPOSab","NK","pDCs","CD8aPOSabT_NK")
                      
saveRDS(overlapped_genes_result, file = "cell.type.DE.DAP.overlapped.genes.result.rds")#11
overlapped_genes_result <- readRDS("cell.type.DE.DAP.overlapped.genes.result.rds")#11
