
#5. identify differentially acessible peaks

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


# integrated <- readRDS("integrated.withRNA.assay.annotated.integration.scRNA.celltype2.4.2.30.after.filter.new.idents.rds")#resol 2.4, +-2000, 2.30, script 4.2
integrated <- readRDS("integrated.withRNA.assay.annotated.integration.scRNA.celltype2.4.2.30.after.filter.new.idents.RNA2NRA.wo.unknown.rds")#resol 2.4, +-2000, 2.30, script 4.2

# Idents(integrated)
# 11 Levels: B CD8abPOSab CD2negGD CD4posab NK Monocytes ... pDCs

DefaultAssay(integrated) <- 'peaks'

# Identifies differentially expressed genes between two groups of cells using a Wilcoxon Rank Sum test (default)
DAP.pig.2.4.2.30.findallmarkers.pos <- FindAllMarkers(integrated,
                         logfc.threshold = 0.25,
                         min.pct = 0.2,
                         only.pos = TRUE)
#assay = "peaks")
DAP.pig.2.4.2.30.findallmarkers.pos <- subset(DAP.pig.2.4.2.30.findallmarkers.pos, p_val_adj < 0.05)
DAP.pig.2.4.2.30.findallmarkers.pos <- DAP.pig.2.4.2.30.findallmarkers.pos[order(DAP.pig.2.4.2.30.findallmarkers.pos$p_val_adj),]

saveRDS(DAP.pig.2.4.2.30.findallmarkers.pos, file = "DAP.pig.2.4.2.30.findallmarkers.pos.integration.with.scRNA.RNA2RNA.rds")#17238     7
DAP.pig.2.4.2.30.findallmarkers.pos <- readRDS("DAP.pig.2.4.2.30.findallmarkers.pos.integration.with.scRNA.RNA2RNA.rds")#17238     7
write.table(DAP.pig.2.4.2.30.findallmarkers.pos, file = "DAP.pig.2.4.2.30.findallmarkers.pos.integration.with.scRNA.RNA2RNA.all.cell.type.txt", sep = "\t",row.names = FALSE,col.names=TRUE,quote = FALSE )

length(unique(DAP.pig.2.4.2.30.findallmarkers.pos$gene))#11872

library(biomaRt)
convertSscrofaGeneList <- function(x){
  require("biomaRt")
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  Sscrofa = useMart("ensembl", dataset = "sscrofa_gene_ensembl")
  genesV2 = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = x , mart = Sscrofa, attributesL = c("hgnc_symbol"), martL = human,verbose = F, uniqueRows=T)
  humanx <- unique(genesV2[, 2])
  # Print the first 6 genes found to the screen
  print(head(humanx))
  return(humanx)
}

unique(DAP.pig.2.4.2.30.findallmarkers.pos$cluster)
# [1] B             CD2negGD      CD4posab      Monocytes     NK           
# [6] CD8abPOSab    CD8aPOSabT_NK cDCs          ASC           pDCs         
# [11] CD2posGD     
# 11 Levels: B CD2negGD CD4posab Monocytes NK CD8abPOSab ... pDCs

all.cell.types.list <- list("cDCs", "Monocytes", "CD4posab","CD2negGD","B","ASC","CD2posGD","CD8abPOSab","NK","pDCs","CD8aPOSabT_NK")

DAP.sub.all.cell.type <- list()
closest_genes_result <- list()

for (i in 1:length(all.cell.types.list))
  # for (i in 1:1)
{
  cell_type <- all.cell.types.list[[i]]
  DAP.sub <- subset(DAP.pig.2.4.2.30.findallmarkers.pos, cluster ==  all.cell.types.list[[i]]) #
  open_DAP_sub <- DAP.sub[,"gene"] #character
  write.table(DAP.sub, file = paste0(all.cell.types.list[[i]], ".celltype.DAP.only.pos.integration.RNA.txt"), sep = "\t",row.names = FALSE,col.names=TRUE,quote = FALSE )
  
  # open_DAP_sub <- DAP.sub[DAP.sub$avg_log2FC > 0.25,"gene"]
  
  closest_genes <- ClosestFeature(integrated, regions = open_DAP_sub)#
  # closest_genes <- subset(closest_genes,distance = 0)
  
  only.gene.names <- closest_genes$gene_name
  # df <- read.table("./integrated.2.30.2.4.closest_genes.and.DAP.findallmarkers.by.cell.type.integration.scRNA/Monocytes.closest_genes.findallmarkers.integration.scRNA.txt", header = FALSE)
  only.gene.names <- na.omit(only.gene.names)
  
  write.table(only.gene.names, file = paste0(all.cell.types.list[[i]], ".closest_genes.findallmarkers.integration.scRNA.txt"), sep = "\t",row.names = FALSE,col.names=FALSE,quote = FALSE )
  
  # humangenes <- convertSscrofaGeneList(only.gene.names)
  # write.table(humangenes, file = paste0(all.cell.types.list[[i]], ".closest_genes.human.genes.findallmarkers.integration.scRNA.txt"), sep = "\t",row.names = FALSE,col.names=FALSE,quote = FALSE )
  #input for HOMER
  write.table(closest_genes$query_region, file = paste0(all.cell.types.list[[i]], ".DAP.findallmarkers.integration.scRNA.txt"), sep = "\t",row.names = FALSE,col.names=FALSE,quote = FALSE )
  
  closest_genes_result[[i]] <- closest_genes
  
  names(DAP.sub) <- list("p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj", "cluster",  "gene")
  
  DAP.sub.all.cell.type[[i]] <- DAP.sub
  # print(paste0("The number of DAP for ", cell_type,  " is : ", dim(closest_genes)[1],"; the number of closest pig genes: ", length(only.gene.names), ";The number of human genes is:", length(humangenes)))
  print(paste0("The number of DAP for ", cell_type,  " is : ", dim(closest_genes)[1],"; the number of closest pig genes: ", length(only.gene.names)))
  
}
names(closest_genes_result) <- list("cDCs", "Monocytes", "CD4posab","CD2negGD","B","ASC","CD2posGD","CD8abPOSab","NK","pDCs","CD8aPOSabT_NK")
names(DAP.sub.all.cell.type) <- list("cDCs", "Monocytes", "CD4posab","CD2negGD","B","ASC","CD2posGD","CD8abPOSab","NK","pDCs","CD8aPOSabT_NK")

saveRDS(closest_genes_result, file = "DAP.pig.2.4.2.30.findallmarkers.pos.integration.with.scRNA.closest.genes.result.RNA2RNA.avg_log2FC.0.25.no.restriction.distance.rds")#11
closest_genes_result <- readRDS("DAP.pig.2.4.2.30.findallmarkers.pos.integration.with.scRNA.closest.genes.result.RNA2RNA.avg_log2FC.0.25.no.restriction.distance.rds")#11

saveRDS(DAP.sub.all.cell.type, file = "DAP.pig.2.4.2.30.findallmarkers.pos.integration.with.scRNA.DAP.sub.RNA2RNA.avg_log2FC.0.25.no.restriction.distance.rds")#11
DAP.sub.all.cell.type <- readRDS("DAP.pig.2.4.2.30.findallmarkers.pos.integration.with.scRNA.DAP.sub.RNA2RNA.avg_log2FC.0.25.no.restriction.distance.rds")#11

list.files(path='/work/abg/pyang19/opt/pig.6798.6800.PBMC.Satija.pipeline.result/pbmc.1x.2x.cellrangeratac1.2.0.wd/strict.filter.new.annotation/strict.filter.new.annotation.common.peaks/commonpeaks.110444.new.directory/GO.input.each.B.cluster.and.Cell.type',pattern='.closest_genes.human.genes.findallmarkers.integration.scRNA.txt')
file <- "data.R"

if (file.exists(file)) {
  unlink(file)
  cat("The file is deleted")
}

celltype.DAP <- list.files(path='/work/abg/pyang19/opt/pig.6798.6800.PBMC.Satija.pipeline.result/pbmc.1x.2x.cellrangeratac1.2.0.wd/strict.filter.new.annotation/strict.filter.new.annotation.common.peaks/commonpeaks.110444.new.directory/',pattern='.celltype.DAP.only.pos.integration.RNA.txt')

folder <- "closest.pig.genes.cell.type"

if (file.exists(folder)) {
  
  cat("The folder already exists")
  
} else {
  
  dir.create(folder)
}


DefaultAssay(integrated) <- 'peaks'
levels(integrated) <- c('CD2negGD',
                        'CD2posGD',
                        'NK',
                        'CD8aPOSabT_NK',
                        'CD8abPOSab',
                        'CD4posab',
                        'ASC',
                        'B',
                        'pDCs',
                        'cDCs',
                        'Monocytes') # # this order will be the order from up to bottom tracks of plot 

