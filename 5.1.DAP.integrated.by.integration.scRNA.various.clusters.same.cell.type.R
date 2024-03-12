
#5. identify differentially acessible peaks

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

################  DAPs for different clusters of B cells ####################################  ####################################   

integrated <- readRDS("integrated.withRNA.assay.annotated.integration.scRNA.celltype2.4.2.30.after.filter.new.idents.RNA2NRA.wo.unknownrds")#resol 2.4, +-2000, 2.30, script 4.2

# Idents(integrated)
# 11 Levels: B CD8abPOSab CD2negGD CD4posab NK Monocytes ... pDCs

DefaultAssay(integrated) <- 'peaks'

Bcells <- subset(x = integrated, idents = "B")
# An object of class Seurat 
# 130479 features across 4993 samples within 2 assays 
# Active assay: peaks (110444 features, 110272 variable features)
# 1 other assay present: RNA
# 3 dimensional reductions calculated: integrated_lsi, umap, lsi

Idents(Bcells) <- Bcells$seurat_clusters#Levels: Levels: 0 3 4 8 11 18 20

DefaultAssay(Bcells)#peaks

DAP.pig.2.4.2.30.findallmarkers.pos.B <- FindAllMarkers(Bcells,
                                                      logfc.threshold = 0.25,
                                                      min.pct = 0.2,
                                                      only.pos = TRUE)

DAP.pig.2.4.2.30.findallmarkers.pos.B <- subset(DAP.pig.2.4.2.30.findallmarkers.pos.B, p_val_adj < 0.05)
DAP.pig.2.4.2.30.findallmarkers.pos.B <- DAP.pig.2.4.2.30.findallmarkers.pos.B[order(DAP.pig.2.4.2.30.findallmarkers.pos.B$p_val_adj),]#1681    7

saveRDS(DAP.pig.2.4.2.30.findallmarkers.pos.B, file = "DAP.pig.2.4.2.30.findallmarkers.pos.B.integration.with.scRNA.RNA2RNA.rds")#1681    7
DAP.pig.2.4.2.30.findallmarkers.pos.B <- readRDS("DAP.pig.2.4.2.30.findallmarkers.pos.B.integration.with.scRNA.RNA2RNA.rds")# 1681    7

length(unique(DAP.pig.2.4.2.30.findallmarkers.pos.B$gene))#1399

ranges.show <- StringToGRanges("4-77405901-77407598")
p <- CoveragePlot(
  object = Bcells,
  region = "4-77405901-77407598",#in cluster 4
  extend.upstream = 20000,
  extend.downstream = 20000,
  show.bulk = TRUE,
  # ranges = ranges.show,
  region.highlight = ranges.show
)#+ ggtitle('integrated.2.4.2.30.annotated.withscRNA.B.cells.4-77405901-77407598.cluster4')+ NoLegend()
p
ggsave("integrated.2.4.2.30.annotated.withscRNA.4.77405901.77407598.cluster4.pdf", path = "./plot.manuscript.final")


ranges.show <- StringToGRanges("16-68768150-68769183")
p <- CoveragePlot(
  object = Bcells,
  region = "16-68768150-68769183",#in cluster 0
  extend.upstream = 20000,
  extend.downstream = 20000,
  show.bulk = TRUE,
  # ranges = ranges.show,
  region.highlight = ranges.show)
p
ggsave("integrated.2.4.2.30.annotated.withscRNA.16.68768150.68769183.cluster0.pdf", path = "./plot.manuscript.final")

ranges.show <- StringToGRanges("1-100659543-100661154")
p <- CoveragePlot(
  object = Bcells,
  region = "1-100659543-100661154",#in cluster 0
  extend.upstream = 20000,
  extend.downstream = 20000,
  show.bulk = TRUE,
  # ranges = ranges.show,
  region.highlight = ranges.show)
p
ggsave("integrated.2.4.2.30.annotated.withscRNA.1.100659543.100661154.cluster0.pdf", path = "./plot.manuscript.final")


ranges.show <- StringToGRanges("3-32320185-32322810")
p <- CoveragePlot(
  object = Bcells,
  region = "3-32320185-32322810",#in cluster 8
  extend.upstream = 40000,
  extend.downstream = 20000,
  show.bulk = TRUE,
  # ranges = ranges.show,
  region.highlight = ranges.show
)#+ ggtitle('integrated.2.4.2.30.annotated.withscRNA.B.cells.3-32320185-32322810.cluster6')+ NoLegend()
p
ggsave("integrated.2.4.2.30.annotated.withscRNA.3.32320185.32322810.cluster8.pdf", path = "./plot.manuscript.final")


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

unique(DAP.pig.2.4.2.30.findallmarkers.pos.B$cluster)
# [1] 8  0  4  18 3  20 11
# Levels: 0 3 4 8 11 18 20

table(DAP.pig.2.4.2.30.findallmarkers.pos.B$cluster)
# 0   3   4   8  11  18  20 
# 487  53 244 474  51 187 185 

all.cell.types.list <- list("0", "3","4","8","11","18","20")

closest_genes_result <- list()
for (i in 1:length(all.cell.types.list))
{
  cell_type <- all.cell.types.list[[i]]
  DAP.sub <- subset(DAP.pig.2.4.2.30.findallmarkers.pos.B, cluster ==  all.cell.types.list[[i]]) #
  open_DAP_sub <- DAP.sub[,"gene"]
  # open_DAP_sub <- DAP.sub[DAP.sub$avg_log2FC > 0.25,"gene"]
  
  closest_genes <- ClosestFeature(integrated, regions = open_DAP_sub)#
  # closest_genes <- subset(closest_genes,distance = 0)
  
  only.gene.names <- closest_genes$gene_name
  # df <- read.table("./integrated.2.30.2.4.closest_genes.and.DAP.findallmarkers.by.cell.type.integration.scRNA/Monocytes.closest_genes.findallmarkers.integration.scRNA.txt", header = FALSE)
  only.gene.names <- na.omit(only.gene.names)
  
  write.table(only.gene.names, file = paste0("Bcellcluster", all.cell.types.list[[i]], ".closest_genes.findallmarkers.integration.scRNA.txt"), sep = "\t",row.names = FALSE,col.names=FALSE,quote = FALSE )
  
  # humangenes <- convertSscrofaGeneList(only.gene.names)
  # write.table(humangenes, file = paste0("B.cell.cluster",all.cell.types.list[[i]], ".closest_genes.human.genes.findallmarkers.integration.scRNA.txt"), sep = "\t",row.names = FALSE,col.names=FALSE,quote = FALSE )
  #input for HOMER
  write.table(closest_genes$query_region, file = paste0("Bcellcluster", all.cell.types.list[[i]], ".DAP.findallmarkers.integration.scRNA.txt"), sep = "\t",row.names = FALSE,col.names=FALSE,quote = FALSE )
  
  closest_genes_result[[i]] <- closest_genes
  # print(paste0("The number of DAP for cluster ", cell_type,  " is : ", dim(closest_genes)[1],"; the number of closest pig genes: ", length(only.gene.names), ";The number of human genes is:", length(humangenes)))
  print(paste0("The number of DAP for ", cell_type,  " is : ", dim(closest_genes)[1],"; the number of closest pig genes: ", length(only.gene.names)))
}
names(closest_genes_result) <- list("0", "3","4","8","11","18","20")

saveRDS(closest_genes_result, file = "DAP.pig.2.4.2.30.findallmarkers.pos.integration.with.scRNA.each.B.cluster.closest.genes.result.RNA2RNA.avg_log2FC.0.25.no.restriction.distance.rds")#
closest_genes_result <- readRDS("DAP.pig.2.4.2.30.findallmarkers.pos.integration.with.scRNA.each.B.cluster.closest.genes.result.RNA2RNA.avg_log2FC.0.25.no.restriction.distance.rds")#

