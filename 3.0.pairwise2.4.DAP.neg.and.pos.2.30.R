
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
set.seed(1234)

integrated <- readRDS("integrated.four.with.cluster.2.4.2.30.rds")#resolution 2.4,2:30

clusters <- unique(Idents(integrated))
pairwise <- combn(clusters, 2)
# dim(pairwise)
# [1]   2 435

p1 <- pairwise[1,]; 
p2 <- pairwise[2,]
comps1 <- data.frame(p1, p2);
# [1] 435   2
colnames(comps1) <- c('pop1', 'pop2')

# head(comps1)
# pop1 pop2
# 1   23   16
# 2   23    2
# 3   23    8
# 4   23   12
# 5   23    6
# 6   23    7

comps2 <- data.frame(p2, p1);
colnames(comps2) <- c('pop1', 'pop2')#> dim(comps2),[1] 435   2
comps <- rbind(comps1, comps2)#[1] 870   2

numDAP.df <- as.data.frame(matrix(nrow=length(clusters),ncol=length(clusters)))#[1] 30 30
names(numDAP.df) <- seq(0,length(clusters) - 1,1)
rownames(numDAP.df) <- seq(0,length(clusters) - 1,1)

results <- list()

for(j in 1:nrow(comps)) {#
  print(paste0(j, " th  pairwise of ",nrow(comps), " of integrated.2 is being processed: "))#nrow(comps1) = 435 in total since C(30,2) = 435
  markers <- FindMarkers(integrated, ident.1 = comps[j,1], ident.2 = comps[j,2], assay = "peaks")#  only.pos = FALSE,so results include both pos and neg
  markers$gene <- rownames(markers)
  markers$pop1 <- paste(comps[j,1])
  markers$pop2 <- paste(comps[j,2])
  markers$comparison <- paste(markers$pop1, markers$pop2, sep = 'v')
  cluster1 <- as.numeric(comps[j,1])#as.numeric will increase cluster number by 1, for example comps1[1,1] is cluster 23, but as.numeric(comps1[1,1]) gives you 24
  cluster2 <- as.numeric(comps[j,2])#as.numeric will increase cluster number by 1, for example comps1[1,2] is cluster 16, but as.numeric(comps1[1,2]) gives you 17
  numDAP.df[cluster1, cluster2] <- dim(markers)[1]
  results[[j]] <- markers
  print(paste0(j, " th  pairwise of ",nrow(comps), " of integrated.2.4 is done!, " ))
}

pwAll <- do.call(rbind, results)
pwAll <- pwAll[order(pwAll$comp, pwAll$p_val_adj),]
write.table(pwAll, file= "integrated.2.4ClustersPairwiseDE_LFC25FDR5_CellExpress20%inAtleast1Cluster.neg.and.pos.2.30.no.filtering.txt", sep="\t", quote = FALSE)

sig <- subset(pwAll, (p_val_adj <= 0.05) & (pct.1 >= 0.2 | pct.2 >= 0.2))
write.table(sig, file= "integrated.2.4ClustersPairwiseDE_LFC25FDR5_CellExpress20%inAtleast1Cluster.neg.and.pos.2.30.no.logFC.txt", sep="\t", quote = FALSE)
write.table(numDAP.df, file= "numDAP.df.integrated.2.4ClustersPairwiseDE_LFC25FDR5_CellExpress20%inAtleast1Cluster.neg.and.pos.2.30.no.logFC.txt", sep="\t", quote = FALSE)
