# ##4. annotate cell types based on the cis-elements
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
library(biomaRt)
set.seed(1234)

integrated <- readRDS("integrated.four.with.cluster.2.4.2.30.rds")#resolution 2.4,2:30, 3.integration.last.strict.non.batch.common.peaks.new.directory.2.30.R
# An object of class Seurat 
# 110444 features across 17207 samples within 1 assay 
# Active assay: peaks (110444 features, 110272 variable features)
# 3 dimensional reductions calculated: integrated_lsi, umap, lsi
DefaultAssay(integrated)#peaks

# table(Idents(integrated))
# 0    1    2    3    4    5    6    7    8    9   10   11   12   13   14   15 
# 1065 1009  974  927  839  838  793  785  729  718  682  640  584  574  563  461 
# 16   17   18   19   20   21   22   23   24   25   26   27   28   29   30   31 
# 450  428  419  404  392  376  334  315  263  239  215  215  197  159  146  142 
# 32   33   34 
# 141  105   86 


#### findallmarker only pos, no distance restriction ###
DE.pig.resol2.4.2.30.findallmarker.onlypos.p0.05 <- readRDS("DE.pig.resol2.4.2.30.findallmarker.onlypos.p0.05.rds")#59275     7, from 3.findallmarkers2.4.DAP.2.30.R

length(unique(DE.pig.resol2.4.2.30.findallmarker.onlypos.p0.05$gene))#20070, 20070/110444 = 0.1817211

all.genomic.regions <- DE.pig.resol2.4.2.30.findallmarker.onlypos.p0.05$gene#59275
closest.genes.all.genomic.regions.integrated.2.30.2.4.findallmarker.onlypos.p0.05 <- ClosestFeature(integrated, regions = all.genomic.regions)#59275    24
DE.pig.resol2.4.2.30.findallmarker.onlypos.p0.05.2 <- DE.pig.resol2.4.2.30.findallmarker.onlypos.p0.05
DE.pig.resol2.4.2.30.findallmarker.onlypos.p0.05.2$closest_genes <- closest.genes.all.genomic.regions.integrated.2.30.2.4.findallmarker.onlypos.p0.05$gene_name
DE.pig.resol2.4.2.30.findallmarker.onlypos.p0.05.2$distance <- closest.genes.all.genomic.regions.integrated.2.30.2.4.findallmarker.onlypos.p0.05$distance# 59275     9
# DE.pig.resol2.4.2.30.findallmarker.onlypos.p0.05.2 <- subset(DE.pig.resol2.4.2.30.findallmarker.onlypos.p0.05.2,distance <= 5000)#     9
DE.pig.resol2.4.2.30.findallmarker.onlypos.p0.05.2 <- DE.pig.resol2.4.2.30.findallmarker.onlypos.p0.05.2[!(is.na(DE.pig.resol2.4.2.30.findallmarker.onlypos.p0.05.2$closest_genes) | DE.pig.resol2.4.2.30.findallmarker.onlypos.p0.05.2$closest_genes==""), ]# 43913     9
DE.pig.resol2.4.2.30.findallmarker.onlypos.p0.05.2 <- DE.pig.resol2.4.2.30.findallmarker.onlypos.p0.05.2[order(DE.pig.resol2.4.2.30.findallmarker.onlypos.p0.05.2$cluster),]

saveRDS(DE.pig.resol2.4.2.30.findallmarker.onlypos.p0.05.2, file = "DE.pig.resol2.4.2.30.findallmarker.onlypos.p0.05.2.with.gene.name.no.dis.rds")#43913     9, NO distance col restriction,no NA
DE.pig.resol2.4.2.30.findallmarker.onlypos.p0.05.2 <- readRDS("DE.pig.resol2.4.2.30.findallmarker.onlypos.p0.05.2.with.gene.name.no.dis.rds")#43913     9, NO distance col restriction ,no NA

# head(DE.pig.resol2.4.2.30.findallmarker.onlypos.p0.05.2)
# p_val avg_log2FC pct.1 pct.2 p_val_adj cluster
# 1-272633753-272634907     0   2.797062 0.204 0.014         0       0
# 15-17315776-17316494      0   2.588860 0.242 0.018         0       0
# 16-68768150-68769183      0   2.466042 0.516 0.042         0       0
# 17-48039835-48040939      0   2.402292 0.350 0.032         0       0
# 7-9835264-9836154         0   2.340606 0.291 0.028         0       0
# 15-56312183-56312852      0   2.265342 0.274 0.026         0       0
# gene closest_genes distance
# 1-272633753-272634907 1-272633753-272634907        SPACA9        0
# 15-17315776-17316494   15-17315776-17316494       TMEM163        0
# 16-68768150-68769183   16-68768150-68769183       GALNT10        0
# 17-48039835-48040939   17-48039835-48040939         ACOT8        0
# 7-9835264-9836154         7-9835264-9836154         GFOD1        0
# 15-56312183-56312852   15-56312183-56312852        MFHAS1    17091

# some of the gene markers used in scRNA paper have alternative names, check whether we have DAP near their alternative gene names
subset(DE.pig.resol2.4.2.30.findallmarker.onlypos.p0.05.2,closest_genes == "DAP10")#0
subset(DE.pig.resol2.4.2.30.findallmarker.onlypos.p0.05.2,closest_genes == "NKG2D")#0
subset(DE.pig.resol2.4.2.30.findallmarker.onlypos.p0.05.2,closest_genes == "CD161")#0
subset(DE.pig.resol2.4.2.30.findallmarker.onlypos.p0.05.2,closest_genes == "DAP12")#0

### subset the cluster DAPs near cell type markers used in scRNA paper
##some of the gene markers used in scRNA paper have alternative names, 
#we did not use alternative gene names here since we knew that altenative genes were not near DAP either based on code above
all.gene.marker.list <- list("CD19", "CD79A", "PAX5", "IRF4","PRDM1",
                             "CD14", "CD163", "NLRP3","TLR4","CSF1R","CD86","SIRPA",
                             "FLT3","FCER1A","SLA-DRB1","SLA-DRA",
                             "TCF4","IRF8","CD4","CD8B","XBP1","CLEC12A","CD93",
                             "CD3E","TRDC","CD4","CD8A","CD8B","CD2","PRF1",
                             "KLRB1","KLRK1","CD5","CD6","TYROBP","HCST")
length(all.gene.marker.list)
# [1] 36
all.result <- list()
for (i in 1:length(all.gene.marker.list))
{
  result <-  subset(DE.pig.resol2.4.2.30.findallmarker.onlypos.p0.05.2,closest_genes == all.gene.marker.list[[i]])#
  all.result[[i]] <- result
}

pwAll <- do.call(rbind, all.result)
pwAll <- pwAll[order(pwAll$cluster, pwAll$p_val_adj),]#max(pwAll$distance): 29427, dim: 434   9
pwAll <- pwAll[order(pwAll$closest_genes),]#max(pwAll$distance): 29427, dim: 434   9
pwAll <- pwAll[order(pwAll$cluster),]#max(pwAll$distance): 29427, dim: 434   9

saveRDS(pwAll, file = "pwAll.nodistance.col.findallmarker.rds")#max(pwAll$distance): 29427, dim: 434   9
pwAll <- readRDS("pwAll.nodistance.col.findallmarker.rds")#max(pwAll$distance): 29427, dim: 434   9

pwAllsub <- pwAll[,c("cluster", "gene", "closest_genes", "distance", "avg_log2FC","p_val_adj")]# 434   6
saveRDS(pwAllsub, file = "DE.pig.resol2.4.2.30.findallmarker.onlypos.p0.05.2.with.gene.name.no.distance.restriction.rds")#434   6
pwAllsub <- readRDS("DE.pig.resol2.4.2.30.findallmarker.onlypos.p0.05.2.with.gene.name.no.distance.restriction.rds")#434   6

write.table(pwAllsub,file = "DAP.pig.resol2.4.2.30.findallmarker.onlypos.p0.05.with.36gene.marker.no.distance.col.sub.txt", row.names = FALSE, sep="\t", quote = FALSE)# there are 35 clusters in total 
#this transfered table is used for cell type annotation only using scATAC-seq dataset

length(unique(pwAllsub$gene))# check whether all 36 genes are near to at least one DAP

############ make a dotplot visualizing the chromatin accessibility at DAP below ######################## 
############ make a dotplot visualizing the chromatin accessibility at DAP below ######################## 

levels(integrated) <- c('5', '7', '9',# Monocytes
                        '0','3', '4', '8', '11', '18', '20', #B # these clusters will be at the bottom of the dotplot actually
                        '27',#ASC
                        '26', # cDCs
                        '29', # pDCs
                        '2', '6', '12', '17', '34',#CD4posab
                        '13', '19', '23','24', '30', '32', #CD8abPOSab
                        '21', #CD4posCD8pos
                        '10',# NK
                        '25','33',#CD8aPOSabT_NK
                        '22', #CD2posGD
                        '1', '14', '16', '28',#CD2negGD
                        '15','31' # unknown
)

cluster.DAP <- subset(pwAllsub,cluster == "0")
cluster.DAP <- subset(pwAllsub,closest_genes == "CD19")

# in the below
# the DAPs below include all DAPs in DotPlot in script 4.2 
p <- DotPlot(integrated, features = c('2-151108982-151113383','2-151125625-151130033', '2-142344198-142344875', '13-138488083-138490791','13-138454263-138457531',# Monocytes, CSF1R, CSF1R, CD14, CD86, CD86
                                      '1-237781345-237790424','3-18591541-18594649', # B,PAX5,CD19
                                      '1-72319331-72321667','1-104711299-104712024','1-104762287-104764023', # ASC, PRDM1, TCF4, TCF4 (TCF was expressed in ASC in scRNA though TCF4 was not described as marker for ASC in scRNA paper.)
                                      '11-5414316-5415903',# cDCs pDCs FLT3 
                                      '7-24889962-24896331','7-24911081-24915245','7-24925141-24926675',#cDCs, SLA-DRB1, SLA-DRB1, SLA-DRB1
                                      '14-46014828-46017493','6-3042818-3044559','6-3046076-3047969','5-63912663-63913826', # pDCs, XBP1, IRF8,IRF8,CD4 (this CD4DAP 5-63912663-63913826 is not open in CD4posab cells) (another CD4 DAP5-63915302-63918975 and a CD8B DAP 3-57990004-57991024 are decribed below)
                                      '9-45618596-45622334', # T, CD3E
                                      '5-63915302-63918975',# CD4posab, CD4
                                      '3-57967427-57970939','3-57990004-57991024', # CD8abPOSab, CD8B CD8B
                                      '3-57998628-58000477',#  CD8abPOSab, CD8A TSS,CD8A ('3-58033618-58035802' is open in all CD4 and CD8 cells)
                                      '14-73519353-73521831','5-61640161-61642007', # NK PRF1,  KLRK1 (excluded "CD8A":'3-58027313-58029236')
                                      '7-76544172-76549080'#, #CD2posGD  NO CD2, TRDC region is same as that of TRDC CD2negGD ##('2-10732450-10737689', '2-10782339-10787079' , CD6, CD6)
                                      ),# SLA-DRB1, SLA-DRB1 were expressed in GD in scRNA though not described as gene marker
             cols = c('yellow', 'red')) + RotatedAxis()+ theme(axis.text.x = element_text(size = 5)) 
p
ggsave("DotPlot.cluster.DAP.pdf",path = "./plot.manuscript.final/findallmarkers.celltype.annotation" )

############ make a dotplot visualizing the chromatin accessibility at DAP above ######################## 
############ make a dotplot visualizing the chromatin accessibility at DAP above ######################## 

dir.create("./plot.manuscript.final/findallmarkers.celltype.annotation")

p <- FeaturePlot(
  object = integrated,
  features = "7-76544172-76549080",
  pt.size = 0.03,
  label = FALSE,
  label.size = 1, ncol = 2, combine = TRUE, coord.fixed = FALSE,
  cols = c("lightgrey", "red")) + ggtitle("TRDC:7-76544172-76549080")
p
ggsave(paste0("cluster1.TRDC.pdf"),path = "./plot.manuscript.final/findallmarkers.celltype.annotation" )


p <- FeaturePlot(
  object = integrated,
  features = "13-138462336-138465245",
  pt.size = 0.03,
  label = FALSE,
  label.size = 1, ncol = 2, combine = TRUE, coord.fixed = FALSE,
  cols = c("lightgrey", "red")) + ggtitle("CD86:13-138462336-138465245")
p
ggsave(paste0("cluster4.CD86.pdf"),path = "./plot.manuscript.final/findallmarkers.celltype.annotation" )

p <- FeaturePlot(
  object = integrated,
  features = "4-103932010-103932891",
  pt.size = 0.03,
  label = FALSE,
  label.size = 1, ncol = 2, combine = TRUE, coord.fixed = FALSE,
  cols = c("lightgrey", "red")) + ggtitle("CD2:4-103932010-103932891")
p
ggsave(paste0("cluster5.CD2.pdf"),path = "./plot.manuscript.final/findallmarkers.celltype.annotation" )


#cluster5 11-5402036-5403846, FLT3
p <- FeaturePlot(
  object = integrated,
  features = "11-5402036-5403846",
  pt.size = 0.03,
  label = FALSE,
  label.size = 1, ncol = 2, combine = TRUE, coord.fixed = FALSE,
  cols = c("lightgrey", "red")) + ggtitle("FLT3:11-5402036-5403846")
p
ggsave(paste0("cluster5.FLT3.pdf"),path = "./plot.manuscript.final/findallmarkers.celltype.annotation" )

#cluster9,3-18596888-18598580, CD19
p <- FeaturePlot(
  object = integrated,
  features = "3-18596888-18598580",
  pt.size = 0.03,
  label = FALSE,
  label.size = 1, ncol = 2, combine = TRUE, coord.fixed = FALSE,
  cols = c("lightgrey", "red")) + ggtitle("CD19:3-18596888-18598580")
p
ggsave(paste0("cluster9.CD19.pdf"),path = "./plot.manuscript.final/findallmarkers.celltype.annotation" )

#9	11-5402036-5403846	FLT3
p <- FeaturePlot(
  object = integrated,
  features = "11-5402036-5403846",
  pt.size = 0.03,
  label = FALSE,
  label.size = 1, ncol = 2, combine = TRUE, coord.fixed = FALSE,
  cols = c("lightgrey", "red")) + ggtitle("FLT3:11-5402036-5403846")
p
ggsave(paste0("cluster9.FLT3.pdf"),path = "./plot.manuscript.final/findallmarkers.celltype.annotation" )

p <- FeaturePlot(
  object = integrated,
  features = "7-147186-147773",
  pt.size = 0.03,
  label = FALSE,
  label.size = 1, ncol = 2, combine = TRUE, coord.fixed = FALSE,
  cols = c("lightgrey", "red")) + ggtitle("IRF4:7-147186-147773")
p
ggsave(paste0("cluster9.IRF4.pdf"),path = "./plot.manuscript.final/findallmarkers.celltype.annotation" )

p <- FeaturePlot(
  object = integrated,
  features = "9-45631519-45635877",
  pt.size = 0.03,
  label = FALSE,
  label.size = 1, ncol = 2, combine = TRUE, coord.fixed = FALSE,
  cols = c("lightgrey", "red")) + ggtitle("CD3E:9-45631519-45635877")
p
ggsave(paste0("cluster10.CD3E.pdf"),path = "./plot.manuscript.final/findallmarkers.celltype.annotation" )


sub.genes <-  subset(pwAllsub,closest_genes == "CD19")#

p <- FeaturePlot(
  object = integrated,
  features = "3-18591541-18594649",
  pt.size = 0.03,
  label = FALSE,
  label.size = 1, ncol = 2, combine = TRUE, coord.fixed = FALSE,
  cols = c("lightgrey", "red")) + ggtitle("CD19:3-18591541-18594649")
p
ggsave(paste0("FeaturePlot.DAP.only.scATAC.CD19.pdf"),path = "./plot.manuscript.final/findallmarkers.celltype.annotation" )

sub.genes <-  subset(pwAllsub,closest_genes == "CSF1R")#

p <- FeaturePlot(
  object = integrated,
  features = "2-151099970-151108375",
  pt.size = 0.03,
  label = FALSE,
  label.size = 1, ncol = 2, combine = TRUE, coord.fixed = FALSE,
  cols = c("lightgrey", "red")) + ggtitle("CSF1R:2-151099970-151108375")
p
ggsave(paste0("FeaturePlot.DAP.only.scATAC.CSF1R.pdf"),path = "./plot.manuscript.final/findallmarkers.celltype.annotation" )

p <- FeaturePlot(
  object = integrated,
  features = "2-151108982-151113383",
  pt.size = 0.03,
  label = FALSE,
  label.size = 1, ncol = 2, combine = TRUE, coord.fixed = FALSE,
  cols = c("lightgrey", "red")) + ggtitle("CSF1R:2-151108982-151113383")
p
ggsave(paste0("FeaturePlot.DAP.only.scATAC.CSF1R.2.pdf"),path = "./plot.manuscript.final/findallmarkers.celltype.annotation" )

p <- FeaturePlot(
  object = integrated,
  features = "2-151125625-151130033",
  pt.size = 0.03,
  label = FALSE,
  label.size = 1, ncol = 2, combine = TRUE, coord.fixed = FALSE,
  cols = c("lightgrey", "red")) + ggtitle("CSF1R:2-151125625-151130033")
p
ggsave(paste0("FeaturePlot.DAP.only.scATAC.CSF1R.3.pdf"),path = "./plot.manuscript.final/findallmarkers.celltype.annotation" )

sub.genes <-  subset(pwAllsub,closest_genes == "CD14")#7,9
p <- FeaturePlot(
  object = integrated,
  features = "2-142344198-142344875",
  pt.size = 0.03,
  label = FALSE,
  label.size = 1, ncol = 2, combine = TRUE, coord.fixed = FALSE,
  cols = c("lightgrey", "red")) + ggtitle("CD14:2-142344198-142344875")
p
ggsave(paste0("FeaturePlot.DAP.only.scATAC.CD14.3.pdf"),path = "./plot.manuscript.final/findallmarkers.celltype.annotation" )


sub.genes <-  subset(pwAllsub,closest_genes == "CD4")#
p <- FeaturePlot(
  object = integrated,
  features = "5-63915302-63918975",
  pt.size = 0.03,
  label = FALSE,
  label.size = 1, ncol = 2, combine = TRUE, coord.fixed = FALSE,
  cols = c("lightgrey", "red")) + ggtitle("CD4:5-63915302-63918975")
p
ggsave(paste0("FeaturePlot.DAP.only.scATAC.CD4.pdf"),path = "./plot.manuscript.final/findallmarkers.celltype.annotation" )


sub.genes <-  subset(pwAllsub,closest_genes == "CD3E")#22 5

p <- FeaturePlot(
  object = integrated,
  features = "9-45618596-45622334",
  pt.size = 0.03,
  label = FALSE,
  label.size = 1, ncol = 2, combine = TRUE, coord.fixed = FALSE,
  cols = c("lightgrey", "red")) + ggtitle("CD3E:9-45618596-45622334")
p
ggsave(paste0("FeaturePlot.DAP.only.scATAC.CD3E.pdf"),path = "./plot.manuscript.final/findallmarkers.celltype.annotation" )


##pairwise find markers only pos

DAP.pig.2.4.2.30.findmarker.onlypos.p0.05 <- read.table("integrated.2.4ClustersPairwiseDE_LFC25FDR5_CellExpress20%inAtleast1Cluster.only.pos.2.30.no.logFC.txt")#1908425       9
# length(unique(DAP.pig.2.4.2.30.findmarker.onlypos.p0.05$gene))# 28471
all.genomic.regions <- DAP.pig.2.4.2.30.findmarker.onlypos.p0.05$gene
closest.genes.findmarker.only.pos <- ClosestFeature(integrated, regions = all.genomic.regions)#1908425      24
DAP.pig.2.4.2.30.findmarker.onlypos.p0.05$closest_genes <- closest.genes.findmarker.only.pos$gene_name
DAP.pig.2.4.2.30.findmarker.onlypos.p0.05$distance <- closest.genes.findmarker.only.pos$distance
DAP.pig.2.4.2.30.findmarker.onlypos.p0.05 <- subset(DAP.pig.2.4.2.30.findmarker.onlypos.p0.05,distance <= 5000)#1598336      11
DAP.pig.2.4.2.30.findmarker.onlypos.p0.05 <- DAP.pig.2.4.2.30.findmarker.onlypos.p0.05[!(is.na(DAP.pig.2.4.2.30.findmarker.onlypos.p0.05$closest_genes) | DAP.pig.2.4.2.30.findmarker.onlypos.p0.05$closest_genes==""), ]# 1265005      11
DAP.pig.2.4.2.30.findmarker.onlypos.p0.05 <- DAP.pig.2.4.2.30.findmarker.onlypos.p0.05[order(DAP.pig.2.4.2.30.findmarker.onlypos.p0.05$pop1),]#1265005      11

saveRDS(DAP.pig.2.4.2.30.findmarker.onlypos.p0.05, file = "DAP.pig.2.4.2.30.findmarker.onlypos.pairwise.p0.05.with.gene.name.rds")#, has distance col, 1265005      11
DAP.pig.2.4.2.30.findmarker.onlypos.p0.05 <- readRDS("DAP.pig.2.4.2.30.findmarker.onlypos.pairwise.p0.05.with.gene.name.rds")#, has distance col, 1265005      11

all.gene.marker.list <- list("CD19", "CD79A", "PAX5", "IRF4","PRDM1",
                             "CD14", "CD163", "NLRP3","TLR4","CSF1R","CD86","SIRPA",
                             "FLT3","FCER1A","SLA-DRB1","SLA-DRA",
                             "TCF4","IRF8","CD4","CD8B","XBP1","CLEC12A","CD93",
                             "CD3E","TRDC","CD4","CD8A","CD8B","CD2","PRF1",
                             "KLRB1","KLRK1","CD5","CD6","TYROBP","HCST")
length(all.gene.marker.list)
# [1] 36
all.result <- list()

for (i in 1:length(all.gene.marker.list))
{
  result <-  subset(DAP.pig.2.4.2.30.findmarker.onlypos.p0.05,closest_genes == all.gene.marker.list[[i]])#
  all.result[[i]] <- result
}

pwAll <- do.call(rbind, all.result)#10754    11
pwAll <- pwAll[order(pwAll$pop1,pwAll$pop2, pwAll$p_val_adj),]#10754    11


saveRDS(pwAll, file = "pwAll.withdistance.col.findmarker.pairwise.only.pos.rds")#max(pwAll$distance): 2602, 
pwAll <- readRDS("pwAll.withdistance.col.findmarker.pairwise.only.pos.rds")#max(pwAll$distance): 2602, 

pwAllsub <- pwAll[,c("gene","pop1", "pop2", "closest_genes", "distance", "avg_log2FC")]# 
saveRDS(pwAllsub, file = "DE.pig.resol2.4.2.30.findmarker.onlypos.p0.05.2.with.gene.name.dis.sub.rds")#
pwAllsub <- readRDS("DE.pig.resol2.4.2.30.findmarker.onlypos.p0.05.2.with.gene.name.dis.sub.rds")#

write.csv(x=pwAllsub,file = "DAP.pig.resol2.4.2.30.findmarker.onlypos.p0.05.with.36gene.marker.with.distance.col.sub.csv")
