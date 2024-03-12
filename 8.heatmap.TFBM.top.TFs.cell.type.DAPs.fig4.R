#7.3. heatmap TFBM result
list.files()
# TFs, p values, Log P values
celltype.top.TFs <- read.table("TFBM.celltype.result.top.20.TFs.each.cell.type.summary.txt", header = FALSE)


celltype.all.TFs.path <- c("path.../TFBM.celltype.output.all.11.tables/ASC.knownResults.txt",
                           "path.../TFBM.celltype.output.all.11.tables/B.knownResults.txt",
                           "path.../CD2negGD.knownResults.txt",
                           "path.../CD2posGD.knownResults.txt",
                           "path.../CD4posab.knownResults.txt",
                           "path.../CD8abPOSab.knownResults.txt",
                           "path.../CD8aPOSabT_NK.knownResults.txt",
                           "path.../cDCs.knownResults.txt",
                           "path.../Monocytes.knownResults.txt",
                           "path.../NK.knownResults.txt",
                           "path.../pDCs.knownResults.txt")

celltype.all.TFs <- read.table("TFBM.celltype.result.all.TFs.each.cell.type.summary.txt", header = FALSE)
dim(celltype.all.TFs)#[1] 428  33

names(celltype.top.TFs) <- c("TF.ASC",
                             "p.ASC",
                             "log.p.ASC",
                             "TF.B",
                             "p.B",
                             "log.p.B",
                             "TF.CD2negGD",
                             "p.CD2negGD",
                             "log.p.CD2negGD",
                             "TF.CD2posGD",
                             "p.CD2posGD",
                             "log.p.CD2posGD",
                             "TF.CD4posab",
                             "p.CD4posab",
                             "log.p.CD4posab",
                             "TF.CD8abPOSab",
                             "p.CD8abPOSab",
                             "log.p.CD8abPOSab",
                             "TF.CD8aPOSabT_NK",
                             "p.CD8aPOSabT_NK",
                             "log.p.CD8aPOSabT_NK",
                             "TF.cDCs",
                             "p.cDCs",
                             "log.p.cDCs",
                             "TF.Monocytes",
                             "p.Monocytes",
                             "log.p.Monocytes",
                             "TF.NK",
                             "p.NK",
                             "log.p.NK",
                             "TF.pDCs",
                             "p.pDCs",
                             "log.p.pDCs") # it is a factor

# dim(celltype.top.TFs)# 20 33

# celltype.all.TFs[,"log.p.pDCs"]

TFs.col <- seq(1, 33, 3)#1  4  7 10 13 16 19 22 25 28 31

celltype.top.TFs[,31]
# [1] "SpiB(ETS)/OCILY3-SPIB-ChIP-Seq(GSE56857)/Homer"                "PU.1:IRF8(ETS:IRF)/pDC-Irf8-ChIP-Seq(GSE66899)/Homer"         

library(stringr) #str_split

all.top20.TFs.all.celltypes <- "" #character
for (i in TFs.col)
{
  print(paste0("i is:" , i))
  celltype.top.TFs[,i] <- sub("\\(.*", "", celltype.top.TFs[,i])  
  all.top20.TFs.all.celltypes <- unique(append(celltype.top.TFs[,i], all.top20.TFs.all.celltypes))
}

all.top20.TFs.all.celltypes
position.empty <- length(all.top20.TFs.all.celltypes)
all.top20.TFs.all.celltypes <- all.top20.TFs.all.celltypes[-position.empty]
all.top20.TFs.all.celltypes# does not include the mepty string when this varialbe was created

write.table(celltype.top.TFs,file = "celltype.top.20.TFs.total.69.txt", row.names = FALSE, sep="\t", quote = FALSE)

TF.df <- data.frame(matrix(NA, nrow = 11, ncol = 69)) #11 cell types and 69 TFs
#"cDCs", "Monocytes", "CD4posab","CD2negGD","B","ASC","CD2posGD","CD8abPOSab","NK","pDCs","CD8aPOSabT_NK"

rownames(TF.df) <- c("ASC", "B","CD2negGD","CD2posGD","CD4posab","CD8abPOSab","CD8aPOSabT_NK","cDCs","Monocytes","NK","pDCs")
names(TF.df) <- all.top20.TFs.all.celltypes

for (i in 1:11)# the order of cell types () in TF.df is the same as that of  celltype.TFBM.known.result
{
  for (j in 1:69)
  {
    celltype.TFBM.known.result <- read.table(celltype.all.TFs.path[i])
    celltype.TFBM.known.result <- celltype.TFBM.known.result[,1:5]# motif.name, Consensus, p value, log p, q value
    celltype.TFBM.known.result <- celltype.TFBM.known.result[,-2]#remove Consensus column
    celltype.TFBM.known.result[,1] <- sub("\\(.*", "", celltype.TFBM.known.result[,1])  # only keep TF name
    names(celltype.TFBM.known.result) <- c("TF.name", "p.value", "Log.p", "q.value")
    
    sub <- subset(celltype.TFBM.known.result, TF.name == names(TF.df)[j])#50  2
    avg.q <- min(sub$q.value) # there are 2 rows when i =1, J =52

    TF.df[i,j] <- avg.q
  }
}


data <- as.matrix(t(TF.df)) # row is TF and column is cell type

library(RColorBrewer)
library(pheatmap)
library(ggplot2)


data <- ifelse( data == 0, 1e-10, data)

p <- pheatmap(-log10(data), 
              cluster_rows=T, 
              cluster_cols=T,
              margins = c(5, 5),
              treeheight_row = 150,
              treeheight_col = 80,
              labels_row= NULL, 
              fontsize_row= 30,
              fontsize_col= 40,
              cellwidth = 100,
              cellheight = 30,
              clustering_distance_rows = "euclidean",#euclidean, cannot clustering using correlation: In cor(t(mat)) : the standard deviation is zero
              clustering_distance_cols = "euclidean",#euclidean, cannot clustering using correlation: In cor(t(mat)) : the standard deviation is zero
              legend = T,
              display_numbers = matrix(ifelse(data < 0.05, "*", ""), nrow = nrow(data),ncol= ncol(data)),
              fontsize_number = 25,
              show_rownames = T, # show TFs names
              show_colnames = T, # show cell type names
              annotation_legend = TRUE,
              # number_color = "white",
              number_color = "green",
              scale = c("none"),
              angle_col = "45",
              # col=(brewer.pal(9,"Blues")))#if not rev (reverse), then higher number, the more significant
              col=(brewer.pal(9,"Reds")))#if not rev (reverse), then higher number, the more significant

ggsave(p, filename = "TFBM.TF.heatmap/celltype/pheatmap.q.value.cluster.by.euclidean.cell.type.scaled.none.neglog10.red.green.pdf", width = 30, height = 40)
