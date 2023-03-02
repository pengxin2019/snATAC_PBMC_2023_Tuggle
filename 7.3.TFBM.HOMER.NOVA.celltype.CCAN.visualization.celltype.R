

##################  visualization of TFBM of  cell type CCAN peaks  result  ###############################################################       
##################  visualization of TFBM of  cell type CCAN peaks  result  ############################################################### 

############  count the significant TF in  cell type CCANs ############  
############  count the significant TF in  cell type CCANs ############  

# celltype.ccan.peaks.list <- c("ASC", "B","CD2posGD","CD4posab","CD8abPOSab","cDCs","Monocytes","NK","pDCs") # it is a factor
33+14+8+13+26+42+11+28+69=244
all.path <- list.files(path='/work/abg/pyang19/scATAC_PBMC_analysis/TFBM_HOMER/input.output.files/output_files/celltype.ccan.peaks/B/TFBM.B.33.ccan.peak.300.bg.output',pattern='300.bg.output')

num.CCAN.peaks <- length(all.path)#14

library(stringr) #str_split
all.unique.sig.TFs.celltype <- "" #character
unique.significant.TF.num <- 0
hub.gene.with.sig.TF <- ""#character
num.unique.hub.gene.with.sig.TF <- 0

# for (i in 1:1)
for (i in 1:num.CCAN.peaks)
{
  hub.gene.name <- sub(".*gene", "", all.path[i]) # Extract characters after pattern,# [1] "B.peak.1-104812482-104818412in.ccan27.gene.TCF4.k.140.31peaks.corralated.0.05.avg.coaccess.0.184216517108376.300.bg.output"--> .TCF4.k.140.31peaks.corralated.0.05.avg.coaccess.0.184216517108376.300.bg.output"
  hub.gene.name <- sub('.', '', hub.gene.name)#remove the first dot: TCF4.k.140.31peaks.corralated.0.05.avg.coaccess.0.184216517108376.300.bg.output
  hub.gene.name <- gsub("\\..*","",hub.gene.name)# remove the substring after the dot
  
  if (file.exists(paste0(all.path[i],"/knownResults.txt")))
  {
    TFBM.CCAN.output <- read.delim(paste0(all.path[i],"/knownResults.txt"), header=TRUE, sep = "\t")
    names(TFBM.CCAN.output) <- c("TF.name","Consensus", "p", "logp", "q.Benjamini","num.Target.Sequences.with.Motif", "percent.of.Targets.Sequences.with.Motif", "num.Background.Sequences.with.Motif", "percent.of.Background.Sequences.with.Motit")
    TFBM.CCAN.output[,"TF.name"] <- sub("\\(.*", "", TFBM.CCAN.output[,"TF.name"])  
    sub.significant <- subset(TFBM.CCAN.output,q.Benjamini <= 0.1)
    
    
    if(dim(sub.significant)[1] != 0)
    {
      hub.gene.with.sig.TF <- append(hub.gene.name, hub.gene.with.sig.TF)
      num.unique.hub.gene.with.sig.TF <- num.unique.hub.gene.with.sig.TF+1
      all.unique.sig.TFs.celltype <- unique(append(sub.significant[,"TF.name"], all.unique.sig.TFs.celltype))
      total.num.unique.TF.sig.celltype <- length(all.unique.sig.TFs.celltype) -1
      
      print(paste0("i: ", i, "; ", hub.gene.name, ": ", sub.significant$TF.name, "; q values: ",sub.significant$q.Benjamini))
    }
  }else 
  {
    i= i+1
  }
}


TF.sig.q.cell.type.ccan.TFBM.summary <- data.frame(matrix(1.0000, nrow = num.unique.hub.gene.with.sig.TF, ncol = total.num.unique.TF.sig.celltype)) #rownames are hub genes; columns are significiant TFs
names(TF.sig.q.cell.type.ccan.TFBM.summary) <- all.unique.sig.TFs.celltype[1:total.num.unique.TF.sig.celltype]
rownames(TF.sig.q.cell.type.ccan.TFBM.summary) <- hub.gene.with.sig.TF[1:num.unique.hub.gene.with.sig.TF]
TF.sig.q.cell.type.ccan.TFBM.summary


############### fill the dataframe: significant cells in this TF.sig.q.cell.type.ccan.TFBM.summary ############### 
############### fill the dataframe: significant cells in this TF.sig.q.cell.type.ccan.TFBM.summary ############### 

for (i in 1:num.CCAN.peaks)
{
  hub.gene.name <- sub(".*gene", "", all.path[i]) # Extract characters after pattern,# [1] "B.peak.1-104812482-104818412in.ccan27.gene.TCF4.k.140.31peaks.corralated.0.05.avg.coaccess.0.184216517108376.300.bg.output"--> .TCF4.k.140.31peaks.corralated.0.05.avg.coaccess.0.184216517108376.300.bg.output"
  hub.gene.name <- sub('.', '', hub.gene.name)#remove the first dot: TCF4.k.140.31peaks.corralated.0.05.avg.coaccess.0.184216517108376.300.bg.output
  hub.gene.name <- gsub("\\..*","",hub.gene.name)# remove the substring after the dot
  
  if (file.exists(paste0(all.path[i],"/knownResults.txt")))
  {
    TFBM.CCAN.output <- read.delim(paste0(all.path[i],"/knownResults.txt"), header=TRUE, sep = "\t")
    names(TFBM.CCAN.output) <- c("TF.name","Consensus", "p", "logp", "q.Benjamini","num.Target.Sequences.with.Motif", "percent.of.Targets.Sequences.with.Motif", "num.Background.Sequences.with.Motif", "percent.of.Background.Sequences.with.Motit")
    TFBM.CCAN.output[,"TF.name"] <- sub("\\(.*", "", TFBM.CCAN.output[,"TF.name"])  
    sub.significant <- subset(TFBM.CCAN.output,q.Benjamini <= 0.1)
    
    if(dim(sub.significant)[1] != 0)
    {
      print(paste0("i: ", i, "; ", hub.gene.name, ": ", sub.significant$TF.name, "; q values: ",sub.significant$q.Benjamini))
      
      for (j in 1:dim(TF.sig.q.cell.type.ccan.TFBM.summary)[1])# all sig TF in this cell type
      {
        if (hub.gene.name == rownames(TF.sig.q.cell.type.ccan.TFBM.summary)[j]) 
        {
          for (m in 1: dim(sub.significant)[1])# row of sub.significant is TF
          {
            for (n in 1: dim(TF.sig.q.cell.type.ccan.TFBM.summary)[2])# colm of TF.sig.q.cell.type.ccan.TFBM.summary is TF
            {
              if (sub.significant[m, "TF.name"] == names(TF.sig.q.cell.type.ccan.TFBM.summary)[n])# if the TF match
              {
                TF.sig.q.cell.type.ccan.TFBM.summary[j,n] = sub.significant[m, "q.Benjamini"]
              }
            }
          }
        }
      }
    }
  }else
  {
    i = i+1
  }
}

write.table(TF.sig.q.cell.type.ccan.TFBM.summary,file = "Monocytescell.celltype.14CCANS.TFBM.plot.2.hub.gene.5TFs.summary.txt", sep="\t", quote = FALSE)


############### make heatmap on local laptop cell type ############################################################  
############### make heatmap on local laptop cell type ############################################################  

#LOCAL LAPTOP
setwd("/Users/15612770360163.com/Desktop/scATAC/scATAC.submission2/scATAC.Satija.cellranger.1.2.batch1.batch2/newest.scATAC.PBMC.files.manuscript.final/plot.manuscript.final/TFBM.celltype.ccan.peaks.plot/Monocytes")

celltype.ccan.TFBM.summary <- read.table("Monocytescell.celltype.14CCANS.TFBM.plot.2.hub.gene.5TFs.summary.txt", header=TRUE, sep = "")

# data <- as.matrix(TF.sig.q.cell.type.ccan.TFBM.summary) # 
data <- as.matrix(celltype.ccan.TFBM.summary) # 
data <- ifelse( data == 0, 1e-10, data)

library(RColorBrewer)
library(pheatmap)
library(ggplot2)
p <- pheatmap(-log10(data), 
              cluster_rows=T, 
              cluster_cols=T,
              labels_row= NULL, 
              clustering_distance_rows = "correlation",#euclidean
              clustering_distance_cols = "correlation",#euclidean
              legend = T,
              display_numbers = matrix(ifelse(data < 0.1, "*", ""), nrow = nrow(data),ncol= ncol(data)),
              show_rownames = T, # show TFs names
              show_colnames = T, # show cell type names
              annotation_legend = TRUE,
              number_color = "white",
              scale = c("none"),
              col=(brewer.pal(9,"Reds")))#if not rev (reverse), then higher number, the more significant

# ggsave(p, filename = "Monocytes.celltypeCCAN.TFBM.14ccans.2hubgenes.5TFs.-log10.0.1.euclidean.pdf")
ggsave(p, filename = "Monocytes.celltypeCCAN.TFBM.14ccans.2hubgenes.5TFs.-log10.0.1.correlation.pdf")

############### make heatmap on local laptop ALL cell type ############################################################  
############### make heatmap on local laptop ALL cell type ############################################################  


path <- c("/Users/15612770360163.com/Desktop/scATAC/scATAC.submission2/scATAC.Satija.cellranger.1.2.batch1.batch2/newest.scATAC.PBMC.files.manuscript.final/plot.manuscript.final/TFBM.celltype.ccan.peaks.plot/ASC",
          "/Users/15612770360163.com/Desktop/scATAC/scATAC.submission2/scATAC.Satija.cellranger.1.2.batch1.batch2/newest.scATAC.PBMC.files.manuscript.final/plot.manuscript.final/TFBM.celltype.ccan.peaks.plot/B",
          "/Users/15612770360163.com/Desktop/scATAC/scATAC.submission2/scATAC.Satija.cellranger.1.2.batch1.batch2/newest.scATAC.PBMC.files.manuscript.final/plot.manuscript.final/TFBM.celltype.ccan.peaks.plot/CD2posGD",
          "/Users/15612770360163.com/Desktop/scATAC/scATAC.submission2/scATAC.Satija.cellranger.1.2.batch1.batch2/newest.scATAC.PBMC.files.manuscript.final/plot.manuscript.final/TFBM.celltype.ccan.peaks.plot/CD8abPOSab",
          "/Users/15612770360163.com/Desktop/scATAC/scATAC.submission2/scATAC.Satija.cellranger.1.2.batch1.batch2/newest.scATAC.PBMC.files.manuscript.final/plot.manuscript.final/TFBM.celltype.ccan.peaks.plot/cDCs",
          "/Users/15612770360163.com/Desktop/scATAC/scATAC.submission2/scATAC.Satija.cellranger.1.2.batch1.batch2/newest.scATAC.PBMC.files.manuscript.final/plot.manuscript.final/TFBM.celltype.ccan.peaks.plot/Monocytes",
          "/Users/15612770360163.com/Desktop/scATAC/scATAC.submission2/scATAC.Satija.cellranger.1.2.batch1.batch2/newest.scATAC.PBMC.files.manuscript.final/plot.manuscript.final/TFBM.celltype.ccan.peaks.plot/NK",
          "/Users/15612770360163.com/Desktop/scATAC/scATAC.submission2/scATAC.Satija.cellranger.1.2.batch1.batch2/newest.scATAC.PBMC.files.manuscript.final/plot.manuscript.final/TFBM.celltype.ccan.peaks.plot/pDCs")

all.hub.genes <- ""
all.TFs <- ""
cell.type <- c("ASC", "B", "CD2posGD", "CD8abPOSab", "cDCs", "Monocytes","NK", "pDCs")

for (i in 1:length(path))
{
  file <- list.files(path =path[i],  pattern='summary.txt')
  result.one.celltype <-  read.table(paste0(path[i],"/", file))
  
  all.hub.genes <- append(all.hub.genes, rownames(result.one.celltype))
  all.TFs <- append(all.TFs, names(result.one.celltype))
  print(rownames(result.one.celltype))
}


all.hub.genes <- all.hub.genes[-1]# remove "" non-unique values when setting 'row.names': ‘IGSF8’, ‘PTPRE’ 
all.TFs <- all.TFs[-1] # remove ""

length(all.hub.genes)#45
length(unique(all.hub.genes))#43

length(all.TFs)#70
length(unique(all.TFs))#41

unique.all.hub.genes <- unique(all.hub.genes)
unique.all.TFs <- unique(all.TFs)

all.celltype.TF.hub.gene.summary <- data.frame(matrix(1.0000, nrow = length(all.hub.genes), ncol = length(unique.all.TFs))) #rownames are hub genes; columns are significiant TFs
names(all.celltype.TF.hub.gene.summary) <- unique.all.TFs
rownames(all.celltype.TF.hub.gene.summary) <- c("SHCBP1", "IGSF8.ASC", "ARID3A", "RPN1",  "SIT1", #
                                                "POU2AF1", "NAPSA", "ZBTB32", "NCF4", "LYN", "HCK", "TCF4", #
                                                "CST7", #
                                                "ARL4C",  #
                                                "CFP","SERPINB1", "PLBD1", "ZNF385A", "SLC25A24", "CFD", "CORO1C", "PTPRE.cDCs", "GRN", 
                                                "HK3", "ADGRE5", 
                                                "ADGRG1", "CD244", "GBP2", "GNLY",  "CD8A", "ITGAM",  "PRF1","PTPRE.NK",
                                                "TNF","TYROBP","MCOLN2","SERHL2","EPS8","IGSF8.pDCs","ATP1B1","TRAM1","CYB561A3","COBLL1","SLC25A38","RHBDF2" )
all.celltype.TF.hub.gene.summary
cell.type <- c("ASC", "B", "CD2posGD", "CD8abPOSab", "cDCs", "Monocytes","NK", "pDCs")


for (i in 1:length(path))
{
  file <- list.files(path =path[i],  pattern='summary.txt')
  result.one.celltype <-  read.table(paste0(path[i],"/", file))
  
  for (m in 1:dim(result.one.celltype)[1]) # rows: hub genes
  {
    for (n in 1:dim(result.one.celltype)[2])# col TFs
    {
      for (j in 1:dim(all.celltype.TF.hub.gene.summary)[1])# rownames of all.celltype.TF.hub.gene.summary
      {
        for (k in 1: dim(all.celltype.TF.hub.gene.summary)[2])# # col TFs of all.celltype.TF.hub.gene.summary
        {
          if (rownames(result.one.celltype)[m] == rownames(all.celltype.TF.hub.gene.summary)[j] & names(result.one.celltype)[n] == names(all.celltype.TF.hub.gene.summary)[k])
          {
            all.celltype.TF.hub.gene.summary[j,k] = result.one.celltype[m,n]
          }
        }
      }
    }
  }
}

all.celltype.TF.hub.gene.summary

all.celltype.TF.hub.gene.summary[2,] # cell type i= 1, p(hub IGSF8, CTCF ) =  0.0004,  p(hub IGSF8, BORIS ) =  0.0004,and cell type 8, common hub genes IGSF8
# cell type i= 5 and cell type 7, common hub genes PTPRE
rownames(all.celltype.TF.hub.gene.summary)
names(all.celltype.TF.hub.gene.summary)


for (l in 1:length(path))
{
  file <- list.files(path =path[l],  pattern='summary.txt')
  result.one.celltype <-  read.table(paste0(path[l],"/", file))
  
  rownames(result.one.celltype)
  names(result.one.celltype)
}

all.celltype.TF.hub.gene.summary["IGSF8.ASC","CTCF"] = 0.0004
all.celltype.TF.hub.gene.summary["IGSF8.ASC","BORIS"] = 0.0004

all.celltype.TF.hub.gene.summary["IGSF8.pDCs","CTCF"] = 0.0240
all.celltype.TF.hub.gene.summary["IGSF8.pDCs","SpiB"] = 0.024

all.celltype.TF.hub.gene.summary["PTPRE.cDCs","BORIS"] = 0.0552
all.celltype.TF.hub.gene.summary["PTPRE.NK","CTCF"] = 0.0913
  
dim(all.celltype.TF.hub.gene.summary)#45 41

write.table(all.celltype.TF.hub.gene.summary,file = "all.celltype.hub.gene.TFs.summary.new.txt", sep="\t", quote = FALSE)

data <- as.matrix(t(all.celltype.TF.hub.gene.summary)) # row is TF and column is cell type

data <- ifelse( data == 0, 1e-10, data)
library(RColorBrewer)
library(pheatmap)
library(ggplot2)


### make larger legend
p <- pheatmap(-log10(data), 
              cluster_rows=T, 
              cluster_cols=F,
              margins = c(5, 5),
              treeheight_row = 100,
              # treeheight_col = 80,
              labels_row= NULL, 
              fontsize_row= 20,
              fontsize_col= 20,
              cellwidth = 20,
              cellheight = 20,
              clustering_distance_rows = "euclidean",
              # clustering_distance_cols = "correlation",
              legend = T,
              display_numbers = matrix(ifelse(data < 0.1, "*", ""), nrow = nrow(data),ncol= ncol(data)),
              fontsize_number = 25,
              show_rownames = T, # show TFs names
              show_colnames = T, # show cell type names
              annotation_legend = TRUE,
              # number_color = "white",
              number_color = "red",
              scale = c("none"),
              angle_col = "45",
              # col=(brewer.pal(9,"Reds")))#if not rev (reverse), then higher number, the more significant
              col=(brewer.pal(9,"Blues")))#if not rev (reverse), then higher number, the more significant

# ggsave(p, filename = "/Users/15612770360163.com/Desktop/scATAC/scATAC.submission2/scATAC.Satija.cellranger.1.2.batch1.batch2/newest.scATAC.PBMC.files.manuscript.final/plot.manuscript.final/TFBM.celltype.ccan.peaks.plot/all.cell.type.41.unique.sig.TF.45nonunique.hub.gene.correlation.celltypenocluster.0.1.pdf", width = 45, height = 45)
ggsave(p, filename = "/Users/15612770360163.com/Desktop/scATAC/scATAC.submission2/scATAC.Satija.cellranger.1.2.batch1.batch2/newest.scATAC.PBMC.files.manuscript.final/plot.manuscript.final/TFBM.celltype.ccan.peaks.plot/all.cell.type.41.unique.sig.TF.45nonunique.hub.gene.correlation.celltypenocluster.0.1.blue.large.big.legend.pdf", width = 49, height = 49)

