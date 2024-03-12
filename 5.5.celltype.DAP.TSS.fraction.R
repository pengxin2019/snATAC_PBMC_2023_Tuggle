# check whether cell type DAP includes TSS of gene or not
library(org.Ss.eg.db)
library(ChIPseeker) #ChIPseeker v1.26.2 
library(GenomicRanges)
library(org.Sscrofa.eg.db) # created by  Penny using custome script  0.generate.OrgDb.geneannotation.R

gtf <- rtracklayer::import('path.../Sus_scrofa.Sscrofa11.1.102.gtf')
library(GenomicFeatures)
TxDb <- makeTxDbFromGRanges(gtf)
library(biomaRt)

## check how many DAP includes TSS
B.DAP <- "B.DAP.findallmarkers.integration.scRNA.txt"
ASC.DAP <- "ASC.DAP.findallmarkers.integration.scRNA.txt"
CD2negGD.DAP <- "CD2negGD.DAP.findallmarkers.integration.scRNA.txt"
CD2posGD.DAP <- "CD2posGD.DAP.findallmarkers.integration.scRNA.txt"
CD4posab.DAP <- "CD4posab.DAP.findallmarkers.integration.scRNA.txt"
CD8abPOSab.DAP <- "CD8abPOSab.DAP.findallmarkers.integration.scRNA.txt"
CD8aPOSabT_NK.DAP <- "CD8aPOSabT_NK.DAP.findallmarkers.integration.scRNA.txt"
cDCs.DAP <- "cDCs.DAP.findallmarkers.integration.scRNA.txt"
Monocytes.DAP <- "Monocytes.DAP.findallmarkers.integration.scRNA.txt"
NK.DAP <- "NK.DAP.findallmarkers.integration.scRNA.txt"
pDCs.DAP <- "pDCs.DAP.findallmarkers.integration.scRNA.txt"

celltype.DAPs <- c(B.DAP, ASC.DAP, CD2negGD.DAP, CD2posGD.DAP,CD4posab.DAP, CD8abPOSab.DAP, CD8aPOSabT_NK.DAP, cDCs.DAP, Monocytes.DAP, NK.DAP, pDCs.DAP)
names(celltype.DAPs) <- c("B.DAP", "ASC.DAP", "CD2negGD.DAP", "CD2posGD.DAP","CD4posab.DAP", "CD8abPOSab.DAP", "CD8aPOSabT_NK.DAP", "cDCs.DAP","Monocytes.DAP", "NK.DAP", "pDCs.DAP")

library(stringr)
library(GenomicRanges)

for (i in 1:11)
{
  DAP.input <- read.table(celltype.DAPs[[i]])
  DAP.list <- str_split(DAP.input$V1, "-") # split it into 3 columns, stringr
  num.DAP <- length(DAP.list)
  
  #generate an empty df, #row = # of DAP in this cell type, 3 cols: chrom, start, end
  DAP.df <- data.frame(matrix(nrow=num.DAP, ncol=3))
  
  #fill each cell of this DAP.df
  for (j in 1:num.DAP)
  {
    for (k in 1:3)
    {
      DAP.df[j,k] <- DAP.list[[j]][k]
    }
  }
  names(DAP.df) <- c("seqnames", "start", "end")
  
  DAP.gr <- makeGRangesFromDataFrame(DAP.df, keep.extra.columns=TRUE)#GenomicRanges
  
  peakAnno <- annotatePeak(DAP.gr, tssRegion=c(-3000, 3000),
                           TxDb=TxDb, annoDb="org.Sscrofa.eg.db")
  
  peakAnno.df <- data.frame(peakAnno)
  
  subsetTSS <- subset(peakAnno.df, distanceToTSS == "0")
  
  num.unique.gene <- length(unique(subsetTSS$geneId))# geneId for example: ENSSSCG00000005743
  
  subsetTSS$DAP.in.TSS <- paste(subsetTSS$seqnames, subsetTSS$start, subsetTSS$end, sep="-") # output DAP list in which all peaks covers TSS of a gene
  DAP.only.in.TSS <- as.data.frame(subsetTSS$DAP.in.TSS) # output DAP list in which all peaks covers TSS of a gene
  names(DAP.only.in.TSS) <- "peaks"
  
  write.table(DAP.only.in.TSS, file = paste0(names(celltype.DAPs)[i], ".DAP.only.in.TSS.txt"), sep = "\t",row.names = FALSE,col.names=FALSE,quote = FALSE )
  
  print(paste0(names(celltype.DAPs)[i], ": number of unique DAP in TSS is: ",length(unique(DAP.only.in.TSS$peaks)) ," ; number of unique gene (geneID) who has a DAP in TSS is ",  num.unique.gene))
}

