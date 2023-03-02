
##pig scATAC 2.4 2:30
# cd /work/abg/pyang19/opt/pig.6798.6800.PBMC.Satija.pipeline.result/pbmc.1x.2x.cellrangeratac1.2.0.wd/strict.filter.new.annotation/strict.filter.new.annotation.common.peaks/commonpeaks.110444.new.directory
setwd("/work/abg/pyang19/opt/pig.6798.6800.PBMC.Satija.pipeline.result/pbmc.1x.2x.cellrangeratac1.2.0.wd/strict.filter.new.annotation/strict.filter.new.annotation.common.peaks/commonpeaks.110444.new.directory")

library(Signac)
library(Seurat)
library(SeuratWrappers)
library(ggplot2)
library(patchwork)
library(cicero)
library(GenomicRanges) #makeGRangesFromDataFrame
library(stringr) #str_split
library(BSgenome.Sscrofa.UCSC.susScr11)
set.seed(1234)

integrated <- readRDS("integrated.withRNA.assay.annotated.integration.scRNA.celltype2.4.2.30.after.filter.new.idents.RNA2NRA.wo.unknown.rds")#resol 2.4, +-2000, 2.30
DefaultAssay(integrated) <- 'peaks'

#### ## ## ##  cell type specific co-accessibility regions

all.cell.types.list <- c("ASC", "B","CD2negGD","CD2posGD","CD4posab","CD8abPOSab","CD8aPOSabT_NK","cDCs","Monocytes","NK","pDCs") # it is a factor
celltype.DAP.onlyin.TSS.path <- c("/work/abg/pyang19/opt/pig.6798.6800.PBMC.Satija.pipeline.result/pbmc.1x.2x.cellrangeratac1.2.0.wd/ArchR_wd_Nova/ASC.DAP.DAP.only.in.TSS.txt",
                                  "/work/abg/pyang19/opt/pig.6798.6800.PBMC.Satija.pipeline.result/pbmc.1x.2x.cellrangeratac1.2.0.wd/ArchR_wd_Nova/B.DAP.DAP.only.in.TSS.txt",
                                  "/work/abg/pyang19/opt/pig.6798.6800.PBMC.Satija.pipeline.result/pbmc.1x.2x.cellrangeratac1.2.0.wd/ArchR_wd_Nova/CD2negGD.DAP.DAP.only.in.TSS.txt",
                                  "/work/abg/pyang19/opt/pig.6798.6800.PBMC.Satija.pipeline.result/pbmc.1x.2x.cellrangeratac1.2.0.wd/ArchR_wd_Nova/CD2posGD.DAP.DAP.only.in.TSS.txt",
                                  "/work/abg/pyang19/opt/pig.6798.6800.PBMC.Satija.pipeline.result/pbmc.1x.2x.cellrangeratac1.2.0.wd/ArchR_wd_Nova/CD4posab.DAP.DAP.only.in.TSS.txt",
                                  "/work/abg/pyang19/opt/pig.6798.6800.PBMC.Satija.pipeline.result/pbmc.1x.2x.cellrangeratac1.2.0.wd/ArchR_wd_Nova/CD8abPOSab.DAP.DAP.only.in.TSS.txt",
                                  "/work/abg/pyang19/opt/pig.6798.6800.PBMC.Satija.pipeline.result/pbmc.1x.2x.cellrangeratac1.2.0.wd/ArchR_wd_Nova/CD8aPOSabT_NK.DAP.DAP.only.in.TSS.txt",
                                  "/work/abg/pyang19/opt/pig.6798.6800.PBMC.Satija.pipeline.result/pbmc.1x.2x.cellrangeratac1.2.0.wd/ArchR_wd_Nova/cDCs.DAP.DAP.only.in.TSS.txt",
                                  "/work/abg/pyang19/opt/pig.6798.6800.PBMC.Satija.pipeline.result/pbmc.1x.2x.cellrangeratac1.2.0.wd/ArchR_wd_Nova/Monocytes.DAP.DAP.only.in.TSS.txt",
                                  "/work/abg/pyang19/opt/pig.6798.6800.PBMC.Satija.pipeline.result/pbmc.1x.2x.cellrangeratac1.2.0.wd/ArchR_wd_Nova/NK.DAP.DAP.only.in.TSS.txt",
                                  "/work/abg/pyang19/opt/pig.6798.6800.PBMC.Satija.pipeline.result/pbmc.1x.2x.cellrangeratac1.2.0.wd/ArchR_wd_Nova/pDCs.DAP.DAP.only.in.TSS.txt")

all.DAP.files <- list.files(path='/work/abg/pyang19/opt/pig.6798.6800.PBMC.Satija.pipeline.result/pbmc.1x.2x.cellrangeratac1.2.0.wd/strict.filter.new.annotation/strict.filter.new.annotation.common.peaks/commonpeaks.110444.new.directory',pattern='.DAP.findallmarkers.integration.scRNA.txt')

all.DE.files <- list.files(path='/work/abg/pyang19/opt/pig.6798.6800.PBMC.Satija.pipeline.result/pbmc.1x.2x.cellrangeratac1.2.0.wd/strict.filter.new.annotation/strict.filter.new.annotation.common.peaks/commonpeaks.110444.new.directory',pattern='.scRNA.DE.gene.df.')

tmpdir <- "CCAN.all.files.all.celltype.batch"
if (!dir.exists(tmpdir))
{
  dir.create(tmpdir)
}

tmpdir <- "CCAN.all.files.all.celltype.batch/CCAN.all.txt.file.DAP.and.peaks.in.same.ccan"
if (!dir.exists(tmpdir))
{
  dir.create(tmpdir)
}

for (i in 1:length(all.cell.types.list)) #1:11 cell types
{
  cell.type.ccan.results <- list()
  cell.type.ccan.results.score <- list()
  
  celltype <- all.cell.types.list[i]
  num.ccan.all.chrom.TSS.DAP.DEG <- 0
  
  #create a folder to contain all the txt files for a specific cell type
  tmpdir <- paste0("CCAN.all.files.all.celltype.batch/CCAN.all.txt.file.DAP.and.peaks.in.same.ccan/",celltype)
  
  if (!dir.exists(tmpdir))
  {
    dir.create(tmpdir)
  }

  tmpdir <- paste0("CCAN.all.files.all.celltype.batch/CCAN.all.txt.file.DAP.and.peaks.in.same.ccan/",celltype, "/gene.name.known.DEG.TSS.DAP.102")
  
  if (!dir.exists(tmpdir))
  {
    dir.create(tmpdir)
  }
  
  sub_seurat <- subset(x = integrated, idents = all.cell.types.list[i])
  sub.cds <- as.cell_data_set(x = sub_seurat) #convert to CellDataSet format
  sub.cicero <- make_cicero_cds(sub.cds, reduced_coordinates = reducedDims(sub.cds)$UMAP)
  DAP <- read.table(all.DAP.files[i])
  
  DEG <- read.table(all.DE.files[i])
  DEG <- as.data.frame(DEG[-1,])
  names(DEG) <- "gene"
  
  each.DAP.onlyin.TSS <- read.table(celltype.DAP.onlyin.TSS.path[i])

  each.DAP.onlyin.TSS.list <- str_split(each.DAP.onlyin.TSS$V1, "-") # split it into 3 columns, stringr
  num.each.DAP.onlyin.TSS <- length(each.DAP.onlyin.TSS.list)
  #generate an empty df, #row = # of DAP in this cell type, 3 cols: chrom, start, end
  each.DAP.onlyin.TSS.df <- data.frame(matrix(nrow = num.each.DAP.onlyin.TSS, ncol=3))
  #fill each cell of this DAP.df
  for (n in 1:num.each.DAP.onlyin.TSS)
  {
    for (m in 1:3)
    {
      each.DAP.onlyin.TSS.df[n,m] <- each.DAP.onlyin.TSS.list[[n]][m]
    }
  }
  names(each.DAP.onlyin.TSS.df) <- c("seqnames", "start", "end")
  each.DAP.onlyin.TSS.gr <- makeGRangesFromDataFrame(each.DAP.onlyin.TSS.df, keep.extra.columns=TRUE)# this each.DAP.in.TSS.gr will be used for ClosestFeature functioin
  closest_genes_each.DAP.in.TSS <- ClosestFeature(integrated, regions = each.DAP.onlyin.TSS.gr)#
  common.genes.in.DAP.DE <- unique(intersect(closest_genes_each.DAP.in.TSS$gene_name,DEG$gene))
  write.table(common.genes.in.DAP.DE,file = paste0("CCAN.all.files.all.celltype.batch/CCAN.all.txt.file.DAP.and.peaks.in.same.ccan/", celltype,"/common.",length(common.genes.in.DAP.DE), ".genes.onlyTSS.DAP.DEG.txt"))
  
  for (j in 1:19)#seqlengths(BSgenome.Sscrofa.UCSC.susScr11)[19] is chrX, seqlengths(BSgenome.Sscrofa.UCSC.susScr11)[20] is chrY
  {
    all.ccan.one.chrom.TSS.DAP.DEG <- "" #character
    num.ccan.one.chrom.TSS.DAP.DEG <- 0
    
    genome <- seqlengths(BSgenome.Sscrofa.UCSC.susScr11)[j]#chrom i
    genome.df <- data.frame("chr" = names(genome), "length" = genome)#two columns, the first is the chromosome name (ex. "chr1") and the second is the chromosome length in base pairs
    cicero_cons <- run_cicero(sub.cicero, genomic_coords = genome.df, sample_num = 100)
    cicero_cons.new <- cicero_cons[which(cicero_cons$Peak1 %in% DAP$V1 | cicero_cons$Peak2 %in% DAP$V1),]#13808     3, 13808/2=6904 < 6988 below, so some of rows in cicero_cons are: both peak1 and peak2 are DAP,also there are some dupilcates
    cicero_cons.new <- cicero_cons.new[!is.na(cicero_cons.new$coaccess),]
    
    # cicero_cons.new <- subset(cicero_cons.new, coaccess > 0.01)
    cicero_cons.new <- subset(cicero_cons.new, coaccess > 0.05)
    cicero_cons.new <- cicero_cons.new[order(cicero_cons.new$coaccess),]#3494    3
    ccans.DAP <- generate_ccans(cicero_cons.new)#Post process cicero co-accessibility scores to extract modules of sites that are co-accessible.
    
    for (k in 1:num.each.DAP.onlyin.TSS) # for each peak in each.DAP.onlyin.TSS/ each.DAP.onlyin.TSS.df/each.DAP.onlyin.TSS.gr
    {
      ccan.info <- subset(ccans.DAP, Peak == each.DAP.onlyin.TSS[k,])# DAP of PAX5,  ccan.info is a df
      # if there is 0 row in the ccan.info, it means the kth peak is not assigned to a ccan, if is.na(closest_genes_each.DAP.in.TSS$gene_name[k]) is true, its because of difference between 102 and 105,
      #if !(closest_genes_each.DAP.in.TSS$gene_name[k] %in% common.genes.in.DAP.DE) is true, the nearest genes is not DE
      if (dim(ccan.info)[1] == 0 | is.na(closest_genes_each.DAP.in.TSS$gene_name[k]) | !(closest_genes_each.DAP.in.TSS$gene_name[k] %in% common.genes.in.DAP.DE)) 
      {
        k = k+1
      }else #  it means the kth TSSDAP peak is  assigned to a ccan, it has a gene name, and this gene is DEG
      {
        all.ccan.one.chrom.TSS.DAP.DEG  <- unique(append(as.character(ccan.info$CCAN),all.ccan.one.chrom.TSS.DAP.DEG)) #in case different each.DAP.onlyin.TSS[k,] are assigned in the same CCAN
        cicero_cons.new.correlated.with.one.peak <- subset(cicero_cons.new, Peak1 == each.DAP.onlyin.TSS[k,] | Peak2 == each.DAP.onlyin.TSS[k,])# find all peaks correlated with each.DAP.onlyin.TSS[k,]
        one.peak.correlated.all <- as.data.frame(unique(append(cicero_cons.new.correlated.with.one.peak$Peak1, cicero_cons.new.correlated.with.one.peak$Peak2)))
        # dim(one.peak.correlated.all)#72  1, all these peaks are correlated to 1-237781345-237790424
        names(one.peak.correlated.all) <- "peaks"
           
        ccan.num <- ccan.info$CCAN  #numeric
        
        ccan <- subset(ccans.DAP, CCAN == ccan.num)#50  2
               # # include the TSS DAP itself:
        same.ccan.correlated.with.one.peak <- as.data.frame(intersect(one.peak.correlated.all$peaks,ccan$Peak))#32 1, in ccan6 which include the TSS of PAX5,  all DAPs that are correlated with this TSS
        # dim(same.ccan.correlated.with.one.peak)#32
        names(same.ccan.correlated.with.one.peak) <- "peaks"
   
        cell.type.ccan.results[[k]] <- same.ccan.correlated.with.one.peak 
        names(cell.type.ccan.results)[[k]] <-  closest_genes_each.DAP.in.TSS$gene_name[k]
        cell.type.ccan.results.score.df <-  subset(cicero_cons.new.correlated.with.one.peak, Peak1 %in% same.ccan.correlated.with.one.peak$peaks)
        cell.type.ccan.results.score[[k]] <- mean(cell.type.ccan.results.score.df$coaccess)
        write.table(same.ccan.correlated.with.one.peak$peaks,file = paste0("CCAN.all.files.all.celltype.batch/CCAN.all.txt.file.DAP.and.peaks.in.same.ccan/", celltype, "/gene.name.known.DEG.TSS.DAP.102/",celltype,".peak.", each.DAP.onlyin.TSS[k,], "in.ccan",ccan.num,".gene.",closest_genes_each.DAP.in.TSS$gene_name[k],".k.",k,".", dim(same.ccan.correlated.with.one.peak)[1], "peaks.corralated.0.05.txt"), row.names = FALSE, sep="\t", quote = FALSE)#coaccess > 0.05
      }
      }#this is the for k loop
    
    num.ccan.one.chrom.TSS.DAP.DEG <- length(all.ccan.one.chrom.TSS.DAP.DEG) - 1 # "" is excluded 
    num.ccan.all.chrom.TSS.DAP.DEG <- num.ccan.all.chrom.TSS.DAP.DEG + num.ccan.one.chrom.TSS.DAP.DEG
    all.ccan.one.chrom.TSS.DAP.DEG
    num.ccan.one.chrom.TSS.DAP.DEG
    num.ccan.all.chrom.TSS.DAP.DEG
    # print(paste0(all.cell.types.list[i], ", chrom: ", j, " ;all.ccan.one.chrom.TSS.DAP.DEG: ", all.ccan.one.chrom.TSS.DAP.DEG, " ;num.ccan.one.chrom.TSS.DAP.DEG: ", num.ccan.one.chrom.TSS.DAP.DEG, ";num.ccan.all.chrom.TSS.DAP.DEG: ", num.ccan.all.chrom.TSS.DAP.DEG))
  }#this is the chromo loop
  
  cell.type.all.ccan.peaks <- do.call(rbind, cell.type.ccan.results)#262
  cell.type.all.ccan.peaks <- as.data.frame(unique(cell.type.all.ccan.peaks[order(cell.type.all.ccan.peaks$peaks),])) #unique(pwAll[order(pwAll$peaks),]) is a character
  names(cell.type.all.ccan.peaks) <- "peaks"
  
  #
  saveRDS(cell.type.ccan.results, file= paste0("CCAN.all.files.all.celltype.batch/CCAN.all.txt.file.DAP.and.peaks.in.same.ccan/",celltype, "/cell.type.ccan.results.DEG.TSSDAP.list.rds"))
  saveRDS(cell.type.ccan.results.score, file= paste0("CCAN.all.files.all.celltype.batch/CCAN.all.txt.file.DAP.and.peaks.in.same.ccan/", celltype, "/cell.type.ccan.results.score.DEG.TSSDAP.list.rds"))
  saveRDS(cell.type.all.ccan.peaks, file= paste0("CCAN.all.files.all.celltype.batch/CCAN.all.txt.file.DAP.and.peaks.in.same.ccan/",celltype, "/", celltype, ".cell.type.all.ccan.peaks.DEG.TSSDAP.rds"))
  write.table(cell.type.all.ccan.peaks, file = paste0("CCAN.all.files.all.celltype.batch/CCAN.all.txt.file.DAP.and.peaks.in.same.ccan/",celltype,"/", celltype, ".cell.type.all.ccan.peaks.DEG.TSSDAP.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE)
  
  print(paste0(all.cell.types.list[i], " :number of ccans near promoter TSS DAPs with DEG on all chrom is: ", num.ccan.all.chrom.TSS.DAP.DEG, " num of DEG with promoter TSS is: ", length(common.genes.in.DAP.DE), ", num of unique peaks in all DEG TSSDAP CCAN is: ", dim(cell.type.all.ccan.peaks)[1]))
} #this is the cell type loop

# write.table(same.ccan.correlated.with.one.peak$peaks,file = paste0("CCAN.all.files.all.celltype.batch/CCAN.all.txt.file.DAP.and.peaks.in.same.ccan/", celltype, "/gene.name.known/peak.", each.DAP.onlyin.TSS[k,], "in.ccan",ccan.num,".gene.",closest_genes_each.DAP.in.TSS$gene_name[k],".k.",k,".", dim(same.ccan.correlated.with.one.peak)[1], "peaks.corralated.txt"), row.names = FALSE, sep="\t", quote = FALSE)#coaccess > 0.05
links <- ConnectionsToLinks(conns = cicero_cons.new.correlated.with.one.peak, ccans = ccan) # it does not matter to use either ccan22.correlated.TSS.PAX5 or ccan22 here
Links(sub_seurat) <- links


