#2.2 comparasion with bulk ATAC
# load the common peaks between scATAC peaks and bulk ATAC peaks created by Ryan

# 1. download it from Box folder uploaded by Ryan and upload it to Nova account

library(Signac)
library(Seurat)
library(SeuratWrappers)

load("PBMC_singleCell_bulk_ATAC_shared_peaks.Rdata")
# shared_sc_peaks
# GRanges object with 86494 ranges and 0 metadata columns:

# generate the count table in shared_sc_peaks for each cell type
integrated <- readRDS("integrated.withRNA.assay.annotated.integration.scRNA.celltype2.4.2.30.after.filter.new.idents.RNA2NRA.wo.unknown.rds")#resol 2.4, +-2000, 2.30, script 4.2
DefaultAssay(integrated) <- 'peaks'

all.cell.types.list <- list("cDCs", "Monocytes", "CD4posab","CD2negGD","B","ASC","CD2posGD","CD8abPOSab","NK","pDCs","CD8aPOSabT_NK")
counts.shared.peaks.sc.bulk.in.sc.df <- data.frame(matrix(ncol = 11, nrow = 86494))
names(counts.shared.peaks.sc.bulk.in.sc.df) <- c("cDCs", "Monocytes", "CD4posab","CD2negGD","B","ASC","CD2posGD","CD8abPOSab","NK","pDCs","CD8aPOSabT_NK")

for (i in 1:length(all.cell.types.list))
{
  celltype <- subset(x = integrated, idents = all.cell.types.list[[i]])
  
    counts.shared.peaks.sc.bulk <- FeatureMatrix(
    fragments = Fragments(celltype),
    features = granges(shared_sc_peaks),  #A GRanges object containing a set of genomic intervals. These will form the rows of the matrix, with each entry recording the number of unique reads falling in the genomic region for each cell.
    cells = Cells(integrated)
  )
  counts.shared.peaks.sc.bulk.in.sc.df[,i] <- as.data.frame(rowSums(as.data.frame(counts.shared.peaks.sc.bulk)))
  
  print(paste0(all.cell.types.list[[i]], "  is done! "))
}

save(counts.shared.peaks.sc.bulk.in.sc.df, file = "counts.shared.peaks.sc.bulk.in.sc.df.celltype.rowSum.RData")# rowSums
# load("counts.shared.peaks.sc.bulk.in.sc.df.celltype.RData")#rowSums


