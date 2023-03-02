
setwd("/work/abg/pyang19/opt/pig.6798.6800.PBMC.Satija.pipeline.result/pbmc.1x.2x.cellrangeratac1.2.0.wd/strict.filter.new.annotation/strict.filter.new.annotation.common.peaks/commonpeaks.110444.new.directory")
library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(ggplot2)
library(patchwork)
library(dplyr)#top_n
library(rtracklayer)
library(GenomicRanges)
set.seed(1234)

#3. scATAC datasets integration and clustering

common.merged.peaks.6798.6800 <- readRDS("/work/abg/pyang19/opt/pig.6798.6800.PBMC.Satija.pipeline.result/pbmc.1x.2x.cellrangeratac1.2.0.wd/common.merged.peaks.6798.6800.filtered.rds") #110444

gtf <- rtracklayer::import('/work/abg/pyang19/opt/pig.6798.6800.PBMC.Satija.pipeline.result/pbmc.1x.2x.cellrangeratac1.2.0.wd/Sus_scrofa.Sscrofa11.1.102.gtf')
#
fragpath.6798.1x <- "/work/abg/pyang19/opt/cellranger.count.input.and.result/PBMC6798_with_motif_count_result_new/outs/fragments.tsv.gz"
fragpath.6800.1x <- "/work/abg/pyang19/opt/cellranger.count.input.and.result/PBMC6800_with_motif_count_result_new/outs/fragments.tsv.gz"
fragpath.6798.2x <- "/work/abg/pyang19/opt/pig.6798.6800.PBMC.Satija.pipeline.result/pbmc.6798.2x.cellranger.count.result.outs/outs/fragments.tsv.gz"
fragpath.6800.2x <- "/work/abg/pyang19/opt/pig.6798.6800.PBMC.Satija.pipeline.result/pbmc.6800.2x.cellranger.count.result.outs/outs/fragments.tsv.gz"
fragcounts.6798.1x <- CountFragments(fragments = fragpath.6798.1x)
fragcounts.6800.1x <- CountFragments(fragments = fragpath.6800.1x)
fragcounts.6798.2x <- CountFragments(fragments = fragpath.6798.2x)
fragcounts.6800.2x <- CountFragments(fragments = fragpath.6800.2x)

atac.cells.6798.1x <- fragcounts.6798.1x[fragcounts.6798.1x$frequency_count > 2000, "CB"]
atac.cells.6800.1x <- fragcounts.6800.1x[fragcounts.6800.1x$frequency_count > 2000, "CB"]
atac.cells.6798.2x <- fragcounts.6798.2x[fragcounts.6798.2x$frequency_count > 2000, "CB"]
atac.cells.6800.2x <- fragcounts.6800.2x[fragcounts.6800.2x$frequency_count > 2000, "CB"]

# create the fragment object
atac.frags.6798.1x <- CreateFragmentObject(path = fragpath.6798.1x, cells = atac.cells.6798.1x)
atac.frags.6800.1x <- CreateFragmentObject(path = fragpath.6800.1x, cells = atac.cells.6800.1x)
atac.frags.6798.2x <- CreateFragmentObject(path = fragpath.6798.2x, cells = atac.cells.6798.2x)
atac.frags.6800.2x <- CreateFragmentObject(path = fragpath.6800.2x, cells = atac.cells.6800.2x)

# # quantify combined.peaks.6798 peaks in the each scATAC-seq dataset
# # argument features : A GRanges object containing a set of genomic intervals.
# # These will form the rows of the matrix, with each entry recording the number of unique reads falling in the genomic region for each cell.
counts.6798.1x.new <- FeatureMatrix(
  fragments = atac.frags.6798.1x,
  features = granges(common.merged.peaks.6798.6800),  #A GRanges object containing a set of genomic intervals. These will form the rows of the matrix, with each entry recording the number of unique reads falling in the genomic region for each cell.
  cells = atac.cells.6798.1x
)
saveRDS(counts.6798.1x.new,file = "counts.6798.1x.new.pbmc.atac.6798.1x.rds")# R 4.0.5
counts.6798.1x.new <- readRDS("counts.6798.1x.new.pbmc.atac.6798.1x.rds")# R 4.0.5

str(counts.6798.1x.new)#'dgCMatrix'
dim(counts.6798.1x.new)#110444   1875


counts.6800.1x.new <- FeatureMatrix(
  fragments = atac.frags.6800.1x,
  features = granges(common.merged.peaks.6798.6800),  #A GRanges object containing a set of genomic intervals. These will form the rows of the matrix, with each entry recording the number of unique reads falling in the genomic region for each cell.
  cells = atac.cells.6800.1x
)
saveRDS(counts.6800.1x.new,file = "counts.6800.1x.new.pbmc.atac.6800.1x.rds")# R 4.0.5
counts.6800.1x.new <- readRDS("counts.6800.1x.new.pbmc.atac.6800.1x.rds")# R 4.0.5

counts.6798.2x.new <- FeatureMatrix(
  fragments = atac.frags.6798.2x,
  features = granges(common.merged.peaks.6798.6800),  #A GRanges object containing a set of genomic intervals. These will form the rows of the matrix, with each entry recording the number of unique reads falling in the genomic region for each cell.
  cells = atac.cells.6798.2x
)

saveRDS(counts.6798.2x.new,file = "counts.6798.2x.new.pbmc.atac.6798.2x.rds")# R 4.0.5
counts.6798.2x.new <- readRDS("counts.6798.2x.new.pbmc.atac.6798.2x.rds")# R 4.0.5

counts.6800.2x.new <- FeatureMatrix(
  fragments = atac.frags.6800.2x,
  features = granges(common.merged.peaks.6798.6800),  #A GRanges object containing a set of genomic intervals. These will form the rows of the matrix, with each entry recording the number of unique reads falling in the genomic region for each cell.
  cells = atac.cells.6800.2x
)
saveRDS(counts.6800.2x.new,file = "counts.6800.2x.new.pbmc.atac.6800.2x.rds")# R 4.0.5
counts.6800.2x.new <- readRDS("counts.6800.2x.new.pbmc.atac.6800.2x.rds")# R 4.0.5

# # create new object for each dataset
atac.assay.6798.1x <- CreateChromatinAssay(
  counts = counts.6798.1x.new,
  min.features = 1000, #Include cells where at least this many features are detected.
  fragments = atac.frags.6798.1x
)

atac.assay.6800.1x <- CreateChromatinAssay(
  counts = counts.6800.1x.new,
  min.features = 1000, #Include cells where at least this many features are detected.
  fragments = atac.frags.6800.1x
)

atac.assay.6798.2x <- CreateChromatinAssay(
  counts = counts.6798.2x.new,
  min.features = 1000, #Include cells where at least this many features are detected.
  fragments = atac.frags.6798.2x
)

atac.assay.6800.2x <- CreateChromatinAssay(
  counts = counts.6800.2x.new,
  min.features = 1000, #Include cells where at least this many features are detected.
  fragments = atac.frags.6800.2x
)

pbmc.atac.6798.1x <- CreateSeuratObject(counts = atac.assay.6798.1x, assay = "peaks")
pbmc.atac.6800.1x <- CreateSeuratObject(counts = atac.assay.6800.1x, assay = "peaks")
pbmc.atac.6798.2x <- CreateSeuratObject(counts = atac.assay.6798.2x, assay = "peaks")
pbmc.atac.6800.2x <- CreateSeuratObject(counts = atac.assay.6800.2x, assay = "peaks")

#
Annotation(pbmc.atac.6798.1x) <- gtf
Annotation(pbmc.atac.6800.1x) <- gtf
Annotation(pbmc.atac.6798.2x) <- gtf
Annotation(pbmc.atac.6800.2x) <- gtf

saveRDS(pbmc.atac.6798.1x,file = "annotated.before.new.doublet.filtered.pbmc.atac.6798.1x.rds")# R 4.0.5
saveRDS(pbmc.atac.6800.1x, file = "annotated.before.new.doublet.filtered.pbmc.atac.6800.1x.rds")# R 4.0.5
saveRDS(pbmc.atac.6798.2x, file = "annotated.before.new.doublet.filtered.pbmc.atac.6798.2x.rds")# R 4.0.5
saveRDS(pbmc.atac.6800.2x, file = "annotated.before.new.doublet.filtered.pbmc.atac.6800.2x.rds")# R 4.0.5

pbmc.atac.6798.1x <- readRDS("annotated.before.new.doublet.filtered.pbmc.atac.6798.1x.rds")# R 4.0.5
pbmc.atac.6800.1x <- readRDS("annotated.before.new.doublet.filtered.pbmc.atac.6800.1x.rds")# R 4.0.5
pbmc.atac.6798.2x <- readRDS("annotated.before.new.doublet.filtered.pbmc.atac.6798.2x.rds")# R 4.0.5
pbmc.atac.6800.2x <- readRDS("annotated.before.new.doublet.filtered.pbmc.atac.6800.2x.rds")# R 4.0.5

##import the cells are predicted as non-doublet for each dataset processed by gg, only keep nuclei that are not doublet

cells.use.6798.1x <- read.table("/work/abg/pyang19/opt/pig.6798.6800.PBMC.Satija.pipeline.result/pbmc.1x.2x.cellrangeratac1.2.0.wd/ArchR/cellsSample.6798.1x.after.Dboulet.removal.txt", header = FALSE, row.names = NULL)
cells.use.6800.1x <- read.table("/work/abg/pyang19/opt/pig.6798.6800.PBMC.Satija.pipeline.result/pbmc.1x.2x.cellrangeratac1.2.0.wd/ArchR/cellsSample.6800.1x.after.Dboulet.removal.txt", header = FALSE, row.names = NULL)
cells.use.6798.2x <- read.table("/work/abg/pyang19/opt/pig.6798.6800.PBMC.Satija.pipeline.result/pbmc.1x.2x.cellrangeratac1.2.0.wd/ArchR/cellsSample.6798.2x.after.Dboulet.removal.txt", header = FALSE, row.names = NULL)
cells.use.6800.2x <- read.table("/work/abg/pyang19/opt/pig.6798.6800.PBMC.Satija.pipeline.result/pbmc.1x.2x.cellrangeratac1.2.0.wd/ArchR/cellsSample.6800.2x.after.Dboulet.removal.txt", header = FALSE, row.names = NULL)

cells.use.6798.1x <- cells.use.6798.1x[-1,] # subset the first row which is a col name
cells.use.6800.1x <- cells.use.6800.1x[-1,] # subset the first row which is a col name
cells.use.6798.2x <- cells.use.6798.2x[-1,] # subset the first row which is a col name
cells.use.6800.2x <- cells.use.6800.2x[-1,] # subset the first row which is a col name

pbmc.atac.6798.1x <- subset(pbmc.atac.6798.1x, cells = cells.use.6798.1x)
pbmc.atac.6800.1x <- subset(pbmc.atac.6800.1x, cells = cells.use.6800.1x)
pbmc.atac.6798.2x <- subset(pbmc.atac.6798.2x, cells = cells.use.6798.2x)
pbmc.atac.6800.2x <- subset(pbmc.atac.6800.2x, cells = cells.use.6800.2x)


saveRDS(pbmc.atac.6798.1x,file = "annotated.doublet.filtered.pbmc.atac.6798.1x.rds")
saveRDS(pbmc.atac.6800.1x, file = "annotated.doublet.filtered.pbmc.atac.6800.1x.rds")
saveRDS(pbmc.atac.6798.2x, file = "annotated.doublet.filtered.pbmc.atac.6798.2x.rds")
saveRDS(pbmc.atac.6800.2x, file = "annotated.doublet.filtered.pbmc.atac.6800.2x.rds")

pbmc.atac.6798.1x <- readRDS("annotated.doublet.filtered.pbmc.atac.6798.1x.rds")
pbmc.atac.6800.1x <- readRDS("annotated.doublet.filtered.pbmc.atac.6800.1x.rds")
pbmc.atac.6798.2x <- readRDS("annotated.doublet.filtered.pbmc.atac.6798.2x.rds")
pbmc.atac.6800.2x <- readRDS("annotated.doublet.filtered.pbmc.atac.6800.2x.rds")

plots.new.doublet.removal <- "plot.manuscript.final"
if (!dir.exists(plots.new.doublet.removal))
{
  dir.create(plots.new.doublet.removal)
}


p <- CoveragePlot(
  object = pbmc.atac.6798.1x,
  region = "1-21403-24616"
) + ggtitle('6798.1x:1.21403-24616')+ NoLegend() + FontSize(x.title = 20, y.title = 20)

p
ggsave("6798.1x.21403.24616.pdf", path = "./plot.manuscript.final")


p <- CoveragePlot(
  object = pbmc.atac.6800.1x,
  region = "1-21403-24616"
) + ggtitle('6800.1x:1-21403-24616')+ NoLegend() + FontSize(x.title = 20, y.title = 20)
p
ggsave("6800.1x.1.21403.24616.pdf", path = "./plot.manuscript.final")


p <- CoveragePlot(
  object = pbmc.atac.6798.2x,
  region = "1-21403-24616"
) + ggtitle('6798.2x:1-21403-24616')+ NoLegend() + FontSize(x.title = 20, y.title = 20)
p
ggsave("6798.2x.21403.24616.pdf", path = "./plot.manuscript.final")

p <- CoveragePlot(
  object = pbmc.atac.6800.2x,
  region = "1-21403-24616"
) + ggtitle('6800.2x:1-21403-24616')+ NoLegend() + FontSize(x.title = 20, y.title = 20)
p
ggsave("6800.2x.21403-24616.pdf", path = "./plot.manuscript.final")

#Computing QC Metrics
# compute nucleosome signal score per cell and compute TSS enrichment score per cell

pbmc.atac.6798.1x <- NucleosomeSignal(object = pbmc.atac.6798.1x)%>% TSSEnrichment(fast = FALSE) 
pbmc.atac.6800.1x <- NucleosomeSignal(object = pbmc.atac.6800.1x)%>% TSSEnrichment(fast = FALSE) 
pbmc.atac.6798.2x <- NucleosomeSignal(object = pbmc.atac.6798.2x)%>% TSSEnrichment(fast = FALSE) 
pbmc.atac.6800.2x <- NucleosomeSignal(object = pbmc.atac.6800.2x)%>% TSSEnrichment(fast = FALSE) 


pbmc.atac.6798.1x <- subset(
  x = pbmc.atac.6798.1x,
  subset = nucleosome_signal < 4 &
    TSS.enrichment > 2 &
    nCount_peaks > 2000 &
    nCount_peaks < 30000
)
pbmc.atac.6800.1x <- subset(
  x = pbmc.atac.6800.1x,
  subset = nucleosome_signal < 4 &
    TSS.enrichment > 2 &
    nCount_peaks > 2000 &
    nCount_peaks < 30000
)

pbmc.atac.6798.2x <- subset(
  x = pbmc.atac.6798.2x,
  subset = nucleosome_signal < 4 &
    TSS.enrichment > 2 &
    nCount_peaks > 2000 &
    nCount_peaks < 30000
)

pbmc.atac.6800.2x <- subset(
  x = pbmc.atac.6800.2x,
  subset = nucleosome_signal < 4 &
    TSS.enrichment > 2 &
    nCount_peaks > 2000 &
    nCount_peaks < 30000
)


# compute LSI
pbmc.atac.6798.1x <- RunTFIDF(pbmc.atac.6798.1x, min.cutoff = 'q0') %>% FindTopFeatures(min.cutoff = 10) %>% RunSVD()
pbmc.atac.6800.1x <- RunTFIDF(pbmc.atac.6800.1x, min.cutoff = 'q0') %>% FindTopFeatures(min.cutoff = 10) %>% RunSVD()
pbmc.atac.6798.2x <- RunTFIDF(pbmc.atac.6798.2x, min.cutoff = 'q0') %>% FindTopFeatures(min.cutoff = 10) %>% RunSVD()
pbmc.atac.6800.2x <- RunTFIDF(pbmc.atac.6800.2x, min.cutoff = 'q0') %>% FindTopFeatures(min.cutoff = 10) %>% RunSVD()


# first add dataset-identifying metadata
pbmc.atac.6798.1x$dataset <- "6798.1x"
pbmc.atac.6800.1x$dataset <- "6800.1x"
pbmc.atac.6798.2x$dataset <- "6798.2x"
pbmc.atac.6800.2x$dataset <- "6800.2x"

saveRDS(pbmc.atac.6798.1x,file = "pbmc.atac.6798.1x.rds")
saveRDS(pbmc.atac.6800.1x, file = "pbmc.atac.6800.1x.rds")
saveRDS(pbmc.atac.6798.2x, file = "pbmc.atac.6798.2x.rds")
saveRDS(pbmc.atac.6800.2x, file = "pbmc.atac.6800.2x.rds")

pbmc.atac.6798.1x <- readRDS("pbmc.atac.6798.1x.rds")
pbmc.atac.6800.1x <- readRDS("pbmc.atac.6800.1x.rds")
pbmc.atac.6798.2x <- readRDS("pbmc.atac.6798.2x.rds")
pbmc.atac.6800.2x <- readRDS("pbmc.atac.6800.2x.rds")

# merge all four datasets, 
pbmc.combined <- merge(pbmc.atac.6798.1x, y = list(pbmc.atac.6800.1x,pbmc.atac.6798.2x, pbmc.atac.6800.2x)) %>% RunTFIDF(min.cutoff = 'q0') %>% FindTopFeatures(min.cutoff = 10)%>% RunSVD()

saveRDS(pbmc.combined,file = "pbmc.combined.four.rds")
pbmc.combined <- readRDS("pbmc.combined.four.rds")

pbmc.combined
# An object of class Seurat 
# 110444 features across 17207 samples within 1 assay 
# Active assay: peaks (110444 features, 110272 variable features)
# 1 dimensional reduction calculated: lsi

pbmc.combined <- RunUMAP(pbmc.combined, reduction = "lsi", dims = 2:30)
p <- DimPlot(pbmc.combined, group.by = "dataset", label = FALSE)
p
ggsave("pbmc.combined.by.dataset.2.30.pdf", path = "./plot.manuscript.final")

pbmc.combined <- RunUMAP(pbmc.combined, reduction = "lsi")

#Integration
# find integration anchors
integration.anchors <- FindIntegrationAnchors(
  object.list = list(pbmc.atac.6798.1x, pbmc.atac.6800.1x,pbmc.atac.6798.2x,pbmc.atac.6800.2x),
  anchor.features = granges(common.merged.peaks.6798.6800),
  reduction = "rlsi",
  dims = 2:30
)

# integrate LSI embeddings
integrated <- IntegrateEmbeddings(
  anchorset = integration.anchors,
  reductions = pbmc.combined[["lsi"]],
  new.reduction.name = "integrated_lsi",
  dims.to.integrate = 1:30 #1:30 of pbmc.combined[["lsi"]]
)
#  No standard deviation info stored for integrated_lsi

integrated <- RunUMAP(integrated, reduction = "integrated_lsi", dims = 2:30)

saveRDS(integrated, file = "integrated.four.rds")
integrated <- readRDS("integrated.four.rds")

# An object of class Seurat 
# 110444 features across 17207 samples within 1 assay 
# Active assay: peaks (110444 features, 0 variable features)
# 2 dimensional reductions calculated: integrated_lsi, umap

#check whether integration is good enough
Idents(integrated) <- integrated$dataset
integrated.6798.1x <- subset(x = integrated, idents = "6798.1x")
integrated.6800.1x <- subset(x = integrated, idents = "6800.1x")
integrated.6798.2x <- subset(x = integrated, idents = "6798.2x")
integrated.6800.2x <- subset(x = integrated, idents = "6800.2x")

p <- DimPlot(object = integrated.6798.1x, label = FALSE) #+ ggtitle('6798.1x')
p
ggsave("integrated.6798.1x.2.30.pdf", path = "./plot.manuscript.final")

p <- DimPlot(object = integrated.6800.1x, label = FALSE) #+ ggtitle('6800.1x')
p
ggsave("integrated.6800.1x.2.30.pdf", path = "./plot.manuscript.final")

p <- DimPlot(object = integrated.6798.2x, label = FALSE) #+ ggtitle('6798.2x')
p
ggsave("integrated.6798.2x.2.30.pdf", path = "./plot.manuscript.final")

p <- DimPlot(object = integrated.6800.2x, label = FALSE)# + ggtitle('6800.2x')
p
ggsave("integrated.6800.2x.2.30.pdf", path = "./plot.manuscript.final")


###clusters
#Normalization and linear dimensional reduction
integrated <- RunTFIDF(integrated, min.cutoff = 'q0') %>% FindTopFeatures(min.cutoff = 10) %>% RunSVD()
integrated <- RunUMAP(integrated, reduction = "integrated_lsi", dims = 2:30)##integrated_lsi is the reduciton which has removed batch effect

p <- DimPlot(integrated, group.by = "dataset", label = FALSE)
# p2+ ggtitle("integrated.2.30.integrated_lsi")
p
ggsave("integrated.four2.30.pdf", path = "./plot.manuscript.final")

integrated.ElbowPlot <- ElbowPlot(integrated, ndims = 50, reduction = "lsi")
integrated.ElbowPlot
ggsave("integrated.ElbowPlot.lsi.after.RunTFIDF.pdf", path = "./plot.manuscript.final")

p3 <- DepthCor(integrated)
p3
ggsave("integrated.no.doublet.cor.2.30.pdf", path = "./plot.manuscript.final")

integrated <- FindNeighbors(object = integrated, reduction = 'integrated_lsi', dims = 2:30)
integrated <- FindClusters(object = integrated, verbose = TRUE, algorithm = 3, resolution = 2.4) #resolution is the best resolution based on the optomization.R result, default resolution is 0.8

unique(integrated$seurat_clusters)#35 Levels: 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 ... 34

p <- DimPlot(object = integrated ,group.by = "peaks_snn_res.2.4", label = TRUE) + NoLegend()# + ggtitle("integrated.resolution2.4.2:30")
p
ggsave("integrated.resolution2.4.2.30.doublet.pdf", path = "./plot.manuscript.final")

saveRDS(integrated, file = "integrated.four.with.cluster.2.4.2.30.rds")#resolution 2.4,2:30
integrated <- readRDS("integrated.four.with.cluster.2.4.2.30.rds")#resolution 2.4,2:30

clusters <- unique(integrated$seurat_clusters)
sub.cluster <- list()

for (i in clusters)
{
  sub.cluster <- subset(integrated, idents = i)
  sub.cells <- Cells(sub.cluster)
  p <- DimPlot(object = integrated, label = TRUE, cells.highlight = sub.cells ) + NoLegend() + ggtitle(paste0("cluster: ", i))
  p
  ggsave(paste0("integrated.2.4.2.30.highlight.cluster.", i, ".pdf"), path = "./plot.manuscript.final")
}

Idents(integrated) <- integrated$seurat_clusters

#install.packages('ape')
library(ape)
integrated.tree <- BuildClusterTree(object = integrated)

PlotClusterTree(integrated.tree, edge.width = 3)  # plot tree with node labels

data.tree <- Tool(object = integrated, 
                  slot = "BuildClusterTree") # pull the tree
ape::plot.phylo(x = data.tree, 
                direction = "downwards", # plot the tree without node labels
                edge.width = 1.5)

SampleTotalCells <- prop.table(table(integrated$dataset)) # What percent of total cells are from each sample?
SamplePercents <- prop.table(table(Idents(integrated),integrated$dataset), margin = 1) # What percent of cells from each cluster belong to each sample?
SamplePercents <- rbind(SamplePercents, SampleTotalCells) # add row of overall percentages to table
# SamplePercents
rowSums(SamplePercents) # make sure all are equal to 1
SamplePercents <- t(SamplePercents)# transpose the table

rownames(SamplePercents)
#[1] "6798.1x" "6798.2x" "6800.1x" "6800.2x"


par(mfrow=c(1, 1), mar=c(5, 5, 4, 8)) #mfrow = c(1, 1)) # Create a 1 x 1 plotting matrix
barplot(SamplePercents, # create stacked bar plot
                     col = rev(c('dodgerblue', 'yellow','orangered', 'orange')),
                     legend = rownames(SamplePercents), #67981x or 68001x or 67982x or 68002x
                     xlab = "cluster number",
                     ylab = "Frequency within cluster",
                     sub= "integrated.2.4.2.30",
                     las = 2, #rotate x axis label
                     cex.names=0.5,
                     border = NA,
                     space = 0.05,
                     legend.text = TRUE,
                     args.legend = list(x = "topright", bty = "n", inset=c(-0.15, 0)))


saveRDS(integrated, file = "integrated.four.with.cluster.2.4.2.30.rds")#resolution 2.4,2:30
integrated <- readRDS("integrated.four.with.cluster.2.4.2.30.rds")#resolution 2.4,2:30

Idents(integrated)#35 Levels: 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 ... 34

Idents(integrated) <- integrated$orig.ident
Idents(integrated)#Levels: SeuratProject

hist(integrated$TSS.enrichment)


p <- VlnPlot(
  object = integrated,
  features = c('TSS.enrichment',  'nucleosome_signal'),
  pt.size = 0.1,
  ncol = 5
)
p
ggsave("integrated.2.4.2.30.TSS.enrichment.nucleosome_signal.pdf", path = "./plot.manuscript.final")
