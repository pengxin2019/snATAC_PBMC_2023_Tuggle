
# on a computation node
library(ArchR) # sessionInfo(), ArchR_1.0.1   
library(BSgenome)

#https://www.archrproject.com/bookdown/getting-set-up.html 15.1:
library(BSgenome.Sscrofa.UCSC.susScr11)
#creating an annotation object
genomeAnnotation <- createGenomeAnnotation(genome = BSgenome.Sscrofa.UCSC.susScr11)

#gene annotation
library(TxDb.Sscrofa.UCSC.susScr11.refGene)
library(org.Ss.eg.db)

geneAnnotation <- createGeneAnnotation(TxDb = TxDb.Sscrofa.UCSC.susScr11.refGene, OrgDb = org.Ss.eg.db)

addArchRThreads(threads = 16) 

frag.6798.1x <- 'fragment.files.4.PBMC.datasets/6798.1x.fragments.tsv.gz'
frag.6800.1x <- 'fragment.files.4.PBMC.datasets/6800.1x.fragments.tsv.gz'
frag.6798.2x <- 'fragment.files.4.PBMC.datasets/6798.2x.fragments.tsv.gz'
frag.6800.2x <- 'fragment.files.4.PBMC.datasets/6800.2x.fragments.tsv.gz'

inputFiles <- c(frag.6798.1x, frag.6800.1x, frag.6798.2x, frag.6800.2x)

names(inputFiles) <- c("pbmc.6798.1x", "pbmc.6800.1x", "pbmc.6798.2x", "pbmc.6800.2x")

#Creation of Arrow files will create a folder in the current working directory called “QualityControl”
#which will contain 2 plots associated with each of your samples. 
ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  sampleNames = names(inputFiles),
  minTSS = 4, #Dont set this too high because you can always increase later
  minFrags = 1000, 
  addTileMat = TRUE,
  genomeAnnotation = genomeAnnotation,
  geneAnnotation = geneAnnotation
)

saveRDS(ArrowFiles,file = "ArrowFiles.4.datasets.rds")
ArrowFiles <- readRDS("ArrowFiles.4.datasets.rds")

#In ArchR, doublet removal is performed in a single step using addDoubletScores().
#This adds the infered doublet scores to each Arrow file and will take approximately 2-5 minutes per sample of the tutorial data.
#doublet scores are stored in cellColData.

doubScores <- addDoubletScores(
  input = ArrowFiles,
  k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
  knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search.
  LSIMethod = 1
)

#Adding doublet scores will create plots in the “QualityControl” directory.

saveRDS(doubScores,file = "doubScores.4.datasets.rds")
doubScores <- readRDS("doubScores.4.datasets.rds")

proj <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = "proj.outs",
  geneAnnotation = geneAnnotation,
  genomeAnnotation = genomeAnnotation,
  copyArrows = TRUE #This is recommened so that you maintain an unaltered copy for later usage.
)

#check how much memory is used to store the ArchRProject in memory within R
paste0("Memory Size = ", round(object.size(proj) / 10^6, 3), " MB")

saveRDS(proj,file = "proj.4.datasets.rds")
proj <- readRDS("proj.4.datasets.rds")

#We can subset the project to keep all cells corresponding to a specific sample:
idxSample.6798.1x <- BiocGenerics::which(proj$Sample %in% "pbmc.6798.1x")
cellsSample.6798.1x <- proj$cellNames[idxSample.6798.1x]
proj[cellsSample.6798.1x, ]

idxSample.6800.1x <- BiocGenerics::which(proj$Sample %in% "pbmc.6800.1x")
cellsSample.6800.1x <- proj$cellNames[idxSample.6800.1x]
proj[cellsSample.6800.1x, ]


idxSample.6798.2x <- BiocGenerics::which(proj$Sample %in% "pbmc.6798.2x")
cellsSample.6798.2x <- proj$cellNames[idxSample.6798.2x]
proj[cellsSample.6798.2x, ]

idxSample.6800.2x <- BiocGenerics::which(proj$Sample %in% "pbmc.6800.2x")
cellsSample.6800.2x <- proj$cellNames[idxSample.6800.2x]


#remove these predicted doublets using filterDoublets()
proj.doublet.filtered <- filterDoublets(ArchRProj = proj)


saveRDS(proj.doublet.filtered,file = "proj.4.datasets.before.LSI.rds")
proj.doublet.filtered <- readRDS("proj.4.datasets.before.LSI.rds")

df.filtered <- getCellColData(proj.doublet.filtered, select = "DoubletScore")


#6798 1x
idxSample.6798.1x.after.Dboulet.removal <- BiocGenerics::which(proj.doublet.filtered$Sample %in% "pbmc.6798.1x")
cellsSample.6798.1x.after.Dboulet.removal <- proj.doublet.filtered$cellNames[idxSample.6798.1x.after.Dboulet.removal]

cellsSample.6798.1x.after.Dboulet.removal <- as.data.frame(cellsSample.6798.1x.after.Dboulet.removal)

# install.packages("stringr")
library(stringr)

cellsSample.6798.1x.after.Dboulet.removal$cellsSample.6798.1x.after.Dboulet.removal <- str_replace(cellsSample.6798.1x.after.Dboulet.removal$cellsSample.6798.1x.after.Dboulet.removal,"pbmc.6798.1x#","") # replace the first matched pattern in each string  #/W for non word characters


write.table(cellsSample.6798.1x.after.Dboulet.removal,file = "ArchR/cellsSample.6798.1x.after.Dboulet.removal.txt", row.names = FALSE, sep="\t", quote = FALSE)

##

#6800 1x
idxSample.6800.1x.after.Dboulet.removal <- BiocGenerics::which(proj.doublet.filtered$Sample %in% "pbmc.6800.1x")
cellsSample.6800.1x.after.Dboulet.removal <- proj.doublet.filtered$cellNames[idxSample.6800.1x.after.Dboulet.removal]

cellsSample.6800.1x.after.Dboulet.removal <- as.data.frame(cellsSample.6800.1x.after.Dboulet.removal)

cellsSample.6800.1x.after.Dboulet.removal$cellsSample.6800.1x.after.Dboulet.removal <- str_replace(cellsSample.6800.1x.after.Dboulet.removal$cellsSample.6800.1x.after.Dboulet.removal,"pbmc.6800.1x#","") # replace the first matched pattern in each string  #/W for non word characters

write.table(cellsSample.6800.1x.after.Dboulet.removal, file = "ArchR/cellsSample.6800.1x.after.Dboulet.removal.txt", row.names = FALSE, sep="\t", quote = FALSE)

#6798 2x
idxSample.6798.2x.after.Dboulet.removal <- BiocGenerics::which(proj.doublet.filtered$Sample %in% "pbmc.6798.2x")
cellsSample.6798.2x.after.Dboulet.removal <- proj.doublet.filtered$cellNames[idxSample.6798.2x.after.Dboulet.removal]

cellsSample.6798.2x.after.Dboulet.removal <- as.data.frame(cellsSample.6798.2x.after.Dboulet.removal)
cellsSample.6798.2x.after.Dboulet.removal$cellsSample.6798.2x.after.Dboulet.removal <- str_replace(cellsSample.6798.2x.after.Dboulet.removal$cellsSample.6798.2x.after.Dboulet.removal,"pbmc.6798.2x#","") # replace the first matched pattern in each string  #/W for non word characters

write.table(cellsSample.6798.2x.after.Dboulet.removal, file = "ArchR/cellsSample.6798.2x.after.Dboulet.removal.txt", row.names = FALSE, sep="\t", quote = FALSE)

#6800 2x
idxSample.6800.2x.after.Dboulet.removal <- BiocGenerics::which(proj.doublet.filtered$Sample %in% "pbmc.6800.2x")
cellsSample.6800.2x.after.Dboulet.removal <- proj.doublet.filtered$cellNames[idxSample.6800.2x.after.Dboulet.removal]

cellsSample.6800.2x.after.Dboulet.removal <- as.data.frame(cellsSample.6800.2x.after.Dboulet.removal)
cellsSample.6800.2x.after.Dboulet.removal$cellsSample.6800.2x.after.Dboulet.removal <- str_replace(cellsSample.6800.2x.after.Dboulet.removal$cellsSample.6800.2x.after.Dboulet.removal,"pbmc.6800.2x#","") # replace the first matched pattern in each string  #/W for non word characters


write.table(cellsSample.6800.2x.after.Dboulet.removal, file = "ArchR/cellsSample.6800.2x.after.Dboulet.removal.txt", row.names = FALSE, sep="\t", quote = FALSE)


