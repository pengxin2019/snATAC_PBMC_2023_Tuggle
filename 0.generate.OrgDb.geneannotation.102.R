
########## 2. generate  org.Sscrofa.eg.db for makeTxDbFromBiomart which only works on R 4.1.1 ########## 

library(AnnotationForge)

# gene_information <- read.table("all.genes.unique.stable.geneID.transcriptID.gene.name.mart_export.txt",sep="\t",header = FALSE,na.strings = c("NA", "NULL", ""))# version 105
gene_information <- read.table("all.genes.unique.stable.geneID.transcriptID.gene.name.mart_export.v.102.txt",sep="\t",header = FALSE,na.strings = c("NA", "NULL", ""))# version 102

names(gene_information)
length(unique(gene_information$V3))#15247, 

subset(gene_information, transcript_id == 'ENSSSCT00000059450', )
subset(gene_information, transcript_id == 'ENSSSCT00000038057', )
subset(gene_information, transcript_id == 'ENSSSCT00000088364', )
subset(gene_information, transcript_id == 'ENSSSCT00000066563', )
subset(gene_information, transcript_id == 'ENSSSCT00000084623', )

gene_information <- gene_information[-1,]

gene_information <- subset(gene_information, !is.na(gene_information$V3))

fSym <- unique(gene_information[, c(1,3)])
colnames(fSym) <- c("GID", "SYMBOL")

ensembl_trans <- unique(gene_information[, c(1:2)])
colnames(ensembl_trans) <- c("GID", "ENSEMBLTRANS")

ensembl <- unique(gene_information[, c(1,1)])
colnames(ensembl) <- c("GID", "ENSEMBL")


#tmpdir <- tempdir()
tmpdir <- "test.102"
if (!dir.exists(tmpdir))
{
  dir.create(tmpdir)
}

makeOrgPackage(gene_info = fSym, 
               ensembl_trans = ensembl_trans,
               ensembl = ensembl,
               version = "0.1",
               maintainer = "Some One <so@someplace.org>",
               author = "Some One <so@someplace.org>",
               outputDir = tmpdir,
               tax_id= "9823",
               genus= "Sus",
               species= "scrofa",
               goTable=NULL)
# Creating package in test/org.Sscrofa.eg.db 

list.files(path = "./test.102")
# [1] "org.Sscrofa.eg.db"
install.packages(file.path(tmpdir, "org.Sscrofa.eg.db"), 
                 type = "source", repos=NULL)

library("org.Sscrofa.eg.db")


###### 3. generate geneAnnotation whose seqnames contains "chr" which matches well with  genomeAnnotation BSgenome.Sscrofa.UCSC.susScr11
library(GenomicFeatures)
TxDb <- makeTxDbFromGRanges(gtf)#102

sscrofaEnsembl <- makeTxDbFromBiomart(biomart="ENSEMBL_MART_ENSEMBL",dataset="sscrofa_gene_ensembl")#this script only works in R4.1.1, does not work in R 4.0.3

geneAnnotation <- createGeneAnnotation(TxDb = sscrofaEnsembl, OrgDb = org.Sscrofa.eg.db) # this geneAnnotation does not have "chr" which makes it unmatch with all chromosomes in  BSgenome.Sscrofa.UCSC.susScr11

  
#add "chr" to the geneAnnotation object
genes <- geneAnnotation$genes
genes.df <- as.data.frame(genes)
genes.df$seqnames <- paste0("chr",genes.df$seqnames)
genes.gr <- makeGRangesFromDataFrame(genes.df, keep.extra.columns=TRUE)

#exons 
exons <- geneAnnotation$exons
exons.df <- as.data.frame(exons)
exons.df$seqnames <- paste0("chr",exons.df$seqnames)
exons.gr <- makeGRangesFromDataFrame(exons.df, keep.extra.columns=TRUE)

#TSS
TSS <- geneAnnotation$TSS
TSS.df <- as.data.frame(TSS)
TSS.df$seqnames <- paste0("chr",TSS.df$seqnames)
TSS.gr <- makeGRangesFromDataFrame(TSS.df, keep.extra.columns=TRUE)

geneAnnotation <- createGeneAnnotation(
  TSS = TSS.gr, 
  exons = exons.gr, 
  genes = genes.gr
)
