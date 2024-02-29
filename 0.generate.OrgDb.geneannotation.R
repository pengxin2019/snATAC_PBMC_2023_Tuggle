
########## 2. generate  org.Sscrofa.eg.db for makeTxDbFromBiomart which only works on R 4.1.1 ########## 
library(AnnotationForge)

gene_information <- read.table("all.genes.unique.stable.geneID.transcriptID.gene.name.mart_export.txt",sep="\t",header = FALSE,na.strings = c("NA", "NULL", ""))#  

gene_information <- gene_information[-1,]

gene_information <- subset(gene_information, !is.na(gene_information$V3))

fSym <- unique(gene_information[, c(1,3)])
colnames(fSym) <- c("GID", "SYMBOL")

ensembl_trans <- unique(gene_information[, c(1:2)])
colnames(ensembl_trans) <- c("GID", "ENSEMBLTRANS")

ensembl <- unique(gene_information[, c(1,1)])
colnames(ensembl) <- c("GID", "ENSEMBL")

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

install.packages(file.path(tmpdir, "org.Sscrofa.eg.db"), 
                 type = "source", repos=NULL)

library("org.Sscrofa.eg.db")

###### 3. generate geneAnnotation whose seqnames contains "chr" which matches well with  genomeAnnotation BSgenome.Sscrofa.UCSC.susScr11
library(GenomicFeatures)
TxDb <- makeTxDbFromGRanges(gtf)#

sscrofaEnsembl <- makeTxDbFromBiomart(biomart="ENSEMBL_MART_ENSEMBL",dataset="sscrofa_gene_ensembl")#this script  works in R4.1.1, does not work in R 4.0.3
geneAnnotation <- createGeneAnnotation(TxDb = sscrofaEnsembl, OrgDb = org.Sscrofa.eg.db) # this geneAnnotation does not have "chr" which makes it unmatch with all chromosomes in  BSgenome.Sscrofa.UCSC.susScr11

# > geneAnnotation
# List of length 3
# names(3): genes exons TSS

# geneAnnotation$genes
# GRanges object with 31908 ranges and 2 metadata columns:
#   seqnames        ranges strand |            gene_id
# <Rle>     <IRanges>  <Rle> |        <character>
#   [1]              1        1-3782      + | ENSSSCG00000048769
# [2]              1    9815-22152      - | ENSSSCG00000037372
# [3]              1   23782-40033      + | ENSSSCG00000027257
# [4]              1  96128-186749      - | ENSSSCG00000029697
# [5]              1 108537-114011      - | ENSSSCG00000049216
# ...            ...           ...    ... .                ...
# [31904] AEMK02000703.1 161357-184831      - | ENSSSCG00000045011
# [31905] AEMK02000704.1    5510-38276      - | ENSSSCG00000035691
# [31906] AEMK02000704.1   61618-62544      + | ENSSSCG00000033102
# [31907] FPKY02000006.1   52921-54362      + | ENSSSCG00000049744
# [31908] FPKY02000009.1   19082-22644      - | ENSSSCG00000044056
# symbol
# <character>
#   [1] NA_ENSSSCG00000048769
# [2] NA_ENSSSCG00000037372
# [3]                 PSMB1
# [4]               FAM120B
# [5] NA_ENSSSCG00000049216
# ...                   ...
# [31904] NA_ENSSSCG00000045011
# [31905] NA_ENSSSCG00000035691
# [31906] NA_ENSSSCG00000033102
# [31907] NA_ENSSSCG00000049744
# [31908] NA_ENSSSCG00000044056
# -------
  
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
# > TSS.gr
# GRanges object with 53507 ranges and 2 metadata columns:
#   seqnames    ranges strand |     tx_id            tx_name
# <Rle> <IRanges>  <Rle> | <integer>        <character>
#   [1]              chr1         1      + |         1 ENSSSCT00000066540
# [2]              chr1       792      + |         2 ENSSSCT00000078085
# [3]              chr1     23782      + |         3 ENSSSCT00000058189
# [4]              chr1     23853      + |         4 ENSSSCT00000023500
# [5]              chr1    199231      + |         5 ENSSSCT00000062318

geneAnnotation <- createGeneAnnotation(
  TSS = TSS.gr, 
  exons = exons.gr, 
  genes = genes.gr
)
# > geneAnnotation$genes
# GRanges object with 31908 ranges and 2 metadata columns:
#   seqnames        ranges strand |            gene_id
# <Rle>     <IRanges>  <Rle> |        <character>
#   [1]              chr1        1-3782      + | ENSSSCG00000048769
# [2]              chr1    9815-22152      - | ENSSSCG00000037372
# [3]              chr1   23782-40033      + | ENSSSCG00000027257
# [4]              chr1  96128-186749      - | ENSSSCG00000029697
# [5]              chr1 108537-114011      - | ENSSSCG00000049216

# > geneAnnotation$exons
# GRanges object with 646004 ranges and 2 metadata columns:
#   seqnames      ranges strand |            gene_id
# <Rle>   <IRanges>  <Rle> |        <character>
#   [1]              chr1      1-2465      + | ENSSSCG00000048769
# [2]              chr1     792-961      + | ENSSSCG00000048769
# [3]              chr1   2371-2465      + | ENSSSCG00000048769
# [4]              chr1   3118-3780      + | ENSSSCG00000048769
# [5]              chr1   3118-3782      + | ENSSSCG00000048769

# > geneAnnotation$TSS
# GRanges object with 53507 ranges and 2 metadata columns:
#   seqnames    ranges strand |     tx_id            tx_name
# <Rle> <IRanges>  <Rle> | <integer>        <character>
#   [1]              chr1         1      + |         1 ENSSSCT00000066540
# [2]              chr1       792      + |         2 ENSSSCT00000078085
# [3]              chr1     23782      + |         3 ENSSSCT00000058189
# [4]              chr1     23853      + |         4 ENSSSCT00000023500
# [5]              chr1    199231      + |         5 ENSSSCT00000062318

saveRDS(geneAnnotation,file = "geneAnnotation.with.chr.rds") ## this geneAnnotation has "chr" which makes it match with all chromosomes in  BSgenome.Sscrofa.UCSC.susScr11 well
geneAnnotation <- readRDS("geneAnnotation.with.chr.rds") ## this geneAnnotation has "chr" which makes it match with all chromosomes in  BSgenome.Sscrofa.UCSC.susScr11 well
