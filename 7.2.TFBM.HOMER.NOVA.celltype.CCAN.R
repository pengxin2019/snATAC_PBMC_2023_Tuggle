
setwd("/work/abg/pyang19/scATAC_PBMC_analysis/TFBM_HOMER/input.output.files/pre_input_files/cell.type.ccan.peaks/gene.name.known.DEG.TSS.DAP.102/B")

all.pre.input.files <- list.files(path='/work/abg/pyang19/scATAC_PBMC_analysis/TFBM_HOMER/input.output.files/pre_input_files/cell.type.ccan.peaks/gene.name.known.DEG.TSS.DAP.102/B/',pattern='.txt')

all.cell.types.list <- c("ASC", "B","CD2negGD","CD2posGD","CD4posab","CD8abPOSab","CD8aPOSabT_NK","cDCs","Monocytes","NK","pDCs") # it is a factor

cd /work/abg/pyang19/scATAC_PBMC_analysis/TFBM_HOMER/input.output.files/input_files/celltype.ccan.peaks/Monocytes
  
for (i in 1:length(all.cell.types.list)) #1:11 cell types
{
  celltype <- all.cell.types.list[i]
  
  tmpdir <- "/work/abg/pyang19/scATAC_PBMC_analysis/TFBM_HOMER/input.output.files/input_files/celltype.ccan.peaks"
  if (!dir.exists(tmpdir))
  {
    dir.create(tmpdir)
  }
  
  tmpdir <- paste0("/work/abg/pyang19/scATAC_PBMC_analysis/TFBM_HOMER/input.output.files/input_files/celltype.ccan.peaks/",celltype)
  if (!dir.exists(tmpdir))
  {
    dir.create(tmpdir)
  }
  
  tmpdir <- "/work/abg/pyang19/scATAC_PBMC_analysis/TFBM_HOMER/input.output.files/output_files/celltype.ccan.peaks"
  if (!dir.exists(tmpdir))
  {
    dir.create(tmpdir)
  }
  
  tmpdir <- paste0("/work/abg/pyang19/scATAC_PBMC_analysis/TFBM_HOMER/input.output.files/output_files/celltype.ccan.peaks/",celltype)
  if (!dir.exists(tmpdir))
  {
    dir.create(tmpdir)
  }
}

library(stringr)
for (i in 1:length(all.pre.input.files))
{
  input <- read.table(all.pre.input.files[i])#names is "V1" now
  input <- as.data.frame(input[-1,])
  names(input) <- "V1"
  input2 <- str_split_fixed(input$V1, "-", 3)# split a vector of strings into a matrix of substrings 
  
  input2 <- as.data.frame(input2)
  input2$strand <- "+"
  input2$ID <- input$V1
  names(input2) <- c("chromosome","start","end","strand","peak_ID")
  input2 <- input2[,c(5,1,2,3,4)]
  write.table(input2, file = (paste0("/work/abg/pyang19/scATAC_PBMC_analysis/TFBM_HOMER/input.output.files/input_files/celltype.ccan.peaks/B/",substr(all.pre.input.files[i],1,nchar(all.pre.input.files[i])-3), "input.HOMER.txt")), sep = "\t",row.names = FALSE,col.names=FALSE,quote = FALSE )
}

#############  RUN HOMER ###################################################################################################################################################################################################                  
#############  RUN HOMER ###################################################################################################################################################################################################                  

#############  RUN HOMER: B cells ##############################################################################      done

cd /work/abg/pyang19/scATAC_PBMC_analysis/TFBM_HOMER/input.output.files/input_files/celltype.ccan.peaks/B
ls *input.HOMER.txt | wc -l
#33
for file in *input.HOMER.txt
do
findMotifsGenome.pl "$file" $HOMER_HOME/data/genomes/susScr11_ensembl_102/Sus_scrofa.Sscrofa11.1.dna.toplevel.fa $HOMER_HOME/input.output.files/output_files/celltype.ccan.peaks/B/TFBM.B.33.ccan.peak.300.bg.output/"${file%.input.HOMER.txt}.300.bg.output"  -size given -mask -mset vertebrates -N 300
done

