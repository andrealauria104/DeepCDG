# R script to download and prepare gene expression data using TCGAbiolinks

# 0.1 Args Parsing ====
args <- commandArgs(trailingOnly = T)

parserr <- 'Incorrect argument parsing. Please specify: \n
1 - TCGA project (Cancer type). \n
2 - Data directory. \n
3 - Format for storing prepared data (if needed: Rdata, csv or both) \n '

if(length(args)<2){
  stop(message(parserr))
}

# Define parsing arguments
cancer_type = as.character(args[1])
dir = as.character(args[2])

if(length(args)<3){
  format = 'none'
} else{
  format = as.character(args[3])
}

# Directory control
dir1 = dirname(dir)
dir2 = basename(dir)

dir = paste0(dir1, dir2)
rm(list=c("dir1","dir2"))

# 0.2 Resources ====

library(TCGAbiolinks)
library(SummarizedExperiment)
library(reshape2)

# 1.0 Download Gene Expression data ====

# 1.1 Query GDC
query.exp.hg38 <- GDCquery( project = cancer_type
                                 ,data.category = "Transcriptome Profiling"
                                 ,data.type = "Gene Expression Quantification"
                                 ,workflow.type = "HTSeq - FPKM"
                                 ,legacy = FALSE)

GDC = paste0(dir,"/GDCdata")
dir.create(GDC, recursive = T)

tryCatch(GDCdownload(query.exp.hg38, directory=GDC), error = function(e) GDCdownload(query.exp.hg38, method = "client", directory=GDC))

if(format!='none'){
  
  # 1.2 Prepare Data 
  FPKM <- GDCprepare(query.exp.hg38, directory=GDC)
  FPKM <- assay(FPKM, "HTSeq - FPKM" )
  FPKM <- melt(FPKM, varnames = c("ensembl_gene_id","Tumor_Sample_Barcode"))
  
  both = F
  if(format=="both"){
    both=T
  }
  
  if(format=="Rdata" | both==T){
    # 1.3 Save Rdata
    Rdata_dir = paste0(dir,"/Rdata")
    dir.create(Rdata_dir, recursive = T)
    
    cat("Saving expression data as .Rdata \n")
    save(FPKM, file = paste0(Rdata_dir,"/",cancer_type,"_Transcriptome_Profiling.Rdata"))
  } else if(format=="csv" | both==T){
    # 1.3 Save CSV
    Rdata_dir = paste0(dir,"/CSV")
    dir.create(Rdata_dir, recursive = T)
    
    cat("Saving expression data as .csv \n")
    write.csv2(FPKM, file = paste0(Rdata_dir,"/",cancer_type,"_Transcriptome_Profiling.csv"))
  }
  
}




