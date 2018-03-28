# R script to download and prepare DNA methylation data using TCGAbiolinks
# 0. Args Parsing ====
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
dir = paste0(dir1,"/",dir2)
rm(list=c("dir1","dir2"))
# utils
source("R/Utils.R")

if(!dir.exists("DownloadEnv")){
  cat("Creating download environment ... \n ")
  dir.create("DownloadEnv")
}
setwd("DownloadEnv")

# 1.0 Download Methylation Data ====
query_met.hg38 <- GDCquery(project = cancer_type
                              , data.category = "DNA Methylation"
                              , platform = "Illumina Human Methylation 450"
                              , legacy = F )

GDC = paste0(dir,"/GDCdata")
dir.create(GDC, recursive = T)

GDCdownload(query_met.hg38, directory=GDC)

if(format!="none"){
  # 1.1 Prepare Data
  met.hg38 <- GDCprepare(query_met.hg38, directory=GDC)
  
  if(format=="csv"){
    met <- assay(met.hg38) 
  }
  # 1.2 Save Data
  store(format = format, dir = dir, cancer_type = cancer_type, data = met, data_category = "Methylation")  
}


