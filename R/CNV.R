# R script to download and prepare Copy Number Variation data using TCGAbiolinks
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

# 1. Download CNV data ====
GDC = paste0(dir,"/GDCdata")
dir.create(GDC, recursive = T)

# Include germline mutations
query.cnv_hg38 <- GDCquery(project = cancer_type
                        , data.category = "Copy Number Variation"
                        , data.type = "Copy Number Segment"
                        , legacy = F)
                  
GDCdownload(query.cnv_hg38, directory = GDC)

# Mask germline mutations
query.cnv.mask_hg38 <- GDCquery(project = cancer_type
                           , data.category = "Copy Number Variation"
                           , data.type = "Masked Copy Number Segment"
                           , legacy = F)

GDCdownload(query.cnv.mask_hg38, directory = GDC)

if(format!="none"){
  # 1.1 Prepare Data
  cnv.hg38 <- GDCprepare(query.cnv_hg38, directory=GDC)
  cnv.mask.hg38 <- GDCprepare(query.cnv.mask_hg38, directory=GDC)
  # 1.2 Save Data
  store(format = format, dir = dir, cancer_type = cancer_type, data = cnv.hg38, data_category = "CNV")  
  store(format = format, dir = dir, cancer_type = cancer_type, data = cnv.mask.hg38, data_category = "Masked_CNV")  
}

