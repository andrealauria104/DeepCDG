# R script to download and prepare Copy Number Variation data using TCGAbiolinks
source("R/Utils.R")
if(!dir.exists("DownloadEnv")){
  cat("Creating download environment ... \n ")
  dir.create("DownloadEnv")
}
setwd("DownloadEnv")

# 1. Download CNV data ====
GDC = paste0(dir,"/GDCdata")
dir.create(GDC, recursive = T)

# Include germiline mutations
query.cnv_hg38 <- GDCquery(project = cancer_type
                        , data.category = "Copy Number Variation"
                        , data.type = "Copy Number Segment"
                        , legacy = F)
                  
GDCdownload(query.cnv_hg38, directory = GDC)

# Mask germiline mutations
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

