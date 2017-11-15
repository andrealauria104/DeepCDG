# R script to download and prepare DNA methylation data using TCGAbiolinks
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

GDCdownload(query_met.hg38, directory=GDC)

if(format!="none"){
  # 1.1 Prepare Data
  met.hg38 <- GDCprepare(query_met.hg38, directory=GDC)
  # 1.2 Save Data
  store(format = format, dir = dir, cancer_type = cancer_type, data = met.hg38, data_category = "Methylation")  
}


