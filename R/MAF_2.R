# R script to download and prepare mutation data (MAF) using TCGAbiolinks
source("R/Utils.R")
setwd("DownloadEnv")

# 1.0 Download MAF - 4 different pipelines ====
if(format=='none'){
  stop(message("Please specify data format for storing prepared data."))
}

GDC = paste0(dir,"/GDCdata")
dir.create(GDC, recursive = T)

ct = sapply(strsplit(cancer_type, "[[:punct:]]"), "[[", 2)

# Download from GDC 
muse.maf <- GDCquery_Maf(ct, pipelines = "muse", directory=GDC)
varscan2.maf <- GDCquery_Maf(ct, pipelines = "varscan2", directory=GDC)
somaticsniper.maf <- GDCquery_Maf(ct, pipelines = "somaticsniper", directory=GDC)
mutect.maf <- GDCquery_Maf(ct, pipelines = "mutect", directory=GDC)

# Save files
datas = ls()[grep("maf", ls())]
for(i in datas){
  store(format = format, dir = dir, cancer_type = cancer_type, data = get(i), data_category = i)  
}
