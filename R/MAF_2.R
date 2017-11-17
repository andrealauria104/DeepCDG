# R script to download and prepare mutation data (MAF) using TCGAbiolinks
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
