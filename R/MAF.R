# R script to download and prepare mutation data (MAF) using TCGAbiolinks

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
dir = paste0(dir1,"/",dir2)
rm(list=c("dir1","dir2"))

# 0.2 Resources ====

library(TCGAbiolinks)
library(SummarizedExperiment)
library(reshape2)

# 1.0 Download MAF - 4 different pipelines ====

ct = sapply(strsplit(cancer_type, "[[:punct:]]"), "[[", 2)

muse.maf <- GDCquery_Maf(ct, pipelines = "muse")
varscan2.maf <- GDCquery_Maf(ct, pipelines = "varscan2")
somaticsniper.maf <- GDCquery_Maf(ct, pipelines = "somaticsniper")
mutect.maf <- GDCquery_Maf(ct, pipelines = "mutect")


