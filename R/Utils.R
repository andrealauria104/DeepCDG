# Utils for GDC data download and preparation 

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

# 1.0 Functions ====

store <- function(format, dir, data, cancer_type, data_category) {
  both = F
  if(format=="both"){
    both=T
  }
  i=1
  j=0
  while(i==1 & j<2){
    i=0
    if(format=="Rdata" | both==T){
      # 1.3 Save Rdata
      Rdata_dir = paste0(dir,"/Rdata")
      dir.create(Rdata_dir, recursive = T)
      
      cat("Saving expression data as .Rdata \n")
      save(data, file = paste0(Rdata_dir,"/",cancer_type,"_", data_category, ".Rdata"))
    } else if(format=="csv"){
      # 1.3 Save CSV
      Rdata_dir = paste0(dir,"/CSV")
      dir.create(Rdata_dir, recursive = T)
      
      cat("Saving expression data as .csv \n")
      write.csv2(data, file = paste0(Rdata_dir,"/",cancer_type,"_", data_category, ".csv"))
    } else{
      stop(message("Invalid format for storing prepared data."))
    }
    if(both){
      i=1
      format="csv"
    }
    j=j+1
  }
}