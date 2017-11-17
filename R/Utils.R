# Utils for GDC data download and preparation 

# 1. Resources ====

library(TCGAbiolinks)
library(SummarizedExperiment)
library(reshape2)
library(biomaRt)

# 2. Functions ====

store <- function(format, dir, data, cancer_type, data_category) {
  both = F
  if(format=="both"){
    both=T
    format="Rdata"
  }
  i=1
  j=0
  while(i==1 & j<2){
    i=0
    if(format=="Rdata"){
      # 1.3 Save Rdata
      Rdata_dir = paste0(dir,"/Rdata")
      dir.create(Rdata_dir, recursive = T)
      
      cat(paste0("Saving ", data_category," as .Rdata \n"))
      save(data, file = paste0(Rdata_dir,"/",cancer_type,"_", data_category, ".Rdata"))
    } else if(format=="csv"){
      # 1.3 Save CSV
      Rdata_dir = paste0(dir,"/CSV")
      dir.create(Rdata_dir, recursive = T)
      
      cat(paste0("Saving ", data_category," as .csv \n"))
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