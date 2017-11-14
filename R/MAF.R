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

GDC = paste0(dir,"/GDCdata")
dir.create(GDC, recursive = T)

ct = sapply(strsplit(cancer_type, "[[:punct:]]"), "[[", 2)

muse.maf <- GDCquery_Maf(ct, pipelines = "muse", directory=GDC)
varscan2.maf <- GDCquery_Maf(ct, pipelines = "varscan2", directory=GDC)
somaticsniper.maf <- GDCquery_Maf(ct, pipelines = "somaticsniper", directory=GDC)
mutect.maf <- GDCquery_Maf(ct, pipelines = "mutect", directory=GDC)

if(format!='none'){

  both = F
  if(format=="both"){
    both=T
  }
  
  i=1
  j=0
  while(i==1 & j<2){
    i=0
    if(format=="Rdata" | both==T){
      # 1.1 Save Rdata
      Rdata_dir = paste0(dir,"/Rdata")
      dir.create(Rdata_dir, recursive = T)
      
      cat("Saving MAF as .Rdata \n")
      save(muse.maf, file = paste0(Rdata_dir,"/",cancer_type,"_muse.MAF.Rdata"))
      save(varscan2.maf, file = paste0(Rdata_dir,"/",cancer_type,"_varscan2.MAF.Rdata"))
      save(somaticsniper.maf, file = paste0(Rdata_dir,"/",cancer_type,"_somaticsniper.MAF.Rdata"))
      save(mutect.maf, file = paste0(Rdata_dir,"/",cancer_type,"_mutect.MAF.Rdata"))
    } else if(format=="csv"){
      # 1.2 Save CSV
      Rdata_dir = paste0(dir,"/CSV")
      dir.create(Rdata_dir, recursive = T)
      
      cat("Saving MAF as .csv \n")
      write.csv2(muse.maf, file = paste0(Rdata_dir,"/",cancer_type,"_muse.MAF.csv"))
      write.csv2(varscan2.maf, file = paste0(Rdata_dir,"/",cancer_type,"_varscan2.MAF.csv"))
      write.csv2(somaticsniper.maf, file = paste0(Rdata_dir,"/",cancer_type,"_somaticsniper.MAF.csv"))
      write.csv2(mutect.maf, file = paste0(Rdata_dir,"/",cancer_type,"_mutect.MAF.csv"))
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



