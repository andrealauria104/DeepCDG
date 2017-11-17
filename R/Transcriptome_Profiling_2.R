# R script to download and prepare gene expression data using TCGAbiolinks
source("R/Utils.R")
if(!dir.exists("DownloadEnv")){
  cat("Creating download environment ... \n ")
  dir.create("DownloadEnv")
}
setwd("DownloadEnv")

# 0. Args Parsing ====
args <- commandArgs(trailingOnly = T)

parserr <- 'Incorrect argument parsing. Please specify: \n
1 - TCGA project (Cancer type). \n
2 - Normalization. \n
3 - Data directory. \n
4 - Format for storing prepared data (if needed: Rdata, csv or both) \n '

if(length(args)<3){
  stop(message(parserr))
}

# Define parsing arguments
cancer_type = as.character(args[1])
my.norm = as.character(args[2])
dir = as.character(args[3])

if(length(args)<3){
  format = 'none'
} else{
  format = as.character(args[4])
}

# Directory control
dir1 = dirname(dir)
dir2 = basename(dir)
dir = paste0(dir1,"/",dir2)
rm(list=c("dir1","dir2"))

# 1.0 Download Gene Expression data ====
GDC = paste0(dir,"/GDCdata")
dir.create(GDC, recursive = T)

if(my.norm=="FPKM"){
  # 1.1 Query GDC
  query.exp.hg38 <- GDCquery( project = cancer_type
                              ,data.category = "Transcriptome Profiling"
                              ,data.type = "Gene Expression Quantification"
                              ,workflow.type = "HTSeq - FPKM"
                              ,legacy = FALSE)
   # 1.2 Download
  tryCatch(GDCdownload(query.exp.hg38, directory=GDC), error = function(e) GDCdownload(query.exp.hg38, method = "client", directory=GDC))
} else if(my.norm=="COUNTS"){
  # 1.1 Query GDC
  query.exp.hg38 <- GDCquery( project = cancer_type
                              ,data.category = "Transcriptome Profiling"
                              ,data.type = "Gene Expression Quantification"
                              ,workflow.type = "HTSeq - Counts"
                              ,legacy = FALSE)
  # 1.2 Download
  tryCatch(GDCdownload(query.exp.hg38, directory=GDC), error = function(e) GDCdownload(query.exp.hg38, method = "client", directory=GDC))
}

if(format!='none'){
  # 1.3 Prepare Data 
  cat("Preparing data ... \n")
  rnasq <- GDCprepare(query.exp.hg38, directory=GDC)
  
  if(my.norm=="FPKM"){
    expr <- assay(rnasq, "HTSeq - FPKM" )
  } else if(my.norm=="COUNTS"){
    expr <- assay(rnasq, "HTSeq - Counts" )
  }
  expr <- melt(expr, varnames = c("ensembl_gene_id","Tumor_Sample_Barcode"))
  
  cat("Adding information from biomart ... \n")
  mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
  # query biomart
  gene_symb <- getBM(attributes = c("ensembl_gene_id"
                                    ,'hgnc_symbol'
                                    ,'external_gene_name'
                                    ,'entrezgene'
                                    ,'gene_biotype')
                     ,mart = mart)
  
  expr$symbol=NA
  expr$symbol=gene_symb$external_gene_name[match(expr$ensembl_gene_id,gene_symb$ensembl_gene_id)]
  expr$gene_type=NA
  expr$gene_type=gene_symb$gene_biotype[match(expr$ensembl_gene_id,gene_symb$ensembl_gene_id)]
  
  # Sample stratification 
  table.code <- c("TP", "TR", "TB", "TRBM", "TAP", "TM",
                  "TAM", "THOC", "TBM", "NB", "NT", "NBC", "NEBV", "NBM",
                  "CELLC", "TRB", "CELL", "XP", "XCL")
  names(table.code) <- c("01", "02", "03", "04", "05", "06", "07",
                         "08", "09", "10", "11", "12", "13", "14", "20", "40",
                         "50", "60", "61")
  
  string <- substr(expr$Tumor_Sample_Barcode , 14, 15)
  expr$sample=table.code[string]
  expr$sample=factor(expr$sample, levels=c("NT","TP","TM"))
  
  if(my.norm=="FPKM"){
    store(format = format, dir = dir, cancer_type = cancer_type, data = expr, data_category = "Transcriptome_Profiling_FPKM")    
  } else if(my.norm=="COUNTS"){
    store(format = format, dir = dir, cancer_type = cancer_type, data = expr, data_category = "Transcriptome_Profiling_COUNTS")
  }
  
}
  





