# R script to download and prepare gene expression data using TCGAbiolinks
source("R/Utils.R")
if(!dir.exists("DownloadEnv")){
  cat("Creating download environment ... \n ")
  dir.create("DownloadEnv")
}
setwd("DownloadEnv")

# 1.0 Download Gene Expression data ====

# 1.1 Query GDC
query.exp.hg38 <- GDCquery( project = cancer_type
                                 ,data.category = "Transcriptome Profiling"
                                 ,data.type = "Gene Expression Quantification"
                                 ,workflow.type = "HTSeq - FPKM"
                                 ,legacy = FALSE)

GDC = paste0(dir,"/GDCdata")
dir.create(GDC, recursive = T)

tryCatch(GDCdownload(query.exp.hg38, directory=GDC), error = function(e) GDCdownload(query.exp.hg38, method = "client", directory=GDC))

if(format!='none'){
  # 1.2 Prepare Data 
  cat("Preparing data ... \n")
  FPKM <- GDCprepare(query.exp.hg38, directory=GDC)
  FPKM <- assay(FPKM, "HTSeq - FPKM" )
  FPKM <- melt(FPKM, varnames = c("ensembl_gene_id","Tumor_Sample_Barcode"))
  
  cat("Adding information from biomart ... \n")
  mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
  
  # query biomart
  gene_symb <- getBM(attributes = c("ensembl_gene_id"
                                    ,'hgnc_symbol'
                                    ,'external_gene_name'
                                    ,'entrezgene'
                                    ,'gene_biotype')
                     ,mart = mart)
  
  FPKM$symbol=NA
  FPKM$symbol=gene_symb$external_gene_name[match(FPKM$ensembl_gene_id,gene_symb$ensembl_gene_id)]
  FPKM$gene_type=NA
  FPKM$gene_type=gene_symb$gene_biotype[match(FPKM$ensembl_gene_id,gene_symb$ensembl_gene_id)]
  
  # Sample stratification 
  table.code <- c("TP", "TR", "TB", "TRBM", "TAP", "TM",
                  "TAM", "THOC", "TBM", "NB", "NT", "NBC", "NEBV", "NBM",
                  "CELLC", "TRB", "CELL", "XP", "XCL")
  names(table.code) <- c("01", "02", "03", "04", "05", "06", "07",
                         "08", "09", "10", "11", "12", "13", "14", "20", "40",
                         "50", "60", "61")
  
  string <- substr(FPKM$Tumor_Sample_Barcode , 14, 15)
  FPKM$sample=table.code[string]
  FPKM$sample=factor(FPKM$sample, levels=c("NT","TP","TM"))
  
  store(format = format, dir = dir, cancer_type = cancer_type, data = FPKM, data_category = "Transcriptome_Profiling")  
}
  





