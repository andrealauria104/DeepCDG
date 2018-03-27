# Prepare CNV data
# 0. Resources ====
library(GenomicRanges)
library(plyr)

CNVDIR <- "/sto1/epigenetics/DB/GDC/171110_TCGAbiolinks/CSV"
ANNOT  <- "/sto1/epigenetics/DB/GENCODE/v27/ucsc_gene_longest_transcript.txt"
OUTDIR <- "/sto1/epigenetics/DB/GDC/171110_TCGAbiolinks/Rdata"

table.code <- c("TP", "TR", "TB", "TRBM", "TAP", "TM",
                "TAM", "THOC", "TBM", "NB", "NT", "NBC", "NEBV", "NBM",
                "CELLC", "TRB", "CELL", "XP", "XCL")
names(table.code) <- c("01", "02", "03", "04", "05", "06", "07",
                       "08", "09", "10", "11", "12", "13", "14", "20", "40",
                       "50", "60", "61")

project <- as.character(commandArgs(trailingOnly = T)[1])

# 1. Read CNV Data ====
cat("Reading CNV data ... \n")
CNVDATA <- list.files(CNVDIR, pattern = paste0(project,"_CNV"), 
                      full.names = T)
cnv <- read.csv2(CNVDATA, 
                 header = T,
                 stringsAsFactors = F)
cnv$X <- NULL
colnames(cnv)[1] <- "Barcode"
cnv$sample <- table.code[substr(cnv$Barcode, 14,15)]

cnv <- subset(cnv, sample=="TP")
cnv <- subset(cnv, abs(Segment_Mean)>0.3 & Segment_Mean<1.5)

cnv$assignedCNV <- ifelse(cnv$Segment_Mean>0.3, "Gain","Loss")
cnv$tcn =  cnv$value = round((2^cnv$Segment_Mean)*2)

# 2. Read gene annotation data =====
cat("Reading gene annotation ... \n")
annot        <- read.delim2(ANNOT, 
                            skip = 1    , 
                            stringsAsFactors = F)
annot$geneId <- sapply(strsplit(annot$geneId, 
                                "[.]"), 
                       "[[", 1)
annot <- subset(annot, transcriptType=="protein_coding")

gene_ids <- annot[,c("geneId","chr","start","end")]
gene_ids <- gene_ids[order(gsub("chr","",gene_ids[,2]),
                           gene_ids[,3]),]

# 3. Find overlaps ====
cat("Find gene/segment overlaps ... \n")
find_geneSegment_overlaps <- function(x, y){
  
  ix <- with(x,GRanges(seqnames = chr, IRanges(start = start, end = end)))  
  iy <- with(y,GRanges(seqnames = paste0('chr', y$Chromosome), IRanges(start = Start, end = End)))
  # countOverlaps(ix,iy)
  o=findOverlaps(ix,iy)
  cnv_overlaps <- cbind(x[queryHits(o),], y[subjectHits(o),])   
  return(cnv_overlaps)
}

cnv_overlaps <- find_geneSegment_overlaps(x = gene_ids, 
                                          y = cnv)
# 4. Compute stats ====
cat("Compute stats ... \n")
df <- unique(cnv_overlaps[,c("geneId","Barcode","assignedCNV", "value")])  

# By gene
cnv_genes <- ddply(df, .(geneId), mutate,
                   Gain      = table(assignedCNV)["Gain"],
                   Loss      = table(assignedCNV)["Loss"]
)
nsamples <- length(unique(cnv_genes$Barcode))
cnv_genes <- ddply(cnv_genes, .(geneId), mutate,
                   totCNV = Gain + Loss,
                   freqGain = Gain/nsamples,
                   freqLoss = Loss/nsamples
)

# By sample
cnv_samples <- ddply(df, .(Barcode), mutate,
                     Gain      = table(assignedCNV)["Gain"],
                     Loss      = table(assignedCNV)["Loss"]
)
cnv_samples <- ddply(cnv_samples, .(Barcode), mutate,
                     totCNV = Gain + Loss
)
cnv_samples <- ddply(cnv_samples, .(Barcode), mutate,
                     freqGain = Gain/totCNV,
                     freqLoss = Loss/totCNV)

# 5. Save files ====
cat("Saving files ...\n")
save(df, cnv_genes, cnv_samples, 
     file=paste0(OUTDIR,"/",
                     project,"_",
                     "CNV_gene_level.Rdata")
      )

cat(" ... all done .\n")
