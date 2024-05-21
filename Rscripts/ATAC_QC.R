#!/usr/bin/env Rscript

library(BiocManager)
BiocManager::install(c("ATACseqQC", "ChIPpeakAnno", "MotifDb", "GenomicAlignments",
                       "BSgenome.Hsapiens.UCSC.hg38", "TxDb.Hsapiens.UCSC.hg38.knownGene",
                       "phastCons100way.UCSC.hg38"))
library(ATACseqQC)

# Define the function to perform QC
do_QC <- function(bamfile, output_png) {
  bamfile.labels <- gsub(".NOmt.filtered.sorted.bam", "", basename(bamfile))
  
  # Open a PNG device for graphics output
  png(output_png)
  
  # Generate fragment size distribution plot
  fragSize <- fragSizeDist(bamfile, bamfile.labels)
  
  # Close the PNG device to save the plot
  dev.off()
}

# Get command-line arguments
args <- commandArgs(trailingOnly = TRUE)
input_bam <- args[1]
output_png <- args[2]

# Perform QC
do_QC(input_bam, output_png)

#test locally:
#Rscript ./ATAC_QC.R results/ATACrep3.NOmt.aligned.filtered.bam results/ATACrep3_FragDistPlot.png
