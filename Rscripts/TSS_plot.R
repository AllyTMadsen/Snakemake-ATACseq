#!/usr/bin/env Rscript

library(BiocManager)
BiocManager::install(c("ATACseqQC", "ChIPpeakAnno", "MotifDb", "GenomicAlignments",
                       "BSgenome.Hsapiens.UCSC.hg38", "TxDb.Hsapiens.UCSC.hg38.knownGene",
                       "phastCons100way.UCSC.hg38"))
library(ATACseqQC)

# Define the function to make TSS plot
do_TSS <- function(bamfile, output_png) {
    ## bamfile tags to be read in
    possibleTag <- combn(LETTERS, 2)
    possibleTag <- c(paste0(possibleTag[1, ], possibleTag[2, ]),
                 paste0(possibleTag[2, ], possibleTag[1, ]))

    library(Rsamtools)
    bamTop100 <- scanBam(BamFile(bamfile, yieldSize = 100),
                        param = ScanBamParam(tag=unlist(possibleTag)))[[1]]$tag
    tags <- names(bamTop100)[lengths(bamTop100)>0]

    ## files will be output into outPath
    outPath <- "splited"
    #dir.create(outPath)
    ## shift the coordinates of 5'ends of alignments in the bam file
    library(TxDb.Hsapiens.UCSC.hg38.knownGene)
    ## if you don't have an available TxDb, please refer
    ## GenomicFeatures::makeTxDbFromGFF to create one from gff3 or gtf file.
    seqinformation <- seqinfo(TxDb.Hsapiens.UCSC.hg38.knownGene)
    which <- as(seqinformation, "GRanges")
    gal <- readBamFile(bamfile, tag=tags, which=which, asMates=TRUE, bigFile=TRUE)
    shiftedBamfile <- file.path(outPath, "shifted.bam")
    gal1 <- shiftGAlignmentsList(gal, outbam=shiftedBamfile)

    library(TxDb.Hsapiens.UCSC.hg38.knownGene)
    txs <- transcripts(TxDb.Hsapiens.UCSC.hg38.knownGene)
    nfr <- NFRscore(gal1, txs)
    tsse <- TSSEscore(gal1, txs)

    library(BSgenome.Hsapiens.UCSC.hg38)
    library(phastCons100way.UCSC.hg38)
    ## run program for chromosome 1 only
    txs <- txs[seqnames(txs)]
    genome <- Hsapiens
    ## split the reads into NucleosomeFree, mononucleosome, 
    ## dinucleosome and trinucleosome.
    ## and save the binned alignments into bam files.
    objs <- splitGAlignmentsByCut(gal1, txs=txs, outPath = outPath)

    TSS <- promoters(txs, upstream=0, downstream=1)
    TSS <- unique(TSS)
    ## estimate the library size for normalization
    (librarySize <- estLibSize(bamfiles))

    ## calculate the signals around TSSs.
    NTILE <- 101
    dws <- ups <- 1010
    sigs <- enrichedFragments(gal=objs[c("NucleosomeFree", 
                                     "mononucleosome",
                                     "dinucleosome",
                                     "trinucleosome")], 
                            TSS=TSS,
                            librarySize=librarySize,
                            TSS.filter=0.5,
                            n.tile = NTILE,
                            upstream = ups,
                            downstream = dws)

    ## get signals normalized for nucleosome-free and nucleosome-bound regions.
    out <- featureAlignedDistribution(sigs, 
                                    reCenterPeaks(TSS, width=ups+dws),
                                    zeroAt=.5, n.tile=NTILE, type="l", 
                                    ylab="Averaged coverage")

    ## rescale the nucleosome-free and nucleosome signals to 0~1
    range01 <- function(x){(x-min(x))/(max(x)-min(x))}
    out <- apply(out, 2, range01)

    png(output_png, width = 800, height = 600)

    matplot(out, type="l", xaxt="n", 
            xlab="Position (bp)", 
            ylab="Fraction of signal")
    axis(1, at=seq(0, 100, by=10)+1, 
        labels=c("-1K", seq(-800, 800, by=200), "1K"), las=2)
    abline(v=seq(0, 100, by=10)+1, lty=2, col="gray")

    dev.off()
}


# Get command-line arguments
args <- commandArgs(trailingOnly = TRUE)
input_bam <- args[1]
output_png <- args[2]

# Perform QC
do_TSS(input_bam, output_png)

#Rscript ./TSS_plot.R results/ATACrep3.NOmt.aligned.sorted.bam results/ATACrep3_TSSplot.png