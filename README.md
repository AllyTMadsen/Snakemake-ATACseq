# bf528-individual-project-AllyTMadsen

<h2>Project Proposal: ATACseq Analysis</h2>

<h4>Methods:</h4>

        Quality control (QC) was conducted on all four samples using FastQC version 0.12.1-0 with default parameters. The Nextera PE adapter sequences were selected for trimming based on the FastQC results. Adapter trimming was then performed using Trimmomatic version 0.39 with the command "trimmomatic PE ILLUMINACLIP:NexteraPE-PE.fa:2:30:10" and a minimum length parameter of 30. Subsequently, QC was performed on the trimmed reads using FastQC version 0.12.1-0, followed by compilation of a final QC report for all samples using MultiQC version 1.20.

        Bowtie2 version 2.5.3 was utilized to build an index for the provided reference genome (hg38), utilizing default parameters. Alignment of forward and reverse paired reads for each sample to the reference genome (hg38) was then conducted using Bowtie2 with parameters "--very-sensitive -X2000 -k 10", allowing a maximum fragment length of 2000 for valid paired-end alignments, and retaining the top 10 distinct, valid alignments for each read.

        Following alignment, the resulting BAM files were indexed and sorted using Samtools version 1.19.2. Samtools idxstats was performed on the indexed BAM files to quantify reads mapped to the mitochondrial chromosome, and any alignments to the mitochondrial chromosome were subsequently removed using Samtools view version 1.19.2 and command line arguments.

        After filtering, the BAM files were indexed and sorted again using Samtools version 1.19.2, and Samtools idxstats was rerun on the filtered indexed BAM files to quantify reads removed due to mapping to the mitochondrial chromosome. To address bias induced by the tagmentation process, reads were shifted using Deeptools version 3.5.4 with the method "alignmentSieve" and the flag "--ATACshift", which is equivalent to --shift 4 -5 5 -4.

        A quality control analysis plotting fragment distribution sizes for the samples was conducted using ATACSeqQC 1.26.0, generating an R script (version 4.1.1) to perform QC checks per sample. Peak calling was performed with MACS3 version 3.0.1 on the paired-end data using the "callPeaks" function with a minimum false discovery rate (FDR) cutoff value of 0.01.

        The intersection between replicates was identified with Bedtools version 2.31.1, retaining peaks with at least 50% overlap in both replicates (-f 0.5 and -r) as specified by Encode's data processing standards for ATACseq data. Reproducible peaks were filtered using Bedtools version 2.31.1 to remove peaks falling into blacklisted regions. HOMER version 4.11 annotatePeaks and findMotifsGenome was used to annotate peaks and identify motifs in the reproducable peaks set using the human reference genome hg38, the -go flag was added to also perform gsea analysis while annotating peaks and motifs were limited to a region size of 200.

        Gene set enrichment analysis was primarily performed using DAVID, and later compared to the results from HOMER. An R script (version 4.1.1) was generated to reformat the Entrez IDs for the reproducable peaks annotated as "promotor TSS" for compatibility with DAVID's web-based interface. This selection was made as to only perform gene set enrichment on peaks found in the TSS site.

        Deeptools version 3.5.4 bamCoverage, computeMatrix, and plotProfile were utilized to create a TSS signal plot displaying the NFR and NBR regions centered around the TSS site, with NBR regions defined to have a minimum fragment length of 100bp, and the NFR region included all other fragments. Reference point was specified in ComputeMatrix to center the data around the TSS, and --perGroup was utilized to plot both signals on one figure. Finally, a pie chart was created using base R functions to display the proportions of regions with accessible chromatin called as a peak.


<h4>Questions to Address:</h4>

        Quality of the sequencing reads and the alignment statistics:
            Are there any concerning aspects of the quality control of your sequencing reads?
            Are there any concerning aspects of the quality control related to alignment?
            Based on all of your quality control, will you exclude any samples from further analysis?

        Calculate how many alignments were generated from each sample in total and how many alignments were against the mitochondrial chromosome:
            Report the total number of alignments per sample
            Report the number of alignments against the mitochondrial genome

        After performing peak calling analysis, generating a set of reproducible peaks and filtering peaks from blacklisted regions:
            How many peaks are present in each of the replicates?
            How many peaks are present in your set of reproducible peaks? What strategy did you use to determine “reproducible” peaks?
            How many peaks remain after filtering out peaks overlapping blacklisted regions?

        After performing motif analysis and gene enrichment on the peak annotations:
            Briefly discuss the main results of both of these analyses
            What can chromatin accessibility let us infer biologically?

<h4>Deliverables:</h4>

        1. Produce a fragment length distribution plot for each of the samples

        2. Produce a table of how many alignments for each sample before and after filtering alignments falling on the mitochondrial chromosome

        3. Create a signal coverage plot centered on the TSS (plotProfile) for the nucleosome-free regions (NFR) and the nucleosome-bound regions (NBR)
            -You may consider fragments (<100bp) to be those from the NFR and the rest as the NBR.
   
        4. A table containing the number of peaks called in each replicate, and the number of reproducible peaks

        5. A single BED file containing the reproducible peaks you determined from the experiment.

        6. Perform motif finding on your reproducible peaks
            -Create a single table / figure with the most interesting results

        7. Perform a gene enrichment analysis on the annotated peaks using a well-validated gene enrichment tool
            -Create a single table / figure with the most interesting results

        8. Produce a figure that displays the proportions of regions that appear to have accessible chromatin called as a peak (Promoter, Intergenic, Intron, Exon, TTS, etc.)

