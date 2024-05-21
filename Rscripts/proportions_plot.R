library(RColorBrewer)

pastel_colors <- brewer.pal(8, "Set3")

anno = read.table("results/annotated_peaks.txt", sep="\t", header=T, quote="")

pietable = table(unlist(lapply(strsplit(as.character(anno$Annotation), " \\("),"[[",1))) #take the part before first ' ('

newtable = table(c("3' UTR","5' UTR","exon","Intergenic","intron","promoter-TSS","TTS"))
newtable[names(newtable)] = 0 #reset everything to 0
newtable[names(pietable)] = pietable

names(newtable) = paste(names(newtable), "(", round(newtable/sum(newtable)*100), "%, ", newtable, ")", sep="")

png("mypeaks.png", width=800, height=600)
pie(newtable, main="Proportions of Regions that have Accessible Chromatin as a Peak", col=pastel_colors)
dev.off()

