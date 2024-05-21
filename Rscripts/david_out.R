
anno = read.table("results/annotated_peaks.txt", sep="\t", header=T, quote="", fill=TRUE)
anno <- anno[complete.cases(anno$Annotation), ]

print(colnames(anno))

#get "Gene.Name" col where where col "Annotation" == "promotor-TSS" and send to a new file:
#promoter_TSS_data <- subset(anno, Annotation == "promoter-TSS", select = "Gene.Name")

annotations <- unlist(lapply(strsplit(as.character(anno$Annotation), " \\("),"[[",1))
print(tail(annotations))
print(class(anno$Annotation))
print(sum(annotations == "promoter-TSS"))

promoter_TSS_data <- anno$Entrez.ID[annotations == "promoter-TSS"]

write.table(promoter_TSS_data, file = "promoter_TSS_genes.txt", row.names = FALSE, col.names = FALSE)