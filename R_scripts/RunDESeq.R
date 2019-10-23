suppressPackageStartupMessages({library(optparse)})

option_list = list(
make_option(c("-i", "--infile"), type="character", default="./", help="Input file", metavar="character"),
make_option(c("-d", "--design_sheet"), type="character", default="./", help="Design sheet", metavar="character"),
make_option(c("-o", "--outpath"), type="character", default="./", help="Path to output files [default= %default]", metavar="character"))

opt_parser = OptionParser(option_list=option_list);
options = parse_args(opt_parser);

if (is.null(options$infile)) {
	print ("Error: No input file specified.")
	stop()
} else {
	cat ("Input coverage file: ", options$infile, "\n")
	in_file <- options$infile
}


if (is.null(options$design_sheet)) {
	print ("Error: No design_sheet specified.")
	stop()
} else {
	cat ("Design_sheet: ", options$design_sheet, "\n")
	d_sheet <- options$design_sheet
}


if (is.null(options$outpath)) {
	print ("Warning: No path to output files specified. Set as default: current directory")
} else {
	cat ("Path to output files: ", options$outpath, "\n")
	out_path <- options$outpath
}



suppressPackageStartupMessages({library(DESeq2)})


cov_data <- read.table(in_file, sep="\t", row.names = 1, header=T)
design_table <- read.table(d_sheet, sep="\t", row.names = 1, header=T)
conditions <- as.character(design_table[,1])

exp_cond <- data.frame(row.names=colnames(cov_data), condition = conditions)
DESeq_dataset <- DESeqDataSetFromMatrix(countData=cov_data, colData=exp_cond, design=~condition)
print ("Start DESeq2")
DESeq_norm <- DESeq(DESeq_dataset)

#pdf(paste(out_path, "Dispersion.pdf", sep=""))
#plotDispEsts(DESeq_norm)
#dev.off()

rlog_tf <- rlog(DESeq_norm)
distances <- dist(t(assay(rlog_tf)))
distance_matrix <- as.matrix(distances)
rownames(distance_matrix) <- colnames(distance_matrix)

suppressPackageStartupMessages({library(RColorBrewer)})
suppressPackageStartupMessages({library(gplots)})

palette <-colorRampPalette(brewer.pal(9,"GnBu"))(100)
hierarchy <- hclust(distances)

print ("Print Dendrogram_and_heatmap.")
pdf(paste(out_path, "Dendrogram_and_heatmap.pdf", sep=""))
heatmap.2(distance_matrix,Rowv=as.dendrogram(hierarchy), symm=TRUE, trace="none", col=rev(palette))
dev.off()
print ("Done.")

print ("Print PCA plot.")
pdf(paste(out_path, "PCA.pdf", sep=""))
plotPCA(rlog_tf, intgroup=c("condition"))
dev.off()
print ("Done.")

unique_conditions <- unique(conditions)

print ("Write constrast table.")
combinations <- combn(unique_conditions,2)
num_combinations <- ncol(combinations)
for (i in 1:num_combinations){
	sample1 <- combinations[,i][1]
	sample2 <- combinations[,i][2]
	result <- results(DESeq_norm, contrast=c("condition", sample1, sample2))
	write.table(result, file=paste(out_path, "Contrast_", sample1, "_vs_", sample2, ".txt", sep= ""), sep="\t")
}
print ("Done.")

DESeq_estSF <- estimateSizeFactors(DESeq_norm)
Norm_counts <- counts(DESeq_estSF, normalized=TRUE)
print ("Write normalized count table.")
write.table(Norm_counts, paste(out_path,"Normalized_expression.txt", sep=""), sep="\t")
print ("Done.")