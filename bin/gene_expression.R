#!/usr/bin/env Rscript

library(biomaRt)
library(DESeq2)
library(dplyr)
library(tximport)
library(RColorBrewer)
library(gplots)
library(tidyverse)

args = commandArgs(trailingOnly=TRUE)
annotation    <- args[1] # experiment.csv
condition     <- args[2] # treatment
kallisto_dir  <- args[3] # kallisto

###############
# Preprocessing
###############

# 1. Load the abundance files
sample_ids <- dir(file.path(kallisto_dir))
sample_ids <- sub("^kallisto_", "", sample_ids)
abundances <- file.path(kallisto_dir, paste0('kallisto_', sample_ids), 'abundance.h5')

# 2. Load the annotation (sample table with metadata)
sampleTable <- read.csv(annotation, row.names = 1)

# 3. Set the condition of interest
# Select/extract the one condition we want to analyse, eg the treatment was used on each sample:
condition <- unlist(sampleTable[condition])

###################################
# Differential Expression Analysis
###################################

# Load the gene info
mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",
  dataset = "hsapiens_gene_ensembl",
  host = "may2015.archive.ensembl.org")
ttg <- biomaRt::getBM(
  attributes = c("ensembl_transcript_id", "transcript_version",
  "ensembl_gene_id", "external_gene_name", "description",
  "transcript_biotype"),
  mart = mart)
ttg <- dplyr::select(ttg, TXNAME = ensembl_transcript_id,
  GENEID = ensembl_gene_id)

# Create the DEseq dataset object using the abundances and the sample data
txi_kallisto <- tximport(abundances, type='kallisto', tx2gene=ttg)
dds <- DESeqDataSetFromTximport(txi=txi_kallisto, colData=sampleTable, design= as.formula(paste("~", condition)))

# Run the Differential Expression
dds <- DESeq(dds)

# Dispersion plot
# Look at the dispersion of read counts across genes in all samples
png("dispersion_plot.png")
plotDispEsts(dds, main="Dispersion plot")

# Log Transformation & Histogram
# Apply log transformation to the read counts and will try to cluster the samples according the log transform counts
# Regularized log transformation for clustering/heatmaps, etc
rld <- rlogTransformation(dds)
png("dispersion_plot.png")
hist(assay(rld))

# Sample Distance Matrix
(mycols <- brewer.pal(8, "Dark2")[1:length(unique(condition))])
sampleDists <- as.matrix(dist(t(assay(rld))))

# Calculating distances across samples from log transformed counts to plot the distance matrix and colour the samples according to the condition
png("sample_distance_matrix.png")
heatmap.2(as.matrix(sampleDists), key=F, trace="none",
          col=colorpanel(100, "black", "white"),
          ColSideColors=mycols[condition], RowSideColors=mycols[condition],
          margin=c(10, 10), main="Sample Distance Matrix")

# PCA Biplot
rld_pca <- function (rld, intgroup = "condition", ntop = 500, colors=NULL, legendpos="bottomleft", main="PCA Biplot", textcx=1, ...) {
  require(genefilter)
  require(calibrate)
  require(RColorBrewer)
  rv = rowVars(assay(rld))
  select = order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  pca = prcomp(t(assay(rld)[select, ]))
  fac = factor(apply(as.data.frame(colData(rld)[, intgroup, drop = FALSE]), 1, paste, collapse = " : "))
  if (is.null(colors)) {
    if (nlevels(fac) >= 3) {
      colors = brewer.pal(nlevels(fac), "Paired")
    }   else {
      colors = c("black", "red")
    }
  }
  pc1var <- round(summary(pca)$importance[2,1]*100, digits=1)
  pc2var <- round(summary(pca)$importance[2,2]*100, digits=1)
  pc1lab <- paste0("PC1 (",as.character(pc1var),"%)")
  pc2lab <- paste0("PC1 (",as.character(pc2var),"%)")
  plot(PC2~PC1, data=as.data.frame(pca$x), bg=colors[fac], pch=21, xlab=pc1lab, ylab=pc2lab, main=main, ...)
  with(as.data.frame(pca$x), textxy(PC1, PC2, labs=rownames(as.data.frame(pca$x)), cex=textcx))
  legend(legendpos, legend=levels(fac), col=colors, pch=20)
}
# PCA visualization to see how the variation separates samples, but also how mixed
png("pca_biplot.png")
rld_pca(rld, colors=mycols, intgroup=params$condition, xlim=c(-75, 35))

# Results table of top genes by adjusted p-value
# Filter and to extract genes based on certain thresholds. Apply filter in order to filter for p-values and adjusted p-values on top of the results table.
res <- results(dds)
# Order by adjusted p-value
res <- res[order(res$padj), ]
# Merge with normalized count data
resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
names(resdata)[1] <- "Gene"

# Histogram of p-value
# Plot for quality control on top of the results we got that we can explore, the histogram of significance levels
png("histogram_of_p_values.png")
hist(res$pvalue, breaks=50, col="grey")

# MA Plot
maplot <- function (res, thresh=0.05, labelsig=TRUE, textcx=1, ...) {
  with(res, plot(baseMean, log2FoldChange, pch=20, cex=.5, log="x", ...))
  with(subset(res, padj<thresh), points(baseMean, log2FoldChange, col="red", pch=20, cex=1.5))
  if (labelsig) {
    require(calibrate)
    with(subset(res, padj<thresh), textxy(baseMean, log2FoldChange, labs=Gene, cex=textcx, col=2))
  }
}
# MA plot, where we can see how the log fold changes compare against the mean averages across genes
png("ma_plot.png")
maplot(resdata, main="MA Plot")

# Volcano plot
volcanoplot <- function (res, lfcthresh=2, sigthresh=0.05, main="Volcano Plot", legendpos="bottomright", labelsig=TRUE, textcx=1, ...) {
  with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main=main, ...))
  with(subset(res, padj<sigthresh ), points(log2FoldChange, -log10(pvalue), pch=20, col="red", ...))
  with(subset(res, abs(log2FoldChange)>lfcthresh), points(log2FoldChange, -log10(pvalue), pch=20, col="orange", ...))
  with(subset(res, padj<sigthresh & abs(log2FoldChange)>lfcthresh), points(log2FoldChange, -log10(pvalue), pch=20, col="green", ...))
  if (labelsig) {
    require(calibrate)
    with(subset(res, padj<sigthresh & abs(log2FoldChange)>lfcthresh), textxy(log2FoldChange, -log10(pvalue), labs=Gene, cex=textcx, ...))
  }
  legend(legendpos, xjust=1, yjust=1, legend=c(paste("FDR<",sigthresh,sep=""), paste("|LogFC|>",lfcthresh,sep=""), "both"), pch=20, col=c("red","orange","green"))
}
png("volcano_plot.png")
volcanoplot(resdata, lfcthresh=2, sigthresh=0.01, textcx=.8, xlim=c(-2.3, 2))

#######################################################
# Write the differential expression results to files
#######################################################
# Save plots
dev.off()
# Write the DESeq2 results
res <- results(dds, tidy = TRUE)
readr::write_csv(res, path="deseq_results.csv")	