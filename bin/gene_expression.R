#!/usr/bin/env Rscript

# args = commandArgs(trailingOnly=TRUE)
# <- args[1]


library(DESeq2)

sampleTable <- data.frame(condition = factor(rep(c("A", "B"), each = 3)))
rownames(sampleTable) <- colnames(txi$counts)
dds <- DESeqDataSetFromTximport(txi, 
                                sampleTable, 
                                ~condition)



dds <- DESeqDataSetFromMatrix(countData=countdata, 
                              colData=sampleTable, 
                              design= as.formula(paste("~", params$condition)))
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design= ~ batch + condition)
dds <- DESeq(dds)
resultsNames(dds) # lists the coefficients
res <- results(dds, name="condition_trt_vs_untrt")
# or to shrink log fold changes association with condition:
res <- lfcShrink(dds, coef="condition_trt_vs_untrt", type="apeglm")