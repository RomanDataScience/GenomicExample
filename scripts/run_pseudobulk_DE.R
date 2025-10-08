#!/usr/bin/env Rscript
library(DESeq2)
counts <- read.csv("results/tables/pseudobulk_counts.csv", row.names=1)
meta <- read.csv("results/tables/pseudobulk_metadata.csv", row.names=1)
dds <- DESeqDataSetFromMatrix(counts, meta, design = ~ diagnosis + age + sex)
dds <- DESeq(dds)
res <- results(dds, contrast=c("diagnosis","AD","Control"))
write.csv(as.data.frame(res), "results/tables/deseq2_results.csv")
cat("✅ DESeq2 results saved.\n")
