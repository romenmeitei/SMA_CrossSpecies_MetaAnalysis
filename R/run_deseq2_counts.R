#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(DESeq2))
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) stop("Usage: run_deseq2_counts.R COUNTS_TSV SAMPLE_MAP ACCESSION OUTPUT_CSV")
counts_path <- args[[1]]
sample_map_path <- args[[2]]
accession <- args[[3]]
output <- args[[4]]

count_table <- read.delim(counts_path, check.names = FALSE, row.names = 1)
metadata <- read.csv(sample_map_path, row.names = 1, stringsAsFactors = FALSE)
missing <- setdiff(rownames(metadata), colnames(count_table))
if (length(missing)) stop("Counts missing annotated samples: ", paste(missing, collapse = ", "))
count_table <- count_table[, rownames(metadata), drop = FALSE]
if (any(abs(count_table - round(count_table)) > 1e-8, na.rm = TRUE)) stop("DESeq2 requires integer counts")
metadata$condition <- factor(metadata$condition, levels = c("control", "disease"))
dds <- DESeqDataSetFromMatrix(round(as.matrix(count_table)), metadata, design = ~ condition)
keep <- rowSums(DESeq2::counts(dds) >= 10) >= 2
dds <- DESeq(dds[keep, ], quiet = TRUE)
qc_dir <- file.path(dirname(dirname(output)), "qc")
dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)
vsd <- vst(dds, blind = FALSE)
pdf(file.path(qc_dir, paste0(accession, "_QC.pdf")), width = 8, height = 6)
pca <- prcomp(t(assay(vsd)))
colors <- c(control = "#26619C", disease = "#B22222")
plot(pca$x[, 1], pca$x[, 2], col = colors[as.character(metadata$condition)], pch = 19,
     xlab = "PC1", ylab = "PC2", main = paste(accession, "variance-stabilized PCA"))
text(pca$x[, 1], pca$x[, 2], labels = rownames(metadata), pos = 3, cex = 0.7)
heatmap(cor(assay(vsd)), symm = TRUE, main = paste(accession, "sample correlation"))
dev.off()
result <- as.data.frame(results(dds, contrast = c("condition", "disease", "control")))
result$gene_id <- rownames(result)
result$gene_symbol <- NA_character_
result$effect <- result$log2FoldChange
result$se <- result$lfcSE
result$z <- result$stat
result$p_value <- result$pvalue
result$fdr <- result$padj
result$accession <- accession
result$species <- "Mus_musculus"
result$tissue <- "spinal_cord"
result$role <- "mouse_neural_validation"
result$eligibility <- "eligible"
result$contrast <- "disease-control"
result$effect_unit <- "log2_fold_change"
columns <- c("accession", "gene_id", "gene_symbol", "effect", "se", "z", "p_value", "fdr", "species", "tissue", "role", "eligibility", "contrast", "effect_unit")
dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
write.csv(result[, columns], output, row.names = FALSE)
writeLines(capture.output(sessionInfo()), paste0(output, ".sessionInfo.txt"))
