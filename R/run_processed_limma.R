#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(limma))
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 7) stop("Usage: run_processed_limma.R MATRIX_TSV SAMPLE_MAP ACCESSION SPECIES TISSUE ELIGIBILITY OUTPUT_CSV")
matrix_path <- args[[1]]
sample_map_path <- args[[2]]
accession <- args[[3]]
species <- args[[4]]
tissue <- args[[5]]
eligibility <- args[[6]]
output <- args[[7]]

input <- read.delim(matrix_path, check.names = FALSE, stringsAsFactors = FALSE)
if (!"gene_id" %in% names(input)) stop("Matrix requires a gene_id column")
gene_symbol <- if ("gene_symbol" %in% names(input)) input$gene_symbol else rep(NA_character_, nrow(input))
metadata <- read.csv(sample_map_path, stringsAsFactors = FALSE)
missing <- setdiff(metadata$sample_id, names(input))
if (length(missing)) stop("Matrix missing samples: ", paste(missing, collapse = ", "))
expression <- as.matrix(input[, metadata$sample_id, drop = FALSE])
rownames(expression) <- input$gene_id
storage.mode(expression) <- "double"
if (unname(quantile(expression, 0.99, na.rm = TRUE)) > 100) expression <- log2(expression + 1)
metadata$condition <- factor(metadata$condition, levels = c("control", "disease"))
qc_dir <- file.path(dirname(dirname(output)), "qc")
dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)
pdf(file.path(qc_dir, paste0(accession, "_QC.pdf")), width = 8, height = 6)
boxplot(expression, las = 2, outline = FALSE, main = paste(accession, "expression distributions"))
pca_expression <- expression[rowSums(is.finite(expression)) >= 2, , drop = FALSE]
for (i in seq_len(nrow(pca_expression))) pca_expression[i, !is.finite(pca_expression[i, ])] <- mean(pca_expression[i, ], na.rm = TRUE)
pca <- prcomp(t(pca_expression), scale. = FALSE)
colors <- c(control = "#26619C", disease = "#B22222")
plot(pca$x[, 1], pca$x[, 2], col = colors[as.character(metadata$condition)], pch = 19,
     xlab = "PC1", ylab = "PC2", main = paste(accession, "PCA"))
text(pca$x[, 1], pca$x[, 2], labels = metadata$sample_id, pos = 3, cex = 0.7)
heatmap(cor(expression, use = "pairwise.complete.obs"), symm = TRUE, main = paste(accession, "sample correlation"))
dev.off()
design <- model.matrix(~0 + condition, metadata)
colnames(design) <- levels(metadata$condition)
if ("biological_unit" %in% names(metadata) && anyDuplicated(metadata$biological_unit)) {
  correlation <- duplicateCorrelation(expression, design, block = metadata$biological_unit)
  fit <- lmFit(expression, design, block = metadata$biological_unit, correlation = correlation$consensus)
} else {
  fit <- lmFit(expression, design)
}
contrast <- makeContrasts(`disease-control` = disease - control, levels = design)
fit <- eBayes(contrasts.fit(fit, contrast), robust = TRUE, trend = TRUE)
table <- topTable(fit, number = Inf, sort.by = "none")
table$gene_id <- rownames(table)
table$gene_symbol <- gene_symbol[match(table$gene_id, input$gene_id)]
table$effect <- table$logFC
table$z <- table$t
table$se <- abs(ifelse(table$t == 0, NA_real_, table$logFC / table$t))
table$p_value <- table$P.Value
table$fdr <- table$adj.P.Val
table$accession <- accession
table$species <- species
table$tissue <- tissue
table$role <- if (species == "Homo_sapiens") "human_neural_replication" else "mouse_neural_discovery"
table$eligibility <- eligibility
table$contrast <- "disease-control"
table$effect_unit <- "log2_fold_change"
columns <- c("accession", "gene_id", "gene_symbol", "effect", "se", "z", "p_value", "fdr", "species", "tissue", "role", "eligibility", "contrast", "effect_unit")
dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
write.csv(table[, columns], output, row.names = FALSE)
writeLines(capture.output(sessionInfo()), paste0(output, ".sessionInfo.txt"))
