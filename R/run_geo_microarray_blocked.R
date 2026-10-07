#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(GEOquery)
  library(limma)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) stop("Usage: run_geo_microarray_blocked.R ACCESSION SAMPLE_MAP OUTPUT_DIR")
accession <- args[[1]]
sample_map_path <- args[[2]]
output_dir <- args[[3]]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

sample_map <- read.csv(sample_map_path, stringsAsFactors = FALSE)
required <- c("sample_id", "condition", "mouse_id")
if (!all(required %in% names(sample_map))) stop("Sample map requires: ", paste(required, collapse = ", "))

sets <- getGEO(accession, GSEMatrix = TRUE, AnnotGPL = TRUE)
if (!is.list(sets)) sets <- list(sets)
eset <- sets[[which.max(vapply(sets, ncol, integer(1)))]]
expression <- exprs(eset)
sample_ids <- sub("\\..*$", "", colnames(expression))
match_index <- match(sample_map$sample_id, sample_ids)
if (anyNA(match_index)) stop("Series matrix is missing samples: ", paste(sample_map$sample_id[is.na(match_index)], collapse = ", "))
expression <- expression[, match_index, drop = FALSE]
colnames(expression) <- sample_map$sample_id

q99 <- unname(quantile(expression, 0.99, na.rm = TRUE))
if (q99 > 100) {
  # Retain zero-valued microarray observations and avoid introducing
  # missing values before model fitting.
  expression <- log2(pmax(expression, 0) + 1)
}

# Remove only probes that still contain non-finite values or have
# no expression variance across the selected samples.
finite_probe <- apply(
  expression,
  1,
  function(values) all(is.finite(values))
)
expression <- expression[
  finite_probe,
  ,
  drop = FALSE
]

variable_probe <- apply(
  expression,
  1,
  function(values) sd(values) > 0
)
expression <- expression[
  variable_probe,
  ,
  drop = FALSE
]

if (nrow(expression) == 0) {
  stop("No finite, variable probes remained after preprocessing")
}

expression <- normalizeBetweenArrays(
  expression,
  method = "quantile"
)

qc_dir <- file.path(dirname(output_dir), "qc")
dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)
pdf(file.path(qc_dir, paste0(accession, "_QC.pdf")), width = 9, height = 7)
boxplot(expression, las = 2, outline = FALSE, main = paste(accession, "normalized distributions"), ylab = "log2 expression")
pca_expression <- expression[rowSums(is.finite(expression)) >= 2, , drop = FALSE]
for (i in seq_len(nrow(pca_expression))) pca_expression[i, !is.finite(pca_expression[i, ])] <- mean(pca_expression[i, ], na.rm = TRUE)
pca <- prcomp(t(pca_expression), scale. = FALSE)
condition_colors <- c(control = "#26619C", disease = "#B22222", treated = "#2E8B57")
plot(pca$x[, 1], pca$x[, 2], col = condition_colors[as.character(sample_map$condition)], pch = 19,
     xlab = "PC1", ylab = "PC2", main = paste(accession, "sample PCA"))
text(pca$x[, 1], pca$x[, 2], labels = sample_map$mouse_id, pos = 3, cex = 0.65)
legend("topright", legend = names(condition_colors), col = condition_colors, pch = 19, bty = "n")
heatmap(cor(expression, use = "pairwise.complete.obs"), symm = TRUE, main = paste(accession, "sample correlation"))
dev.off()

annotation <- fData(eset)
annotation <- annotation[match(rownames(expression), rownames(annotation)), , drop = FALSE]
symbol_candidates <- c("Gene Symbol", "GENE_SYMBOL", "gene_assignment", "Symbol", "SYMBOL")
symbol_column <- symbol_candidates[symbol_candidates %in% names(annotation)][1]
if (is.na(symbol_column)) {
  warning("No gene-symbol annotation found; retaining probe identifiers")
  symbols <- rownames(expression)
} else {
  symbols <- as.character(annotation[[symbol_column]])
  symbols <- sub(" ///.*$", "", symbols)
  symbols <- sub(" //.*$", "", symbols)
  symbols <- trimws(symbols)
  symbols[is.na(symbols) | symbols == ""] <- rownames(expression)[is.na(symbols) | symbols == ""]
}

# Select one pre-model probe per gene by expression variability, never by p-value.
probe_iqr <- apply(expression, 1, IQR, na.rm = TRUE)
order_index <- order(symbols, -probe_iqr, rownames(expression))
keep <- order_index[!duplicated(symbols[order_index])]
expression <- expression[keep, , drop = FALSE]
symbols <- symbols[keep]

sample_map$condition <- factor(sample_map$condition, levels = c("control", "disease", "treated"))
design <- model.matrix(~0 + condition, sample_map)
colnames(design) <- levels(sample_map$condition)
corfit <- duplicateCorrelation(expression, design, block = sample_map$mouse_id)
fit <- lmFit(expression, design, block = sample_map$mouse_id, correlation = corfit$consensus)
contrasts <- makeContrasts(`disease-control` = disease - control, `treated-disease` = treated - disease, levels = design)
fit <- eBayes(contrasts.fit(fit, contrasts), robust = TRUE, trend = TRUE)

format_table <- function(fit_object, coefficient, suffix) {
  table <- topTable(fit_object, coef = coefficient, number = Inf, sort.by = "none")
  table$gene_id <- rownames(table)
  table$gene_symbol <- symbols[match(rownames(table), rownames(expression))]
  table$effect <- table$logFC
  table$z <- table$t
  table$se <- abs(ifelse(table$t == 0, NA_real_, table$logFC / table$t))
  table$p_value <- table$P.Value
  table$fdr <- table$adj.P.Val
  table$accession <- accession
  table$species <- "Mus_musculus"
  table$tissue <- "skeletal_muscle"
  table$role <- "mouse_muscle_treatment_reversal"
  table$eligibility <- "eligible"
  table$contrast <- coefficient
  table$effect_unit <- "log2_fold_change"
  columns <- c("accession", "gene_id", "gene_symbol", "effect", "se", "z", "p_value", "fdr", "species", "tissue", "role", "eligibility", "contrast", "effect_unit")
  write.csv(table[, columns], file.path(output_dir, paste0(accession, "_", suffix, "_primary.csv")), row.names = FALSE)
}
format_table(fit, "disease-control", "disease-control")
format_table(fit, "treated-disease", "treated-disease")

# Prespecified sensitivity: aggregate fibers to the mouse, then fit mice as samples.
mouse_ids <- unique(sample_map$mouse_id)
mouse_expression <- vapply(mouse_ids, function(id) rowMeans(expression[, sample_map$mouse_id == id, drop = FALSE], na.rm = TRUE), numeric(nrow(expression)))
rownames(mouse_expression) <- rownames(expression)
mouse_condition <- vapply(mouse_ids, function(id) as.character(unique(sample_map$condition[sample_map$mouse_id == id])[1]), character(1))
mouse_condition <- factor(mouse_condition, levels = levels(sample_map$condition))
mouse_design <- model.matrix(~0 + mouse_condition)
colnames(mouse_design) <- levels(mouse_condition)
mouse_fit <- eBayes(contrasts.fit(lmFit(mouse_expression, mouse_design), contrasts), robust = TRUE, trend = TRUE)

format_mouse <- function(fit_object, coefficient, suffix) {
  table <- topTable(fit_object, coef = coefficient, number = Inf, sort.by = "none")
  table$gene_id <- rownames(table)
  table$gene_symbol <- symbols[match(rownames(table), rownames(expression))]
  table$effect <- table$logFC
  table$z <- table$t
  table$se <- abs(ifelse(table$t == 0, NA_real_, table$logFC / table$t))
  table$p_value <- table$P.Value
  table$fdr <- table$adj.P.Val
  table$accession <- accession
  table$species <- "Mus_musculus"
  table$tissue <- "skeletal_muscle"
  table$role <- "mouse_muscle_treatment_reversal"
  table$eligibility <- "eligible"
  table$contrast <- coefficient
  table$effect_unit <- "log2_fold_change"
  columns <- c("accession", "gene_id", "gene_symbol", "effect", "se", "z", "p_value", "fdr", "species", "tissue", "role", "eligibility", "contrast", "effect_unit")
  write.csv(table[, columns], file.path(output_dir, paste0(accession, "_", suffix, "_mouse_level.csv")), row.names = FALSE)
}
format_mouse(mouse_fit, "disease-control", "disease-control")
format_mouse(mouse_fit, "treated-disease", "treated-disease")

writeLines(capture.output(list(log2_trigger_q99 = q99, duplicate_correlation = corfit$consensus, session = sessionInfo())), file.path(output_dir, paste0(accession, "_model_audit.txt")))
