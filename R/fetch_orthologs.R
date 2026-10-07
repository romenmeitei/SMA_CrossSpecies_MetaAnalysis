#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
output <- if (length(args) >= 1) args[[1]] else "resources/ensembl_human_mouse_one_to_one.tsv"
release <- if (length(args) >= 2) as.integer(args[[2]]) else 113L
if (!requireNamespace("biomaRt", quietly = TRUE)) stop("Install biomaRt with environment/install_bioc.R")

human <- biomaRt::useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl", version = release)
mapping <- biomaRt::getBM(
  attributes = c(
    "ensembl_gene_id", "external_gene_name",
    "mmusculus_homolog_ensembl_gene", "mmusculus_homolog_associated_gene_name",
    "mmusculus_homolog_orthology_type"
  ),
  mart = human,
  uniqueRows = TRUE
)
names(mapping) <- c("human_ensembl_gene_id", "human_symbol", "mouse_ensembl_gene_id", "mouse_symbol", "orthology_type")
mapping <- mapping[mapping$orthology_type == "ortholog_one2one", , drop = FALSE]
mapping$ensembl_release <- release
dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
write.table(mapping, output, sep = "\t", row.names = FALSE, quote = FALSE)
message("Wrote ", nrow(mapping), " one-to-one orthologs to ", output)
