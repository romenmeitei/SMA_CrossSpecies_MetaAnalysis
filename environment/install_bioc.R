options(repos = c(CRAN = "https://cloud.r-project.org"))
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
cran <- c("optparse", "metafor", "data.table", "readxl", "jsonlite")
missing_cran <- cran[!vapply(cran, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_cran)) install.packages(missing_cran)
bioc <- c("limma", "edgeR", "DESeq2", "GEOquery", "biomaRt", "sva", "org.Hs.eg.db", "org.Mm.eg.db")
missing_bioc <- bioc[!vapply(bioc, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_bioc)) BiocManager::install(missing_bioc, ask = FALSE, update = FALSE)
writeLines(capture.output(sessionInfo()), "results/qc/R_sessionInfo.txt")

