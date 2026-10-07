suppressPackageStartupMessages({
  library(tximport)
  library(DESeq2)
})

args <- commandArgs(trailingOnly=TRUE)

if (length(args) != 12) {
  stop(
    "Expected 12 arguments: quant_root tx2gene annotation ",
    "accession species tissue role eligibility output ",
    "min_count qc_dir status_file"
  )
}

quant_root     <- args[[1]]
tx2gene_path   <- args[[2]]
annotation_path <- args[[3]]
accession      <- args[[4]]
species        <- args[[5]]
tissue         <- args[[6]]
role           <- args[[7]]
eligibility    <- args[[8]]
output         <- args[[9]]
minimum_count  <- as.integer(args[[10]])
qc_dir         <- args[[11]]
status_file    <- args[[12]]

dir.create(
  dirname(output),
  recursive=TRUE,
  showWarnings=FALSE
)

dir.create(
  qc_dir,
  recursive=TRUE,
  showWarnings=FALSE
)

ann <- read.csv(
  annotation_path,
  stringsAsFactors=FALSE,
  check.names=FALSE
)

required <- c(
  "sample_id",
  "condition",
  "biological_unit"
)

if (!all(required %in% names(ann))) {
  stop(
    "Annotation columns missing: ",
    paste(
      setdiff(required, names(ann)),
      collapse=", "
    )
  )
}

if (anyDuplicated(ann$sample_id)) {
  stop("Duplicate sample_id values in annotation")
}

if (!all(
  ann$condition %in% c("control", "disease")
)) {
  stop("Unexpected condition values")
}

unit_condition_count <- tapply(
  ann$condition,
  ann$biological_unit,
  function(x) length(unique(x))
)

if (any(unit_condition_count != 1)) {
  stop(
    "A biological unit was assigned to multiple conditions"
  )
}

files <- file.path(
  quant_root,
  accession,
  ann$sample_id,
  "quant.sf"
)

names(files) <- ann$sample_id

if (any(!file.exists(files))) {
  stop(
    "Missing quant.sf: ",
    paste(
      files[!file.exists(files)],
      collapse=", "
    )
  )
}

tx2gene_full <- read.delim(
  tx2gene_path,
  stringsAsFactors=FALSE,
  check.names=FALSE
)

required_tx <- c(
  "transcript_id",
  "gene_id",
  "gene_symbol"
)

if (!all(required_tx %in% names(tx2gene_full))) {
  stop(
    "tx2gene columns missing: ",
    paste(
      setdiff(required_tx, names(tx2gene_full)),
      collapse=", "
    )
  )
}

tx2gene <- unique(
  tx2gene_full[
    ,
    c("transcript_id", "gene_id")
  ]
)

message(
  accession,
  ": importing ",
  length(files),
  " Salmon quantifications"
)

txi <- tximport(
  files,
  type="salmon",
  tx2gene=tx2gene,
  countsFromAbundance="no",
  ignoreTxVersion=FALSE,
  dropInfReps=TRUE
)

rownames(ann) <- ann$sample_id

ann$condition <- relevel(
  factor(ann$condition),
  ref="control"
)

dds <- DESeqDataSetFromTximport(
  txi,
  colData=ann,
  design=~condition
)

# Collapse sequencing/technical replicates belonging
# to the same independent biological unit.
dds <- collapseReplicates(
  dds,
  groupby=as.character(
    dds$biological_unit
  ),
  run=rownames(colData(dds)),
  renameCols=TRUE
)

dds$condition <- relevel(
  factor(dds$condition),
  ref="control"
)

sample_table <- as.data.frame(
  colData(dds)
)

sample_table$biological_unit <- rownames(
  sample_table
)

sample_table$collapsed_run_count <- vapply(
  rownames(sample_table),
  function(unit) {
    sum(ann$biological_unit == unit)
  },
  integer(1)
)

write.csv(
  sample_table,
  file.path(
    qc_dir,
    paste0(
      accession,
      "_collapsed_sample_table.csv"
    )
  ),
  row.names=FALSE
)

group_sizes <- table(dds$condition)

if (
  !all(
    c("control", "disease")
    %in% names(group_sizes)
  )
) {
  stop(
    "Both control and disease groups are required"
  )
}

minimum_group_n <- min(group_sizes)

keep <- (
  rowSums(
    counts(dds) >= minimum_count
  )
  >= minimum_group_n
)

filter_audit <- data.frame(
  accession=accession,
  total_genes_before=nrow(dds),
  retained_genes=sum(keep),
  removed_genes=sum(!keep),
  minimum_count=minimum_count,
  minimum_group_n=minimum_group_n,
  stringsAsFactors=FALSE
)

write.csv(
  filter_audit,
  file.path(
    qc_dir,
    paste0(
      accession,
      "_gene_filter_audit.csv"
    )
  ),
  row.names=FALSE
)

dds <- dds[keep, ]

if (nrow(dds) == 0) {
  stop("No genes remained after count filtering")
}

gene_map <- tx2gene_full[
  !duplicated(tx2gene_full$gene_id),
  c("gene_id", "gene_symbol")
]

if (minimum_group_n < 2) {
  # One biological unit per condition after collapse:
  # report descriptive normalized log2 ratios only.
  dds <- estimateSizeFactors(dds)

  normalized <- counts(
    dds,
    normalized=TRUE
  )

  control_cols <- which(
    dds$condition == "control"
  )

  disease_cols <- which(
    dds$condition == "disease"
  )

  effect <- log2(
    (
      rowMeans(
        normalized[
          ,
          disease_cols,
          drop=FALSE
        ]
      ) + 0.5
    )
    /
    (
      rowMeans(
        normalized[
          ,
          control_cols,
          drop=FALSE
        ]
      ) + 0.5
    )
  )

  out <- data.frame(
    accession=accession,
    gene_id=rownames(dds),
    effect=effect,
    se=NA_real_,
    z=NA_real_,
    p_value=NA_real_,
    fdr=NA_real_,
    species=species,
    tissue=tissue,
    role=role,
    eligibility="descriptive_only_confounded",
    contrast="disease-control",
    effect_unit=(
      "DESeq2_normalized_descriptive_log2_ratio"
    ),
    base_mean=rowMeans(normalized),
    stringsAsFactors=FALSE
  )

  out <- merge(
    out,
    gene_map,
    by="gene_id",
    all.x=TRUE,
    sort=FALSE
  )

  out <- out[
    order(
      -abs(out$effect),
      na.last=TRUE
    ),
  ]

  write.csv(
    out,
    output,
    row.names=FALSE
  )

  write.csv(
    as.data.frame(normalized),
    file.path(
      qc_dir,
      paste0(
        accession,
        "_normalized_counts.csv"
      )
    )
  )

  writeLines(
    c(
      "status=DESCRIPTIVE_ONLY",
      "reason=one biological unit per condition",
      paste0(
        "control_biological_units=",
        group_sizes[["control"]]
      ),
      paste0(
        "disease_biological_units=",
        group_sizes[["disease"]]
      ),
      paste0(
        "reported_genes=",
        nrow(out)
      )
    ),
    status_file
  )

  capture.output(
    sessionInfo(),
    file=file.path(
      qc_dir,
      paste0(
        accession,
        "_R_sessionInfo.txt"
      )
    )
  )

  quit(
    save="no",
    status=0
  )
}

dds <- DESeq(
  dds,
  test="Wald",
  quiet=TRUE
)

res <- results(
  dds,
  contrast=c(
    "condition",
    "disease",
    "control"
  ),
  alpha=0.05,
  independentFiltering=TRUE
)

out <- data.frame(
  accession=accession,
  gene_id=rownames(res),
  effect=res$log2FoldChange,
  se=res$lfcSE,
  z=res$stat,
  p_value=res$pvalue,
  fdr=res$padj,
  species=species,
  tissue=tissue,
  role=role,
  eligibility=eligibility,
  contrast="disease-control",
  effect_unit="DESeq2_log2_fold_change",
  base_mean=res$baseMean,
  stringsAsFactors=FALSE
)

out <- merge(
  out,
  gene_map,
  by="gene_id",
  all.x=TRUE,
  sort=FALSE
)

out <- out[
  order(
    out$fdr,
    out$p_value,
    na.last=TRUE
  ),
]

write.csv(
  out,
  output,
  row.names=FALSE
)

normalized <- counts(
  dds,
  normalized=TRUE
)

write.csv(
  as.data.frame(normalized),
  file.path(
    qc_dir,
    paste0(
      accession,
      "_normalized_counts.csv"
    )
  )
)

png(
  file.path(
    qc_dir,
    paste0(
      accession,
      "_MA_plot.png"
    )
  ),
  width=1800,
  height=1500,
  res=200
)

plotMA(
  res,
  ylim=c(-6, 6),
  alpha=0.05
)

dev.off()

if (
  ncol(dds) >= 3
  && nrow(dds) >= 100
) {
  transformed <- tryCatch(
    {
      vst(
        dds,
        blind=FALSE,
        nsub=min(
          1000,
          nrow(dds)
        )
      )
    },
    error=function(e) {
      varianceStabilizingTransformation(
        dds,
        blind=FALSE
      )
    }
  )

  pca <- plotPCA(
    transformed,
    intgroup="condition",
    returnData=TRUE
  )

  percentVar <- round(
    100 * attr(
      pca,
      "percentVar"
    )
  )

  png(
    file.path(
      qc_dir,
      paste0(
        accession,
        "_PCA.png"
      )
    ),
    width=1800,
    height=1500,
    res=200
  )

  condition_factor <- factor(
    pca$condition
  )

  plot(
    pca$PC1,
    pca$PC2,
    col=as.integer(condition_factor),
    pch=19,
    xlab=paste0(
      "PC1: ",
      percentVar[1],
      "% variance"
    ),
    ylab=paste0(
      "PC2: ",
      percentVar[2],
      "% variance"
    ),
    main=accession
  )

  text(
    pca$PC1,
    pca$PC2,
    labels=rownames(pca),
    pos=3,
    cex=0.7
  )

  legend(
    "topright",
    legend=levels(condition_factor),
    col=seq_along(
      levels(condition_factor)
    ),
    pch=19
  )

  dev.off()
}

writeLines(
  c(
    "status=INFERENTIAL_DESEQ2",
    paste0(
      "control_biological_units=",
      group_sizes[["control"]]
    ),
    paste0(
      "disease_biological_units=",
      group_sizes[["disease"]]
    ),
    paste0(
      "tested_genes=",
      nrow(out)
    ),
    paste0(
      "fdr_le_0.05=",
      sum(
        out$fdr <= 0.05,
        na.rm=TRUE
      )
    )
  ),
  status_file
)

capture.output(
  sessionInfo(),
  file=file.path(
    qc_dir,
    paste0(
      accession,
      "_R_sessionInfo.txt"
    )
  )
)
