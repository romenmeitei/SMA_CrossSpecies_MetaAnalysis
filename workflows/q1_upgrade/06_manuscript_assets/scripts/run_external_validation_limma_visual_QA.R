args <- commandArgs(trailingOnly=TRUE)

if (length(args) != 4) {
  stop(
    "Usage: Rscript run_external_validation_limma.R ",
    "<prepared_root> <audit_root> ",
    "<result_root> <figure_root>"
  )
}

prepared_root <- normalizePath(args[[1]])
audit_root <- normalizePath(args[[2]])
result_root <- normalizePath(args[[3]])
figure_root <- normalizePath(args[[4]])

dir.create(
  result_root,
  recursive=TRUE,
  showWarnings=FALSE
)

dir.create(
  figure_root,
  recursive=TRUE,
  showWarnings=FALSE
)

if (!requireNamespace("limma", quietly=TRUE)) {
  stop("The R package 'limma' is required.")
}

if (!requireNamespace("statmod", quietly=TRUE)) {
  stop(
    "The R package 'statmod' is required for ",
    "robust empirical Bayes estimation."
  )
}

suppressPackageStartupMessages(
  library(limma)
)

`%+%` <- function(left, right) {
  paste0(left, right)
}


extract_result <- function(
  fit,
  coefficient,
  accession,
  contrast_label,
  model_label,
  effect_unit
) {
  table <- topTable(
    fit,
    coef=coefficient,
    number=Inf,
    sort.by="none"
  )

  coefficient_index <- if (
    is.character(coefficient)
  ) {
    match(
      coefficient,
      colnames(fit$coefficients)
    )
  } else {
    coefficient
  }

  if (is.na(coefficient_index)) {
    stop(
      "Coefficient not found: ",
      coefficient
    )
  }

  moderated_se <- (
    fit$stdev.unscaled[, coefficient_index]
    * sqrt(fit$s2.post)
  )

  output <- data.frame(
    accession=accession,
    gene_symbol=rownames(table),
    contrast=contrast_label,
    model=model_label,
    effect=table$logFC,
    se=moderated_se,
    moderated_t=table$t,
    p_value=table$P.Value,
    fdr=p.adjust(
      table$P.Value,
      method="BH"
    ),
    average_expression=table$AveExpr,
    effect_unit=effect_unit,
    stringsAsFactors=FALSE
  )

  output <- output[
    order(
      output$p_value,
      output$gene_symbol
    ),
    ,
    drop=FALSE
  ]

  rownames(output) <- NULL
  output
}


write_pca <- function(
  expression,
  annotations,
  accession,
  output_prefix
) {
  annotation_order <- match(
    colnames(expression),
    annotations$sample_id
  )

  if (any(is.na(annotation_order))) {
    stop(
      accession,
      ": PCA annotation mismatch."
    )
  }

  annotation <- annotations[
    annotation_order,
    ,
    drop=FALSE
  ]

  pca <- prcomp(
    t(expression),
    center=TRUE,
    scale.=FALSE
  )

  variance <- (
    100
    * pca$sdev^2
    / sum(pca$sdev^2)
  )

  coordinates <- data.frame(
    sample_id=rownames(pca$x),
    condition=annotation$condition,
    PC1=pca$x[, 1],
    PC2=pca$x[, 2],
    PC1_variance_percent=variance[[1]],
    PC2_variance_percent=variance[[2]],
    stringsAsFactors=FALSE
  )

  write.csv(
    coordinates,
    paste0(output_prefix, "_PCA_coordinates.csv"),
    row.names=FALSE
  )

  conditions <- unique(
    as.character(annotation$condition)
  )

  palette <- c(
    control="#2166AC",
    disease="#B2182B",
    corrected="#1B7837"
  )

  colors <- palette[
    as.character(annotation$condition)
  ]

  draw_plot <- function() {
    finite_x <- coordinates$PC1[
      is.finite(coordinates$PC1)
    ]

    finite_y <- coordinates$PC2[
      is.finite(coordinates$PC2)
    ]

    x_range <- range(finite_x)
    y_range <- range(finite_y)

    x_padding <- max(
      diff(x_range) * 0.14,
      1
    )

    y_padding <- max(
      diff(y_range) * 0.14,
      1
    )

    x_limits <- c(
      x_range[[1]] - x_padding,
      x_range[[2]] + x_padding
    )

    y_limits <- c(
      y_range[[1]] - y_padding,
      y_range[[2]] + y_padding
    )

    x_midpoint <- mean(x_range)

    label_positions <- ifelse(
      coordinates$PC1 >= x_midpoint,
      2,
      4
    )

    old_par <- par(
      no.readonly=TRUE
    )

    on.exit(
      par(old_par),
      add=TRUE
    )

    par(
      mar=c(5.2, 5.2, 4.5, 9.0),
      xpd=NA
    )

    plot(
      coordinates$PC1,
      coordinates$PC2,
      pch=19,
      cex=1.5,
      col=colors,
      xlim=x_limits,
      ylim=y_limits,
      xlab=sprintf(
        "PC1 (%.1f%%)",
        variance[[1]]
      ),
      ylab=sprintf(
        "PC2 (%.1f%%)",
        variance[[2]]
      ),
      main=paste0(
        accession,
        " external-validation PCA"
      )
    )

    text(
      coordinates$PC1,
      coordinates$PC2,
      labels=coordinates$sample_id,
      pos=label_positions,
      offset=0.55,
      cex=0.72
    )

    legend(
      "topright",
      inset=c(-0.34, 0),
      legend=conditions,
      col=palette[conditions],
      pch=19,
      bty="n",
      xpd=NA
    )
  }

  png(
    paste0(output_prefix, "_PCA.png"),
    width=1800,
    height=1500,
    res=220
  )
  draw_plot()
  dev.off()

  pdf(
    paste0(output_prefix, "_PCA.pdf"),
    width=8,
    height=7
  )
  draw_plot()
  dev.off()

  coordinates
}


summarize_result <- function(result) {
  data.frame(
    accession=result$accession[[1]],
    contrast=result$contrast[[1]],
    model=result$model[[1]],
    tested_genes=nrow(result),
    finite_effects=sum(
      is.finite(result$effect)
    ),
    nominal_p_le_0.05=sum(
      result$p_value <= 0.05,
      na.rm=TRUE
    ),
    fdr_le_0.05=sum(
      result$fdr <= 0.05,
      na.rm=TRUE
    ),
    minimum_p=min(
      result$p_value,
      na.rm=TRUE
    ),
    minimum_fdr=min(
      result$fdr,
      na.rm=TRUE
    ),
    stringsAsFactors=FALSE
  )
}


# ==================================================================
# GSE51735: in-vivo mouse motor neurons, n=2 versus n=2
# ==================================================================

rpkm_51735 <- read.delim(
  file.path(
    prepared_root,
    "GSE51735_motor_neuron_RPKM.tsv"
  ),
  row.names=1,
  check.names=FALSE,
  stringsAsFactors=FALSE
)

rpkm_51735 <- as.matrix(rpkm_51735)
storage.mode(rpkm_51735) <- "double"

annotation_51735 <- read.csv(
  file.path(
    audit_root,
    "GSE51735_sample_map.csv"
  ),
  stringsAsFactors=FALSE
)

annotation_51735 <- annotation_51735[
  tolower(
    as.character(
      annotation_51735$
        include_primary_validation
    )
  ) == "true",
  ,
  drop=FALSE
]

sample_order_51735 <- match(
  colnames(rpkm_51735),
  annotation_51735$sample_id
)

if (any(is.na(sample_order_51735))) {
  stop("GSE51735 annotation mismatch.")
}

annotation_51735 <- annotation_51735[
  sample_order_51735,
  ,
  drop=FALSE
]

keep_51735 <- (
  rowSums(rpkm_51735 >= 1) >= 2
)

expression_51735 <- log2(
  rpkm_51735[keep_51735, , drop=FALSE]
  + 0.5
)

condition_51735 <- factor(
  annotation_51735$condition,
  levels=c(
    "control",
    "disease"
  )
)

design_51735 <- model.matrix(
  ~0 + condition_51735
)

colnames(design_51735) <- levels(
  condition_51735
)

contrast_51735 <- makeContrasts(
  disease_control=disease-control,
  levels=design_51735
)

# Main deposited-expression model.
fit_51735 <- lmFit(
  expression_51735,
  design_51735
)

fit_51735 <- contrasts.fit(
  fit_51735,
  contrast_51735
)

fit_51735 <- eBayes(
  fit_51735,
  trend=TRUE,
  robust=TRUE
)

result_51735 <- extract_result(
  fit=fit_51735,
  coefficient="disease_control",
  accession="GSE51735",
  contrast_label="disease-control",
  model_label="deposited_log2_RPKM",
  effect_unit="log2_RPKM_difference"
)

write.csv(
  result_51735,
  file.path(
    result_root,
    "GSE51735_motor_neuron_"
    %+%
    "disease-control_limma.csv"
  ),
  row.names=FALSE
)

# Quantile-normalization sensitivity model.
expression_51735_quantile <- (
  normalizeBetweenArrays(
    expression_51735,
    method="quantile"
  )
)

fit_51735_quantile <- lmFit(
  expression_51735_quantile,
  design_51735
)

fit_51735_quantile <- contrasts.fit(
  fit_51735_quantile,
  contrast_51735
)

fit_51735_quantile <- eBayes(
  fit_51735_quantile,
  trend=TRUE,
  robust=TRUE
)

result_51735_quantile <- extract_result(
  fit=fit_51735_quantile,
  coefficient="disease_control",
  accession="GSE51735",
  contrast_label="disease-control",
  model_label="quantile_normalized_sensitivity",
  effect_unit="normalized_log2_RPKM_difference"
)

write.csv(
  result_51735_quantile,
  file.path(
    result_root,
    "GSE51735_motor_neuron_"
    %+%
    "disease-control_"
    %+%
    "quantile_sensitivity.csv"
  ),
  row.names=FALSE
)

pca_51735 <- write_pca(
  expression=expression_51735,
  annotations=annotation_51735,
  accession="GSE51735",
  output_prefix=file.path(
    figure_root,
    "GSE51735_motor_neuron"
  )
)


# ==================================================================
# GSE27205: human iPSC-derived motor neurons, three groups
# ==================================================================

expression_27205 <- read.delim(
  file.path(
    prepared_root,
    "GSE27205_gene_expression_"
    %+%
    "untransformed.tsv"
  ),
  row.names=1,
  check.names=FALSE,
  stringsAsFactors=FALSE
)

expression_27205 <- as.matrix(
  expression_27205
)

storage.mode(expression_27205) <- "double"

# Matrix values are already log2-scaled.
if (
  max(expression_27205, na.rm=TRUE) >= 100
) {
  stop(
    "GSE27205 no longer appears log2-scaled; ",
    "do not continue."
  )
}

annotation_27205 <- read.csv(
  file.path(
    audit_root,
    "GSE27205_sample_map.csv"
  ),
  stringsAsFactors=FALSE
)

sample_order_27205 <- match(
  colnames(expression_27205),
  annotation_27205$sample_id
)

if (any(is.na(sample_order_27205))) {
  stop("GSE27205 annotation mismatch.")
}

annotation_27205 <- annotation_27205[
  sample_order_27205,
  ,
  drop=FALSE
]

condition_27205 <- factor(
  annotation_27205$condition,
  levels=c(
    "control",
    "disease",
    "corrected"
  )
)

design_27205 <- model.matrix(
  ~0 + condition_27205
)

colnames(design_27205) <- levels(
  condition_27205
)

contrasts_27205 <- makeContrasts(
  disease_control=disease-control,
  corrected_disease=corrected-disease,
  levels=design_27205
)

fit_27205 <- lmFit(
  expression_27205,
  design_27205
)

fit_27205 <- contrasts.fit(
  fit_27205,
  contrasts_27205
)

fit_27205 <- eBayes(
  fit_27205,
  trend=TRUE,
  robust=TRUE
)

result_27205_disease <- extract_result(
  fit=fit_27205,
  coefficient="disease_control",
  accession="GSE27205",
  contrast_label="disease-control",
  model_label="deposited_log2_expression",
  effect_unit="log2_expression_difference"
)

result_27205_correction <- extract_result(
  fit=fit_27205,
  coefficient="corrected_disease",
  accession="GSE27205",
  contrast_label="corrected-disease",
  model_label="deposited_log2_expression",
  effect_unit="log2_expression_difference"
)

write.csv(
  result_27205_disease,
  file.path(
    result_root,
    "GSE27205_disease-control_limma.csv"
  ),
  row.names=FALSE
)

write.csv(
  result_27205_correction,
  file.path(
    result_root,
    "GSE27205_corrected-disease_limma.csv"
  ),
  row.names=FALSE
)

pca_27205 <- write_pca(
  expression=expression_27205,
  annotations=annotation_27205,
  accession="GSE27205",
  output_prefix=file.path(
    figure_root,
    "GSE27205_motor_neuron"
  )
)


# ==================================================================
# Cross-model summaries
# ==================================================================

summary_table <- rbind(
  summarize_result(result_51735),
  summarize_result(
    result_51735_quantile
  ),
  summarize_result(
    result_27205_disease
  ),
  summarize_result(
    result_27205_correction
  )
)

write.csv(
  summary_table,
  file.path(
    audit_root,
    "external_validation_limma_summary.csv"
  ),
  row.names=FALSE
)

priority_symbols <- c(
  "FAM107A",
  "UGP2"
)

priority_results <- rbind(
  result_51735[
    toupper(result_51735$gene_symbol)
      %in% priority_symbols,
    ,
    drop=FALSE
  ],
  result_51735_quantile[
    toupper(
      result_51735_quantile$gene_symbol
    ) %in% priority_symbols,
    ,
    drop=FALSE
  ],
  result_27205_disease[
    toupper(
      result_27205_disease$gene_symbol
    ) %in% priority_symbols,
    ,
    drop=FALSE
  ],
  result_27205_correction[
    toupper(
      result_27205_correction$gene_symbol
    ) %in% priority_symbols,
    ,
    drop=FALSE
  ]
)

priority_results$normalized_symbol <- (
  toupper(priority_results$gene_symbol)
)

priority_results <- priority_results[
  order(
    priority_results$normalized_symbol,
    priority_results$accession,
    priority_results$contrast,
    priority_results$model
  ),
  ,
  drop=FALSE
]

write.csv(
  priority_results,
  file.path(
    result_root,
    "FAM107A_UGP2_external_"
    %+%
    "validation_results.csv"
  ),
  row.names=FALSE
)

capture.output(
  sessionInfo(),
  file=file.path(
    audit_root,
    "external_validation_R_sessionInfo.txt"
  )
)

print(
  paste(rep("=", 100), collapse="")
)
print("EXTERNAL-VALIDATION MODEL SUMMARY")
print(
  paste(rep("=", 100), collapse="")
)
print(summary_table)

print("")
print(
  paste(rep("=", 100), collapse="")
)
print("FAM107A AND UGP2")
print(
  paste(rep("=", 100), collapse="")
)

print(
  priority_results[
    ,
    c(
      "accession",
      "gene_symbol",
      "contrast",
      "model",
      "effect",
      "moderated_t",
      "p_value",
      "fdr"
    )
  ]
)

print("")
cat(
  "
PASS: external-validation limma models completed.
"
)
