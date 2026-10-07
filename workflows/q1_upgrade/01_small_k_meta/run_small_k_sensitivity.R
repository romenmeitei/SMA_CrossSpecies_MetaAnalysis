args <- commandArgs(trailingOnly=TRUE)

if (length(args) != 3) {
  stop(
    "Usage: Rscript run_small_k_sensitivity.R ",
    "<project_root> <input_csv> <output_root>"
  )
}

project_root <- normalizePath(
  args[[1]],
  mustWork=TRUE
)

input_csv <- normalizePath(
  args[[2]],
  mustWork=TRUE
)

output_root <- normalizePath(
  args[[3]],
  mustWork=TRUE
)

if (!requireNamespace("metafor", quietly=TRUE)) {
  stop("The R package 'metafor' is required.")
}

effects <- read.csv(
  input_csv,
  stringsAsFactors=FALSE,
  check.names=FALSE
)

required <- c(
  "ortholog_id",
  "gene_symbol",
  "accession",
  "species",
  "effect",
  "se"
)

missing <- setdiff(
  required,
  colnames(effects)
)

if (length(missing) > 0) {
  stop(
    "Study-effect table is missing: ",
    paste(missing, collapse=", ")
  )
}

effects$effect <- as.numeric(effects$effect)
effects$se <- as.numeric(effects$se)

effects <- effects[
  is.finite(effects$effect)
  & is.finite(effects$se)
  & effects$se > 0,
]

effects <- effects[
  !duplicated(
    effects[, c(
      "ortholog_id",
      "accession"
    )]
  ),
]

first_nonmissing <- function(values) {
  values <- as.character(values)

  values <- values[
    !is.na(values)
    & nzchar(values)
    & values != "NA"
  ]

  if (length(values) == 0) {
    return(NA_character_)
  }

  values[[1]]
}

fit_one_gene <- function(frame) {
  frame <- frame[
    order(frame$accession),
  ]

  k <- nrow(frame)
  n_species <- length(
    unique(frame$species)
  )

  if (k < 2 || n_species < 2) {
    return(list(
      row=NULL,
      error="Insufficient studies or species"
    ))
  }

  yi <- frame$effect
  sei <- frame$se

  result <- tryCatch(
    {
      fixed <- suppressWarnings(
        metafor::rma.uni(
          yi=yi,
          sei=sei,
          method="FE",
          test="z"
        )
      )

      reml_z <- suppressWarnings(
        metafor::rma.uni(
          yi=yi,
          sei=sei,
          method="REML",
          test="z",
          control=list(
            maxiter=1000,
            stepadj=0.5
          )
        )
      )

      hartung_knapp <- suppressWarnings(
        metafor::rma.uni(
          yi=yi,
          sei=sei,
          method="REML",
          test="knha",
          control=list(
            maxiter=1000,
            stepadj=0.5
          )
        )
      )

      estimate_reml <- as.numeric(
        reml_z$b
      )

      se_reml <- as.numeric(
        reml_z$se
      )

      se_hk <- as.numeric(
        hartung_knapp$se
      )

      # Modified Hartung–Knapp safeguard:
      # never permit an HK standard error smaller
      # than the conventional REML standard error.
      se_mkh <- max(
        se_reml,
        se_hk,
        na.rm=TRUE
      )

      df_mkh <- k - 1
      statistic_mkh <- (
        estimate_reml
        / se_mkh
      )

      p_mkh <- 2 * pt(
        abs(statistic_mkh),
        df=df_mkh,
        lower.tail=FALSE
      )

      critical_mkh <- qt(
        0.975,
        df=df_mkh
      )

      positive_fraction <- mean(
        yi > 0
      )

      negative_fraction <- mean(
        yi < 0
      )

      sign_concordance <- max(
        positive_fraction,
        negative_fraction
      )

      row <- data.frame(
        ortholog_id=as.character(
          frame$ortholog_id[[1]]
        ),
        gene_symbol=first_nonmissing(
          frame$gene_symbol
        ),
        k=k,
        n_species=n_species,
        datasets=paste(
          sort(unique(frame$accession)),
          collapse=";"
        ),
        species=paste(
          sort(unique(frame$species)),
          collapse=";"
        ),
        sign_concordance=sign_concordance,

        estimate_fixed=as.numeric(
          fixed$b
        ),
        se_fixed=as.numeric(
          fixed$se
        ),
        ci_low_fixed=as.numeric(
          fixed$ci.lb
        ),
        ci_high_fixed=as.numeric(
          fixed$ci.ub
        ),
        z_fixed=as.numeric(
          fixed$zval
        ),
        p_fixed=as.numeric(
          fixed$pval
        ),

        estimate_reml=estimate_reml,
        se_reml=se_reml,
        ci_low_reml=as.numeric(
          reml_z$ci.lb
        ),
        ci_high_reml=as.numeric(
          reml_z$ci.ub
        ),
        z_reml=as.numeric(
          reml_z$zval
        ),
        p_reml=as.numeric(
          reml_z$pval
        ),

        se_hk=se_hk,
        ci_low_hk=as.numeric(
          hartung_knapp$ci.lb
        ),
        ci_high_hk=as.numeric(
          hartung_knapp$ci.ub
        ),
        statistic_hk=as.numeric(
          hartung_knapp$zval
        ),
        p_hk=as.numeric(
          hartung_knapp$pval
        ),

        se_mkh=se_mkh,
        ci_low_mkh=(
          estimate_reml
          - critical_mkh * se_mkh
        ),
        ci_high_mkh=(
          estimate_reml
          + critical_mkh * se_mkh
        ),
        statistic_mkh=statistic_mkh,
        df_mkh=df_mkh,
        p_mkh=p_mkh,

        hk_scale_ratio=(
          se_hk / se_reml
        ),
        mkh_safeguard_used=(
          se_hk < se_reml
        ),

        tau2=as.numeric(
          reml_z$tau2
        ),
        i2=as.numeric(
          reml_z$I2
        ),
        q=as.numeric(
          reml_z$QE
        ),
        q_p_value=as.numeric(
          reml_z$QEp
        ),
        stringsAsFactors=FALSE
      )

      list(
        row=row,
        error=NA_character_
      )
    },
    error=function(error) {
      list(
        row=NULL,
        error=conditionMessage(error)
      )
    }
  )

  result
}

groups <- split(
  effects,
  effects$ortholog_id
)

detected_cores <- parallel::detectCores()

if (
  is.na(detected_cores)
  || detected_cores < 2
) {
  detected_cores <- 2
}

workers <- min(
  8L,
  detected_cores - 1L
)

cat(
  "Ortholog groups:",
  length(groups),
  "\n"
)

cat(
  "Parallel workers:",
  workers,
  "\n"
)

fits <- parallel::mclapply(
  groups,
  fit_one_gene,
  mc.cores=workers,
  mc.preschedule=TRUE
)

successful_rows <- lapply(
  fits,
  function(item) item$row
)

successful_rows <- successful_rows[
  !vapply(
    successful_rows,
    is.null,
    logical(1)
  )
]

results <- do.call(
  rbind,
  successful_rows
)

failures <- data.frame(
  ortholog_id=names(groups),
  error=vapply(
    fits,
    function(item) {
      if (is.null(item$error)) {
        return(NA_character_)
      }

      item$error
    },
    character(1)
  ),
  stringsAsFactors=FALSE
)

failures <- failures[
  !is.na(failures$error),
]

results$fdr_fixed <- p.adjust(
  results$p_fixed,
  method="BH"
)

results$fdr_reml <- p.adjust(
  results$p_reml,
  method="BH"
)

results$fdr_hk <- p.adjust(
  results$p_hk,
  method="BH"
)

results$fdr_mkh <- p.adjust(
  results$p_mkh,
  method="BH"
)

results$significant_existing_style <- (
  results$fdr_reml <= 0.05
)

results$significant_hk <- (
  results$fdr_hk <= 0.05
)

results$significant_mkh <- (
  results$fdr_mkh <= 0.05
)

results$significant_all_four <- (
  results$fdr_fixed <= 0.05
  & results$fdr_reml <= 0.05
  & results$fdr_hk <= 0.05
  & results$fdr_mkh <= 0.05
)

results <- results[
  order(
    results$fdr_mkh,
    results$p_mkh,
    results$ortholog_id
  ),
]

row.names(results) <- NULL

current_meta <- read.csv(
  file.path(
    project_root,
    "results",
    "meta",
    "gene_meta_neural_joint.csv"
  ),
  stringsAsFactors=FALSE,
  check.names=FALSE
)

current_columns <- current_meta[, c(
  "ortholog_id",
  "estimate",
  "se",
  "z",
  "p_value",
  "fdr",
  "tau2",
  "i2",
  "k"
)]

colnames(current_columns) <- c(
  "ortholog_id",
  "existing_estimate",
  "existing_se",
  "existing_z",
  "existing_p_value",
  "existing_fdr",
  "existing_tau2",
  "existing_i2",
  "existing_k"
)

comparison <- merge(
  results,
  current_columns,
  by="ortholog_id",
  all=TRUE
)

comparison$absolute_estimate_difference <- abs(
  comparison$estimate_reml
  - comparison$existing_estimate
)

comparison$absolute_se_difference <- abs(
  comparison$se_reml
  - comparison$existing_se
)

summarize_method <- function(
  method,
  p_values,
  fdr_values,
  k_values
) {
  data.frame(
    method=method,
    genes_tested=sum(
      is.finite(p_values)
    ),
    fdr_le_0.05=sum(
      fdr_values <= 0.05,
      na.rm=TRUE
    ),
    fdr_le_0.10=sum(
      fdr_values <= 0.10,
      na.rm=TRUE
    ),
    k2_fdr_le_0.05=sum(
      k_values == 2
      & fdr_values <= 0.05,
      na.rm=TRUE
    ),
    k3_fdr_le_0.05=sum(
      k_values == 3
      & fdr_values <= 0.05,
      na.rm=TRUE
    ),
    minimum_p=min(
      p_values,
      na.rm=TRUE
    ),
    minimum_fdr=min(
      fdr_values,
      na.rm=TRUE
    ),
    stringsAsFactors=FALSE
  )
}

method_summary <- rbind(
  summarize_method(
    "existing_Python_REML",
    current_meta$p_value,
    current_meta$fdr,
    current_meta$k
  ),
  summarize_method(
    "fixed_effect",
    results$p_fixed,
    results$fdr_fixed,
    results$k
  ),
  summarize_method(
    "REML_normal_z",
    results$p_reml,
    results$fdr_reml,
    results$k
  ),
  summarize_method(
    "Hartung_Knapp",
    results$p_hk,
    results$fdr_hk,
    results$k
  ),
  summarize_method(
    "modified_Hartung_Knapp",
    results$p_mkh,
    results$fdr_mkh,
    results$k
  )
)

robustness <- read.csv(
  file.path(
    project_root,
    "results",
    "meta",
    "raw_Ensembl113_robustness",
    "all_genes_robustness_matrix.csv"
  ),
  stringsAsFactors=FALSE,
  check.names=FALSE
)

as_flag <- function(values) {
  tolower(
    as.character(values)
  ) %in% c(
    "true",
    "t",
    "1",
    "yes"
  )
}

tier3 <- as_flag(
  robustness$tier3_GSE56284_independent
)

tier4 <- as_flag(
  robustness$tier4_all_available_LOO_core
)

priority_symbols <- unique(c(
  "FAM107A",
  "UGP2",
  robustness$gene_symbol[tier3],
  robustness$gene_symbol[tier4]
))

priority_symbols <- priority_symbols[
  !is.na(priority_symbols)
  & nzchar(priority_symbols)
]

priority <- results[
  results$gene_symbol %in% priority_symbols,
]

priority <- merge(
  priority,
  robustness[, c(
    "ortholog_id",
    "tier1_full_discovery",
    "tier2_concordant_low_heterogeneity",
    "tier3_GSE56284_independent",
    "tier4_all_available_LOO_core",
    "fdr_without_GSE56284",
    "fdr_without_GSE209926"
  )],
  by="ortholog_id",
  all.x=TRUE
)

priority <- priority[
  order(
    priority$gene_symbol
  ),
]

agreement_summary <- data.frame(
  criterion=c(
    "REML_FDR05",
    "HK_FDR05",
    "modified_HK_FDR05",
    "all_four_methods_FDR05",
    "modified_HK_FDR05_and_sign_concordance_1",
    "modified_HK_FDR05_I2_below_75_and_concordance_1"
  ),
  genes=c(
    sum(results$fdr_reml <= 0.05),
    sum(results$fdr_hk <= 0.05),
    sum(results$fdr_mkh <= 0.05),
    sum(results$significant_all_four),
    sum(
      results$fdr_mkh <= 0.05
      & results$sign_concordance == 1
    ),
    sum(
      results$fdr_mkh <= 0.05
      & results$i2 < 75
      & results$sign_concordance == 1
    )
  )
)

write.csv(
  results,
  file.path(
    output_root,
    "small_k_meta_all_genes.csv"
  ),
  row.names=FALSE
)

write.csv(
  method_summary,
  file.path(
    output_root,
    "small_k_method_summary.csv"
  ),
  row.names=FALSE
)

write.csv(
  comparison,
  file.path(
    output_root,
    "small_k_vs_existing_meta.csv"
  ),
  row.names=FALSE
)

write.csv(
  priority,
  file.path(
    output_root,
    "priority_gene_small_k_sensitivity.csv"
  ),
  row.names=FALSE
)

write.csv(
  agreement_summary,
  file.path(
    output_root,
    "small_k_agreement_summary.csv"
  ),
  row.names=FALSE
)

write.csv(
  failures,
  file.path(
    output_root,
    "small_k_fit_failures.csv"
  ),
  row.names=FALSE
)

capture.output(
  sessionInfo(),
  file=file.path(
    output_root,
    "small_k_R_sessionInfo.txt"
  )
)

cat("\n")
cat(paste(rep("=", 100), collapse=""))
cat("\nSMALL-k METHOD SUMMARY\n")
cat(paste(rep("=", 100), collapse=""))
cat("\n")

print(
  method_summary,
  row.names=FALSE
)

cat("\n")
cat(paste(rep("=", 100), collapse=""))
cat("\nMETHOD-AGREEMENT SUMMARY\n")
cat(paste(rep("=", 100), collapse=""))
cat("\n")

print(
  agreement_summary,
  row.names=FALSE
)

cat("\n")
cat(paste(rep("=", 100), collapse=""))
cat("\nPRIORITY-GENE SENSITIVITY\n")
cat(paste(rep("=", 100), collapse=""))
cat("\n")

priority_print <- priority[, c(
  "ortholog_id",
  "gene_symbol",
  "k",
  "estimate_reml",
  "i2",
  "sign_concordance",
  "fdr_reml",
  "fdr_hk",
  "fdr_mkh",
  "mkh_safeguard_used",
  "tier3_GSE56284_independent",
  "tier4_all_available_LOO_core",
  "fdr_without_GSE56284",
  "fdr_without_GSE209926"
)]

print(
  priority_print,
  row.names=FALSE
)

cat("\nFit failures:", nrow(failures), "\n")

cat(
  "Maximum absolute REML estimate difference ",
  "versus existing implementation: ",
  max(
    comparison$absolute_estimate_difference,
    na.rm=TRUE
  ),
  "\n",
  sep=""
)

cat(
  "Median absolute REML estimate difference: ",
  median(
    comparison$absolute_estimate_difference,
    na.rm=TRUE
  ),
  "\n",
  sep=""
)

if (nrow(results) != 13693) {
  stop(
    "Expected 13,693 fitted orthologs; obtained ",
    nrow(results)
  )
}

if (nrow(failures) != 0) {
  stop(
    "One or more models failed. Review ",
    "small_k_fit_failures.csv."
  )
}

if (anyDuplicated(results$ortholog_id)) {
  stop("Duplicated ortholog IDs in small-k results.")
}

cat(
  "\nPASS: all 13,693 orthologs completed ",
  "the small-k sensitivity analysis.\n",
  sep=""
)
