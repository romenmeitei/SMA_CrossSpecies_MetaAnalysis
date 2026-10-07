args <- commandArgs(trailingOnly=TRUE)

if (length(args) != 4) {
    stop(
        "Expected arguments: study_effects, tier3, ",
        "corrected_meta, output_root"
    )
}

study_path <- normalizePath(args[[1]])
tier3_path <- normalizePath(args[[2]])
meta_path <- normalizePath(args[[3]])
output_root <- normalizePath(
    args[[4]],
    mustWork=FALSE
)

if (!requireNamespace("metafor", quietly=TRUE)) {
    stop("The R package 'metafor' is required.")
}

dir.create(
    output_root,
    recursive=TRUE,
    showWarnings=FALSE
)

pdf_root <- file.path(
    output_root,
    "individual_pdf"
)

png_root <- file.path(
    output_root,
    "individual_png"
)

audit_root <- file.path(
    output_root,
    "audit"
)

dir.create(pdf_root, recursive=TRUE, showWarnings=FALSE)
dir.create(png_root, recursive=TRUE, showWarnings=FALSE)
dir.create(audit_root, recursive=TRUE, showWarnings=FALSE)

studies <- read.csv(
    study_path,
    stringsAsFactors=FALSE,
    check.names=FALSE
)

tier3 <- read.csv(
    tier3_path,
    stringsAsFactors=FALSE,
    check.names=FALSE
)

meta <- read.csv(
    meta_path,
    stringsAsFactors=FALSE,
    check.names=FALSE
)

required_study <- c(
    "accession",
    "ortholog_id",
    "effect",
    "se",
    "species",
    "tissue"
)

required_meta <- c(
    "ortholog_id",
    "estimate",
    "se",
    "p_value",
    "fdr",
    "tau2",
    "i2",
    "sign_concordance"
)

missing_study <- setdiff(
    required_study,
    colnames(studies)
)

missing_meta <- setdiff(
    required_meta,
    colnames(meta)
)

if (length(missing_study) > 0) {
    stop(
        "Study-effect table is missing: ",
        paste(missing_study, collapse=", ")
    )
}

if (length(missing_meta) > 0) {
    stop(
        "Corrected meta table is missing: ",
        paste(missing_meta, collapse=", ")
    )
}

# Standardize gene-symbol fields.
if (!"gene_symbol" %in% colnames(studies)) {
    if ("ortholog_symbol" %in% colnames(studies)) {
        studies$gene_symbol <- studies$ortholog_symbol
    } else {
        studies$gene_symbol <- NA_character_
    }
}

if (!"gene_symbol" %in% colnames(meta)) {
    if ("ortholog_symbol" %in% colnames(meta)) {
        meta$gene_symbol <- meta$ortholog_symbol
    } else {
        meta$gene_symbol <- NA_character_
    }
}

numeric_study <- c(
    "effect",
    "se",
    "z",
    "p_value"
)

for (column in intersect(
    numeric_study,
    colnames(studies)
)) {
    studies[[column]] <- suppressWarnings(
        as.numeric(studies[[column]])
    )
}

numeric_meta <- c(
    "estimate",
    "se",
    "z",
    "p_value",
    "fdr",
    "tau2",
    "i2",
    "sign_concordance",
    "ci_low",
    "ci_high"
)

for (column in intersect(
    numeric_meta,
    colnames(meta)
)) {
    meta[[column]] <- suppressWarnings(
        as.numeric(meta[[column]])
    )
}

if ("stratum" %in% colnames(meta)) {
    meta <- meta[
        meta$stratum == "neural_joint",
        ,
        drop=FALSE
    ]
}

if (anyDuplicated(meta$ortholog_id)) {
    stop(
        "Duplicate ortholog IDs occur in the corrected ",
        "neural meta-analysis."
    )
}

tier3_ids <- unique(
    as.character(tier3$ortholog_id)
)

ugp2_id <- meta$ortholog_id[
    toupper(meta$gene_symbol) == "UGP2"
]

if (length(ugp2_id) != 1) {
    stop(
        "Expected exactly one UGP2 ortholog; found ",
        length(ugp2_id)
    )
}

selected_ids <- unique(
    c(tier3_ids, ugp2_id)
)

selection <- meta[
    meta$ortholog_id %in% selected_ids,
    ,
    drop=FALSE
]

selection$selection_reason <- ifelse(
    selection$ortholog_id == ugp2_id,
    "Cross-modal splicing candidate; corrected Tier 2",
    "Corrected Tier 3: independent of GSE56284"
)

selection$selection_reason[
    toupper(selection$gene_symbol) == "FAM107A"
] <- paste(
    "Corrected Tier 4:",
    "stable in both cross-species LOO analyses"
)

# Place FAM107A and UGP2 first, followed by the
# remaining genes alphabetically.
priority_order <- ifelse(
    toupper(selection$gene_symbol) == "FAM107A",
    1,
    ifelse(
        toupper(selection$gene_symbol) == "UGP2",
        2,
        3
    )
)

selection <- selection[
    order(
        priority_order,
        toupper(selection$gene_symbol)
    ),
    ,
    drop=FALSE
]

if (nrow(selection) != length(selected_ids)) {
    stop(
        "Some selected orthologs are absent from the ",
        "corrected meta-analysis."
    )
}

write.csv(
    selection,
    file.path(
        audit_root,
        "priority_gene_forest_plot_selection.csv"
    ),
    row.names=FALSE
)

accession_order <- c(
    "GSE108094",
    "GSE209926",
    "GSE56284"
)

format_species <- function(value) {
    value <- gsub("_", " ", value)

    value <- sub(
        "^Homo sapiens$",
        "human",
        value,
        ignore.case=TRUE
    )

    value <- sub(
        "^Mus musculus$",
        "mouse",
        value,
        ignore.case=TRUE
    )

    value
}

format_tissue <- function(value) {
    gsub("_", " ", value)
}

format_p <- function(value) {
    if (!is.finite(value)) {
        return("NA")
    }

    format.pval(
        value,
        digits=3,
        eps=1e-300
    )
}

extract_fit_row <- function(
    fit,
    ortholog_id,
    gene_symbol,
    method,
    selection_reason
) {
    data.frame(
        ortholog_id=ortholog_id,
        gene_symbol=gene_symbol,
        row_type="pooled",
        accession=NA_character_,
        species="cross-species",
        tissue="neural",
        method=method,
        estimate=as.numeric(fit$b[1, 1]),
        se=as.numeric(fit$se),
        ci_low=as.numeric(fit$ci.lb),
        ci_high=as.numeric(fit$ci.ub),
        statistic=as.numeric(fit$zval),
        p_value=as.numeric(fit$pval),
        fdr=NA_real_,
        tau2=as.numeric(fit$tau2),
        i2=as.numeric(fit$I2),
        random_effect_weight_percent=NA_real_,
        selection_reason=selection_reason,
        stringsAsFactors=FALSE
    )
}

all_statistics <- list()
plot_payloads <- list()

for (index in seq_len(nrow(selection))) {
    selected <- selection[index, ]

    ortholog_id <- as.character(
        selected$ortholog_id
    )

    gene_symbol <- as.character(
        selected$gene_symbol
    )

    selection_reason <- as.character(
        selected$selection_reason
    )

    current <- studies[
        studies$ortholog_id == ortholog_id
        & is.finite(studies$effect)
        & is.finite(studies$se)
        & studies$se > 0,
        ,
        drop=FALSE
    ]

    current$accession_order <- match(
        current$accession,
        accession_order
    )

    current <- current[
        order(
            current$accession_order,
            current$accession
        ),
        ,
        drop=FALSE
    ]

    if (nrow(current) < 2) {
        stop(
            gene_symbol,
            " has fewer than two usable studies."
        )
    }

    yi <- current$effect
    sei <- current$se

    fit_normal <- metafor::rma.uni(
        yi=yi,
        sei=sei,
        method="REML",
        test="z"
    )

    fit_hk <- metafor::rma.uni(
        yi=yi,
        sei=sei,
        method="REML",
        test="knha",
        adhoc="none"
    )

    # Modified Hartung-Knapp with the ad-hoc safeguard:
    # the KH-adjusted standard error cannot be smaller
    # than the conventional random-effects standard error.
    fit_mkh <- metafor::rma.uni(
        yi=yi,
        sei=sei,
        method="REML",
        test="adhoc"
    )

    corrected_estimate <- as.numeric(
        selected$estimate
    )

    corrected_tau2 <- as.numeric(
        selected$tau2
    )

    estimate_difference <- abs(
        as.numeric(fit_normal$b[1, 1])
        - corrected_estimate
    )

    tau2_difference <- abs(
        as.numeric(fit_normal$tau2)
        - corrected_tau2
    )

    if (
        estimate_difference > 1e-4
        || tau2_difference > 1e-4
    ) {
        stop(
            "Independent metafor result does not match ",
            "the corrected meta-analysis for ",
            gene_symbol,
            ". Estimate difference=",
            estimate_difference,
            "; tau2 difference=",
            tau2_difference
        )
    }

    weights <- 1 / (
        current$se^2
        + corrected_tau2
    )

    weights <- 100 * weights / sum(weights)

    study_rows <- data.frame(
        ortholog_id=ortholog_id,
        gene_symbol=gene_symbol,
        row_type="study",
        accession=current$accession,
        species=vapply(
            current$species,
            format_species,
            character(1)
        ),
        tissue=vapply(
            current$tissue,
            format_tissue,
            character(1)
        ),
        method="Study estimate",
        estimate=current$effect,
        se=current$se,
        ci_low=current$effect - 1.96 * current$se,
        ci_high=current$effect + 1.96 * current$se,
        statistic=current$effect / current$se,
        p_value=2 * pnorm(
            -abs(current$effect / current$se)
        ),
        fdr=NA_real_,
        tau2=NA_real_,
        i2=NA_real_,
        random_effect_weight_percent=weights,
        selection_reason=selection_reason,
        stringsAsFactors=FALSE
    )

    normal_row <- extract_fit_row(
        fit_normal,
        ortholog_id,
        gene_symbol,
        "REML normal-z",
        selection_reason
    )

    # Use the corrected production values for the
    # normal-z REML row.
    normal_row$estimate <- corrected_estimate
    normal_row$se <- as.numeric(selected$se)

    if (
        "ci_low" %in% colnames(selected)
        && is.finite(selected$ci_low)
    ) {
        normal_row$ci_low <- as.numeric(
            selected$ci_low
        )
    } else {
        normal_row$ci_low <- (
            corrected_estimate
            - 1.96 * normal_row$se
        )
    }

    if (
        "ci_high" %in% colnames(selected)
        && is.finite(selected$ci_high)
    ) {
        normal_row$ci_high <- as.numeric(
            selected$ci_high
        )
    } else {
        normal_row$ci_high <- (
            corrected_estimate
            + 1.96 * normal_row$se
        )
    }

    normal_row$statistic <- as.numeric(
        selected$z
    )

    normal_row$p_value <- as.numeric(
        selected$p_value
    )

    normal_row$fdr <- as.numeric(
        selected$fdr
    )

    normal_row$tau2 <- corrected_tau2
    normal_row$i2 <- as.numeric(selected$i2)

    hk_row <- extract_fit_row(
        fit_hk,
        ortholog_id,
        gene_symbol,
        "Hartung-Knapp",
        selection_reason
    )

    mkh_row <- extract_fit_row(
        fit_mkh,
        ortholog_id,
        gene_symbol,
        "Modified Hartung-Knapp",
        selection_reason
    )

    pooled_rows <- rbind(
        normal_row,
        hk_row,
        mkh_row
    )

    all_statistics[[ortholog_id]] <- rbind(
        study_rows,
        pooled_rows
    )

    plot_payloads[[ortholog_id]] <- list(
        study=study_rows,
        pooled=pooled_rows,
        gene_symbol=gene_symbol,
        ortholog_id=ortholog_id,
        reason=selection_reason,
        k=nrow(study_rows)
    )
}

statistics <- do.call(
    rbind,
    all_statistics
)

rownames(statistics) <- NULL

statistics_path <- file.path(
    output_root,
    "priority_gene_forest_plot_statistics.csv"
)

write.csv(
    statistics,
    statistics_path,
    row.names=FALSE
)

draw_forest <- function(payload) {
    study <- payload$study
    pooled <- payload$pooled

    number_studies <- nrow(study)
    maximum_y <- number_studies + 4

    study_y <- seq(
        maximum_y,
        5,
        by=-1
    )

    pooled_y <- c(3, 2, 1)

    all_low <- c(
        study$ci_low,
        pooled$ci_low,
        0
    )

    all_high <- c(
        study$ci_high,
        pooled$ci_high,
        0
    )

    x_limits <- range(
        c(all_low, all_high),
        finite=TRUE
    )

    padding <- max(
        0.15,
        diff(x_limits) * 0.10
    )

    x_limits <- c(
        x_limits[1] - padding,
        x_limits[2] + padding
    )

    study_labels <- paste0(
        study$accession,
        " | ",
        study$species,
        " ",
        study$tissue,
        " | weight ",
        sprintf(
            "%.1f%%",
            study$random_effect_weight_percent
        )
    )

    pooled_labels <- paste0(
        pooled$method,
        " | p=",
        vapply(
            pooled$p_value,
            format_p,
            character(1)
        )
    )

    row_positions <- c(
        study_y,
        pooled_y
    )

    row_labels <- c(
        study_labels,
        pooled_labels
    )

    par(
        mar=c(5.5, 16, 5, 2),
        las=1
    )

    plot(
        NA,
        xlim=x_limits,
        ylim=c(0.3, maximum_y + 0.7),
        xlab="Disease-control log2 fold change",
        ylab="",
        yaxt="n",
        main=paste0(
            payload$gene_symbol,
            " (",
            payload$ortholog_id,
            ")"
        ),
        cex.main=1.15,
        bty="n"
    )

    axis(
        side=2,
        at=row_positions,
        labels=row_labels,
        tick=FALSE,
        cex.axis=0.72
    )

    abline(
        v=0,
        col="grey45",
        lty=2,
        lwd=1.2
    )

    abline(
        h=4,
        col="grey80",
        lty=3
    )

    for (row in seq_len(nrow(study))) {
        segments(
            study$ci_low[row],
            study_y[row],
            study$ci_high[row],
            study_y[row],
            col="#2455A4",
            lwd=2
        )

        point_size <- 0.9 + (
            study$random_effect_weight_percent[row]
            / max(
                study$random_effect_weight_percent
            )
        )

        points(
            study$estimate[row],
            study_y[row],
            pch=15,
            cex=point_size,
            col="#2455A4"
        )
    }

    pooled_colours <- c(
        "#222222",
        "#C44E52",
        "#7A3E9D"
    )

    for (row in seq_len(nrow(pooled))) {
        y <- pooled_y[row]

        polygon(
            x=c(
                pooled$ci_low[row],
                pooled$estimate[row],
                pooled$ci_high[row],
                pooled$estimate[row]
            ),
            y=c(
                y,
                y + 0.20,
                y,
                y - 0.20
            ),
            col=pooled_colours[row],
            border=pooled_colours[row]
        )
    }

    mtext(
        payload$reason,
        side=3,
        line=0.35,
        cex=0.78,
        col="grey30"
    )

    mtext(
        paste0(
            "Squares: study estimates and 95% CIs; ",
            "diamonds: pooled estimates. ",
            "Positive values indicate higher expression in SMA."
        ),
        side=1,
        line=4.1,
        cex=0.72,
        col="grey30"
    )
}

combined_pdf <- file.path(
    output_root,
    "priority_genes_study_level_forest_plots.pdf"
)

pdf(
    combined_pdf,
    width=12,
    height=7,
    onefile=TRUE
)

for (ortholog_id in selection$ortholog_id) {
    draw_forest(
        plot_payloads[[ortholog_id]]
    )
}

dev.off()

for (ortholog_id in selection$ortholog_id) {
    payload <- plot_payloads[[ortholog_id]]

    safe_symbol <- gsub(
        "[^A-Za-z0-9._-]",
        "_",
        payload$gene_symbol
    )

    file_stem <- paste0(
        safe_symbol,
        "_",
        ortholog_id,
        "_forest_plot"
    )

    pdf(
        file.path(
            pdf_root,
            paste0(file_stem, ".pdf")
        ),
        width=12,
        height=7
    )

    draw_forest(payload)
    dev.off()

    png(
        file.path(
            png_root,
            paste0(file_stem, ".png")
        ),
        width=2400,
        height=1400,
        res=200
    )

    draw_forest(payload)
    dev.off()
}

pooled_summary <- statistics[
    statistics$row_type == "pooled",
    c(
        "ortholog_id",
        "gene_symbol",
        "method",
        "estimate",
        "se",
        "ci_low",
        "ci_high",
        "statistic",
        "p_value",
        "fdr",
        "tau2",
        "i2",
        "selection_reason"
    )
]

write.csv(
    pooled_summary,
    file.path(
        output_root,
        "priority_gene_pooled_method_summary.csv"
    ),
    row.names=FALSE
)

cat(
    "\n",
    paste(rep("=", 108), collapse=""),
    "\n",
    sep=""
)

cat("PRIORITY-GENE FOREST-PLOT SELECTION\n")

cat(
    paste(rep("=", 108), collapse=""),
    "\n",
    sep=""
)

print(
    selection[
        ,
        c(
            "ortholog_id",
            "gene_symbol",
            "estimate",
            "fdr",
            "tau2",
            "i2",
            "selection_reason"
        )
    ],
    row.names=FALSE
)

cat(
    "\n",
    paste(rep("=", 108), collapse=""),
    "\n",
    sep=""
)

cat("POOLED METHOD SUMMARY\n")

cat(
    paste(rep("=", 108), collapse=""),
    "\n",
    sep=""
)

print(
    pooled_summary[
        ,
        c(
            "gene_symbol",
            "method",
            "estimate",
            "ci_low",
            "ci_high",
            "p_value",
            "fdr"
        )
    ],
    row.names=FALSE
)

cat(
    "\nCombined forest-plot PDF: ",
    combined_pdf,
    "\n",
    sep=""
)

cat(
    "Statistics table: ",
    statistics_path,
    "\n",
    sep=""
)

cat(
    "\nPASS: study-level forest plots created ",
    "for all corrected Tier 3 genes plus UGP2.\n"
)
