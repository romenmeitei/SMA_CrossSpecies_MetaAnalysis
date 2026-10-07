from pathlib import Path
from itertools import combinations
from datetime import datetime, timezone
import json
import sys

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from scipy.stats import spearmanr


raw_root = Path(sys.argv[1]).resolve()
project_root = Path(sys.argv[2]).resolve()
output_root = Path(sys.argv[3]).resolve()

figure_root = output_root / "figures"
audit_root = output_root / "audit"

figure_root.mkdir(parents=True, exist_ok=True)
audit_root.mkdir(parents=True, exist_ok=True)

accession = "GSE56284"

de_path = (
    raw_root
    / "differential_expression"
    / "GSE56284_disease-control_raw_Ensembl113.csv"
)

deseq_qc_root = (
    raw_root
    / "qc"
    / "deseq2"
    / accession
)

normalized_path = (
    deseq_qc_root
    / "GSE56284_normalized_counts.csv"
)

sample_table_path = (
    deseq_qc_root
    / "GSE56284_collapsed_sample_table.csv"
)

run_manifest_path = (
    raw_root
    / "audit"
    / "selected_sra_ena_run_manifest.csv"
)

quant_root = (
    raw_root
    / "salmon_quant"
    / accession
)

fastp_root = (
    raw_root
    / "qc"
    / accession
)

tier3_path = (
    project_root
    / "results"
    / "meta"
    / "raw_Ensembl113_robustness"
    / "tier3_GSE56284_independent.csv"
)

required_paths = [
    de_path,
    normalized_path,
    sample_table_path,
    run_manifest_path,
    tier3_path,
]

for path in required_paths:
    if not path.is_file():
        raise FileNotFoundError(path)


# ------------------------------------------------------------------
# Read and standardize the existing results.
# ------------------------------------------------------------------

de = pd.read_csv(de_path)

normalized = pd.read_csv(
    normalized_path,
    index_col=0,
)

sample_table = pd.read_csv(sample_table_path)
run_manifest = pd.read_csv(run_manifest_path)
tier3 = pd.read_csv(tier3_path)

de["gene_id"] = (
    de["gene_id"]
    .astype(str)
    .str.split(".")
    .str[0]
)

normalized.index = (
    normalized.index
    .astype(str)
    .str.split(".")
    .str[0]
)

if "gene_symbol" not in de.columns:
    raise ValueError(
        "The GSE56284 DE table lacks gene_symbol."
    )

for column in [
    "effect",
    "se",
    "z",
    "p_value",
    "fdr",
    "base_mean",
]:
    if column in de.columns:
        de[column] = pd.to_numeric(
            de[column],
            errors="coerce",
        )

sample_table["sample_id"] = (
    sample_table["sample_id"].astype(str)
)

sample_table["biological_unit"] = (
    sample_table["biological_unit"].astype(str)
)

sample_table["condition"] = (
    sample_table["condition"].astype(str)
)

condition_by_unit = (
    sample_table
    .set_index("biological_unit")["condition"]
)

missing_units = sorted(
    set(normalized.columns)
    - set(condition_by_unit.index)
)

if missing_units:
    raise RuntimeError(
        "Normalized-count columns absent from "
        f"sample metadata: {missing_units}"
    )

conditions = condition_by_unit.loc[
    normalized.columns
]

if set(conditions) != {"control", "disease"}:
    raise RuntimeError(
        f"Unexpected conditions: {set(conditions)}"
    )


# ------------------------------------------------------------------
# Salmon and fastp technical QC.
# ------------------------------------------------------------------

technical_records = []

for row in sample_table.itertuples(index=False):
    sample_id = str(row.sample_id)
    biological_unit = str(row.biological_unit)
    condition = str(row.condition)

    meta_path = (
        quant_root
        / sample_id
        / "aux_info"
        / "meta_info.json"
    )

    if not meta_path.is_file():
        raise FileNotFoundError(meta_path)

    salmon = json.loads(meta_path.read_text())

    fastp_files = sorted(
        (fastp_root / sample_id).glob(
            "*.fastp.json"
        )
    )

    if not fastp_files:
        raise FileNotFoundError(
            f"No fastp JSON for {sample_id}"
        )

    before_reads = 0
    after_reads = 0
    before_bases = 0
    after_bases = 0
    weighted_q30_bases = 0.0
    weighted_gc_bases = 0.0

    for fastp_path in fastp_files:
        fastp = json.loads(
            fastp_path.read_text()
        )

        before = (
            fastp.get("summary", {})
            .get("before_filtering", {})
        )

        after = (
            fastp.get("summary", {})
            .get("after_filtering", {})
        )

        before_reads += int(
            before.get("total_reads", 0)
        )

        after_reads += int(
            after.get("total_reads", 0)
        )

        before_bases += int(
            before.get("total_bases", 0)
        )

        current_after_bases = int(
            after.get("total_bases", 0)
        )

        after_bases += current_after_bases

        weighted_q30_bases += (
            float(after.get("q30_rate", np.nan))
            * current_after_bases
        )

        weighted_gc_bases += (
            float(after.get("gc_content", np.nan))
            * current_after_bases
        )

    q30_rate = (
        weighted_q30_bases / after_bases
        if after_bases > 0
        else np.nan
    )

    gc_content = (
        weighted_gc_bases / after_bases
        if after_bases > 0
        else np.nan
    )

    technical_records.append({
        "sample_id": sample_id,
        "biological_unit": biological_unit,
        "condition": condition,
        "fastp_files": len(fastp_files),
        "fastp_reads_before": before_reads,
        "fastp_reads_after": after_reads,
        "fastp_retained_fraction": (
            after_reads / before_reads
            if before_reads > 0
            else np.nan
        ),
        "fastp_q30_rate_after": q30_rate,
        "fastp_gc_content_after": gc_content,
        "salmon_processed_reads":
            salmon.get("num_processed"),
        "salmon_mapped_reads":
            salmon.get("num_mapped"),
        "salmon_mapping_rate":
            salmon.get("percent_mapped"),
        "salmon_version":
            salmon.get("salmon_version"),
        "library_types": ";".join(
            salmon.get("library_types", [])
        ),
    })

technical = pd.DataFrame(technical_records)

manifest_subset = run_manifest.loc[
    run_manifest["sample_id"].isin(
        sample_table["sample_id"]
    )
].copy()

manifest_columns = [
    column
    for column in [
        "sample_id",
        "run_accession",
        "library_strategy",
        "library_source",
        "library_selection",
        "library_layout",
        "instrument_model",
        "scientific_name",
        "read_count",
        "base_count",
        "download_bytes",
        "fastq_ftp",
    ]
    if column in manifest_subset.columns
]

manifest_subset = manifest_subset[
    manifest_columns
].copy()

technical = technical.merge(
    manifest_subset,
    on="sample_id",
    how="left",
    validate="one_to_one",
)

technical.to_csv(
    output_root
    / "GSE56284_sample_technical_QC.csv",
    index=False,
)

technical_group = (
    technical
    .groupby("condition", as_index=False)
    .agg(
        samples=("sample_id", "nunique"),
        mean_processed_reads=(
            "salmon_processed_reads",
            "mean",
        ),
        minimum_mapping_rate=(
            "salmon_mapping_rate",
            "min",
        ),
        mean_mapping_rate=(
            "salmon_mapping_rate",
            "mean",
        ),
        maximum_mapping_rate=(
            "salmon_mapping_rate",
            "max",
        ),
        mean_q30_rate=(
            "fastp_q30_rate_after",
            "mean",
        ),
        mean_gc_content=(
            "fastp_gc_content_after",
            "mean",
        ),
    )
)

technical_group.to_csv(
    audit_root
    / "GSE56284_technical_QC_by_condition.csv",
    index=False,
)


# ------------------------------------------------------------------
# PCA and sample-correlation diagnostics.
# ------------------------------------------------------------------

log_counts = np.log2(
    normalized.astype(float) + 1.0
)

gene_variance = log_counts.var(axis=1)

variable_genes = gene_variance.loc[
    gene_variance > 0
].sort_values(ascending=False)

top_n = min(1000, len(variable_genes))
top_gene_ids = variable_genes.head(top_n).index

top_matrix = (
    log_counts.loc[top_gene_ids]
    .T
    .astype(float)
)

centered = (
    top_matrix
    - top_matrix.mean(axis=0)
)

u, singular_values, vt = np.linalg.svd(
    centered.to_numpy(),
    full_matrices=False,
)

scores = (
    u
    * singular_values
)

variance_explained = (
    singular_values ** 2
    / np.sum(singular_values ** 2)
)

pca = pd.DataFrame({
    "biological_unit": top_matrix.index,
    "condition": conditions.loc[
        top_matrix.index
    ].values,
    "PC1": scores[:, 0],
    "PC2": scores[:, 1],
    "PC1_variance_fraction":
        variance_explained[0],
    "PC2_variance_fraction":
        variance_explained[1],
})

pca = pca.merge(
    sample_table[
        [
            "sample_id",
            "biological_unit",
        ]
    ],
    on="biological_unit",
    how="left",
    validate="one_to_one",
)

pca.to_csv(
    output_root
    / "GSE56284_PCA_coordinates.csv",
    index=False,
)

colours = {
    "control": "#3B6FB6",
    "disease": "#C94C4C",
}

fig, ax = plt.subplots(
    figsize=(7.5, 6.2)
)

for condition, group in pca.groupby(
    "condition",
    sort=False,
):
    ax.scatter(
        group["PC1"],
        group["PC2"],
        s=90,
        color=colours[condition],
        label=condition,
        edgecolor="black",
        linewidth=0.7,
    )

    for row in group.itertuples(index=False):
        ax.annotate(
            row.biological_unit,
            (row.PC1, row.PC2),
            xytext=(5, 5),
            textcoords="offset points",
            fontsize=8,
        )

ax.axhline(0, color="grey", linewidth=0.7)
ax.axvline(0, color="grey", linewidth=0.7)

ax.set_xlabel(
    "PC1 "
    f"({100 * variance_explained[0]:.1f}% variance)"
)

ax.set_ylabel(
    "PC2 "
    f"({100 * variance_explained[1]:.1f}% variance)"
)

ax.set_title(
    "GSE56284: PCA of log2 normalized counts\n"
    f"Top {top_n:,} variable genes"
)

ax.legend(frameon=False)
fig.tight_layout()

for suffix in ["png", "pdf"]:
    fig.savefig(
        figure_root
        / f"GSE56284_PCA_top{top_n}.{suffix}",
        dpi=220,
        bbox_inches="tight",
    )

plt.close(fig)

correlation = log_counts.corr(
    method="pearson"
)

correlation.to_csv(
    output_root
    / "GSE56284_sample_correlation_matrix.csv"
)

fig, ax = plt.subplots(
    figsize=(7.5, 6.5)
)

image = ax.imshow(
    correlation.to_numpy(),
    vmin=0.80,
    vmax=1.00,
    cmap="viridis",
)

ax.set_xticks(
    np.arange(len(correlation.columns))
)

ax.set_xticklabels(
    correlation.columns,
    rotation=45,
    ha="right",
)

ax.set_yticks(
    np.arange(len(correlation.index))
)

ax.set_yticklabels(correlation.index)

for i in range(len(correlation.index)):
    for j in range(len(correlation.columns)):
        ax.text(
            j,
            i,
            f"{correlation.iloc[i, j]:.3f}",
            ha="center",
            va="center",
            fontsize=7,
            color=(
                "white"
                if correlation.iloc[i, j] < 0.90
                else "black"
            ),
        )

fig.colorbar(
    image,
    ax=ax,
    label="Pearson correlation",
)

ax.set_title(
    "GSE56284 sample correlation\n"
    "log2 normalized counts"
)

fig.tight_layout()

for suffix in ["png", "pdf"]:
    fig.savefig(
        figure_root
        / f"GSE56284_sample_correlation.{suffix}",
        dpi=220,
        bbox_inches="tight",
    )

plt.close(fig)


# ------------------------------------------------------------------
# Sample-level outlier diagnostics.
# ------------------------------------------------------------------

outlier_records = []

for sample in normalized.columns:
    condition = conditions[sample]

    same_group = [
        other
        for other in normalized.columns
        if other != sample
        and conditions[other] == condition
    ]

    other_group = [
        other
        for other in normalized.columns
        if conditions[other] != condition
    ]

    same_correlation = correlation.loc[
        sample,
        same_group,
    ]

    other_correlation = correlation.loc[
        sample,
        other_group,
    ]

    sample_vector = top_matrix.loc[
        sample
    ].to_numpy()

    same_centroid = (
        top_matrix.loc[same_group]
        .mean(axis=0)
        .to_numpy()
    )

    other_centroid = (
        top_matrix.loc[other_group]
        .mean(axis=0)
        .to_numpy()
    )

    own_distance = float(
        np.sqrt(
            np.mean(
                (
                    sample_vector
                    - same_centroid
                ) ** 2
            )
        )
    )

    other_distance = float(
        np.sqrt(
            np.mean(
                (
                    sample_vector
                    - other_centroid
                ) ** 2
            )
        )
    )

    outlier_records.append({
        "biological_unit": sample,
        "condition": condition,
        "median_same_condition_correlation":
            float(np.median(same_correlation)),
        "median_other_condition_correlation":
            float(np.median(other_correlation)),
        "correlation_separation_margin": (
            float(np.median(same_correlation))
            - float(np.median(other_correlation))
        ),
        "distance_to_leave_one_out_own_centroid":
            own_distance,
        "distance_to_other_condition_centroid":
            other_distance,
        "closer_to_own_condition": (
            own_distance < other_distance
        ),
    })

outliers = pd.DataFrame(outlier_records)

outliers.to_csv(
    output_root
    / "GSE56284_sample_outlier_diagnostics.csv",
    index=False,
)


# ------------------------------------------------------------------
# Exact three-versus-three global-separation permutation.
# ------------------------------------------------------------------

sample_names = list(top_matrix.index)
observed_disease = set(
    conditions.loc[
        conditions.eq("disease")
    ].index
)

data_matrix = top_matrix.to_numpy()


def separation_statistic(disease_indices):
    disease_indices = sorted(disease_indices)

    control_indices = [
        index
        for index in range(len(sample_names))
        if index not in disease_indices
    ]

    disease_values = data_matrix[
        disease_indices,
        :
    ]

    control_values = data_matrix[
        control_indices,
        :
    ]

    disease_centroid = disease_values.mean(axis=0)
    control_centroid = control_values.mean(axis=0)

    between = np.mean(
        (
            disease_centroid
            - control_centroid
        ) ** 2
    )

    within_disease = np.mean(
        (
            disease_values
            - disease_centroid
        ) ** 2
    )

    within_control = np.mean(
        (
            control_values
            - control_centroid
        ) ** 2
    )

    within = (
        within_disease
        + within_control
    ) / 2.0

    return float(
        between / (within + np.finfo(float).eps)
    )


permutation_records = []

for disease_indices in combinations(
    range(len(sample_names)),
    3,
):
    disease_samples = {
        sample_names[index]
        for index in disease_indices
    }

    permutation_records.append({
        "putative_disease_samples":
            ";".join(sorted(disease_samples)),
        "separation_statistic":
            separation_statistic(
                disease_indices
            ),
        "is_observed_assignment": (
            disease_samples
            == observed_disease
        ),
        "is_observed_or_complement": (
            disease_samples
            == observed_disease
            or disease_samples
            == (
                set(sample_names)
                - observed_disease
            )
        ),
    })

permutations = pd.DataFrame(
    permutation_records
).sort_values(
    "separation_statistic",
    ascending=False,
).reset_index(drop=True)

observed_statistic = float(
    permutations.loc[
        permutations[
            "is_observed_assignment"
        ],
        "separation_statistic",
    ].iloc[0]
)

exact_p = float(
    (
        permutations["separation_statistic"]
        >= observed_statistic - 1e-12
    ).mean()
)

permutations["rank"] = (
    permutations[
        "separation_statistic"
    ].rank(
        method="min",
        ascending=False,
    ).astype(int)
)

permutations.to_csv(
    output_root
    / "GSE56284_exact_label_permutation.csv",
    index=False,
)


# ------------------------------------------------------------------
# Leave-one-sample-out expression-effect stability.
# ------------------------------------------------------------------

full_control = [
    sample
    for sample in normalized.columns
    if conditions[sample] == "control"
]

full_disease = [
    sample
    for sample in normalized.columns
    if conditions[sample] == "disease"
]

full_log_effect = (
    log_counts[full_disease].mean(axis=1)
    - log_counts[full_control].mean(axis=1)
)

loo_records = []

for omitted in normalized.columns:
    retained = [
        sample
        for sample in normalized.columns
        if sample != omitted
    ]

    control_samples = [
        sample
        for sample in retained
        if conditions[sample] == "control"
    ]

    disease_samples = [
        sample
        for sample in retained
        if conditions[sample] == "disease"
    ]

    loo_effect = (
        log_counts[disease_samples].mean(axis=1)
        - log_counts[control_samples].mean(axis=1)
    )

    finite = (
        np.isfinite(full_log_effect)
        & np.isfinite(loo_effect)
    )

    correlation_value = spearmanr(
        full_log_effect.loc[finite],
        loo_effect.loc[finite],
    ).statistic

    nonzero = (
        finite
        & full_log_effect.ne(0)
        & loo_effect.ne(0)
    )

    sign_agreement = float(
        (
            np.sign(full_log_effect.loc[nonzero])
            == np.sign(loo_effect.loc[nonzero])
        ).mean()
    )

    loo_records.append({
        "omitted_biological_unit": omitted,
        "omitted_condition":
            conditions[omitted],
        "remaining_control_units":
            len(control_samples),
        "remaining_disease_units":
            len(disease_samples),
        "effect_spearman_with_full":
            correlation_value,
        "effect_sign_agreement_with_full":
            sign_agreement,
        "median_absolute_effect_change":
            float(
                np.median(
                    np.abs(
                        loo_effect.loc[finite]
                        - full_log_effect.loc[finite]
                    )
                )
            ),
    })

loo_stability = pd.DataFrame(loo_records)

loo_stability.to_csv(
    output_root
    / "GSE56284_leave_one_sample_out_stability.csv",
    index=False,
)


# ------------------------------------------------------------------
# DE distribution and potential inflation diagnostics.
# ------------------------------------------------------------------

finite_fdr = np.isfinite(de["fdr"])
significant = finite_fdr & (de["fdr"] <= 0.05)
finite_p = np.isfinite(de["p_value"])
finite_z = np.isfinite(de["z"])
finite_effect = np.isfinite(de["effect"])

median_chisq = float(
    np.median(
        de.loc[finite_z, "z"] ** 2
    )
)

lambda_gc_descriptive = (
    median_chisq / 0.454936423119572
)

de_summary = pd.DataFrame([{
    "tested_genes": len(de),
    "finite_effects": int(finite_effect.sum()),
    "finite_p_values": int(finite_p.sum()),
    "finite_fdr_values": int(finite_fdr.sum()),
    "fdr_le_0.05": int(significant.sum()),
    "fdr_le_0.05_fraction": float(
        significant.sum() / finite_fdr.sum()
    ),
    "significant_positive_effect": int(
        (
            significant
            & (de["effect"] > 0)
        ).sum()
    ),
    "significant_negative_effect": int(
        (
            significant
            & (de["effect"] < 0)
        ).sum()
    ),
    "nominal_p_lt_0.05": int(
        (
            finite_p
            & (de["p_value"] < 0.05)
        ).sum()
    ),
    "absolute_effect_ge_1": int(
        (
            finite_effect
            & (de["effect"].abs() >= 1)
        ).sum()
    ),
    "absolute_effect_ge_2": int(
        (
            finite_effect
            & (de["effect"].abs() >= 2)
        ).sum()
    ),
    "absolute_z_gt_8": int(
        (
            finite_z
            & (de["z"].abs() > 8)
        ).sum()
    ),
    "absolute_z_gt_20": int(
        (
            finite_z
            & (de["z"].abs() > 20)
        ).sum()
    ),
    "maximum_absolute_z": float(
        de.loc[finite_z, "z"].abs().max()
    ),
    "exact_zero_p_values": int(
        (
            finite_p
            & de["p_value"].eq(0)
        ).sum()
    ),
    "descriptive_lambda_gc": (
        lambda_gc_descriptive
    ),
}])

de_summary.to_csv(
    output_root
    / "GSE56284_DE_distribution_summary.csv",
    index=False,
)

expression_bins = pd.cut(
    de["base_mean"],
    bins=[
        -np.inf,
        1,
        10,
        100,
        1000,
        np.inf,
    ],
    labels=[
        "<=1",
        "1-10",
        "10-100",
        "100-1000",
        ">1000",
    ],
)

expression_bin_summary = (
    de.assign(
        expression_bin=expression_bins,
        significant=significant,
    )
    .groupby(
        "expression_bin",
        observed=False,
        as_index=False,
    )
    .agg(
        genes=("gene_id", "size"),
        significant_genes=(
            "significant",
            "sum",
        ),
        median_absolute_effect=(
            "effect",
            lambda values: float(
                np.nanmedian(
                    np.abs(values)
                )
            ),
        ),
        median_absolute_z=(
            "z",
            lambda values: float(
                np.nanmedian(
                    np.abs(values)
                )
            ),
        ),
    )
)

expression_bin_summary[
    "significant_fraction"
] = (
    expression_bin_summary[
        "significant_genes"
    ]
    / expression_bin_summary["genes"]
)

expression_bin_summary.to_csv(
    output_root
    / "GSE56284_DE_by_expression_bin.csv",
    index=False,
)

fig, axes = plt.subplots(
    1,
    2,
    figsize=(12, 4.8),
)

positive_p = de.loc[
    finite_p
    & (de["p_value"] > 0),
    "p_value",
]

axes[0].hist(
    de.loc[finite_p, "p_value"],
    bins=50,
    color="#4C72B0",
    edgecolor="white",
)

axes[0].set_xlabel("DESeq2 p-value")
axes[0].set_ylabel("Genes")
axes[0].set_title("P-value distribution")

axes[1].hist(
    de.loc[
        finite_effect
        & (de["effect"].abs() <= 10),
        "effect",
    ],
    bins=60,
    color="#C44E52",
    edgecolor="white",
)

axes[1].axvline(
    0,
    color="black",
    linestyle="--",
    linewidth=0.8,
)

axes[1].set_xlabel(
    "DESeq2 log2 fold change"
)

axes[1].set_ylabel("Genes")
axes[1].set_title(
    "Effect distribution\n"
    "(display restricted to |effect| <= 10)"
)

fig.suptitle(
    "GSE56284 differential-expression diagnostics"
)

fig.tight_layout()

for suffix in ["png", "pdf"]:
    fig.savefig(
        figure_root
        / f"GSE56284_DE_distributions.{suffix}",
        dpi=220,
        bbox_inches="tight",
    )

plt.close(fig)


# ------------------------------------------------------------------
# Verify DE orientation against normalized counts.
# ------------------------------------------------------------------

simple_effect = full_log_effect.rename(
    "normalized_count_log_ratio"
).reset_index()

simple_effect.columns = [
    "gene_id",
    "normalized_count_log_ratio",
]

orientation = de[
    [
        "gene_id",
        "effect",
    ]
].merge(
    simple_effect,
    on="gene_id",
    how="inner",
    validate="one_to_one",
)

orientation_finite = (
    np.isfinite(orientation["effect"])
    & np.isfinite(
        orientation[
            "normalized_count_log_ratio"
        ]
    )
)

orientation_spearman = spearmanr(
    orientation.loc[
        orientation_finite,
        "effect",
    ],
    orientation.loc[
        orientation_finite,
        "normalized_count_log_ratio",
    ],
).statistic

orientation_nonzero = (
    orientation_finite
    & orientation["effect"].ne(0)
    & orientation[
        "normalized_count_log_ratio"
    ].ne(0)
)

orientation_sign_agreement = float(
    (
        np.sign(
            orientation.loc[
                orientation_nonzero,
                "effect",
            ]
        )
        == np.sign(
            orientation.loc[
                orientation_nonzero,
                "normalized_count_log_ratio",
            ]
        )
    ).mean()
)

orientation_summary = pd.DataFrame([{
    "genes_compared":
        int(orientation_finite.sum()),
    "effect_spearman":
        orientation_spearman,
    "sign_agreement":
        orientation_sign_agreement,
}])

orientation_summary.to_csv(
    audit_root
    / "GSE56284_DE_orientation_check.csv",
    index=False,
)


# ------------------------------------------------------------------
# Priority genes and originally reported validation markers.
# ------------------------------------------------------------------

priority_symbols = sorted(
    set(
        tier3["gene_symbol"]
        .dropna()
        .astype(str)
    )
    | {"UGP2"}
)

priority_columns = [
    column
    for column in [
        "gene_id",
        "gene_symbol",
        "effect",
        "se",
        "z",
        "p_value",
        "fdr",
        "base_mean",
    ]
    if column in de.columns
]

priority = de.loc[
    de["gene_symbol"]
    .astype(str)
    .str.upper()
    .isin({
        symbol.upper()
        for symbol in priority_symbols
    }),
    priority_columns,
].copy()

priority["absolute_z"] = (
    priority["z"].abs()
)

priority = priority.sort_values(
    ["fdr", "gene_symbol"],
    na_position="last",
)

priority.to_csv(
    output_root
    / "GSE56284_priority_gene_DE_evidence.csv",
    index=False,
)

marker_symbols = {
    "CRABP1",
    "CRABP2",
    "NKX2-2",
    "NKX2.2",
    "SMN1",
}

marker_check = de.loc[
    de["gene_symbol"]
    .astype(str)
    .str.upper()
    .isin(marker_symbols),
    priority_columns,
].copy()

marker_check.to_csv(
    output_root
    / "GSE56284_published_marker_check.csv",
    index=False,
)


# ------------------------------------------------------------------
# Original-publication comparison and design limitation.
# ------------------------------------------------------------------

published_comparison = pd.DataFrame([
    {
        "analysis": (
            "Original publication/Cuffdiff"
        ),
        "feature_unit": "transcripts",
        "increased": 3094,
        "decreased": 6964,
        "total_significant": 10058,
        "note": (
            "Published transcript-level result; "
            "not directly numerically equivalent "
            "to the current Ensembl 113 gene-level "
            "DESeq2 result."
        ),
    },
    {
        "analysis": (
            "Current Salmon-DESeq2 reanalysis"
        ),
        "feature_unit": "genes",
        "increased": int(
            (
                significant
                & (de["effect"] > 0)
            ).sum()
        ),
        "decreased": int(
            (
                significant
                & (de["effect"] < 0)
            ).sum()
        ),
        "total_significant": int(
            significant.sum()
        ),
        "note": (
            "Raw-read reanalysis with Salmon 2.7, "
            "Ensembl release 113, and DESeq2."
        ),
    },
])

published_comparison.to_csv(
    output_root
    / "GSE56284_published_vs_current_DE_counts.csv",
    index=False,
)

design_context = {
    "accession": accession,
    "control_cell_line": (
        "Hb9 normal mouse embryonic stem-cell line"
    ),
    "disease_cell_line": (
        "A2 severe-SMA mouse embryonic stem-cell line"
    ),
    "control_genotype": "wildtype",
    "disease_genotype": "SMN2+/+;mSmn-/-",
    "sequencing_replicates_per_cell_line": 3,
    "sequencing_protocol": (
        "50-cycle single-end Illumina HiSeq 2500"
    ),
    "genotype_confounded_with_cell_line": True,
    "interpretation": (
        "The within-study contrast measures the "
        "difference between the Hb9 control and A2 "
        "SMA cell-line systems. It cannot by itself "
        "separate an SMA genotype effect from other "
        "stable cell-line-specific effects."
    ),
    "recommended_role": (
        "Influential discovery dataset requiring "
        "leave-one-dataset-out and independent "
        "cross-species support."
    ),
    "created_utc": (
        datetime.now(timezone.utc).isoformat()
    ),
}

(
    audit_root
    / "GSE56284_design_context.json"
).write_text(
    json.dumps(
        design_context,
        indent=2,
    )
    + "\n"
)

source_audit = pd.DataFrame([
    {
        "source": "NCBI GEO",
        "identifier": "GSE56284",
        "url": (
            "https://www.ncbi.nlm.nih.gov/geo/"
            "query/acc.cgi?acc=GSE56284"
        ),
        "information": (
            "Study design, n=3/genotype, cell-line "
            "names, sequencing platform and original "
            "DE transcript counts."
        ),
    },
    {
        "source": "Original publication",
        "identifier": "PMID:25191843",
        "url": (
            "https://pubmed.ncbi.nlm.nih.gov/"
            "25191843/"
        ),
        "information": (
            "Triplicates for each cell line and "
            "50-cycle single-end RNA sequencing."
        ),
    },
])

source_audit.to_csv(
    audit_root
    / "GSE56284_external_source_audit.csv",
    index=False,
)


# ------------------------------------------------------------------
# Concise machine-readable audit summary.
# ------------------------------------------------------------------

summary = {
    "samples": len(sample_table),
    "control_units": int(
        conditions.eq("control").sum()
    ),
    "disease_units": int(
        conditions.eq("disease").sum()
    ),
    "minimum_salmon_mapping_rate": float(
        technical[
            "salmon_mapping_rate"
        ].min()
    ),
    "maximum_salmon_mapping_rate": float(
        technical[
            "salmon_mapping_rate"
        ].max()
    ),
    "all_samples_closer_to_own_condition": bool(
        outliers[
            "closer_to_own_condition"
        ].all()
    ),
    "PC1_variance_fraction": float(
        variance_explained[0]
    ),
    "PC2_variance_fraction": float(
        variance_explained[1]
    ),
    "observed_global_separation_statistic":
        observed_statistic,
    "exact_label_permutation_p": exact_p,
    "minimum_LOO_effect_spearman": float(
        loo_stability[
            "effect_spearman_with_full"
        ].min()
    ),
    "minimum_LOO_sign_agreement": float(
        loo_stability[
            "effect_sign_agreement_with_full"
        ].min()
    ),
    "DE_FDR_le_0.05": int(
        significant.sum()
    ),
    "DE_FDR_le_0.05_fraction": float(
        significant.sum() / finite_fdr.sum()
    ),
    "DE_orientation_spearman": float(
        orientation_spearman
    ),
    "DE_orientation_sign_agreement": float(
        orientation_sign_agreement
    ),
    "genotype_confounded_with_cell_line": True,
}

(
    audit_root
    / "GSE56284_investigation_summary.json"
).write_text(
    json.dumps(summary, indent=2) + "\n"
)


# ------------------------------------------------------------------
# Console report.
# ------------------------------------------------------------------

print("\n" + "=" * 108)
print("GSE56284 TECHNICAL QC")
print("=" * 108)
print(technical.to_string(index=False))

print("\n" + "=" * 108)
print("GSE56284 CONDITION-LEVEL TECHNICAL QC")
print("=" * 108)
print(technical_group.to_string(index=False))

print("\n" + "=" * 108)
print("GSE56284 PCA")
print("=" * 108)
print(pca.to_string(index=False))

print("\n" + "=" * 108)
print("GSE56284 SAMPLE OUTLIER DIAGNOSTICS")
print("=" * 108)
print(outliers.to_string(index=False))

print("\n" + "=" * 108)
print("GSE56284 EXACT LABEL PERMUTATION")
print("=" * 108)
print(
    "Observed statistic:",
    observed_statistic,
)
print(
    "Exact permutation p-value:",
    exact_p,
)
print(
    permutations.head(6).to_string(
        index=False
    )
)

print("\n" + "=" * 108)
print("GSE56284 LEAVE-ONE-SAMPLE-OUT STABILITY")
print("=" * 108)
print(loo_stability.to_string(index=False))

print("\n" + "=" * 108)
print("GSE56284 DIFFERENTIAL-EXPRESSION SUMMARY")
print("=" * 108)
print(de_summary.to_string(index=False))

print("\n" + "=" * 108)
print("PUBLISHED VERSUS CURRENT DE COUNTS")
print("=" * 108)
print(
    published_comparison.to_string(
        index=False
    )
)

print("\n" + "=" * 108)
print("PRIORITY-GENE EVIDENCE IN GSE56284")
print("=" * 108)
print(priority.to_string(index=False))

print("\n" + "=" * 108)
print("CRITICAL DESIGN LIMITATION")
print("=" * 108)
print(design_context["interpretation"])

print(
    "\nPASS: focused GSE56284 technical, "
    "sample-level, DE-distribution, influence, "
    "and design audit completed."
)
