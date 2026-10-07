from pathlib import Path
from datetime import datetime, timezone
import hashlib
import json
import shutil
import sys

import matplotlib
matplotlib.use("Agg")

import matplotlib.image as mpimg
import matplotlib.pyplot as plt
from matplotlib.patches import (
    FancyArrowPatch,
    FancyBboxPatch,
)
import numpy as np
import pandas as pd


project = Path(sys.argv[1]).resolve()
q1_root = Path(sys.argv[2]).resolve()
raw_root = Path(sys.argv[3]).resolve()
asset_root = Path(sys.argv[4]).resolve()

main_fig = asset_root / "figures" / "main"
supp_fig = asset_root / "figures" / "supplementary"
source_fig = asset_root / "figures" / "source_panels"
main_tab = asset_root / "tables" / "main"
supp_tab = asset_root / "tables" / "supplementary"
audit_root = asset_root / "audit"

for directory in [
    main_fig,
    supp_fig,
    source_fig,
    main_tab,
    supp_tab,
    audit_root,
]:
    directory.mkdir(
        parents=True,
        exist_ok=True,
    )

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 10,
    "axes.titlesize": 12,
    "axes.labelsize": 10,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
})

NAVY = "#17365D"
BLUE = "#4472C4"
LIGHT_BLUE = "#D9EAF7"
RED = "#C00000"
GREEN = "#548235"
GOLD = "#C69214"
GRAY = "#666666"
LIGHT_GRAY = "#E7E6E6"

generated = []


def require(path):
    path = Path(path)

    if not path.is_file():
        raise FileNotFoundError(path)

    return path


def save_figure(fig, stem, folder):
    png = folder / f"{stem}.png"
    pdf = folder / f"{stem}.pdf"

    fig.savefig(
        png,
        dpi=600,
        bbox_inches="tight",
        facecolor="white",
    )

    fig.savefig(
        pdf,
        bbox_inches="tight",
        facecolor="white",
    )

    plt.close(fig)

    generated.extend([png, pdf])

    return png, pdf


def find_column(frame, candidates):
    normalized = {
        str(column).strip().lower(): column
        for column in frame.columns
    }

    for candidate in candidates:
        key = candidate.lower()

        if key in normalized:
            return normalized[key]

    return None


def add_panel_label(axis, label):
    axis.text(
        -0.07,
        1.04,
        label,
        transform=axis.transAxes,
        fontsize=15,
        fontweight="bold",
        va="bottom",
        ha="left",
    )


def copy_preserving(source, destination):
    source = require(source)
    destination = Path(destination)
    destination.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    shutil.copy2(source, destination)
    generated.append(destination)


# ------------------------------------------------------------------
# Figure 1
# ------------------------------------------------------------------

fig, ax = plt.subplots(
    figsize=(8.2, 10.2)
)

ax.set_xlim(0, 1)
ax.set_ylim(0, 1)
ax.axis("off")

ax.text(
    0.5,
    0.975,
    "Cross-species SMA transcriptomic analysis",
    ha="center",
    va="top",
    fontsize=17,
    fontweight="bold",
    color="black",
)

steps = [
    (
        "Raw-read reprocessing",
        "24 of 24 samples completed\n"
        "Salmon 2.7.0  Ensembl release 113",
        NAVY,
    ),
    (
        "Study-level inference",
        "GSE108094  human motor neuron\n"
        "GSE209926  mouse spinal cord\n"
        "GSE56284  mouse motor neuron\n"
        "GSE69175  descriptive only",
        BLUE,
    ),
    (
        "One-to-one ortholog synthesis",
        "13,693 cross-species orthologs\n"
        "12,519 represented in all three studies",
        "#5B9BD5",
    ),
    (
        "Corrected REML screening",
        "731 Tier 1 genes\n"
        "370 concordant low-heterogeneity genes\n"
        "9 GSE56284-independent genes\n"
        "1 Tier 4 gene  FAM107A",
        GREEN,
    ),
    (
        "Small-k inference",
        "Hartung-Knapp and modified Hartung-Knapp\n"
        "0 transcriptome-wide FDR discoveries",
        GOLD,
    ),
    (
        "Independent validation",
        "GSE51735 and GSE27205\n"
        "Directional support without corrected replication",
        RED,
    ),
]

y_positions = np.linspace(
    0.86,
    0.16,
    len(steps),
)

for index, (
    heading,
    details,
    color,
) in enumerate(steps):
    y = y_positions[index]

    box = FancyBboxPatch(
        (0.13, y - 0.055),
        0.74,
        0.105,
        boxstyle="round,pad=0.012,rounding_size=0.015",
        linewidth=1.5,
        edgecolor=color,
        facecolor="white",
    )

    ax.add_patch(box)

    ax.text(
        0.18,
        y + 0.036,
        heading,
        ha="left",
        va="top",
        fontsize=11.2,
        fontweight="bold",
        color=color,
    )

    ax.text(
        0.18,
        y + 0.006,
        details,
        ha="left",
        va="top",
        fontsize=8.5,
        color="black",
        linespacing=1.12,
    )

    if index < len(steps) - 1:
        arrow = FancyArrowPatch(
            (0.5, y - 0.06),
            (
                0.5,
                y_positions[index + 1] + 0.06,
            ),
            arrowstyle="-|>",
            mutation_scale=16,
            linewidth=1.5,
            color=GRAY,
        )

        ax.add_patch(arrow)

ax.text(
    0.5,
    0.055,
    "Primary conclusion",
    ha="center",
    va="center",
    fontsize=11,
    fontweight="bold",
)

ax.text(
    0.5,
    0.023,
    "FAM107A is the leading internally consistent candidate; "
    "UGP2, splicing and treatment-pathway results remain secondary.",
    ha="center",
    va="center",
    fontsize=9.5,
)

save_figure(
    fig,
    "Figure_1_Study_design_and_evidence_hierarchy",
    main_fig,
)


# ------------------------------------------------------------------
# Figure 2
# ------------------------------------------------------------------

forest_root = (
    q1_root
    / "02_forest_plots"
    / "individual_png"
)

forest_paths = [
    require(
        forest_root
        / "FAM107A_ENSG00000168309_forest_plot.png"
    ),
    require(
        forest_root
        / "UGP2_ENSG00000169764_forest_plot.png"
    ),
]

fig, axes = plt.subplots(
    1,
    2,
    figsize=(14, 6.5),
)

for axis, path, label in zip(
    axes,
    forest_paths,
    ["A", "B"],
):
    image = mpimg.imread(path)
    axis.imshow(image)
    axis.axis("off")
    add_panel_label(axis, label)

fig.suptitle(
    "Study-level effects and small-k meta-analysis",
    fontsize=16,
    fontweight="bold",
    y=0.99,
)

fig.text(
    0.5,
    0.015,
    "Positive estimates indicate higher expression in SMA. "
    "The individual forest plots are reproduced from the "
    "verified small-k analysis.",
    ha="center",
    fontsize=9.5,
)

fig.tight_layout(
    rect=(0, 0.035, 1, 0.96)
)

save_figure(
    fig,
    "Figure_2_FAM107A_UGP2_forest_plots",
    main_fig,
)


# Preserve individual source panels.
for source in forest_paths:
    copy_preserving(
        source,
        source_fig / source.name,
    )

copy_preserving(
    q1_root
    / "02_forest_plots"
    / "priority_genes_study_level_forest_plots.pdf",
    source_fig
    / "priority_genes_study_level_forest_plots.pdf",
)


# ------------------------------------------------------------------
# Figure 3
# ------------------------------------------------------------------

gse_fig_root = (
    q1_root
    / "03_GSE56284_investigation"
    / "figures"
)

gse_panels = [
    (
        require(
            gse_fig_root
            / "GSE56284_PCA_top1000.png"
        ),
        "A",
        "Principal component analysis",
    ),
    (
        require(
            gse_fig_root
            / "GSE56284_sample_correlation.png"
        ),
        "B",
        "Sample correlation",
    ),
    (
        require(
            gse_fig_root
            / "GSE56284_DE_distributions.png"
        ),
        "C",
        "Differential-expression distributions",
    ),
]

fig, axes = plt.subplots(
    1,
    3,
    figsize=(17, 5.8),
)

for axis, (
    path,
    label,
    title,
) in zip(axes, gse_panels):
    image = mpimg.imread(path)
    axis.imshow(image)
    axis.axis("off")
    axis.set_title(
        title,
        fontsize=11,
        fontweight="bold",
        pad=8,
    )
    add_panel_label(axis, label)

fig.suptitle(
    "GSE56284 technical quality and dataset influence",
    fontsize=16,
    fontweight="bold",
)

fig.text(
    0.5,
    0.015,
    "The dataset shows strong internal coherence, but genotype "
    "is perfectly confounded with cell-line identity.",
    ha="center",
    fontsize=9.5,
)

fig.tight_layout(
    rect=(0, 0.04, 1, 0.93)
)

save_figure(
    fig,
    "Figure_3_GSE56284_quality_and_influence",
    main_fig,
)

for path, _, _ in gse_panels:
    copy_preserving(
        path,
        source_fig / path.name,
    )


# ------------------------------------------------------------------
# Figure 4
# ------------------------------------------------------------------

external_specs = [
    (
        "GSE51735",
        "Deposited log2 RPKM",
        "disease-control",
        q1_root
        / "04_external_validation"
        / "results"
        / "GSE51735_motor_neuron_disease-control_limma.csv",
        RED,
    ),
    (
        "GSE51735",
        "Quantile sensitivity",
        "disease-control",
        q1_root
        / "04_external_validation"
        / "results"
        / "GSE51735_motor_neuron_disease-control_quantile_sensitivity.csv",
        "#ED7D31",
    ),
    (
        "GSE27205",
        "Disease minus control",
        "disease-control",
        q1_root
        / "04_external_validation"
        / "results"
        / "GSE27205_disease-control_limma.csv",
        RED,
    ),
    (
        "GSE27205",
        "Corrected minus disease",
        "corrected-disease",
        q1_root
        / "04_external_validation"
        / "results"
        / "GSE27205_corrected-disease_limma.csv",
        GREEN,
    ),
]

external_records = []

for (
    accession,
    model_label,
    contrast,
    path,
    color,
) in external_specs:
    frame = pd.read_csv(require(path))

    symbol_column = find_column(
        frame,
        [
            "gene_symbol",
            "symbol",
            "gene",
        ],
    )

    effect_column = find_column(
        frame,
        [
            "effect",
            "logfc",
            "log2foldchange",
        ],
    )

    t_column = find_column(
        frame,
        [
            "moderated_t",
            "t",
            "stat",
        ],
    )

    p_column = find_column(
        frame,
        [
            "p_value",
            "p.value",
            "pvalue",
        ],
    )

    fdr_column = find_column(
        frame,
        [
            "fdr",
            "adj.p.val",
            "padj",
        ],
    )

    if symbol_column is None or effect_column is None:
        raise RuntimeError(
            f"Could not identify gene/effect columns in {path}"
        )

    normalized_symbols = (
        frame[symbol_column]
        .astype(str)
        .str.upper()
        .str.strip()
    )

    for gene in ["FAM107A", "UGP2"]:
        selected = frame.loc[
            normalized_symbols.eq(gene)
        ].copy()

        if selected.empty:
            continue

        row = selected.iloc[0]

        effect = pd.to_numeric(
            pd.Series([row[effect_column]]),
            errors="coerce",
        ).iloc[0]

        moderated_t = (
            pd.to_numeric(
                pd.Series([row[t_column]]),
                errors="coerce",
            ).iloc[0]
            if t_column
            else np.nan
        )

        if (
            np.isfinite(effect)
            and np.isfinite(moderated_t)
            and abs(moderated_t) > 1e-12
        ):
            moderated_se = abs(
                effect / moderated_t
            )
        else:
            moderated_se = np.nan

        external_records.append({
            "accession": accession,
            "model": model_label,
            "contrast": contrast,
            "gene_symbol": gene,
            "effect": effect,
            "moderated_t": moderated_t,
            "moderated_se": moderated_se,
            "p_value": (
                row[p_column]
                if p_column
                else np.nan
            ),
            "fdr": (
                row[fdr_column]
                if fdr_column
                else np.nan
            ),
            "color": color,
        })

external = pd.DataFrame(external_records)

external.to_csv(
    main_tab
    / "Table_4_FAM107A_UGP2_external_validation.csv",
    index=False,
)

fig, axes = plt.subplots(
    1,
    2,
    figsize=(12.5, 5.7),
    sharex=False,
)

for axis, gene, label in zip(
    axes,
    ["FAM107A", "UGP2"],
    ["A", "B"],
):
    subset = external.loc[
        external["gene_symbol"].eq(gene)
    ].copy()

    subset["display"] = (
        subset["accession"]
        + "  "
        + subset["model"]
    )

    positions = np.arange(len(subset))

    for position, (_, row) in zip(
        positions,
        subset.iterrows(),
    ):
        xerr = (
            1.96 * row["moderated_se"]
            if np.isfinite(row["moderated_se"])
            else None
        )

        axis.errorbar(
            row["effect"],
            position,
            xerr=xerr,
            fmt="o",
            color=row["color"],
            ecolor=row["color"],
            markersize=7,
            capsize=3,
            linewidth=1.5,
        )

        p_value = pd.to_numeric(
            pd.Series([row["p_value"]]),
            errors="coerce",
        ).iloc[0]

        if np.isfinite(p_value):
            axis.text(
                row["effect"],
                position - 0.20,
                f"p = {p_value:.3g}",
                ha="center",
                va="bottom",
                fontsize=8,
                color=GRAY,
            )

    axis.axvline(
        0,
        color="black",
        linewidth=1,
        linestyle="--",
    )

    axis.set_yticks(positions)
    axis.set_yticklabels(
        subset["display"],
        fontsize=8.5,
    )
    axis.set_ylim(
        len(subset) - 0.50,
        -0.50,
    )
    axis.set_xlabel("Expression effect")
    axis.set_title(
        gene,
        fontweight="bold",
    )
    axis.grid(
        axis="x",
        alpha=0.2,
    )
    add_panel_label(axis, label)

fig.suptitle(
    "External directional validation of priority genes",
    fontsize=16,
    fontweight="bold",
)

fig.text(
    0.5,
    0.01,
    "Error bars show approximate 95% intervals derived from "
    "the moderated t statistics. No priority result passed "
    "multiplicity-corrected external validation.",
    ha="center",
    fontsize=9,
)

fig.tight_layout(
    rect=(0, 0.05, 1, 0.93)
)

save_figure(
    fig,
    "Figure_4_External_validation_priority_genes",
    main_fig,
)


# Preserve external PCA figures as supplementary panels.
external_pca_root = (
    q1_root
    / "04_external_validation"
    / "figures"
)

for filename in [
    "GSE51735_motor_neuron_PCA.png",
    "GSE51735_motor_neuron_PCA.pdf",
    "GSE27205_motor_neuron_PCA.png",
    "GSE27205_motor_neuron_PCA.pdf",
]:
    copy_preserving(
        external_pca_root / filename,
        supp_fig / filename,
    )


# ------------------------------------------------------------------
# Figure 5
# ------------------------------------------------------------------

splicing_root = project / "results" / "splicing"
ora_root = (
    project
    / "results"
    / "enrichment"
    / "GSE108094_splicing_ORA"
)

events = pd.read_csv(
    require(
        splicing_root
        / "GSE108094_significant_events.csv"
    )
)

brain_genes = pd.read_csv(
    require(
        ora_root
        / "GSE108094_Brain_Development_splicing_gene_evidence.csv"
    )
)

ugp2 = pd.read_csv(
    require(
        ora_root
        / "UGP2_splicing_event_evidence.csv"
    )
)

ora = pd.read_csv(
    require(
        ora_root
        / "GSE108094_splicing_ORA_all_results_with_global_FDR.csv"
    )
)

delta_column = find_column(
    events,
    ["delta_psi"],
)

if delta_column is None:
    raise RuntimeError(
        "delta_psi is missing from significant events."
    )

delta = pd.to_numeric(
    events[delta_column],
    errors="coerce",
)

higher_count = int((delta > 0).sum())
lower_count = int((delta < 0).sum())

gene_delta_column = find_column(
    brain_genes,
    [
        "largest_effect_delta_psi",
        "delta_psi",
    ],
)

gene_symbol_column = find_column(
    brain_genes,
    [
        "gene_symbol",
        "symbol",
    ],
)

if (
    gene_delta_column is None
    or gene_symbol_column is None
):
    raise RuntimeError(
        "Brain Development evidence columns were not found."
    )

brain_genes[gene_delta_column] = pd.to_numeric(
    brain_genes[gene_delta_column],
    errors="coerce",
)

brain_plot = (
    brain_genes
    .dropna(
        subset=[
            gene_symbol_column,
            gene_delta_column,
        ]
    )
    .sort_values(gene_delta_column)
    .copy()
)

term_column = find_column(
    ora,
    ["term"],
)

global_fdr_column = find_column(
    ora,
    [
        "global_fdr_all_eight_analyses",
        "global_fdr",
    ],
)

enrichment_column = find_column(
    ora,
    ["enrichment_ratio"],
)

brain_row = None

if term_column:
    matches = ora.loc[
        ora[term_column]
        .astype(str)
        .str.contains(
            "Brain Development",
            case=False,
            na=False,
        )
    ]

    if not matches.empty:
        brain_row = matches.iloc[0]

fig = plt.figure(
    figsize=(15, 6.7)
)

grid = fig.add_gridspec(
    1,
    3,
    width_ratios=[0.9, 1.75, 1.0],
    wspace=0.42,
)

ax1 = fig.add_subplot(grid[0, 0])
ax2 = fig.add_subplot(grid[0, 1])
ax3 = fig.add_subplot(grid[0, 2])

# Panel A
bars = ax1.bar(
    ["Higher\nin SMA", "Lower\nin SMA"],
    [higher_count, lower_count],
    color=[BLUE, RED],
    edgecolor="black",
    linewidth=0.6,
)

for bar, value in zip(
    bars,
    [higher_count, lower_count],
):
    ax1.text(
        bar.get_x() + bar.get_width() / 2,
        bar.get_height() + 18,
        f"{value:,}",
        ha="center",
        va="bottom",
        fontweight="bold",
    )

ax1.set_ylabel("Significant skipped-exon events")
ax1.set_title(
    "Event direction",
    fontweight="bold",
)
ax1.set_ylim(
    0,
    max(higher_count, lower_count) * 1.18,
)
add_panel_label(ax1, "A")

# Panel B
y = np.arange(len(brain_plot))
values = brain_plot[
    gene_delta_column
].to_numpy(dtype=float)

colors = np.where(
    values >= 0,
    BLUE,
    RED,
)

for index, value in enumerate(values):
    ax2.plot(
        [0, value],
        [index, index],
        color=colors[index],
        linewidth=1.5,
    )

ax2.scatter(
    values,
    y,
    color=colors,
    s=35,
    edgecolor="black",
    linewidth=0.4,
    zorder=3,
)

ax2.axvline(
    0,
    color="black",
    linewidth=0.9,
)

ax2.set_yticks(y)
ax2.set_yticklabels(
    brain_plot[gene_symbol_column],
    fontsize=8,
)

ax2.set_xlabel("Largest event delta PSI")
ax2.set_title(
    "Brain Development overlap",
    fontweight="bold",
)

if brain_row is not None:
    global_fdr = (
        pd.to_numeric(
            pd.Series([
                brain_row[global_fdr_column]
            ]),
            errors="coerce",
        ).iloc[0]
        if global_fdr_column
        else np.nan
    )

    enrichment_ratio = (
        pd.to_numeric(
            pd.Series([
                brain_row[enrichment_column]
            ]),
            errors="coerce",
        ).iloc[0]
        if enrichment_column
        else np.nan
    )

    annotation = []

    if np.isfinite(enrichment_ratio):
        annotation.append(
            f"enrichment ratio = {enrichment_ratio:.2f}"
        )

    if np.isfinite(global_fdr):
        annotation.append(
            f"global FDR = {global_fdr:.3g}"
        )

    if annotation:
        ax2.text(
            0.02,
            0.01,
            "\n".join(annotation),
            transform=ax2.transAxes,
            ha="left",
            va="bottom",
            fontsize=8.5,
            color=GRAY,
        )

add_panel_label(ax2, "B")

# Panel C
if ugp2.empty:
    raise RuntimeError(
        "UGP2 splicing event table is empty."
    )

ugp2_row = ugp2.iloc[0]

control_column = find_column(
    ugp2,
    ["mean_psi_control"],
)

disease_column = find_column(
    ugp2,
    ["mean_psi_disease"],
)

ugp2_delta_column = find_column(
    ugp2,
    ["delta_psi"],
)

ugp2_fdr_column = find_column(
    ugp2,
    ["fdr"],
)

if control_column is None or disease_column is None:
    raise RuntimeError(
        "UGP2 mean PSI columns are missing."
    )

control_psi = float(ugp2_row[control_column])
disease_psi = float(ugp2_row[disease_column])

bars = ax3.bar(
    ["Control", "SMA"],
    [control_psi, disease_psi],
    color=[BLUE, RED],
    edgecolor="black",
    linewidth=0.7,
    width=0.62,
)

for bar, value in zip(
    bars,
    [control_psi, disease_psi],
):
    ax3.text(
        bar.get_x() + bar.get_width() / 2,
        value + 0.015,
        f"{value:.3f}",
        ha="center",
        va="bottom",
        fontweight="bold",
    )

ugp2_delta = (
    float(ugp2_row[ugp2_delta_column])
    if ugp2_delta_column
    else disease_psi - control_psi
)

ugp2_fdr = (
    float(ugp2_row[ugp2_fdr_column])
    if ugp2_fdr_column
    else np.nan
)

ax3.set_ylabel("Mean exon inclusion PSI")
ax3.set_ylim(
    0,
    max(control_psi, disease_psi) + 0.18,
)

ax3.set_title(
    "UGP2 skipped-exon event",
    fontweight="bold",
)

annotation = f"delta PSI = {ugp2_delta:.3f}"

if np.isfinite(ugp2_fdr):
    annotation += f"\nFDR = {ugp2_fdr:.3g}"

ax3.text(
    0.5,
    0.93,
    annotation,
    transform=ax3.transAxes,
    ha="center",
    va="top",
    fontsize=9,
)

add_panel_label(ax3, "C")

fig.suptitle(
    "GSE108094 skipped-exon changes and neural-development enrichment",
    fontsize=16,
    fontweight="bold",
)

fig.tight_layout(
    rect=(0, 0.02, 1, 0.93)
)

save_figure(
    fig,
    "Figure_5_GSE108094_splicing_and_UGP2",
    main_fig,
)


# ------------------------------------------------------------------
# Supplementary AMPA figure
# ------------------------------------------------------------------

ampa_path = (
    project
    / "results"
    / "enrichment"
    / "GSE207890_pathway_reversal"
    / "GSE207890_nonredundant_pathway_reversal_theme.csv"
)

ampa = pd.read_csv(require(ampa_path))

if ampa.empty:
    raise RuntimeError(
        "AMPA pathway theme table is empty."
    )

row = ampa.iloc[0]

ampa_records = [
    {
        "model": "Primary blocked-fiber",
        "contrast": "Disease minus control",
        "nes": row["disease_primary_NES"],
        "fdr": row["disease_primary_FDR"],
        "color": RED,
    },
    {
        "model": "Primary blocked-fiber",
        "contrast": "Treated minus disease",
        "nes": row["treatment_primary_NES"],
        "fdr": row["treatment_primary_FDR"],
        "color": GREEN,
    },
    {
        "model": "Mouse-level sensitivity",
        "contrast": "Disease minus control",
        "nes": row["disease_mouse_NES"],
        "fdr": row["disease_mouse_FDR"],
        "color": RED,
    },
    {
        "model": "Mouse-level sensitivity",
        "contrast": "Treated minus disease",
        "nes": row["treatment_mouse_NES"],
        "fdr": row["treatment_mouse_FDR"],
        "color": GREEN,
    },
]

ampa_plot = pd.DataFrame(ampa_records)

ampa_plot.to_csv(
    main_tab
    / "Table_5_AMPA_pathway_reversal_summary.csv",
    index=False,
)

fig, ax = plt.subplots(
    figsize=(9.5, 5.6)
)

labels = [
    (
        f'{record["model"]}\n'
        f'{record["contrast"]}'
    )
    for record in ampa_records
]

positions = np.arange(
    len(ampa_records)
)

bars = ax.barh(
    positions,
    ampa_plot["nes"],
    color=ampa_plot["color"],
    edgecolor="black",
    linewidth=0.6,
)

ax.axvline(
    0,
    color="black",
    linewidth=1,
)

ax.set_yticks(positions)
ax.set_yticklabels(
    labels,
    fontsize=9,
)

ax.invert_yaxis()
ax.set_xlabel("Normalized enrichment score")
ax.set_title(
    "AMPA receptor trafficking pathway reversal",
    fontsize=15,
    fontweight="bold",
)

for bar, fdr in zip(
    bars,
    ampa_plot["fdr"],
):
    value = bar.get_width()

    x_position = (
        value - 0.08
        if value >= 0
        else value + 0.08
    )

    alignment = (
        "right"
        if value >= 0
        else "left"
    )

    ax.text(
        x_position,
        bar.get_y() + bar.get_height() / 2,
        f"FDR = {float(fdr):.3g}",
        ha=alignment,
        va="center",
        fontsize=9,
        color="white",
        fontweight="bold",
    )

limit = max(
    abs(ampa_plot["nes"].min()),
    abs(ampa_plot["nes"].max()),
) + 0.65

ax.set_xlim(-limit, limit)
ax.grid(
    axis="x",
    alpha=0.2,
)

fig.text(
    0.5,
    0.02,
    "The treatment contrasts meet the prespecified secondary "
    "FDR < 0.10 threshold; no gene-level reversal passed "
    "the conjunction test.",
    ha="center",
    fontsize=9,
)

fig.tight_layout(
    rect=(0, 0.055, 1, 0.94)
)

save_figure(
    fig,
    "Supplementary_Figure_S1_AMPA_pathway_reversal",
    supp_fig,
)


# ------------------------------------------------------------------
# Main tables
# ------------------------------------------------------------------

study_table = pd.DataFrame([
    {
        "accession": "GSE108094",
        "species": "Homo sapiens",
        "tissue": "motor neuron",
        "control_biological_units": 2,
        "disease_biological_units": 2,
        "tested_genes": 19551,
        "FDR_le_0.05": 383,
        "analysis_role": "eligible discovery",
    },
    {
        "accession": "GSE209926",
        "species": "Mus musculus",
        "tissue": "spinal cord",
        "control_biological_units": 3,
        "disease_biological_units": 3,
        "tested_genes": 16254,
        "FDR_le_0.05": 15,
        "analysis_role": "eligible discovery",
    },
    {
        "accession": "GSE56284",
        "species": "Mus musculus",
        "tissue": "motor neuron",
        "control_biological_units": 3,
        "disease_biological_units": 3,
        "tested_genes": 16884,
        "FDR_le_0.05": 11697,
        "analysis_role": "eligible discovery",
    },
    {
        "accession": "GSE69175",
        "species": "Homo sapiens",
        "tissue": "motor neuron",
        "control_biological_units": 1,
        "disease_biological_units": 1,
        "tested_genes": 14062,
        "FDR_le_0.05": np.nan,
        "analysis_role": "descriptive sensitivity only",
    },
])

study_path = (
    main_tab
    / "Table_1_Study_design_and_DE_summary.csv"
)

study_table.to_csv(
    study_path,
    index=False,
)

generated.append(study_path)

tier_table = pd.DataFrame([
    {
        "tier": "Tier 1",
        "definition": (
            "Corrected full neural REML "
            "meta-analysis FDR <= 0.05"
        ),
        "genes": 731,
    },
    {
        "tier": "Tier 2",
        "definition": (
            "Tier 1 plus all-study direction "
            "agreement and I2 < 75%"
        ),
        "genes": 370,
    },
    {
        "tier": "Tier 3",
        "definition": (
            "Full-model significance, I2 < 75%, "
            "and significance after excluding GSE56284"
        ),
        "genes": 9,
    },
    {
        "tier": "Tier 4",
        "definition": (
            "Tier 2 plus significant concordant "
            "results in both available cross-species "
            "leave-one-dataset-out analyses"
        ),
        "genes": 1,
    },
])

tier_path = (
    main_tab
    / "Table_2_Corrected_REML_evidence_hierarchy.csv"
)

tier_table.to_csv(
    tier_path,
    index=False,
)

generated.append(tier_path)

robust_path = (
    project
    / "results"
    / "meta"
    / "raw_Ensembl113_robustness"
    / "FAM107A_UGP2_corrected_REML_robustness.csv"
)

small_k_path = (
    q1_root
    / "01_small_k_meta"
    / "priority_gene_exact_small_k_statistics.csv"
)

robust = pd.read_csv(require(robust_path))
small_k = pd.read_csv(require(small_k_path))

join_keys = [
    column
    for column in [
        "ortholog_id",
        "gene_symbol",
    ]
    if (
        column in robust.columns
        and column in small_k.columns
    )
]

if not join_keys:
    raise RuntimeError(
        "No common priority-gene join key was found."
    )

priority = robust.merge(
    small_k,
    on=join_keys,
    how="left",
    suffixes=(
        "_corrected_REML",
        "_small_k",
    ),
)

priority_path = (
    main_tab
    / "Table_3_FAM107A_UGP2_primary_evidence.csv"
)

priority.to_csv(
    priority_path,
    index=False,
)

generated.append(priority_path)

secondary_table = pd.DataFrame([
    {
        "finding": "GSE108094 skipped-exon events",
        "result": "1,069 events in 844 genes",
        "evidence_role": "secondary",
    },
    {
        "finding": "Brain Development splicing ORA",
        "result": (
            "20 genes; enrichment ratio 3.12; "
            "global FDR 0.0422"
        ),
        "evidence_role": "secondary",
    },
    {
        "finding": "UGP2 skipped-exon event",
        "result": (
            "delta PSI -0.174; FDR 0.0298"
        ),
        "evidence_role": "secondary cross-modal",
    },
    {
        "finding": "AMPA receptor trafficking",
        "result": (
            "Disease FDR < 0.05 and treatment "
            "reversal FDR < 0.10 in both models"
        ),
        "evidence_role": "exploratory secondary",
    },
])

secondary_path = (
    main_tab
    / "Table_5_Secondary_findings_summary.csv"
)

secondary_table.to_csv(
    secondary_path,
    index=False,
)

generated.append(secondary_path)


# ------------------------------------------------------------------
# Excel workbook
# ------------------------------------------------------------------

workbook_path = (
    main_tab
    / "SMA_manuscript_main_tables.xlsx"
)

with pd.ExcelWriter(
    workbook_path,
    engine="openpyxl",
) as writer:
    study_table.to_excel(
        writer,
        sheet_name="Study summary",
        index=False,
    )

    tier_table.to_excel(
        writer,
        sheet_name="Evidence hierarchy",
        index=False,
    )

    priority.to_excel(
        writer,
        sheet_name="Priority genes",
        index=False,
    )

    external.drop(
        columns=["color"],
        errors="ignore",
    ).to_excel(
        writer,
        sheet_name="External validation",
        index=False,
    )

    secondary_table.to_excel(
        writer,
        sheet_name="Secondary findings",
        index=False,
    )

generated.append(workbook_path)


# ------------------------------------------------------------------
# Supplementary source tables
# ------------------------------------------------------------------

supplementary_sources = {
    (
        q1_root
        / "01_small_k_meta"
        / "small_k_meta_all_genes.csv"
    ): "Supplementary_Table_S1_small_k_all_genes.csv",

    (
        project
        / "results"
        / "meta"
        / "gene_meta_neural_joint.csv"
    ): "Supplementary_Table_S2_corrected_REML_all_genes.csv",

    (
        project
        / "results"
        / "meta"
        / "leave_one_dataset_out.csv"
    ): "Supplementary_Table_S3_leave_one_dataset_out.csv",

    (
        project
        / "results"
        / "meta"
        / "raw_Ensembl113_robustness"
        / "all_genes_robustness_matrix.csv"
    ): "Supplementary_Table_S4_robustness_matrix.csv",

    (
        q1_root
        / "04_external_validation"
        / "results"
        / "priority_gene_external_validation_evidence.csv"
    ): "Supplementary_Table_S5_external_priority_genes.csv",

    (
        splicing_root
        / "GSE108094_significant_events.csv"
    ): "Supplementary_Table_S6_significant_SE_events.csv",

    (
        ora_root
        / "GSE108094_splicing_ORA_all_results_with_global_FDR.csv"
    ): "Supplementary_Table_S7_splicing_ORA.csv",

    ampa_path:
        "Supplementary_Table_S8_AMPA_pathway_theme.csv",
}

for source, filename in supplementary_sources.items():
    copy_preserving(
        source,
        supp_tab / filename,
    )


# ------------------------------------------------------------------
# Figure legends
# ------------------------------------------------------------------

legends = """FIGURE LEGENDS

Figure 1
Cross-species SMA transcriptomic analysis and evidence hierarchy.
Twenty-four samples from four studies underwent harmonized raw-read
reprocessing. Three eligible neural studies contributed to the primary
one-to-one ortholog meta-analysis. Corrected REML screening was followed
by Hartung-Knapp and modified Hartung-Knapp small-k sensitivity
analysis, internal influence testing and external directional
validation.

Figure 2
Study-level effects for FAM107A and UGP2. Points and intervals show
study-level disease-minus-control estimates and the corresponding
meta-analytic summaries reproduced from the verified small-k forest
plots. Positive effects indicate higher expression in SMA.

Figure 3
Technical quality and influence assessment of GSE56284. The panels show
principal-component separation, normalized-expression correlation and
the distribution of study-level differential-expression statistics.
The signal was internally coherent, although genotype was perfectly
confounded with cell-line identity.

Figure 4
External directional validation of FAM107A and UGP2 in GSE51735 and
GSE27205. Positive disease-minus-control values indicate higher
expression in SMA. Negative corrected-minus-disease values indicate
directional reversal after gene correction. None of the priority-gene
tests achieved multiplicity-corrected external replication.

Figure 5
GSE108094 skipped-exon changes. Panel A summarizes significant event
directions. Panel B shows the largest delta-PSI event for the 20 Brain
Development overlap genes. Panel C shows the UGP2 skipped-exon event,
which had lower inclusion in SMA.

Supplementary Figure S1
AMPA receptor trafficking pathway reversal in GSE207890. Disease
contrasts showed negative enrichment, whereas treated-minus-disease
contrasts showed positive enrichment in the primary blocked-fiber and
mouse-level models. Treatment evidence met the prespecified secondary
FDR threshold of 0.10 but not gene-level reversal criteria.
"""

legend_path = (
    asset_root
    / "figures"
    / "figure_legends.txt"
)

legend_path.write_text(
    legends,
    encoding="utf-8",
)

generated.append(legend_path)


# ------------------------------------------------------------------
# Manifest and checksums
# ------------------------------------------------------------------

manifest_records = []

for path in sorted(
    set(
        Path(path).resolve()
        for path in generated
        if Path(path).is_file()
    )
):
    digest = hashlib.sha256(
        path.read_bytes()
    ).hexdigest()

    manifest_records.append({
        "relative_path": str(
            path.relative_to(asset_root)
        ),
        "bytes": path.stat().st_size,
        "sha256": digest,
    })

manifest = pd.DataFrame(
    manifest_records
)

manifest_path = (
    audit_root
    / "generated_manuscript_asset_manifest.csv"
)

manifest.to_csv(
    manifest_path,
    index=False,
)

checksum_path = (
    audit_root
    / "generated_manuscript_assets.sha256"
)

with checksum_path.open(
    "w",
    encoding="utf-8",
) as handle:
    for record in manifest_records:
        handle.write(
            f'{record["sha256"]}  '
            f'{record["relative_path"]}\n'
        )

metadata = {
    "created_utc": (
        datetime.now(timezone.utc)
        .isoformat()
    ),
    "project_root": str(project),
    "q1_root": str(q1_root),
    "raw_root": str(raw_root),
    "asset_root": str(asset_root),
    "generated_files": len(manifest_records),
}

metadata_path = (
    audit_root
    / "manuscript_asset_metadata.json"
)

metadata_path.write_text(
    json.dumps(
        metadata,
        indent=2,
    )
    + "\n",
    encoding="utf-8",
)

print("=" * 100)
print("MANUSCRIPT ASSETS GENERATED")
print("=" * 100)

print("\nMain figures:")

for path in sorted(main_fig.iterdir()):
    if path.is_file():
        print(
            path.name,
            f"{path.stat().st_size / 1024:.1f} KB",
        )

print("\nSupplementary figures:")

for path in sorted(supp_fig.iterdir()):
    if path.is_file():
        print(
            path.name,
            f"{path.stat().st_size / 1024:.1f} KB",
        )

print("\nMain tables:")

for path in sorted(main_tab.iterdir()):
    if path.is_file():
        print(
            path.name,
            f"{path.stat().st_size / 1024:.1f} KB",
        )

print("\nManifest:", manifest_path)
print("Checksums:", checksum_path)

print(
    "\nPASS: publication assets generated "
    "without modifying the frozen release."
)
