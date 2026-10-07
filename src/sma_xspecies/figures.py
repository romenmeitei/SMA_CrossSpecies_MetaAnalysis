from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd


def _style() -> None:
    import matplotlib.pyplot as plt
    import seaborn as sns

    sns.set_theme(style="whitegrid", context="paper")
    plt.rcParams.update({"figure.dpi": 150, "savefig.dpi": 300, "font.family": "DejaVu Sans"})


def plot_meta_summary(table: pd.DataFrame, destination: str | Path, top_n: int = 25) -> None:
    import matplotlib.pyplot as plt

    data = table.dropna(subset=["estimate", "ci_low", "ci_high"]).sort_values(["fdr", "p_value"]).head(top_n).copy()
    if data.empty:
        return
    data = data.sort_values("estimate")
    y = np.arange(len(data))
    fig, ax = plt.subplots(figsize=(7.2, max(3.5, 0.28 * len(data))))
    ax.errorbar(data["estimate"], y, xerr=[data["estimate"] - data["ci_low"], data["ci_high"] - data["estimate"]], fmt="o", color="#26619C", ecolor="#8AA7C4", capsize=2)
    ax.axvline(0, color="black", lw=0.8, ls="--")
    ax.set_yticks(y, data["gene_symbol"].fillna(data["ortholog_id"]))
    ax.set_xlabel("Random-effects pooled log2 fold change")
    ax.set_ylabel("")
    ax.set_title("Top cross-species meta-analysis effects")
    fig.tight_layout()
    fig.savefig(destination, bbox_inches="tight")
    plt.close(fig)


def plot_reversal(table: pd.DataFrame, destination: str | Path) -> None:
    import matplotlib.pyplot as plt

    data = table.dropna(subset=["disease_effect", "treatment_effect"]).copy()
    if data.empty:
        return
    highlighted = data.get("passes_conjunction", pd.Series(False, index=data.index)).astype(bool)
    fig, ax = plt.subplots(figsize=(6, 5.2))
    ax.scatter(data.loc[~highlighted, "disease_effect"], data.loc[~highlighted, "treatment_effect"], s=8, alpha=0.25, color="#7A7A7A", rasterized=True)
    ax.scatter(data.loc[highlighted, "disease_effect"], data.loc[highlighted, "treatment_effect"], s=16, alpha=0.85, color="#B22222", label="passes conjunction")
    ax.axhline(0, color="black", lw=0.7)
    ax.axvline(0, color="black", lw=0.7)
    ax.set_xlabel("Disease − control log2 fold change")
    ax.set_ylabel("Treated − disease log2 fold change")
    ax.set_title("GSE207890 treatment-reversal geometry")
    if highlighted.any():
        ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(destination, bbox_inches="tight")
    plt.close(fig)


def plot_splicing(table: pd.DataFrame, destination: str | Path, fdr: float = 0.05, delta: float = 0.10) -> None:
    import matplotlib.pyplot as plt

    data = table.dropna(subset=["delta_psi", "fdr"]).copy()
    if data.empty:
        return
    data["minus_log10_fdr"] = -np.log10(data["fdr"].clip(lower=np.finfo(float).tiny))
    selected = (data["fdr"] <= fdr) & (data["delta_psi"].abs() >= delta)
    fig, ax = plt.subplots(figsize=(6, 5.2))
    ax.scatter(data.loc[~selected, "delta_psi"], data.loc[~selected, "minus_log10_fdr"], s=8, alpha=0.25, color="#7A7A7A", rasterized=True)
    ax.scatter(data.loc[selected, "delta_psi"], data.loc[selected, "minus_log10_fdr"], s=14, alpha=0.8, color="#5B2C6F")
    ax.axvline(-delta, color="black", lw=0.7, ls="--")
    ax.axvline(delta, color="black", lw=0.7, ls="--")
    ax.axhline(-np.log10(fdr), color="black", lw=0.7, ls="--")
    ax.set_xlabel("ΔPSI")
    ax.set_ylabel("−log10(FDR)")
    ax.set_title("GSE108094 exon-skipping events")
    fig.tight_layout()
    fig.savefig(destination, bbox_inches="tight")
    plt.close(fig)


def generate_available_figures(root: str | Path) -> list[Path]:
    root = Path(root)
    output = root / "figures"
    output.mkdir(parents=True, exist_ok=True)
    _style()
    generated: list[Path] = []
    jobs = [
        (root / "results/meta/gene_meta_by_tissue.csv", output / "Figure_meta_effects.pdf", plot_meta_summary),
        (root / "results/rescue/GSE207890_treatment_reversal_primary.csv", output / "Figure_treatment_reversal.pdf", plot_reversal),
        (root / "results/splicing/GSE108094_all_events.csv", output / "Figure_splicing_volcano.pdf", plot_splicing),
    ]
    for source, destination, function in jobs:
        if source.exists() and source.stat().st_size > 0:
            function(pd.read_csv(source), destination)
            if destination.exists():
                generated.append(destination)
    return generated

