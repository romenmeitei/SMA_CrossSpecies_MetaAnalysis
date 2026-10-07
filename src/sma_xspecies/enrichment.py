from __future__ import annotations

from pathlib import Path

import pandas as pd


def run_preranked(ranking: pd.DataFrame, gene_sets: str | Path, output_dir: str | Path, permutations: int = 1000, seed: int = 20260814) -> pd.DataFrame:
    try:
        import gseapy as gp
    except ImportError as exc:
        raise RuntimeError("Install the enrichment extra: pip install -e '.[enrichment]'") from exc
    required = {"gene_symbol", "score"}
    if missing := required.difference(ranking.columns):
        raise ValueError(f"Ranking is missing {sorted(missing)}")
    ranked = ranking[["gene_symbol", "score"]].dropna().drop_duplicates("gene_symbol").sort_values("score", ascending=False)
    result = gp.prerank(
        rnk=ranked,
        gene_sets=str(gene_sets),
        outdir=str(output_dir),
        permutation_num=permutations,
        seed=seed,
        min_size=10,
        max_size=500,
        verbose=False,
    )
    table = result.res2d.reset_index()
    table.to_csv(Path(output_dir) / "preranked_gsea.csv", index=False)
    return table

