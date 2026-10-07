from __future__ import annotations

import numpy as np
import pandas as pd


def treatment_reversal(
    disease: pd.DataFrame,
    treatment: pd.DataFrame,
    *,
    fdr_threshold: float = 0.10,
    minimum_absolute_z: float = 1.0,
) -> pd.DataFrame:
    keys = [key for key in ("gene_id", "gene_symbol") if key in disease.columns and key in treatment.columns]
    if "gene_id" not in keys:
        raise ValueError("Both contrasts require gene_id")
    left = disease[keys + ["effect", "z", "fdr"]].rename(columns={"effect": "disease_effect", "z": "disease_z", "fdr": "disease_fdr"})
    right = treatment[keys + ["effect", "z", "fdr"]].rename(columns={"effect": "treatment_effect", "z": "treatment_z", "fdr": "treatment_fdr"})
    out = left.merge(right, on=keys, how="inner")
    out["reversal_score"] = -out["disease_z"] * out["treatment_z"]
    out["opposite_direction"] = out["disease_effect"] * out["treatment_effect"] < 0
    out["restoration_fraction"] = -out["treatment_effect"] / out["disease_effect"].replace(0, np.nan)
    out["passes_conjunction"] = (
        out["opposite_direction"]
        & (out["disease_fdr"] <= fdr_threshold)
        & (out["treatment_fdr"] <= fdr_threshold)
        & (out["disease_z"].abs() >= minimum_absolute_z)
        & (out["treatment_z"].abs() >= minimum_absolute_z)
    )
    out["both_contrasts_supported"] = (
        (out["disease_fdr"] <= fdr_threshold)
        & (out["treatment_fdr"] <= fdr_threshold)
        & (out["disease_z"].abs() >= minimum_absolute_z)
        & (out["treatment_z"].abs() >= minimum_absolute_z)
    )

    out["reversal_class"] = np.select(
        [
            (
                out["passes_conjunction"]
                & out["restoration_fraction"].between(
                    0,
                    1.5,
                )
            ),
            (
                out["passes_conjunction"]
                & (
                    out["restoration_fraction"]
                    > 1.5
                )
            ),
            (
                out["both_contrasts_supported"]
                & ~out["opposite_direction"]
            ),
        ],
        [
            "restored",
            "overshoot",
            "reinforced",
        ],
        default="inconclusive",
    )
    return out.sort_values(["passes_conjunction", "reversal_score"], ascending=[False, False])


def compare_primary_and_mouse_level(primary: pd.DataFrame, sensitivity: pd.DataFrame) -> pd.DataFrame:
    merged = primary[["gene_id", "disease_effect", "treatment_effect", "passes_conjunction"]].merge(
        sensitivity[["gene_id", "disease_effect", "treatment_effect", "passes_conjunction"]],
        on="gene_id", suffixes=("_fiber_blocked", "_mouse_level"),
    )
    merged["disease_direction_stable"] = np.sign(merged["disease_effect_fiber_blocked"]) == np.sign(merged["disease_effect_mouse_level"])
    merged["treatment_direction_stable"] = np.sign(merged["treatment_effect_fiber_blocked"]) == np.sign(merged["treatment_effect_mouse_level"])
    merged["robust_reversal"] = merged["passes_conjunction_fiber_blocked"] & merged["passes_conjunction_mouse_level"]
    return merged

