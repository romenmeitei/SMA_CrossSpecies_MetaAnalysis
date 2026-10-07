from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from .ortholog import attach_ortholog_ids, read_one_to_one_mapping
from .stats import random_effects_reml, weighted_stouffer
from .utils import bh_adjust


def read_effect_files(directory: str | Path) -> pd.DataFrame:
    files = sorted(Path(directory).glob("*.csv"))
    if not files:
        raise FileNotFoundError(f"No standardized effect CSV files found in {directory}")
    frames = [pd.read_csv(path) for path in files]
    required = {"accession", "gene_id", "effect", "se", "z", "species", "tissue", "eligibility", "contrast"}
    for path, frame in zip(files, frames):
        missing = required.difference(frame.columns)
        if missing:
            raise ValueError(f"{path.name} is missing {sorted(missing)}")
    return pd.concat(frames, ignore_index=True)


def _finite_median(values: pd.Series) -> float:
    numeric = pd.to_numeric(
        values,
        errors="coerce",
    ).to_numpy(dtype=float)

    numeric = numeric[np.isfinite(numeric)]

    return (
        float(np.median(numeric))
        if len(numeric)
        else np.nan
    )


def collapse_within_study_orthologs(
    data: pd.DataFrame,
) -> pd.DataFrame:
    """Collapse multiple source rows per study and ortholog.

    Median aggregation prevents arbitrary row selection and
    avoids selecting transcripts using statistical significance.
    """

    group_columns = [
        column
        for column in [
            "accession",
            "ortholog_id",
            "ortholog_symbol",
            "species",
            "tissue",
            "eligibility",
            "contrast",
            "stratum",
        ]
        if column in data.columns
    ]

    collapsed = (
        data.groupby(
            group_columns,
            as_index=False,
            dropna=False,
            sort=False,
        )
        .agg(
            effect=("effect", _finite_median),
            z=("z", _finite_median),
            se=("se", _finite_median),
            n_source_rows=("effect", "size"),
        )
    )

    return collapsed


def _meta_group(frame: pd.DataFrame, minimum_studies: int, minimum_species: int) -> dict[str, object] | None:
    frame = frame.drop_duplicates("accession").copy()

    valid_se = (
        np.isfinite(frame["effect"])
        & np.isfinite(frame["se"])
        & (frame["se"] > 0)
    )

    se_frame = frame.loc[valid_se].copy()

    if (
        se_frame["accession"].nunique() >= minimum_studies
        and se_frame["species"].nunique() >= minimum_species
    ):
        result = random_effects_reml(
            se_frame["effect"].to_numpy(),
            se_frame["se"].to_numpy(),
        )
        method = "random_effects_reml"
        payload = result.__dict__.copy()
        analysis_frame = se_frame

    else:
        z_frame = frame.loc[
            np.isfinite(frame["z"])
            & np.isfinite(frame["effect"])
        ].copy()

        if (
            z_frame["accession"].nunique() < minimum_studies
            or z_frame["species"].nunique() < minimum_species
        ):
            return None

        # Numerical-safety guard for the Stouffer fallback.
        #
        # Extremely small study-level p-values may produce very large
        # finite Wald z statistics. In an equal-weight Stouffer analysis,
        # one such value can dominate the cross-study result. Preserve the
        # original standardized z values, but cap only the values supplied
        # to the Stouffer combination.
        stouffer_z_cap = 8.0
        original_z = z_frame["z"].to_numpy(dtype=float)
        capped_z = np.clip(
            original_z,
            -stouffer_z_cap,
            stouffer_z_cap,
        )

        combined_z, p_value = weighted_stouffer(
            capped_z
        )

        method = "weighted_stouffer"
        payload = {
            "k": int(
                z_frame["accession"].nunique()
            ),
            "estimate": np.nan,
            "se": np.nan,
            "ci_low": np.nan,
            "ci_high": np.nan,
            "z": combined_z,
            "p_value": p_value,
            "tau2": np.nan,
            "q": np.nan,
            "q_p_value": np.nan,
            "i2": np.nan,
            "stouffer_z_cap": stouffer_z_cap,
            "maximum_absolute_input_z": float(
                np.max(np.abs(original_z))
            ),
            "number_z_values_capped": int(
                np.sum(np.abs(original_z) > stouffer_z_cap)
            ),
        }
        analysis_frame = z_frame

    payload.update({
        "method": method,
        "n_species": int(
            analysis_frame["species"].nunique()
        ),
        "datasets": ";".join(
            sorted(
                analysis_frame["accession"].unique()
            )
        ),
        "species": ";".join(
            sorted(
                analysis_frame["species"].unique()
            )
        ),
        "sign_concordance": float(
            max(
                (analysis_frame["effect"] > 0).mean(),
                (analysis_frame["effect"] < 0).mean(),
            )
        ),
    })

    return payload


def run_gene_meta(
    effects: pd.DataFrame,
    mapping_path: str | Path,
    *,
    minimum_studies: int = 2,
    strict: bool = True,
    neural_joint: bool = False,
) -> pd.DataFrame:
    data = effects.loc[effects["contrast"] == "disease-control"].copy()
    if strict:
        data = data.loc[data["eligibility"] == "eligible"]
    mapping = read_one_to_one_mapping(mapping_path)
    data = attach_ortholog_ids(data, mapping)
    if neural_joint:
        data = data.loc[data["tissue"].isin(["motor_neuron", "spinal_cord"])].copy()
        data["stratum"] = "neural_joint"
    else:
        data["stratum"] = data["tissue"]

    data = collapse_within_study_orthologs(data)

    records: list[dict[str, object]] = []
    for (ortholog_id, symbol, stratum), group in data.groupby(["ortholog_id", "ortholog_symbol", "stratum"], sort=False):
        result = _meta_group(group, minimum_studies, minimum_species=2)
        if result is not None:
            records.append({"ortholog_id": ortholog_id, "gene_symbol": symbol, "stratum": stratum, **result})
    out = pd.DataFrame(records)
    if not out.empty:
        out["fdr"] = bh_adjust(out["p_value"])
        out = out.sort_values(["fdr", "p_value", "gene_symbol"], na_position="last")
    return out


def leave_one_dataset_out(effects: pd.DataFrame, mapping_path: str | Path, strict: bool = True) -> pd.DataFrame:
    neural_tissues = {"motor_neuron", "spinal_cord"}

    eligible_inputs = effects.loc[
        (effects["contrast"] == "disease-control")
        & effects["tissue"].isin(neural_tissues)
    ].copy()

    if strict:
        eligible_inputs = eligible_inputs.loc[
            eligible_inputs["eligibility"] == "eligible"
        ]

    records = []

    for accession in sorted(
        eligible_inputs["accession"].dropna().unique()
    ):
        subset = effects.loc[
            effects["accession"] != accession
        ]

        result = run_gene_meta(
            subset,
            mapping_path,
            strict=strict,
            neural_joint=True,
        )

        if result.empty:
            continue

        result.insert(
            0,
            "omitted_dataset",
            accession,
        )
        records.append(result)

    return (
        pd.concat(records, ignore_index=True)
        if records
        else pd.DataFrame()
    )


def species_concordance(effects: pd.DataFrame, mapping_path: str | Path, strict: bool = True) -> pd.DataFrame:
    neural_tissues = {"motor_neuron", "spinal_cord"}

    data = effects.loc[
        (effects["contrast"] == "disease-control")
        & effects["tissue"].isin(neural_tissues)
    ].copy()

    if strict:
        data = data.loc[
            data["eligibility"] == "eligible"
        ]

    data = attach_ortholog_ids(
        data,
        read_one_to_one_mapping(mapping_path),
    )

    data = collapse_within_study_orthologs(data)

    means = (
        data.groupby(
            [
                "ortholog_id",
                "ortholog_symbol",
                "species",
            ],
            as_index=False,
        )["z"]
        .mean()
    )

    wide = means.pivot(
        index=["ortholog_id", "ortholog_symbol"],
        columns="species",
        values="z",
    ).dropna()

    required_species = {
        "Homo_sapiens",
        "Mus_musculus",
    }

    if (
        wide.empty
        or not required_species.issubset(wide.columns)
    ):
        return pd.DataFrame(
            [{
                "n_genes": 0,
                "spearman_rho": np.nan,
                "p_value": np.nan,
                "sign_agreement": np.nan,
            }]
        )

    rho, p_value = spearmanr(
        wide["Homo_sapiens"],
        wide["Mus_musculus"],
    )

    sign_agreement = np.mean(
        np.sign(wide["Homo_sapiens"])
        == np.sign(wide["Mus_musculus"])
    )

    return pd.DataFrame(
        [{
            "n_genes": len(wide),
            "spearman_rho": rho,
            "p_value": p_value,
            "sign_agreement": sign_agreement,
        }]
    )


def ortholog_mapping_audit(effects: pd.DataFrame, mapping_path: str | Path, strict: bool = True) -> pd.DataFrame:
    neural_tissues = {"motor_neuron", "spinal_cord"}

    data = effects.loc[
        (effects["contrast"] == "disease-control")
        & effects["tissue"].isin(neural_tissues)
    ].copy()

    if strict:
        data = data.loc[
            data["eligibility"] == "eligible"
        ]
    mapped = attach_ortholog_ids(data, read_one_to_one_mapping(mapping_path))
    total = data.groupby("accession")["gene_id"].nunique().rename("n_input_genes")
    mapped_total = mapped.groupby("accession")["gene_id"].nunique().rename("n_mapped_genes")
    audit = pd.concat([total, mapped_total], axis=1).fillna(0).reset_index()
    audit["n_mapped_genes"] = audit["n_mapped_genes"].astype(int)
    audit["mapping_fraction"] = audit["n_mapped_genes"] / audit["n_input_genes"]
    modes = mapped.groupby(["accession", "mapping_key"]).size().unstack(fill_value=0).reset_index()
    return audit.merge(modes, on="accession", how="left")
