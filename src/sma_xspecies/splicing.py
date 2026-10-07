from __future__ import annotations

from pathlib import Path
import gzip

import numpy as np
import pandas as pd

from .standardize import read_table


def parse_inc_levels(value: object) -> float:
    numbers = pd.to_numeric(pd.Series(str(value).split(",")), errors="coerce")
    return float(numbers.mean()) if numbers.notna().any() else np.nan


def read_rmats_table(
    path: str | Path,
) -> pd.DataFrame:
    """Read an rMATS table containing a metadata preamble."""

    path = Path(path)

    compressed = path.name.lower().endswith(".gz")

    open_function = gzip.open if compressed else open
    open_mode = "rt" if compressed else "r"

    header_index = None

    with open_function(
        path,
        open_mode,
        encoding="utf-8",
        errors="replace",
    ) as handle:
        for line_number, line in enumerate(handle):
            if (
                "GeneID" in line
                and "IncLevelDifference" in line
            ):
                header_index = line_number
                break

    if header_index is None:
        raise ValueError(
            f"Could not locate rMATS header in {path}"
        )

    return pd.read_csv(
        path,
        sep="\t",
        compression="gzip" if compressed else None,
        skiprows=header_index,
        engine="python",
        quotechar='"',
    )


def standardize_rmats(
    path: str | Path,
    accession: str,
    event_type: str = "SE",
) -> pd.DataFrame:
    frame = read_rmats_table(path)

    frame.columns = [
        str(column).strip()
        for column in frame.columns
    ]

    required = {
        "GeneID",
        "FDR",
        "IncLevelDifference",
        "IncLevel_CTL",
        "IncLevel_SMA",
    }

    missing = required.difference(frame.columns)

    if missing:
        raise ValueError(
            f"rMATS file is missing {sorted(missing)}"
        )

    symbol_column = (
        "geneSymbol"
        if "geneSymbol" in frame.columns
        else None
    )

    mean_control = frame[
        "IncLevel_CTL"
    ].map(parse_inc_levels)

    mean_disease = frame[
        "IncLevel_SMA"
    ].map(parse_inc_levels)

    # Project-wide orientation: disease minus control.
    delta_disease_control = (
        mean_disease - mean_control
    )

    reported_delta = pd.to_numeric(
        frame["IncLevelDifference"],
        errors="coerce",
    )

    # Validate that the deposited value is CTL minus SMA.
    comparable = (
        reported_delta.notna()
        & delta_disease_control.notna()
    )

    orientation_error = (
        reported_delta[comparable]
        + delta_disease_control[comparable]
    ).abs()

    if (
        len(orientation_error)
        and orientation_error.quantile(0.99) > 0.01
    ):
        raise ValueError(
            "Deposited IncLevelDifference does not "
            "match CTL-minus-SMA orientation."
        )

    output = pd.DataFrame({
        "accession": accession,
        "event_type": event_type,
        "gene_id": (
            frame["GeneID"]
            .astype("string")
            .str.replace(
                '"',
                "",
                regex=False,
            )
            .str.replace(
                r"\.\d+$",
                "",
                regex=True,
            )
        ),
        "gene_symbol": (
            frame[symbol_column]
            .astype("string")
            .str.replace(
                '"',
                "",
                regex=False,
            )
            .str.strip()
            if symbol_column
            else pd.Series(
                pd.NA,
                index=frame.index,
                dtype="string",
            )
        ),
        "mean_psi_control": mean_control,
        "mean_psi_disease": mean_disease,
        "delta_psi": delta_disease_control,
        "reported_ctl_minus_sma":
            reported_delta,
        "p_value": (
            pd.to_numeric(
                frame["PValue"],
                errors="coerce",
            )
            if "PValue" in frame.columns
            else np.nan
        ),
        "fdr": pd.to_numeric(
            frame["FDR"],
            errors="coerce",
        ),
        "direction_definition":
            "SMA_minus_control",
    })

    coordinate_columns = [
        column
        for column in [
            "chr",
            "strand",
            "exonStart_0base",
            "exonEnd",
            "upstreamES",
            "upstreamEE",
            "downstreamES",
            "downstreamEE",
        ]
        if column in frame.columns
    ]

    if coordinate_columns:
        output["event_id"] = (
            frame[coordinate_columns]
            .astype(str)
            .agg(":".join, axis=1)
        )
    else:
        output["event_id"] = (
            accession
            + ":"
            + event_type
            + ":"
            + output.index.astype(str)
        )

    return output.dropna(
        subset=[
            "gene_id",
            "delta_psi",
            "fdr",
        ]
    )


def significant_events(events: pd.DataFrame, fdr: float = 0.05, minimum_absolute_delta_psi: float = 0.10) -> pd.DataFrame:
    return events.loc[(events["fdr"] <= fdr) & (events["delta_psi"].abs() >= minimum_absolute_delta_psi)].copy()

