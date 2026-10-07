from __future__ import annotations

from pathlib import Path

import pandas as pd

from .config import load_parameters


REQUIRED_COLUMNS = {
    "accession", "enabled", "eligibility", "species", "tissue", "assay", "role",
    "biological_unit", "n_control_units", "n_disease_units", "analysis_engine", "geo_url",
}


def read_manifest(path: str | Path) -> pd.DataFrame:
    frame = pd.read_csv(path, keep_default_na=False)
    missing = REQUIRED_COLUMNS.difference(frame.columns)
    if missing:
        raise ValueError(f"Manifest is missing required columns: {sorted(missing)}")
    frame["enabled"] = frame["enabled"].astype(str).str.lower().map({"true": True, "false": False})
    if frame["enabled"].isna().any():
        raise ValueError("Manifest enabled column must contain only true/false")
    for column in ("n_control_units", "n_disease_units", "n_treated_units"):
        if column in frame:
            frame[column] = pd.to_numeric(frame[column], errors="raise").astype(int)
    return frame


def validate_manifest(manifest_path: str | Path, parameters_path: str | Path) -> pd.DataFrame:
    frame = read_manifest(manifest_path)
    parameters = load_parameters(parameters_path)
    rules = parameters["eligibility"]
    issues: list[dict[str, str]] = []

    def add(severity: str, accession: str, check: str, message: str) -> None:
        issues.append({"severity": severity, "accession": accession, "check": check, "message": message})

    duplicates = frame[frame["accession"].duplicated(keep=False)]
    for accession in duplicates["accession"].unique():
        add("ERROR", accession, "unique_accession", "Accession appears more than once")

    forbidden = set(rules.get("forbidden_accessions", []))
    for row in frame.itertuples(index=False):
        accession = row.accession
        if row.enabled and accession in forbidden:
            add("ERROR", accession, "forbidden_accession", "A known ineligible accession is enabled")
        if row.enabled and not str(row.geo_url).startswith("https://www.ncbi.nlm.nih.gov/geo/"):
            add("ERROR", accession, "geo_provenance", "Enabled GEO dataset lacks an official GEO URL")
        if row.enabled and row.eligibility == "eligible":
            minimum = int(rules["minimum_independent_units_per_group"])
            if row.n_control_units < minimum or row.n_disease_units < minimum:
                add("ERROR", accession, "independent_units", f"Eligible core dataset has fewer than {minimum} units per group")
        if row.assay == "single_cell_rna_seq" and row.analysis_engine != "scrna_localization_only":
            minimum = int(rules["single_cell_minimum_independent_samples_per_group_for_de"])
            if row.n_control_units < minimum or row.n_disease_units < minimum:
                add("ERROR", accession, "single_cell_pseudoreplication", "Cells cannot substitute for independent samples")
        if row.role == "mouse_muscle_treatment_reversal" and row.biological_unit != "mouse":
            add("ERROR", accession, "repeated_measure_unit", "Single fibers must be blocked or aggregated by mouse")
        if row.eligibility in {"sensitivity_only", "localization_only"}:
            add("INFO", accession, "restricted_role", f"Dataset is restricted to {row.eligibility}")
        if not row.enabled:
            add("INFO", accession, "excluded", str(row.design_caveat))

    if not issues:
        issues.append({"severity": "INFO", "accession": "PROJECT", "check": "manifest", "message": "No issues detected"})
    return pd.DataFrame(issues)


def assert_manifest_valid(manifest_path: str | Path, parameters_path: str | Path) -> pd.DataFrame:
    issues = validate_manifest(manifest_path, parameters_path)
    errors = issues[issues["severity"] == "ERROR"]
    if not errors.empty:
        rendered = "\n".join(f"{r.accession}: {r.message}" for r in errors.itertuples())
        raise ValueError(f"Manifest validation failed:\n{rendered}")
    return issues

