from pathlib import Path

from sma_xspecies.manifest import validate_manifest


ROOT = Path(__file__).resolve().parents[1]


def test_manifest_has_no_errors():
    issues = validate_manifest(ROOT / "config/dataset_manifest.csv", ROOT / "config/analysis_parameters.yaml")
    assert not (issues["severity"] == "ERROR").any(), issues.to_string(index=False)


def test_known_wrong_disease_is_disabled():
    import pandas as pd

    manifest = pd.read_csv(ROOT / "config/dataset_manifest.csv")
    row = manifest.loc[manifest["accession"] == "GSE38417"].iloc[0]
    assert not bool(row["enabled"])
    assert row["eligibility"] == "excluded_wrong_disease"

