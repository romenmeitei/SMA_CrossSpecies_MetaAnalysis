# Reproducibility guide

## Software setup

```bash
conda env create -f environment.yml
conda activate sma_downstream
python -m pip install -e ".[test,enrichment]"
```

## Basic verification

```bash
python -m pytest -q tests
python -m sma_xspecies.cli validate --root .
```

## Evidence hierarchy

The primary statistical conclusion is based on Hartung-Knapp and modified
Hartung-Knapp small-k sensitivity analyses. No gene passed transcriptome-wide
FDR control under these procedures.

Corrected REML identified 731 BH-FDR screening genes. These remain screening
evidence because only three studies contributed to the neural meta-analysis.

FAM107A was the sole Tier-4 internally robust candidate. External datasets
supported its direction but did not provide nominal or multiplicity-corrected
statistical replication.

UGP2, Brain Development splicing enrichment and AMPA-receptor trafficking
treatment-response results are explicitly secondary findings.

## Frozen outputs

The complete results, audit records and environment specifications are
distributed in the checksummed full-release archive attached to the GitHub
Release.

Manuscript figures and tables are distributed in the separate post-visual-QA
archive.

Large public datasets should be retrieved from their authoritative GEO/NCBI
locations using the accessions in `config/dataset_manifest.csv`.

The expected release assets and their SHA-256 values are recorded in
`docs/RELEASE_ASSETS.md`.
