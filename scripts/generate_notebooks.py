#!/usr/bin/env python3
"""Generate the Colab master runner and focused stage notebooks."""
from __future__ import annotations

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
NOTEBOOKS = ROOT / "notebooks"


def md(source: str) -> dict:
    return {"cell_type": "markdown", "metadata": {}, "source": source.splitlines(keepends=True)}


def code(source: str) -> dict:
    return {"cell_type": "code", "execution_count": None, "metadata": {}, "outputs": [], "source": source.splitlines(keepends=True)}


def notebook(cells: list[dict]) -> dict:
    return {
        "cells": cells,
        "metadata": {
            "colab": {"provenance": [], "toc_visible": True},
            "kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"},
            "language_info": {"name": "python", "version": "3"},
        },
        "nbformat": 4,
        "nbformat_minor": 5,
    }


MASTER_CELLS = [
    md("""# SMA cross-species meta-analysis — master Colab runner

This notebook runs the rebuilt workflow from audit to publication outputs. It does
not use the discarded conserved-gene or chunked-network results. Read the printed
audit messages: missing or ineligible evidence is surfaced rather than silently
filled in.

**Before running:** download the complete project ZIP. The next cell will ask you
to upload it if the project is not already present in Colab. Start with light
downloads. Raw-read reprocessing is a separate confirmatory profile documented in
`docs/RAW_REPROCESSING.md`."""),
    code("""from pathlib import Path
import zipfile

# Colab-safe setup: locate an existing copy or request the complete project ZIP.
PREFERRED_PROJECT_ROOT = Path('/content/SMA_CrossSpecies_MetaAnalysis')
AUTO_UPLOAD_PROJECT_ZIP = True

def locate_project(search_root=Path('/content')):
    candidates = []
    preferred = PREFERRED_PROJECT_ROOT / 'pyproject.toml'
    if preferred.exists():
        candidates.append(preferred.parent)
    for pyproject in search_root.rglob('pyproject.toml'):
        try:
            if 'name = \"sma-xspecies\"' in pyproject.read_text(encoding='utf-8'):
                candidates.append(pyproject.parent)
        except (OSError, UnicodeDecodeError):
            pass
    unique = sorted(set(candidates), key=lambda path: (len(path.parts), str(path)))
    return unique[0] if unique else None

PROJECT_ROOT = locate_project()

if PROJECT_ROOT is None and AUTO_UPLOAD_PROJECT_ZIP:
    try:
        from google.colab import files
    except ImportError as exc:
        raise FileNotFoundError(
            'Project folder is missing. Unzip SMA_CrossSpecies_MetaAnalysis_v0.1.1.zip '
            'and set PREFERRED_PROJECT_ROOT to that folder.'
        ) from exc

    print('Upload the COMPLETE project ZIP (not only 00_Master_Run.ipynb).')
    uploaded = files.upload()
    zip_names = [name for name in uploaded if name.lower().endswith('.zip')]
    if len(zip_names) != 1:
        raise ValueError(f'Upload exactly one project ZIP; received: {list(uploaded)}')

    zip_path = Path('/content') / Path(zip_names[0]).name
    if not zip_path.exists():
        zip_path.write_bytes(uploaded[zip_names[0]])

    extraction_root = Path('/content').resolve()
    with zipfile.ZipFile(zip_path) as archive:
        for member in archive.infolist():
            destination = (extraction_root / member.filename).resolve()
            if destination != extraction_root and extraction_root not in destination.parents:
                raise ValueError(f'Unsafe path in ZIP: {member.filename}')
        archive.extractall(extraction_root)

    PROJECT_ROOT = locate_project()

if PROJECT_ROOT is None or not (PROJECT_ROOT / 'pyproject.toml').exists():
    raise FileNotFoundError(
        'The SMA project could not be located. Upload the complete project ZIP, '
        'not the notebook file by itself.'
    )

RUN_DOWNLOADS = True
INCLUDE_LARGE_ARCHIVES = False     # GSE207890/GSE208629 raw archives are >200 MB
INSTALL_R_PACKAGES = True
FETCH_ORTHOLOGS = True
RUN_GSE209926_DE = True
RUN_GSE207890_BLOCKED_DE = True
RUN_META_ANALYSIS = True
RUN_TREATMENT_REVERSAL = True
RUN_SPLICING = True

print('Project ready:', PROJECT_ROOT)"""),
    md("## 1. Install and validate the locked design"),
    code("""import os, subprocess, sys
os.chdir(PROJECT_ROOT)
subprocess.run([sys.executable, '-m', 'pip', 'install', '-q', '-e', '.[enrichment,notebook,test]'], check=True)

def run(*parts, check=True):
    command = [str(part) for part in parts]
    print('$', ' '.join(command))
    return subprocess.run(command, check=check)

run('sma-xspecies', 'validate')"""),
    md("""The validation report must contain no errors. Informational restrictions
for GSE69175 and GSE208629 are expected and are part of the design."""),
    md("## 2. Retrieve small public supplementary files and record checksums"),
    code("""if RUN_DOWNLOADS:
    command = ['sma-xspecies', 'download']
    if INCLUDE_LARGE_ARCHIVES:
        command.append('--include-large')
    run(*command, check=False)
else:
    print('Download stage disabled')"""),
    md("## 3. Standardize published study-level summaries"),
    code("""run('sma-xspecies', 'standardize', check=False)
from pathlib import Path
print('\\nCurrent DE inputs:')
for path in sorted((PROJECT_ROOT / 'results/de').glob('*.csv')):
    print(' ', path.name)"""),
    md("## 4. Install Bioconductor and pin one-to-one orthology"),
    code("""if INSTALL_R_PACKAGES:
    run('Rscript', 'environment/install_bioc.R')
if FETCH_ORTHOLOGS:
    run('Rscript', 'R/fetch_orthologs.R', 'resources/ensembl_human_mouse_one_to_one.tsv', '113')"""),
    md("## 4b. Fit any processed expression matrices discovered during standardization"),
    code("""processed_models = [
    ('GSE56284', 'Mus_musculus', 'motor_neuron', 'eligible'),
    ('GSE69175', 'Homo_sapiens', 'motor_neuron', 'sensitivity_only'),
]
for accession, species, tissue, eligibility in processed_models:
    matrix = PROJECT_ROOT / f'data/processed/{accession}_expression.tsv'
    output = PROJECT_ROOT / f'results/de/{accession}_disease-control.csv'
    if matrix.exists() and not output.exists():
        run('Rscript', 'R/run_processed_limma.R', matrix,
            f'config/sample_annotations/{accession}.csv', accession, species,
            tissue, eligibility, output)
    else:
        print(accession, ': published summary already standardized or matrix absent')"""),
    md("## 5. GSE209926 spinal-cord DE from sample-level counts"),
    code("""archive = PROJECT_ROOT / 'data/raw/GSE209926/GSE209926_RAW.tar'
counts = PROJECT_ROOT / 'data/processed/GSE209926_counts.tsv'
if RUN_GSE209926_DE and archive.exists():
    run(sys.executable, 'scripts/prepare_gse209926_counts.py', archive,
        'config/sample_annotations/GSE209926.csv', counts)
    run('Rscript', 'R/run_deseq2_counts.R', counts,
        'config/sample_annotations/GSE209926.csv', 'GSE209926',
        'results/de/GSE209926_disease-control.csv')
else:
    print('Skipped: archive absent or RUN_GSE209926_DE=False')"""),
    md("## 6. GSE207890 blocked myofiber DE and animal-level sensitivity"),
    code("""if RUN_GSE207890_BLOCKED_DE:
    run('Rscript', 'R/run_geo_microarray_blocked.R', 'GSE207890',
        'config/sample_annotations/GSE207890.csv', 'results/de')
else:
    print('Blocked myofiber stage disabled')"""),
    md("## 7. Cross-species, tissue-aware meta-analysis"),
    code("""mapping = PROJECT_ROOT / 'resources/ensembl_human_mouse_one_to_one.tsv'
if RUN_META_ANALYSIS and mapping.exists():
    run('sma-xspecies', 'meta')
    # Prespecified sensitivity including the line-confounded GSE69175 study.
    run('sma-xspecies', 'meta', '--include-sensitivity-only')
else:
    print('Skipped: ortholog mapping absent or RUN_META_ANALYSIS=False')"""),
    md("## 8. Treatment reversal, splicing, and final audit"),
    code("""if RUN_TREATMENT_REVERSAL:
    run('sma-xspecies', 'rescue', check=False)
if RUN_SPLICING:
    run('sma-xspecies', 'splicing', check=False)
run('pytest')
run('sma-xspecies', 'figures', check=False)
run('sma-xspecies', 'report')"""),
    md("""## 9. Review before interpretation

Open `data/audit/manifest_validation.csv`, `results/run_record.json`, the per-study
QC outputs, and `docs/PUBLICATION_GATES.md`. A nonempty table is not automatically
valid evidence: check participating studies, species, tissue stratum, heterogeneity,
and the leave-one-dataset-out result for every prioritized gene.

The notebook intentionally does not run enrichment or a network unless versioned
gene-set and STRING bulk files are provided. This prevents silent dependence on a
changing web service or an incomplete chunked graph."""),
]


STAGES = {
    "01_Dataset_Audit_QC.ipynb": (
        "Dataset audit and QC",
        "Validate eligibility, units of replication, sources, and checksums before any biology is interpreted.",
        "run('sma-xspecies', 'validate')\nrun('sma-xspecies', 'download', check=False)",
    ),
    "02_Differential_Expression.ipynb": (
        "Within-dataset differential expression",
        "Standardize published summaries, run DESeq2 for counts, and fit the blocked myofiber model.",
        "run('sma-xspecies', 'standardize', check=False)\n# Dataset-specific R commands are listed in 00_Master_Run.ipynb.",
    ),
    "03_CrossSpecies_MetaAnalysis.ipynb": (
        "Cross-species meta-analysis",
        "Attach versioned one-to-one orthologs and run tissue-stratified REML or signed-z synthesis.",
        "run('sma-xspecies', 'meta')",
    ),
    "04_Tissue_Treatment_Analysis.ipynb": (
        "Tissue and treatment-reversal analysis",
        "Keep neural evidence separate from muscle and require two independent contrasts for reversal.",
        "run('sma-xspecies', 'rescue')",
    ),
    "05_Differential_Splicing.ipynb": (
        "Event-level differential splicing",
        "Parse rMATS events and threshold FDR plus absolute delta-PSI; no symbol-substring modules.",
        "run('sma-xspecies', 'splicing')",
    ),
    "06_Enrichment_Networks.ipynb": (
        "Versioned enrichment and networks",
        "Use preranked GSEA with a supplied GMT and the complete STRING bulk graph with degree-matched permutations.",
        "print('Provide a versioned GMT and STRING bulk files; see resources/README.md. Network analysis remains off by default.')",
    ),
    "07_Figures_Tables_Validation.ipynb": (
        "Figures, tables, and validation",
        "Run unit tests, write the run record, and evaluate every publication gate.",
        "run('pytest')\nrun('sma-xspecies', 'figures', check=False)\nrun('sma-xspecies', 'report')",
    ),
}


def stage_cells(title: str, description: str, command: str) -> list[dict]:
    return [
        md(f"# {title}\n\n{description}\n\nRun this notebook from the repository root in Colab."),
        code("""from pathlib import Path
import os, subprocess, sys
PROJECT_ROOT = Path('/content/SMA_CrossSpecies_MetaAnalysis')  # edit if needed
if not (PROJECT_ROOT / 'pyproject.toml').exists():
    raise FileNotFoundError('Set PROJECT_ROOT to the unzipped/cloned repository')
os.chdir(PROJECT_ROOT)
subprocess.run([sys.executable, '-m', 'pip', 'install', '-q', '-e', '.[enrichment,notebook,test]'], check=True)
def run(*parts, check=True):
    print('$', ' '.join(map(str, parts)))
    return subprocess.run([str(part) for part in parts], check=check)"""),
        code(command),
        md("Review the machine-readable audit and result tables before proceeding to the next stage."),
    ]


def main() -> None:
    NOTEBOOKS.mkdir(parents=True, exist_ok=True)
    (NOTEBOOKS / "00_Master_Run.ipynb").write_text(json.dumps(notebook(MASTER_CELLS), indent=1) + "\n", encoding="utf-8")
    for filename, (title, description, command) in STAGES.items():
        (NOTEBOOKS / filename).write_text(json.dumps(notebook(stage_cells(title, description, command)), indent=1) + "\n", encoding="utf-8")
    print(f"Generated {1 + len(STAGES)} notebooks in {NOTEBOOKS}")


if __name__ == "__main__":
    main()
