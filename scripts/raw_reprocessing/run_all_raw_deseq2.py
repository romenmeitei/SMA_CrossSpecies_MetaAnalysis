#!/usr/bin/env python3

from pathlib import Path
from datetime import datetime, timezone
import subprocess
import pandas as pd

ROOT = Path(
    "/home/romen/SMA_Cross/work/"
    "SMA_RawRead_Reprocessing"
)

R_SCRIPT = Path(
    "/home/romen/SMA_Cross/"
    "run_raw_deseq2_linux.R"
)

MANIFEST = pd.read_csv(
    ROOT
    / "audit"
    / "selected_sra_ena_run_manifest.csv",
    dtype=str,
).fillna("")

QUANT_ROOT = ROOT / "salmon_quant"
REFERENCE_ROOT = (
    ROOT
    / "references"
    / "ensembl_113"
)
ANNOTATION_ROOT = ROOT / "annotations"
DE_ROOT = ROOT / "differential_expression"
QC_ROOT = ROOT / "qc" / "deseq2"
LOG_ROOT = ROOT / "logs" / "deseq2"

for directory in [
    DE_ROOT,
    QC_ROOT,
    LOG_ROOT,
]:
    directory.mkdir(
        parents=True,
        exist_ok=True,
    )

ORDER = [
    "GSE56284",
    "GSE209926",
    "GSE108094",
    "GSE69175",
]

records = []

for accession in ORDER:
    study = MANIFEST.loc[
        MANIFEST["accession"] == accession
    ].copy()

    if study.empty:
        raise RuntimeError(
            f"No manifest entries for {accession}"
        )

    unique_value = {}

    for field in [
        "reference_key",
        "species",
        "tissue",
        "role",
        "eligibility",
    ]:
        values = sorted(
            set(
                value
                for value in study[field]
                if value
            )
        )

        if len(values) != 1:
            raise RuntimeError(
                f"{accession}: expected one "
                f"{field}, observed {values}"
            )

        unique_value[field] = values[0]

    reference_key = unique_value[
        "reference_key"
    ]

    suffix = (
        "descriptive_only"
        if accession == "GSE69175"
        else "disease-control_raw"
    )

    output = (
        DE_ROOT
        / (
            f"{accession}_{suffix}_"
            "Ensembl113.csv"
        )
    )

    qc_dir = QC_ROOT / accession

    status_file = (
        qc_dir
        / f"{accession}_analysis_status.txt"
    )

    command = [
        "Rscript",
        str(R_SCRIPT),
        str(QUANT_ROOT),
        str(
            REFERENCE_ROOT
            / reference_key
            / "tx2gene.tsv"
        ),
        str(
            ANNOTATION_ROOT
            / f"{accession}.csv"
        ),
        accession,
        unique_value["species"],
        unique_value["tissue"],
        unique_value["role"],
        unique_value["eligibility"],
        str(output),
        "10",
        str(qc_dir),
        str(status_file),
    ]

    log_path = (
        LOG_ROOT
        / f"{accession}.log"
    )

    print(
        "\nRunning:",
        accession,
        flush=True,
    )

    print(
        "$",
        " ".join(command),
        flush=True,
    )

    with open(log_path, "w") as log:
        completed = subprocess.run(
            command,
            text=True,
            stdout=log,
            stderr=subprocess.STDOUT,
        )

    if completed.returncode != 0:
        print(
            log_path.read_text()[-8000:]
        )

        raise RuntimeError(
            f"DESeq2 failed for {accession}"
        )

    if (
        not output.is_file()
        or output.stat().st_size <= 1000
    ):
        raise RuntimeError(
            f"Invalid DE output: {output}"
        )

    if not status_file.is_file():
        raise RuntimeError(
            f"Status file missing: "
            f"{status_file}"
        )

    result = pd.read_csv(output)

    status_lines = (
        status_file
        .read_text()
        .splitlines()
    )

    status = (
        status_lines[0].split(
            "=",
            1,
        )[-1]
    )

    records.append({
        "accession": accession,
        "status": status,
        "output": str(output),
        "rows": len(result),
        "finite_z": int(
            pd.to_numeric(
                result["z"],
                errors="coerce",
            ).notna().sum()
        ),
        "fdr_le_0.05": int(
            (
                pd.to_numeric(
                    result["fdr"],
                    errors="coerce",
                )
                <= 0.05
            ).sum()
        ),
        "completed_utc": datetime.now(
            timezone.utc
        ).isoformat(),
    })

    print(
        status_file.read_text(),
        flush=True,
    )

summary = pd.DataFrame(records)

summary.to_csv(
    ROOT
    / "audit"
    / "raw_deseq2_summary.csv",
    index=False,
)

print("\nDESeq2 summary:")
print(summary.to_string(index=False))
print("\nPASS: all four analyses completed.")
