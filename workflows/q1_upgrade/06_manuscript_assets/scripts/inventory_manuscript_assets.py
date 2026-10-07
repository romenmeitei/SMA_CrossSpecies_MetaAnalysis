from pathlib import Path
import csv
import sys

project = Path(sys.argv[1]).resolve()
q1_root = Path(sys.argv[2]).resolve()
raw_root = Path(sys.argv[3]).resolve()
asset_root = Path(sys.argv[4]).resolve()

audit_root = asset_root / "audit"
audit_root.mkdir(parents=True, exist_ok=True)

search_roots = [
    project / "results",
    project / "data" / "audit",
    q1_root / "01_small_k_meta",
    q1_root / "02_forest_plots",
    q1_root / "03_GSE56284_investigation",
    q1_root / "04_external_validation" / "results",
    q1_root / "04_external_validation" / "figures",
    q1_root / "04_external_validation" / "audit",
    raw_root / "qc" / "deseq2",
    raw_root / "reports" / "multiqc",
    raw_root / "audit",
]

figure_extensions = {
    ".png",
    ".pdf",
    ".svg",
    ".tif",
    ".tiff",
    ".jpg",
    ".jpeg",
}

table_extensions = {
    ".csv",
    ".tsv",
    ".xlsx",
    ".xls",
}

report_extensions = {
    ".html",
    ".txt",
    ".json",
    ".log",
    ".sha256",
}

records = []

for root in search_roots:
    if not root.exists():
        continue

    for path in sorted(root.rglob("*")):
        if not path.is_file():
            continue

        suffix = path.suffix.lower()

        if suffix in figure_extensions:
            kind = "figure"
        elif suffix in table_extensions:
            kind = "table"
        elif suffix in report_extensions:
            kind = "report"
        else:
            continue

        records.append({
            "kind": kind,
            "bytes": path.stat().st_size,
            "size_kb": round(
                path.stat().st_size / 1024,
                2,
            ),
            "filename": path.name,
            "path": str(path),
        })

inventory_path = (
    audit_root
    / "existing_manuscript_asset_inventory.csv"
)

with inventory_path.open(
    "w",
    newline="",
    encoding="utf-8",
) as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=[
            "kind",
            "bytes",
            "size_kb",
            "filename",
            "path",
        ],
    )

    writer.writeheader()
    writer.writerows(records)


def files_matching(directory, words):
    if not directory.exists():
        return []

    matches = []

    for path in directory.rglob("*"):
        if not path.is_file():
            continue

        if path.suffix.lower() not in figure_extensions:
            continue

        lower_name = path.name.lower()

        if any(
            word.lower() in lower_name
            for word in words
        ):
            matches.append(path)

    return sorted(matches)


figure_checks = [
    {
        "figure": "Figure 1",
        "description": (
            "Study workflow and evidence hierarchy"
        ),
        "candidates": files_matching(
            asset_root,
            [
                "workflow",
                "study_design",
                "evidence_hierarchy",
            ],
        ),
    },
    {
        "figure": "Figure 2",
        "description": (
            "FAM107A and UGP2 forest plots"
        ),
        "candidates": files_matching(
            q1_root / "02_forest_plots",
            [
                "fam107a",
                "ugp2",
                "forest",
            ],
        ),
    },
    {
        "figure": "Figure 3",
        "description": (
            "GSE56284 PCA correlation and influence"
        ),
        "candidates": files_matching(
            q1_root
            / "03_GSE56284_investigation",
            [
                "pca",
                "correlation",
                "heatmap",
                "dispersion",
                "influence",
                "stability",
            ],
        ),
    },
    {
        "figure": "Figure 4",
        "description": (
            "External validation effects"
        ),
        "candidates": files_matching(
            q1_root
            / "04_external_validation"
            / "figures",
            [
                "external",
                "validation",
                "fam107a",
                "ugp2",
                "pca",
            ],
        ),
    },
    {
        "figure": "Figure 5",
        "description": (
            "GSE108094 splicing and "
            "Brain Development enrichment"
        ),
        "candidates": files_matching(
            project / "results",
            [
                "splicing",
                "brain_development",
                "ugp2",
                "delta_psi",
            ],
        ),
    },
    {
        "figure": "Supplementary Figure",
        "description": (
            "AMPA receptor trafficking reversal"
        ),
        "candidates": files_matching(
            project
            / "results"
            / "enrichment"
            / "GSE207890_pathway_reversal",
            [
                "ampa",
                "reversal",
                "enrichment",
            ],
        ),
    },
]

status_path = (
    audit_root
    / "planned_figure_status.csv"
)

status_records = []

for check in figure_checks:
    candidates = check["candidates"]

    status_records.append({
        "figure": check["figure"],
        "description": check["description"],
        "status": (
            "GENERATED"
            if candidates
            else "NOT_GENERATED"
        ),
        "candidate_files": ";".join(
            str(path)
            for path in candidates
        ),
    })

with status_path.open(
    "w",
    newline="",
    encoding="utf-8",
) as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=[
            "figure",
            "description",
            "status",
            "candidate_files",
        ],
    )

    writer.writeheader()
    writer.writerows(status_records)

print("=" * 100)
print("MANUSCRIPT ASSET INVENTORY")
print("=" * 100)
print("Figures:", sum(
    record["kind"] == "figure"
    for record in records
))
print("Tables:", sum(
    record["kind"] == "table"
    for record in records
))
print("Reports:", sum(
    record["kind"] == "report"
    for record in records
))

print()
print("=" * 100)
print("PLANNED FIGURE STATUS")
print("=" * 100)

for record in status_records:
    print(
        f'{record["figure"]:22s} '
        f'{record["status"]:15s} '
        f'{record["description"]}'
    )

    if record["candidate_files"]:
        for path in record[
            "candidate_files"
        ].split(";"):
            print("   ", path)

print()
print("Inventory:", inventory_path)
print("Figure status:", status_path)
print()
print(
    "PASS: existing manuscript assets "
    "were inventoried without changing "
    "analysis outputs."
)
