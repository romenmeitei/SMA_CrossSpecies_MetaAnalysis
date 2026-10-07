from __future__ import annotations

from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

from .utils import runtime_versions, write_json


def write_run_record(root: str | Path) -> Path:
    root = Path(root)
    outputs = []
    for path in sorted((root / "results").rglob("*")):
        if path.is_file() and not path.name.startswith(".") and path.name not in {"run_record.json", "RESULT_INDEX.md"}:
            outputs.append({"path": str(path.relative_to(root)), "bytes": path.stat().st_size})
    record = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "versions": runtime_versions(),
        "outputs": outputs,
        "interpretation_guards": [
            "Tissues are not pooled without an explicit joint neural sensitivity model.",
            "GSE69175 is sensitivity-only because line and disease status are confounded.",
            "GSE207890 fibers are repeated measures within mice.",
            "GSE208629 is localization-only because there is one pooled sample per condition.",
            "Orthology is restricted to versioned Ensembl one-to-one relationships.",
        ],
    }
    destination = root / "results" / "run_record.json"
    write_json(record, destination)
    return destination


def write_result_index(root: str | Path) -> Path:
    root = Path(root)
    lines = ["# Analysis output index", "", "This file is generated from the current results directory.", ""]
    files = [path for path in sorted((root / "results").rglob("*")) if path.is_file() and not path.name.startswith(".") and path.name != "RESULT_INDEX.md"]
    for path in files:
        relative = path.relative_to(root)
        lines.append(f"- `{relative}` ({path.stat().st_size:,} bytes)")
    if not files:
        lines.append("- No scientific outputs have been generated yet.")
    destination = root / "results" / "RESULT_INDEX.md"
    destination.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return destination
