from pathlib import Path
from datetime import datetime, timezone
import csv
import hashlib
import json
import sys

root = Path(sys.argv[1]).resolve()
audit_root = root / "audit"

audit_root.mkdir(
    parents=True,
    exist_ok=True,
)

included_roots = [
    root / "figures",
    root / "tables",
]

files = []

for included_root in included_roots:
    if not included_root.is_dir():
        raise FileNotFoundError(included_root)

    files.extend(
        path
        for path in included_root.rglob("*")
        if path.is_file()
    )

files = sorted(set(files))

records = []

for path in files:
    digest = hashlib.sha256()

    with path.open("rb") as handle:
        for block in iter(
            lambda: handle.read(1024 * 1024),
            b"",
        ):
            digest.update(block)

    records.append({
        "relative_path": str(
            path.relative_to(root)
        ),
        "bytes": path.stat().st_size,
        "sha256": digest.hexdigest(),
    })

manifest_path = (
    audit_root
    / "generated_manuscript_asset_manifest.csv"
)

with manifest_path.open(
    "w",
    newline="",
    encoding="utf-8",
) as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=[
            "relative_path",
            "bytes",
            "sha256",
        ],
    )

    writer.writeheader()
    writer.writerows(records)

checksum_path = (
    audit_root
    / "generated_manuscript_assets.sha256"
)

with checksum_path.open(
    "w",
    encoding="utf-8",
) as handle:
    for record in records:
        handle.write(
            f'{record["sha256"]}  '
            f'{record["relative_path"]}\n'
        )

metadata = {
    "finalized_utc": (
        datetime.now(timezone.utc)
        .isoformat()
    ),
    "asset_root": str(root),
    "files_in_manifest": len(records),
    "figures": sum(
        record["relative_path"]
        .startswith("figures/")
        for record in records
    ),
    "tables": sum(
        record["relative_path"]
        .startswith("tables/")
        for record in records
    ),
}

metadata_path = (
    audit_root
    / "finalized_manuscript_asset_metadata.json"
)

metadata_path.write_text(
    json.dumps(
        metadata,
        indent=2,
    )
    + "\n",
    encoding="utf-8",
)

print("=" * 80)
print("FINAL MANUSCRIPT-ASSET MANIFEST")
print("=" * 80)
print(json.dumps(metadata, indent=2))
print()
print("Manifest:", manifest_path)
print("Checksums:", checksum_path)
print()
print("PASS: complete asset manifest regenerated.")
