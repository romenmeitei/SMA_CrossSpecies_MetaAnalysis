from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

import yaml


@dataclass(frozen=True)
class ProjectPaths:
    root: Path

    @property
    def config(self) -> Path:
        return self.root / "config"

    @property
    def data(self) -> Path:
        return self.root / "data"

    @property
    def results(self) -> Path:
        return self.root / "results"

    @property
    def resources(self) -> Path:
        return self.root / "resources"

    @classmethod
    def discover(cls, start: str | Path | None = None) -> "ProjectPaths":
        here = Path(start or Path.cwd()).resolve()
        for candidate in (here, *here.parents):
            if (candidate / "config" / "dataset_manifest.csv").exists():
                return cls(candidate)
        raise FileNotFoundError("Could not find config/dataset_manifest.csv in this directory or its parents")


def load_parameters(path: str | Path) -> dict[str, Any]:
    with Path(path).open("r", encoding="utf-8") as handle:
        params = yaml.safe_load(handle)
    if not isinstance(params, dict):
        raise ValueError("Analysis parameters must be a YAML mapping")
    return params

