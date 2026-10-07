from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

from .config import ProjectPaths, load_parameters
from .io import download_manifest_files, file_inventory
from .manifest import assert_manifest_valid, read_manifest
from .meta import leave_one_dataset_out, ortholog_mapping_audit, read_effect_files, run_gene_meta, species_concordance
from .report import write_result_index, write_run_record
from .rescue import compare_primary_and_mouse_level, treatment_reversal
from .splicing import significant_events, standardize_rmats
from .standardize import prepare_expression_matrix, read_table, standardize_published_summary
from .utils import write_json


def _paths(root: str | None) -> ProjectPaths:
    return ProjectPaths.discover(root)


def cmd_validate(args: argparse.Namespace) -> int:
    paths = _paths(args.root)
    issues = assert_manifest_valid(paths.config / "dataset_manifest.csv", paths.config / "analysis_parameters.yaml")
    destination = paths.data / "audit" / "manifest_validation.csv"
    destination.parent.mkdir(parents=True, exist_ok=True)
    issues.to_csv(destination, index=False)
    print(issues.to_string(index=False))
    print(f"\nWrote {destination}")
    return 0


def cmd_download(args: argparse.Namespace) -> int:
    paths = _paths(args.root)
    manifest = read_manifest(paths.config / "dataset_manifest.csv")
    audit = download_manifest_files(manifest, paths.data / "raw", include_large=args.include_large)
    destination = paths.data / "audit" / "download_audit.csv"
    destination.parent.mkdir(parents=True, exist_ok=True)
    audit.to_csv(destination, index=False)
    file_inventory(paths.data / "raw").to_csv(paths.data / "audit" / "raw_file_inventory.csv", index=False)
    print(audit.to_string(index=False))
    return int((audit.get("status", pd.Series(dtype=str)) == "failed").any())


def cmd_standardize(args: argparse.Namespace) -> int:
    paths = _paths(args.root)
    manifest = read_manifest(paths.config / "dataset_manifest.csv")
    outdir = paths.results / "de"
    outdir.mkdir(parents=True, exist_ok=True)
    records = []
    for _, row in manifest.loc[manifest["enabled"]].iterrows():
        if row["analysis_engine"] != "published_cuffdiff_summary":
            records.append({"accession": row["accession"], "status": "deferred_to_dataset_engine", "message": row["analysis_engine"]})
            continue
        candidates = list((paths.data / "raw" / row["accession"]).glob("*"))
        candidates = [path for path in candidates if "gene_exp" in path.name and "MATS" not in path.name]
        if not candidates:
            records.append({"accession": row["accession"], "status": "missing_input", "message": "Run the download stage"})
            continue
        try:
            table = standardize_published_summary(candidates[0], row)
            destination = outdir / f"{row['accession']}_disease-control.csv"
            table.to_csv(destination, index=False)
            records.append({"accession": row["accession"], "status": "standardized", "message": str(destination)})
        except Exception as summary_error:
            try:
                annotations = pd.read_csv(paths.config / "sample_annotations" / f"{row['accession']}.csv")
                matrix = prepare_expression_matrix(read_table(candidates[0]), annotations)
                destination = paths.data / "processed" / f"{row['accession']}_expression.tsv"
                destination.parent.mkdir(parents=True, exist_ok=True)
                matrix.to_csv(destination, sep="\t", index=False)
                records.append({"accession": row["accession"], "status": "prepared_expression_matrix", "message": str(destination)})
            except Exception as matrix_error:
                records.append({"accession": row["accession"], "status": "failed", "message": f"summary: {summary_error}; matrix: {matrix_error}"})
    audit = pd.DataFrame(records)
    audit.to_csv(paths.data / "audit" / "standardization_audit.csv", index=False)
    print(audit.to_string(index=False))
    return int((audit["status"] == "failed").any())


def cmd_meta(args: argparse.Namespace) -> int:
    paths = _paths(args.root)
    params = load_parameters(paths.config / "analysis_parameters.yaml")
    effects = read_effect_files(paths.results / "de")
    mapping = paths.root / params["orthology"]["mapping_file"]
    if not mapping.exists():
        raise FileNotFoundError(f"Missing {mapping}. Run Rscript R/fetch_orthologs.R first.")
    strict = not args.include_sensitivity_only
    tissue = run_gene_meta(effects, mapping, minimum_studies=params["meta_analysis"]["minimum_studies"], strict=strict, neural_joint=False)
    neural = run_gene_meta(effects, mapping, minimum_studies=params["meta_analysis"]["minimum_studies"], strict=strict, neural_joint=True)
    loo = leave_one_dataset_out(effects, mapping, strict=strict)
    concordance = species_concordance(effects, mapping, strict=strict)
    mapping_audit = ortholog_mapping_audit(effects, mapping, strict=strict)
    outdir = paths.results / "meta"
    outdir.mkdir(parents=True, exist_ok=True)
    suffix = "_with_sensitivity_studies" if not strict else ""
    tissue.to_csv(outdir / f"gene_meta_by_tissue{suffix}.csv", index=False)
    neural.to_csv(outdir / f"gene_meta_neural_joint{suffix}.csv", index=False)
    loo.to_csv(outdir / f"leave_one_dataset_out{suffix}.csv", index=False)
    concordance.to_csv(outdir / f"species_concordance{suffix}.csv", index=False)
    mapping_audit.to_csv(outdir / f"ortholog_mapping_audit{suffix}.csv", index=False)
    print(f"Wrote {len(tissue):,} tissue-stratified and {len(neural):,} neural-joint meta-analysis rows")
    return 0


def _read_contrast(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(path)
    return pd.read_csv(path)


def cmd_rescue(args: argparse.Namespace) -> int:
    paths = _paths(args.root)
    params = load_parameters(paths.config / "analysis_parameters.yaml")["treatment_reversal"]
    source = paths.results / "de"
    disease = _read_contrast(source / "GSE207890_disease-control_primary.csv")
    treated = _read_contrast(source / "GSE207890_treated-disease_primary.csv")
    result = treatment_reversal(disease, treated, fdr_threshold=params["fdr_each_contrast"], minimum_absolute_z=params["minimum_absolute_z"])
    outdir = paths.results / "rescue"
    outdir.mkdir(parents=True, exist_ok=True)
    result.to_csv(outdir / "GSE207890_treatment_reversal_primary.csv", index=False)
    sensitivity_disease = source / "GSE207890_disease-control_mouse_level.csv"
    sensitivity_treated = source / "GSE207890_treated-disease_mouse_level.csv"
    if sensitivity_disease.exists() and sensitivity_treated.exists():
        mouse = treatment_reversal(_read_contrast(sensitivity_disease), _read_contrast(sensitivity_treated), fdr_threshold=params["fdr_each_contrast"], minimum_absolute_z=params["minimum_absolute_z"])
        mouse.to_csv(outdir / "GSE207890_treatment_reversal_mouse_level.csv", index=False)
        compare_primary_and_mouse_level(result, mouse).to_csv(outdir / "GSE207890_reversal_robustness.csv", index=False)
    print(f"Wrote {len(result):,} gene-level reversal assessments")
    return 0


def cmd_splicing(args: argparse.Namespace) -> int:
    paths = _paths(args.root)
    params = load_parameters(paths.config / "analysis_parameters.yaml")["splicing"]
    candidates = list((paths.data / "raw" / "GSE108094").glob("*MATS*"))
    if not candidates:
        raise FileNotFoundError("GSE108094 rMATS file not found; run the download stage")
    events = standardize_rmats(candidates[0], "GSE108094", "SE")
    selected = significant_events(events, params["event_fdr"], params["minimum_absolute_delta_psi"])
    outdir = paths.results / "splicing"
    outdir.mkdir(parents=True, exist_ok=True)
    events.to_csv(outdir / "GSE108094_all_events.csv", index=False)
    selected.to_csv(outdir / "GSE108094_significant_events.csv", index=False)
    print(f"Wrote {len(events):,} events; {len(selected):,} pass the prespecified thresholds")
    return 0


def cmd_report(args: argparse.Namespace) -> int:
    paths = _paths(args.root)
    print(write_run_record(paths.root))
    print(write_result_index(paths.root))
    return 0


def cmd_enrichment(args: argparse.Namespace) -> int:
    paths = _paths(args.root)
    from .enrichment import run_preranked

    ranking = pd.read_csv(args.ranking)
    result = run_preranked(ranking, args.gmt, paths.results / "enrichment", permutations=args.permutations)
    print(f"Wrote {len(result):,} gene-set results")
    return 0


def cmd_network(args: argparse.Namespace) -> int:
    paths = _paths(args.root)
    from .network import degree_matched_edge_test, induced_subgraph_statistics, load_string_graph

    frame = pd.read_csv(args.proteins, sep=None, engine="python")
    column = args.column or frame.columns[0]
    proteins = set(frame[column].dropna().astype(str))
    graph = load_string_graph(args.links, minimum_score=args.minimum_score)
    stats = induced_subgraph_statistics(graph, proteins)
    stats.update(degree_matched_edge_test(graph, proteins, permutations=args.permutations))
    outdir = paths.results / "network"
    outdir.mkdir(parents=True, exist_ok=True)
    write_json(stats, outdir / "induced_network_statistics.json")
    nx_edges = pd.DataFrame(graph.subgraph(proteins).edges(), columns=["protein1", "protein2"])
    nx_edges.to_csv(outdir / "induced_network_edges.csv", index=False)
    print(stats)
    return 0


def cmd_figures(args: argparse.Namespace) -> int:
    paths = _paths(args.root)
    from .figures import generate_available_figures

    generated = generate_available_figures(paths.root)
    print("Generated:", *(path.name for path in generated), sep="\n  ")
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="sma-xspecies", description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    for command, function in (
        ("validate", cmd_validate), ("download", cmd_download), ("standardize", cmd_standardize),
        ("meta", cmd_meta), ("rescue", cmd_rescue), ("splicing", cmd_splicing),
        ("enrichment", cmd_enrichment), ("network", cmd_network), ("figures", cmd_figures), ("report", cmd_report),
    ):
        child = subparsers.add_parser(command)
        child.add_argument("--root", default=None, help="Project root; auto-detected by default")
        child.set_defaults(function=function)
        if command == "download":
            child.add_argument("--include-large", action="store_true", help="Also download >200 MB raw archives")
        if command == "meta":
            child.add_argument("--include-sensitivity-only", action="store_true")
        if command == "enrichment":
            child.add_argument("--ranking", required=True, help="CSV with gene_symbol and score")
            child.add_argument("--gmt", required=True, help="Versioned GMT file")
            child.add_argument("--permutations", type=int, default=1000)
        if command == "network":
            child.add_argument("--proteins", required=True, help="CSV/TSV containing STRING protein IDs")
            child.add_argument("--column", default=None)
            child.add_argument("--links", required=True, help="Complete STRING protein links file")
            child.add_argument("--minimum-score", type=int, default=700)
            child.add_argument("--permutations", type=int, default=1000)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        return int(args.function(args))
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
