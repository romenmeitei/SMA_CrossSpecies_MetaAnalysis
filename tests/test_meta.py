from pathlib import Path

import pandas as pd

from sma_xspecies.meta import ortholog_mapping_audit, run_gene_meta


def test_cross_species_meta_uses_versioned_one_to_one_mapping(tmp_path: Path):
    effects = pd.DataFrame([
        {"accession": "H1", "gene_id": "legacy_human_id", "gene_symbol": "A", "effect": 0.5, "se": 0.1, "z": 5.0, "p_value": 0.01, "fdr": 0.02, "species": "Homo_sapiens", "tissue": "motor_neuron", "role": "x", "eligibility": "eligible", "contrast": "disease-control"},
        {"accession": "M1", "gene_id": "ENSMUSG1", "gene_symbol": "A", "effect": 0.7, "se": 0.2, "z": 3.5, "p_value": 0.01, "fdr": 0.02, "species": "Mus_musculus", "tissue": "motor_neuron", "role": "x", "eligibility": "eligible", "contrast": "disease-control"},
    ])
    mapping_path = tmp_path / "mapping.tsv"
    pd.DataFrame({
        "human_ensembl_gene_id": ["ENSG1"], "human_symbol": ["A"],
        "mouse_ensembl_gene_id": ["ENSMUSG1"], "mouse_symbol": ["A"],
        "orthology_type": ["ortholog_one2one"], "ensembl_release": [113],
    }).to_csv(mapping_path, sep="\t", index=False)
    result = run_gene_meta(effects, mapping_path)
    assert len(result) == 1
    assert result.loc[0, "k"] == 2
    audit = ortholog_mapping_audit(effects, mapping_path)
    assert audit["n_mapped_genes"].sum() == 2

