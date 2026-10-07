import pandas as pd

from sma_xspecies.rescue import treatment_reversal


def test_reversal_requires_opposite_sign_and_both_fdrs():
    disease = pd.DataFrame({
        "gene_id": ["g1", "g2", "g3"], "gene_symbol": ["A", "B", "C"],
        "effect": [1.0, 1.0, -1.0], "z": [3.0, 3.0, -3.0], "fdr": [0.01, 0.01, 0.2],
    })
    treatment = pd.DataFrame({
        "gene_id": ["g1", "g2", "g3"], "gene_symbol": ["A", "B", "C"],
        "effect": [-0.8, 0.8, 0.8], "z": [-2.5, 2.5, 2.5], "fdr": [0.02, 0.02, 0.02],
    })
    out = treatment_reversal(disease, treatment)
    passed = set(out.loc[out["passes_conjunction"], "gene_id"])
    assert passed == {"g1"}
    assert out.loc[out["gene_id"] == "g1", "reversal_class"].item() == "restored"

