import numpy as np
import pandas as pd

from sma_xspecies.standardize import standardize_de_frame


def test_standardization_derives_signed_z_without_se():
    raw = pd.DataFrame({
        "test_id": ["ENSG000001.4", "ENSG000002.1"],
        "gene": ["A", "B"],
        "log2(fold_change)": [1.0, -2.0],
        "p_value": [0.01, 0.001],
        "q_value": [0.02, 0.002],
    })
    out = standardize_de_frame(raw, accession="TEST", species="Homo_sapiens", tissue="motor_neuron", role="test", eligibility="eligible")
    assert out["gene_id"].tolist() == ["ENSG000001", "ENSG000002"]
    assert out.loc[out["gene_symbol"] == "A", "z"].item() > 0
    assert out.loc[out["gene_symbol"] == "B", "z"].item() < 0
    assert np.isfinite(out["z"]).all()

