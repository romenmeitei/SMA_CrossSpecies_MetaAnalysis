import numpy as np

from sma_xspecies.stats import (
    random_effects_reml,
)


def test_reml_finds_boundary_global_minimum():
    # AMPH study effects from the audited analysis.
    effects = np.array([
        -4.164483,
        0.252454,
        0.356822,
    ])

    standard_errors = np.array([
        1.650890,
        0.193298,
        0.110139,
    ])

    result = random_effects_reml(
        effects,
        standard_errors,
    )

    assert result.tau2 <= 1e-8
    assert np.isclose(
        result.estimate,
        0.316187,
        atol=1e-4,
    )


def test_reml_retains_positive_heterogeneity():
    effects = np.array([
        -2.0,
        0.0,
        2.0,
    ])

    standard_errors = np.array([
        0.1,
        0.1,
        0.1,
    ])

    result = random_effects_reml(
        effects,
        standard_errors,
    )

    assert result.tau2 > 0.1
    assert np.isclose(
        result.estimate,
        0.0,
        atol=1e-8,
    )
