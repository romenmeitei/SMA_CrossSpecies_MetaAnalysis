import numpy as np

from sma_xspecies.stats import random_effects_reml, weighted_stouffer


def test_random_effects_direction_and_interval():
    result = random_effects_reml(np.array([0.5, 0.7, 0.3]), np.array([0.1, 0.2, 0.15]))
    assert result.k == 3
    assert 0.3 < result.estimate < 0.7
    assert result.ci_low < result.estimate < result.ci_high
    assert 0 <= result.p_value <= 1
    assert result.tau2 >= 0


def test_weighted_stouffer_preserves_direction():
    z, p = weighted_stouffer(np.array([2.0, 3.0]), np.array([10, 20]))
    assert z > 0
    assert 0 <= p < 0.05

