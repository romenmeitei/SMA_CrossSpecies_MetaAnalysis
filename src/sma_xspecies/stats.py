from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.optimize import minimize_scalar
from scipy.stats import chi2, norm


@dataclass(frozen=True)
class MetaResult:
    k: int
    estimate: float
    se: float
    ci_low: float
    ci_high: float
    z: float
    p_value: float
    tau2: float
    q: float
    q_p_value: float
    i2: float


def random_effects_reml(effects: np.ndarray, ses: np.ndarray) -> MetaResult:
    y = np.asarray(effects, dtype=float)
    s = np.asarray(ses, dtype=float)
    keep = np.isfinite(y) & np.isfinite(s) & (s > 0)
    y, s = y[keep], s[keep]
    if not len(y):
        raise ValueError("No finite effect/SE pairs")
    v = s**2
    fixed_w = 1.0 / v
    fixed_mean = np.sum(fixed_w * y) / np.sum(fixed_w)
    q = float(np.sum(fixed_w * (y - fixed_mean) ** 2))
    q_p = float(chi2.sf(q, max(len(y) - 1, 1))) if len(y) > 1 else np.nan
    i2 = max(0.0, (q - (len(y) - 1)) / q * 100.0) if len(y) > 1 and q > 0 else 0.0

    if len(y) == 1:
        tau2 = 0.0
    else:
        upper = max(
            float(np.var(y, ddof=1) * 20.0),
            float(np.max(v) * 100.0),
            1e-6,
        )

        def reml_nll(tau: float) -> float:
            if not np.isfinite(tau) or tau < 0:
                return np.inf

            total_variance = v + tau

            if (
                np.any(~np.isfinite(total_variance))
                or np.any(total_variance <= 0)
            ):
                return np.inf

            weights = 1.0 / total_variance
            weight_sum = np.sum(weights)

            if (
                not np.isfinite(weight_sum)
                or weight_sum <= 0
            ):
                return np.inf

            mean = np.sum(weights * y) / weight_sum
            residual = np.sum(
                weights * (y - mean) ** 2
            )

            return 0.5 * (
                np.sum(np.log(total_variance))
                + np.log(weight_sum)
                + residual
            )

        # A single bounded optimization is unsafe here:
        # with very small k and highly unequal study SEs,
        # the profile can contain a boundary optimum and
        # a separate interior basin. Search globally on a
        # logarithmic grid, then refine the best interval.
        lower_positive = max(
            upper * 1e-12,
            1e-12,
        )

        tau_grid = np.concatenate([
            np.array([0.0]),
            np.geomspace(
                lower_positive,
                upper,
                num=257,
            ),
        ])

        objective_grid = np.array([
            reml_nll(tau)
            for tau in tau_grid
        ])

        best_index = int(
            np.nanargmin(objective_grid)
        )

        candidate_tau2 = [
            0.0,
            float(tau_grid[best_index]),
            float(upper),
        ]

        if best_index > 0:
            local_lower = float(
                tau_grid[best_index - 1]
            )

            local_upper = float(
                tau_grid[
                    min(
                        len(tau_grid) - 1,
                        best_index + 1,
                    )
                ]
            )

            if local_upper > local_lower:
                local_fit = minimize_scalar(
                    reml_nll,
                    bounds=(
                        local_lower,
                        local_upper,
                    ),
                    method="bounded",
                    options={
                        "xatol": 1e-14,
                        "maxiter": 2000,
                    },
                )

                if (
                    local_fit.success
                    and np.isfinite(local_fit.x)
                ):
                    candidate_tau2.append(
                        float(local_fit.x)
                    )

        candidate_tau2 = [
            max(0.0, candidate)
            for candidate in candidate_tau2
            if (
                np.isfinite(candidate)
                and np.isfinite(
                    reml_nll(candidate)
                )
            )
        ]

        tau2 = float(
            min(
                candidate_tau2,
                key=reml_nll,
            )
        )

    weights = 1.0 / (v + tau2)
    estimate = float(np.sum(weights * y) / np.sum(weights))
    se = float(np.sqrt(1.0 / np.sum(weights)))
    z = estimate / se
    p = float(2.0 * norm.sf(abs(z)))
    return MetaResult(len(y), estimate, se, estimate - 1.96 * se, estimate + 1.96 * se, z, p, tau2, q, q_p, i2)


def weighted_stouffer(z_scores: np.ndarray, sample_sizes: np.ndarray | None = None) -> tuple[float, float]:
    z = np.asarray(z_scores, dtype=float)
    if sample_sizes is None:
        w = np.ones_like(z)
    else:
        n = np.asarray(sample_sizes, dtype=float)
        w = np.sqrt(n)
    keep = np.isfinite(z) & np.isfinite(w) & (w > 0)
    if not keep.any():
        raise ValueError("No finite z scores")
    combined = float(np.sum(w[keep] * z[keep]) / np.sqrt(np.sum(w[keep] ** 2)))
    return combined, float(2.0 * norm.sf(abs(combined)))

