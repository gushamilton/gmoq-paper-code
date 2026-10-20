"""Data-generating models and exact population targets for the simulations.

What: draw(case, n, seed) simulates genotypes and phenotypes under a normal model (genotype
  effects on the mean, a, and on the log-variance, c) or a two-component normal mixture
  (effects on membership, gamma, and on the lower and upper component locations, delta1 and
  delta2). population_truth(case, taus) computes the exact marginal quantiles, the
  population RIF slopes at the 45 quantile levels and their OLS projection onto the model's
  weight columns: the target that Stage 2 estimates. fit_baseline fits the null model to
  one sample. None of this uses the G-MoQ package, so the targets are independent of it.
Paper: Supplementary Figs. S3-S6; Supplementary Tables S1, S2 (studies 3-7)
Input: none
Output: none (imported by the drivers)
Runs: anywhere, imported

Mixture fitting uses the same bounded, ordered parameterisation for sample and population
fits. The population likelihood is an expected log likelihood computed by deterministic
Gauss-Hermite quadrature, repeated at increasing orders; the fit diagnostics are returned
in population_truth rather than silently accepting an unrefined or failed fit.
"""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence
from typing import Any

import numpy as np

try:  # SciPy is part of the validation environment, but fail explicitly if absent.
    from scipy.optimize import minimize
    from scipy.special import ndtr, ndtri
except Exception:  # pragma: no cover - exercised only by a missing environment
    minimize = None
    ndtr = None
    ndtri = None


_LOG_2PI = math.log(2.0 * math.pi)
_MIXTURE_P_LO = 0.02
_MIXTURE_P_HI = 0.98
_MIXTURE_SD_LO = 0.1
_MIXTURE_SD_HI = 10.0
_MIXTURE_MEAN_GAP = 0.05
_ROOT_TOL = 1.0e-10
_DENSITY_FLOOR = 1.0e-12


class TruthNumericalError(ValueError):
    """Raised when a truth or nuisance calculation cannot meet its contract."""


# ---------------------------------------------------------------------------
# Case validation and elementary model functions


def _finite_float(value: Any, name: str) -> float:
    try:
        result = float(value)
    except (TypeError, ValueError) as exc:
        raise TruthNumericalError(f"{name} must be numeric, got {value!r}") from exc
    if not np.isfinite(result):
        raise TruthNumericalError(f"{name} must be finite, got {value!r}")
    return result


def _case(case: Mapping[str, Any]) -> dict[str, Any]:
    if not isinstance(case, Mapping):
        raise TruthNumericalError("case must be a mapping")

    family = str(case.get("family", "")).lower()
    if family not in {"normal", "mixture"}:
        raise TruthNumericalError("case['family'] must be 'normal' or 'mixture'")

    # The non-applicable effects are zero by contract.  The mixture baseline is
    # the fixed validation baseline; explicit values are still accepted so the
    # truth code can be used for bounded follow-up cases.
    out: dict[str, Any] = {
        "id": str(case.get("id", family)),
        "family": family,
        "maf": _finite_float(case.get("maf", 0.30), "maf"),
        "a": _finite_float(case.get("a", 0.0), "a"),
        "c": _finite_float(case.get("c", 0.0), "c"),
        "gamma": _finite_float(case.get("gamma", 0.0), "gamma"),
        "delta1": _finite_float(case.get("delta1", 0.0), "delta1"),
        "delta2": _finite_float(case.get("delta2", 0.0), "delta2"),
        "pi0": _finite_float(case.get("pi0", 0.5), "pi0"),
        "mu1": _finite_float(case.get("mu1", -2.0), "mu1"),
        "mu2": _finite_float(case.get("mu2", 2.0), "mu2"),
        "sd1": _finite_float(case.get("sd1", 1.0), "sd1"),
        "sd2": _finite_float(case.get("sd2", 1.0), "sd2"),
    }
    if not (0.0 < out["maf"] <= 0.5):
        raise TruthNumericalError("maf must be in (0, 0.5]")
    if family == "mixture":
        if not (0.0 < out["pi0"] < 1.0):
            raise TruthNumericalError("pi0 must be strictly between zero and one")
        if out["sd1"] <= 0.0 or out["sd2"] <= 0.0:
            raise TruthNumericalError("mixture component SDs must be positive")

        # The validation design identifies component 1 by its lower baseline
        # mean and requires the alternatives to retain that ordering at each of
        # the three genotype values.  Refuse a label crossing rather than
        # silently relabelling a requested data-generating case.
        for x in (0.0 - 2.0 * out["maf"], 1.0 - 2.0 * out["maf"], 2.0 - 2.0 * out["maf"]):
            m1 = out["mu1"] + out["delta1"] * x
            m2 = out["mu2"] + out["delta2"] * x
            if not m1 + 1.0e-12 < m2:
                raise TruthNumericalError(
                    "mixture component means cross or violate the required lower-label "
                    f"ordering at x={x:g}: {m1:g} >= {m2:g}"
                )
    return out


def _expit(x: np.ndarray | float) -> np.ndarray | float:
    arr = np.asarray(x, dtype=float)
    # This branch avoids overflow for deterministic fits with large trial values.
    result = np.empty_like(arr)
    positive = arr >= 0.0
    result[positive] = 1.0 / (1.0 + np.exp(-arr[positive]))
    e = np.exp(arr[~positive])
    result[~positive] = e / (1.0 + e)
    if np.ndim(x) == 0:
        return float(result)
    return result


def _logit(p: float) -> float:
    return math.log(p / (1.0 - p))


def _normal_cdf(y: np.ndarray | float, mean: np.ndarray | float, sd: np.ndarray | float) -> np.ndarray:
    z = (np.asarray(y, dtype=float) - np.asarray(mean, dtype=float)) / np.asarray(sd, dtype=float)
    if ndtr is not None:
        return np.asarray(ndtr(z), dtype=float)
    # Small fallback for truth evaluation.  Optimisation still requires SciPy.
    flat = np.ravel(z)
    values = [0.5 * (1.0 + math.erf(float(v) / math.sqrt(2.0))) for v in flat]
    return np.asarray(values, dtype=float).reshape(z.shape)


def _normal_pdf(y: np.ndarray | float, mean: np.ndarray | float, sd: np.ndarray | float) -> np.ndarray:
    z = (np.asarray(y, dtype=float) - np.asarray(mean, dtype=float)) / np.asarray(sd, dtype=float)
    return np.exp(-0.5 * z * z) / (np.asarray(sd, dtype=float) * math.sqrt(2.0 * math.pi))


def _normal_ppf(taus: np.ndarray) -> np.ndarray:
    if ndtri is not None:
        return np.asarray(ndtri(taus), dtype=float)
    # statistics.NormalDist is deterministic and adequate for the fixed interior grid.
    from statistics import NormalDist

    return np.asarray([NormalDist().inv_cdf(float(t)) for t in taus], dtype=float)


def _genotypes(c: Mapping[str, Any]) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    p = float(c["maf"])
    g = np.arange(3, dtype=float)
    probs = np.array(
        [(1.0 - p) ** 2, 2.0 * p * (1.0 - p), p**2], dtype=float
    )
    x = g - 2.0 * p
    return g, x, probs


def _conditional_components(c: Mapping[str, Any]) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Return genotype x, component means/scales and P(L=1|g)."""
    _, x, probs = _genotypes(c)
    if c["family"] == "normal":
        means = (float(c["a"]) * x)[:, None]
        scales = np.exp((float(c["c"]) * x / 2.0))[:, None]
        if not np.all(np.isfinite(scales)) or np.any(scales <= 0.0):
            raise TruthNumericalError("normal conditional SDs are not finite and positive")
        return x, means, scales, np.ones(3, dtype=float)

    logit_pi = _logit(float(c["pi0"])) + float(c["gamma"]) * x
    pi = np.asarray(_expit(logit_pi), dtype=float)
    means = np.column_stack(
        [float(c["mu1"]) + float(c["delta1"]) * x,
         float(c["mu2"]) + float(c["delta2"]) * x]
    )
    scales = np.column_stack(
        [np.full(3, float(c["sd1"])), np.full(3, float(c["sd2"]))]
    )
    return x, means, scales, pi


def _marginal_cdf_density(c: Mapping[str, Any], y: float) -> tuple[float, float, np.ndarray]:
    """Full marginal CDF/density and conditional genotype CDFs at scalar y."""
    _, x, probs = _genotypes(c)
    _, means, scales, pi = _conditional_components(c)
    if c["family"] == "normal":
        Fg = _normal_cdf(y, means[:, 0], scales[:, 0])
        fg = _normal_pdf(y, means[:, 0], scales[:, 0])
    else:
        F1 = _normal_cdf(y, means[:, 0], scales[:, 0])
        F2 = _normal_cdf(y, means[:, 1], scales[:, 1])
        f1 = _normal_pdf(y, means[:, 0], scales[:, 0])
        f2 = _normal_pdf(y, means[:, 1], scales[:, 1])
        Fg = pi * F1 + (1.0 - pi) * F2
        fg = pi * f1 + (1.0 - pi) * f2
    F = float(np.dot(probs, Fg))
    f = float(np.dot(probs, fg))
    if not (np.isfinite(F) and np.isfinite(f)):
        raise TruthNumericalError(f"non-finite marginal CDF/density at y={y!r}")
    return F, f, np.asarray(Fg, dtype=float)


def _quantile(c: Mapping[str, Any], tau: float) -> tuple[float, float, float, np.ndarray]:
    if not (0.0 < tau < 1.0) or not np.isfinite(tau):
        raise TruthNumericalError(f"tau must be strictly inside (0, 1), got {tau!r}")
    _, means, scales, _ = _conditional_components(c)
    lo = float(np.min(means - 14.0 * scales))
    hi = float(np.max(means + 14.0 * scales))
    Flo = _marginal_cdf_density(c, lo)[0]
    Fhi = _marginal_cdf_density(c, hi)[0]
    for _ in range(60):
        if Flo <= tau and Fhi >= tau:
            break
        width = max(hi - lo, 1.0)
        if Flo > tau:
            lo -= width
            Flo = _marginal_cdf_density(c, lo)[0]
        if Fhi < tau:
            hi += width
            Fhi = _marginal_cdf_density(c, hi)[0]
    else:
        raise TruthNumericalError(f"could not bracket marginal quantile tau={tau:g}")

    for _ in range(260):
        mid = 0.5 * (lo + hi)
        Fmid = _marginal_cdf_density(c, mid)[0]
        if Fmid < tau:
            lo = mid
        else:
            hi = mid
        if hi - lo <= 2.0e-14 * max(1.0, abs(mid)):
            break
    candidates = [0.5 * (lo + hi), lo, hi]
    q = min(candidates, key=lambda value: abs(_marginal_cdf_density(c, value)[0] - tau))
    F, f, Fg = _marginal_cdf_density(c, q)
    residual = abs(F - tau)
    if not np.isfinite(residual) or residual > _ROOT_TOL:
        raise TruthNumericalError(
            f"marginal quantile root residual {residual:.3g} exceeds {_ROOT_TOL:g} "
            f"at tau={tau:g}"
        )
    return float(q), float(f), float(residual), Fg


# ---------------------------------------------------------------------------
# Drawing and checking weights


def draw(case: Mapping[str, Any], n: int, seed: int) -> dict[str, np.ndarray]:
    """Draw one independent dataset from a validation case.

    The returned dictionary contains only ``g`` and ``y``.  For a mixture the
    latent component draw is deliberately private so estimator code cannot
    accidentally receive oracle memberships.
    """
    c = _case(case)
    if not isinstance(n, (int, np.integer)) or int(n) <= 0:
        raise TruthNumericalError("n must be a positive integer")
    if not isinstance(seed, (int, np.integer)):
        raise TruthNumericalError("seed must be an integer")
    rng = np.random.default_rng(int(seed))
    g = rng.binomial(2, float(c["maf"]), size=int(n)).astype(np.int8)
    x = g.astype(float) - 2.0 * float(c["maf"])
    if c["family"] == "normal":
        scale = np.exp(float(c["c"]) * x / 2.0)
        y = float(c["a"]) * x + rng.normal(size=int(n)) * scale
    else:
        p1 = np.asarray(_expit(_logit(float(c["pi0"])) + float(c["gamma"]) * x))
        latent_one = rng.binomial(1, p1, size=int(n))
        mean = np.where(
            latent_one == 1,
            float(c["mu1"]) + float(c["delta1"]) * x,
            float(c["mu2"]) + float(c["delta2"]) * x,
        )
        sd = np.where(latent_one == 1, float(c["sd1"]), float(c["sd2"]))
        y = mean + rng.normal(size=int(n)) * sd
    y = np.asarray(y, dtype=float)
    if not np.all(np.isfinite(y)):
        raise TruthNumericalError("draw produced non-finite phenotype values")
    return {"g": g, "y": y}


def _baseline_values(baseline: Mapping[str, Any], family: str) -> tuple[float, ...]:
    if not isinstance(baseline, Mapping):
        raise TruthNumericalError("baseline must be a mapping")
    if family == "normal":
        mu = _finite_float(baseline.get("mu"), "baseline.mu")
        sd = _finite_float(baseline.get("sd"), "baseline.sd")
        if sd <= 0.0:
            raise TruthNumericalError("baseline.sd must be positive")
        return mu, sd
    names = ("p1", "mu1", "sd1", "mu2", "sd2")
    vals = tuple(_finite_float(baseline.get(name), f"baseline.{name}") for name in names)
    p1, mu1, sd1, mu2, sd2 = vals
    if not (_MIXTURE_P_LO <= p1 <= _MIXTURE_P_HI):
        raise TruthNumericalError("baseline.p1 must lie in [0.02, 0.98]")
    if not (_MIXTURE_SD_LO <= sd1 <= _MIXTURE_SD_HI and _MIXTURE_SD_LO <= sd2 <= _MIXTURE_SD_HI):
        raise TruthNumericalError("baseline component SDs must lie in [0.1, 10]")
    if not mu1 + _MIXTURE_MEAN_GAP <= mu2:
        raise TruthNumericalError("baseline means must be ordered with gap at least 0.05")
    return vals


def independent_weights(
    taus: Sequence[float],
    q_tau: Sequence[float],
    baseline: Mapping[str, Any],
    family: str,
) -> np.ndarray:
    """Build the independent checking matrix using package conventions.

    For a Normal model the columns are ``(mu, sigma2)`` and the second column
    is ``qnorm(tau)/(2*sd)``.  This intentionally follows the pinned package's
    model-density convention; the supplied ``q_tau`` and ``mu`` do not replace
    the model quantile in that column.

    For a mixture the columns are ``(gamma, mu1, mu2)``.  ``gamma`` is a
    *logit* membership slope, hence the factor ``p1*(1-p1)``.
    """
    family = str(family).lower()
    if family not in {"normal", "mixture"}:
        raise TruthNumericalError("family must be 'normal' or 'mixture'")
    taus_arr = np.asarray(taus, dtype=float)
    q_arr = np.asarray(q_tau, dtype=float)
    if taus_arr.ndim != 1 or q_arr.ndim != 1 or taus_arr.size == 0 or taus_arr.size != q_arr.size:
        raise TruthNumericalError("taus and q_tau must be non-empty one-dimensional arrays of equal length")
    if not np.all(np.isfinite(taus_arr)) or np.any((taus_arr <= 0.0) | (taus_arr >= 1.0)):
        raise TruthNumericalError("taus must be finite and strictly inside (0, 1)")
    if not np.all(np.isfinite(q_arr)):
        raise TruthNumericalError("q_tau must be finite")

    if family == "normal":
        _, sd = _baseline_values(baseline, family)
        return np.column_stack([np.ones(taus_arr.size), _normal_ppf(taus_arr) / (2.0 * sd)])

    p1, mu1, sd1, mu2, sd2 = _baseline_values(baseline, family)
    F1 = _normal_cdf(q_arr, mu1, sd1)
    F2 = _normal_cdf(q_arr, mu2, sd2)
    f1 = _normal_pdf(q_arr, mu1, sd1)
    f2 = _normal_pdf(q_arr, mu2, sd2)
    f = p1 * f1 + (1.0 - p1) * f2
    if np.any(~np.isfinite(f)) or np.any(f <= _DENSITY_FLOOR):
        raise TruthNumericalError("mixture checking weights encounter a non-positive density")
    left = p1 * f1 / f
    right = (1.0 - p1) * f2 / f
    membership = -p1 * (1.0 - p1) * (F1 - F2) / f
    return np.column_stack([membership, left, right])


# ---------------------------------------------------------------------------
# Bounded mixture nuisance fitting


def _softplus(x: float) -> float:
    if x > 35.0:
        return float(x)
    if x < -35.0:
        return float(math.exp(x))
    return float(math.log1p(math.exp(x)))


def _bounded_logit(value: float, lo: float, hi: float) -> float:
    fraction = (float(value) - lo) / (hi - lo)
    fraction = min(max(fraction, 1.0e-8), 1.0 - 1.0e-8)
    return _logit(fraction)


def _decode_mixture(z: np.ndarray) -> dict[str, float]:
    z = np.asarray(z, dtype=float)
    if z.shape != (5,) or not np.all(np.isfinite(z)):
        raise TruthNumericalError("invalid transformed mixture parameter vector")
    p1 = _MIXTURE_P_LO + (_MIXTURE_P_HI - _MIXTURE_P_LO) * float(_expit(z[0]))
    mu1 = float(z[1])
    gap = _MIXTURE_MEAN_GAP + _softplus(float(z[2]))
    sd1 = _MIXTURE_SD_LO + (_MIXTURE_SD_HI - _MIXTURE_SD_LO) * float(_expit(z[3]))
    sd2 = _MIXTURE_SD_LO + (_MIXTURE_SD_HI - _MIXTURE_SD_LO) * float(_expit(z[4]))
    result = {"p1": p1, "mu1": mu1, "sd1": sd1, "mu2": mu1 + gap, "sd2": sd2}
    if not all(np.isfinite(v) for v in result.values()):
        raise TruthNumericalError("non-finite decoded mixture parameters")
    return result


def _encode_mixture(params: Mapping[str, float]) -> np.ndarray:
    p1, mu1, sd1, mu2, sd2 = _baseline_values(params, "mixture")
    return np.array(
        [
            _bounded_logit(p1, _MIXTURE_P_LO, _MIXTURE_P_HI),
            mu1,
            math.log(math.expm1(max(mu2 - mu1 - _MIXTURE_MEAN_GAP, 1.0e-8))),
            _bounded_logit(sd1, _MIXTURE_SD_LO, _MIXTURE_SD_HI),
            _bounded_logit(sd2, _MIXTURE_SD_LO, _MIXTURE_SD_HI),
        ],
        dtype=float,
    )


def _candidate_logpdf(y: np.ndarray, params: Mapping[str, float]) -> np.ndarray:
    p1, mu1, sd1, mu2, sd2 = _baseline_values(params, "mixture")
    y = np.asarray(y, dtype=float)
    z1 = (y - mu1) / sd1
    z2 = (y - mu2) / sd2
    a = math.log(p1) - 0.5 * (z1 * z1 + _LOG_2PI) - math.log(sd1)
    b = math.log1p(-p1) - 0.5 * (z2 * z2 + _LOG_2PI) - math.log(sd2)
    return np.logaddexp(a, b)


def _mixture_starts(mean: float, sd: float, quantiles: Sequence[float], extra: Sequence[Mapping[str, float]] = ()) -> list[dict[str, float]]:
    scale = max(float(sd), 0.25)
    q = list(float(v) for v in quantiles if np.isfinite(v))
    centre = float(mean)
    pairs: list[tuple[float, float]] = [
        (centre - 0.80 * scale, centre + 0.80 * scale),
        (centre - 1.00 * scale, centre + 1.00 * scale),
    ]
    if len(q) >= 4:
        pairs.extend([(q[0], q[-1]), (q[1], q[-2])])
    if len(q) >= 3:
        pairs.append((q[0], q[2]))
    starts: list[dict[str, float]] = []
    for left, right in pairs:
        left, right = sorted((float(left), float(right)))
        if right - left < _MIXTURE_MEAN_GAP:
            right = left + _MIXTURE_MEAN_GAP + 0.25 * scale
        for p1 in (0.25, 0.5, 0.75):
            for sd1, sd2 in ((0.55 * scale, 0.55 * scale), (0.40 * scale, 0.75 * scale)):
                starts.append(
                    {
                        "p1": p1,
                        "mu1": left,
                        "sd1": min(max(sd1, _MIXTURE_SD_LO + 1.0e-5), _MIXTURE_SD_HI - 1.0e-5),
                        "mu2": right,
                        "sd2": min(max(sd2, _MIXTURE_SD_LO + 1.0e-5), _MIXTURE_SD_HI - 1.0e-5),
                    }
                )
    # Put explicitly supplied deterministic starts first (the population fit
    # supplies the generating baseline as one such start).  They are starts,
    # not fixed nuisance values: the objective still selects the maximiser.
    starts = [dict(v) for v in extra] + starts
    # Stable de-duplication keeps the deterministic start count bounded.
    unique: list[dict[str, float]] = []
    seen: set[tuple[float, ...]] = set()
    for start in starts:
        try:
            key = tuple(round(float(start[k]), 10) for k in ("p1", "mu1", "sd1", "mu2", "sd2"))
            _baseline_values(start, "mixture")
        except TruthNumericalError:
            continue
        if key not in seen:
            seen.add(key)
            unique.append(start)
    if not unique:
        raise TruthNumericalError("could not construct a valid deterministic mixture start set")
    return unique[:24]


def _fit_mixture_objective(
    objective: Any,
    starts: Sequence[Mapping[str, float]],
    maxiter: int = 800,
) -> tuple[dict[str, float], dict[str, Any]]:
    if minimize is None:
        raise TruthNumericalError("mixture fitting requires scipy.optimize")
    # Only the lower transformed bounds are required for p and SD.  In
    # particular, do not impose an undocumented upper bound on the ordered
    # component mean gap.
    bounds = [(-12.0, 12.0), (None, None), (-12.0, None), (-12.0, 12.0), (-12.0, 12.0)]
    best: dict[str, float] | None = None
    best_value = math.inf
    converged = 0
    failures: list[str] = []
    for start in starts:
        try:
            z0 = _encode_mixture(start)
            result = minimize(
                objective,
                z0,
                method="L-BFGS-B",
                bounds=bounds,
                options={"maxiter": int(maxiter), "ftol": 1.0e-13, "gtol": 1.0e-8, "maxls": 40},
            )
            value = float(result.fun)
            params = _decode_mixture(result.x)
            if not np.isfinite(value) or not result.success:
                failures.append(str(result.message))
                continue
            converged += 1
            if value < best_value - 1.0e-12:
                best_value = value
                best = params
        except (TruthNumericalError, FloatingPointError, OverflowError, ValueError) as exc:
            failures.append(str(exc))
    if best is None:
        detail = "; ".join(failures[:3])
        raise TruthNumericalError(
            f"all {len(starts)} deterministic mixture starts failed to converge" +
            (f": {detail}" if detail else "")
        )
    return best, {
        "objective": float(best_value),
        "starts": int(len(starts)),
        "converged_starts": int(converged),
        "failed_starts": int(len(starts) - converged),
    }


def fit_baseline(y: Sequence[float], family: str) -> dict[str, float]:
    """Fit the declared marginal nuisance model to a replicate phenotype.

    Normal nuisance values are the sample mean and sample SD (``ddof=1``).
    Mixture fitting is an unlabelled two-Normal maximum likelihood fit with
    ordered means, a minimum gap, bounded mixing proportion and bounded SDs.
    The return value contains no observations or latent labels.
    """
    family = str(family).lower()
    if family not in {"normal", "mixture"}:
        raise TruthNumericalError("family must be 'normal' or 'mixture'")
    values = np.asarray(y, dtype=float)
    if values.ndim != 1 or values.size < 2 or not np.all(np.isfinite(values)):
        raise TruthNumericalError("y must be a finite one-dimensional array with at least two values")
    if family == "normal":
        sd = float(np.std(values, ddof=1))
        if not np.isfinite(sd) or sd <= 0.0:
            raise TruthNumericalError("sample normal SD is not finite and positive")
        return {"mu": float(np.mean(values)), "sd": sd}

    if values.size < 8:
        raise TruthNumericalError("mixture baseline fitting requires at least eight observations")
    mean = float(np.mean(values))
    sd = float(np.std(values, ddof=0))
    if not np.isfinite(sd) or sd <= 0.0:
        raise TruthNumericalError("mixture fitting requires positive sample variation")
    qs = np.quantile(values, [0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95], method="linear")
    starts = _mixture_starts(mean, sd, qs)

    def objective(z: np.ndarray) -> float:
        try:
            params = _decode_mixture(z)
            value = -float(np.mean(_candidate_logpdf(values, params)))
            return value if np.isfinite(value) else 1.0e300
        except (TruthNumericalError, FloatingPointError, ValueError):
            return 1.0e300

    params, _ = _fit_mixture_objective(objective, starts)
    # Revalidate after optimisation, including the declared constraints.
    _baseline_values(params, "mixture")
    return {key: float(params[key]) for key in ("p1", "mu1", "sd1", "mu2", "sd2")}


# ---------------------------------------------------------------------------
# Population expected-likelihood fit and truth projection


def _population_components(c: Mapping[str, Any]) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Unconditional component means, SDs and weights of the full DGP."""
    _, x, probs = _genotypes(c)
    _, means, scales, pi = _conditional_components(c)
    if c["family"] == "normal":
        return means[:, 0].copy(), scales[:, 0].copy(), probs.copy()
    all_means = np.column_stack([means[:, 0], means[:, 1]]).ravel()
    all_scales = np.column_stack([scales[:, 0], scales[:, 1]]).ravel()
    all_weights = np.column_stack([probs * pi, probs * (1.0 - pi)]).ravel()
    return all_means, all_scales, all_weights


def _population_moments(c: Mapping[str, Any]) -> tuple[float, float]:
    means, scales, weights = _population_components(c)
    mean = float(np.dot(weights, means))
    variance = float(np.dot(weights, scales * scales + means * means) - mean * mean)
    if not np.isfinite(variance) or variance <= 0.0:
        raise TruthNumericalError("population marginal variance is not finite and positive")
    return mean, math.sqrt(variance)


def _gauss_hermite_data(c: Mapping[str, Any], order: int) -> tuple[np.ndarray, np.ndarray]:
    if order < 16:
        raise TruthNumericalError("Gauss-Hermite order is too small")
    try:
        nodes, weights = np.polynomial.hermite.hermgauss(int(order))
    except Exception as exc:
        raise TruthNumericalError(f"Gauss-Hermite construction failed at order {order}") from exc
    means, scales, mixture_weights = _population_components(c)
    y = (means[:, None] + math.sqrt(2.0) * scales[:, None] * nodes[None, :]).ravel()
    w = (mixture_weights[:, None] * weights[None, :] / math.sqrt(math.pi)).ravel()
    if not np.all(np.isfinite(y)) or not np.all(np.isfinite(w)) or np.any(w < 0.0):
        raise TruthNumericalError(f"non-finite Gauss-Hermite data at order {order}")
    return y, w


def _population_mixture_fit(c: Mapping[str, Any]) -> tuple[dict[str, float], dict[str, Any]]:
    mean, sd = _population_moments(c)
    quantiles = [_quantile(c, tau)[0] for tau in (0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95)]
    # The generating baseline is merely a deterministic start, never the
    # returned answer unless it maximises the expected likelihood.
    oracle_start = {
        "p1": float(c["pi0"]),
        "mu1": float(c["mu1"]),
        "sd1": float(c["sd1"]),
        "mu2": float(c["mu2"]),
        "sd2": float(c["sd2"]),
    }
    starts = _mixture_starts(mean, sd, quantiles, extra=(oracle_start,))
    previous: dict[str, float] | None = None
    records: list[dict[str, Any]] = []

    for order in (64, 128, 256):
        y_nodes, w_nodes = _gauss_hermite_data(c, order)
        order_starts = list(starts)
        if previous is not None:
            order_starts.insert(0, previous)

        def objective(z: np.ndarray) -> float:
            try:
                params = _decode_mixture(z)
                logpdf = _candidate_logpdf(y_nodes, params)
                value = -float(np.dot(w_nodes, logpdf))
                return value if np.isfinite(value) else 1.0e300
            except (TruthNumericalError, FloatingPointError, ValueError):
                return 1.0e300

        fitted, details = _fit_mixture_objective(objective, order_starts)
        parameter_delta = None
        objective_delta = None
        if previous is not None:
            parameter_delta = max(abs(fitted[k] - previous[k]) for k in ("p1", "mu1", "sd1", "mu2", "sd2"))
            # Re-evaluate the previous answer on this order so the recorded
            # objective change is a quadrature refinement diagnostic, not an
            # optimiser-scale artefact.
            previous_objective = objective(_encode_mixture(previous))
            objective_delta = abs(float(details["objective"]) - float(previous_objective))
        record = dict(details)
        record.update({"order": int(order), "parameter_max_abs_change": parameter_delta,
                       "objective_abs_change": objective_delta})
        records.append(record)
        previous = fitted

    if previous is None:  # defensive; the loop always executes
        raise TruthNumericalError("population mixture fit produced no result")
    last = records[-1]
    refinement_ok = (
        len(records) >= 2
        and last["parameter_max_abs_change"] is not None
        and float(last["parameter_max_abs_change"]) <= 2.0e-4
        and last["objective_abs_change"] is not None
        and float(last["objective_abs_change"]) <= 2.0e-6
    )
    if not refinement_ok:
        raise TruthNumericalError(
            "population mixture expected-log-likelihood fit did not refine stably: "
            f"parameter change={last['parameter_max_abs_change']!r}, "
            f"objective change={last['objective_abs_change']!r}"
        )
    return {key: float(previous[key]) for key in ("p1", "mu1", "sd1", "mu2", "sd2")}, {
        "fit_kind": "expected_log_likelihood",
        "quadrature": "standard-normal Gauss-Hermite via hermgauss",
        "quadrature_orders": [64, 128, 256],
        "refinement": records,
        "refinement_converged": bool(refinement_ok),
    }


def _json_list(values: np.ndarray | Sequence[float]) -> list[float]:
    return [float(v) for v in np.asarray(values, dtype=float).ravel()]


def population_truth(case: Mapping[str, Any], taus: Sequence[float]) -> dict[str, Any]:
    """Compute the fixed population RIF projection target independently.

    ``b_tau`` is the exact three-genotype population RIF regression coefficient
    at each requested quantile.  ``theta`` is the OLS projection of that vector
    onto the independently constructed Stage 2 weight matrix.  The two are not
    replaced by the generating scientific coefficients.
    """
    c = _case(case)
    taus_arr = np.asarray(taus, dtype=float)
    if taus_arr.ndim != 1 or taus_arr.size == 0 or not np.all(np.isfinite(taus_arr)):
        raise TruthNumericalError("taus must be a non-empty finite one-dimensional sequence")
    if np.any((taus_arr <= 0.0) | (taus_arr >= 1.0)):
        raise TruthNumericalError("taus must be strictly inside (0, 1)")
    if np.any(np.diff(taus_arr) <= 0.0):
        raise TruthNumericalError("taus must be strictly increasing")

    _, x, probs = _genotypes(c)
    variance_x = float(np.dot(probs, x * x))
    if not np.isfinite(variance_x) or variance_x <= 0.0:
        raise TruthNumericalError("genotype variance is not finite and positive")

    q_values: list[float] = []
    f_values: list[float] = []
    b_values: list[float] = []
    root_residuals: list[float] = []
    for tau in taus_arr:
        q, f, residual, Fg = _quantile(c, float(tau))
        numerator = float(np.dot(probs * x, Fg))
        b = -numerator / (variance_x * f)
        if not np.isfinite(b) or f <= _DENSITY_FLOOR:
            raise TruthNumericalError(f"invalid population truth at tau={float(tau):g}")
        q_values.append(q)
        f_values.append(f)
        b_values.append(float(b))
        root_residuals.append(residual)

    if c["family"] == "normal":
        _, means, scales, _ = _conditional_components(c)
        baseline_mu = float(np.dot(probs, means[:, 0]))
        baseline_var = float(np.dot(probs, scales[:, 0] ** 2 + means[:, 0] ** 2) - baseline_mu**2)
        if baseline_var <= 0.0 or not np.isfinite(baseline_var):
            raise TruthNumericalError("normal population baseline variance is invalid")
        baseline = {"mu": baseline_mu, "sd": math.sqrt(baseline_var)}
        baseline_diagnostics: dict[str, Any] = {"fit_kind": "population_moments"}
        W = independent_weights(taus_arr, q_values, baseline, "normal")
        generating = {"a": float(c["a"]), "c": float(c["c"])}
        theta_names = ["a", "c"]
    else:
        baseline, baseline_diagnostics = _population_mixture_fit(c)
        W = independent_weights(taus_arr, q_values, baseline, "mixture")
        generating = {
            "gamma": float(c["gamma"]),
            "delta1": float(c["delta1"]),
            "delta2": float(c["delta2"]),
        }
        theta_names = ["gamma", "delta1", "delta2"]

    W = np.asarray(W, dtype=float)
    b_arr = np.asarray(b_values, dtype=float)
    singular_values = np.linalg.svd(W, compute_uv=False)
    rank = int(np.linalg.matrix_rank(W, tol=max(W.shape) * np.finfo(float).eps * singular_values[0]))
    condition = float(singular_values[0] / singular_values[-1]) if singular_values[-1] > 0.0 else math.inf
    if rank < W.shape[1] or not np.isfinite(condition) or condition > 1.0e12:
        raise TruthNumericalError(f"checking weight matrix is rank deficient or ill-conditioned: rank={rank}, condition={condition:g}")
    try:
        theta = np.linalg.solve(W.T @ W, W.T @ b_arr)
    except np.linalg.LinAlgError as exc:
        raise TruthNumericalError("could not solve the population OLS projection") from exc
    if not np.all(np.isfinite(theta)):
        raise TruthNumericalError("population projection target is non-finite")

    F_fit_at_q = []
    for q in q_values:
        if c["family"] == "normal":
            F_fit_at_q.append(float(_normal_cdf(q, baseline["mu"], baseline["sd"])))
        else:
            F1 = float(_normal_cdf(q, baseline["mu1"], baseline["sd1"]))
            F2 = float(_normal_cdf(q, baseline["mu2"], baseline["sd2"]))
            F_fit_at_q.append(float(baseline["p1"] * F1 + (1.0 - baseline["p1"]) * F2))

    diagnostics: dict[str, Any] = {
        "root_residual_max": float(max(root_residuals)),
        "root_residuals": [float(v) for v in root_residuals],
        "density_min": float(min(f_values)),
        "genotype_probabilities": _json_list(probs),
        "genotype_centered": _json_list(x),
        "genotype_variance": float(variance_x),
        "W_rank": rank,
        "W_condition": condition,
        "projection_residual_norm": float(np.linalg.norm(W @ theta - b_arr)),
        "F_fit_at_q_minus_tau": [float(v - t) for v, t in zip(F_fit_at_q, taus_arr)],
        "baseline_fit": baseline_diagnostics,
        "constraints": {
            "mixture_p1": [_MIXTURE_P_LO, _MIXTURE_P_HI],
            "mixture_sd": [_MIXTURE_SD_LO, _MIXTURE_SD_HI],
            "mixture_mean_gap": _MIXTURE_MEAN_GAP,
        },
    }
    return {
        "case_id": c["id"],
        "family": c["family"],
        "taus": _json_list(taus_arr),
        "q_tau": [float(v) for v in q_values],
        "f_tau": [float(v) for v in f_values],
        "b_tau": [float(v) for v in b_values],
        "baseline": {key: float(value) for key, value in baseline.items()},
        "W_independent": [[float(v) for v in row] for row in W],
        "theta": _json_list(theta),
        "theta_names": theta_names,
        "generating": generating,
        "diagnostics": diagnostics,
    }


__all__ = [
    "TruthNumericalError",
    "draw",
    "fit_baseline",
    "independent_weights",
    "population_truth",
]
