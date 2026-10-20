"""Direct maximum-likelihood fits used as comparators.

What: fits the normal and two-component mixture models directly to genotype and phenotype
  by maximum likelihood (SciPy), with standard errors from the observed Hessian. Study 4
  uses fit_mixture_ml: the conditional mixture likelihood refitted at every variant, the
  comparator for G-MoQ's mixture estimates. The fits use only g, y and the minor allele
  frequency; generating effects and latent memberships are never read.
Paper: Supplementary Figs. S4, S5; Supplementary Table S2 (studies 4 and 5)
Input: none (imported by figures/mixture_comparison.sh)
Output: none
Runs: anywhere, imported. The DGLM adapter (fit_dglm) needs dglm.R, which was not used in
  the paper and is not included.
"""

from __future__ import annotations

import csv
import math
import os
import subprocess
import tempfile
import time
from pathlib import Path
from typing import Any, Callable, Mapping, Sequence

import numpy as np
from scipy.optimize import minimize
from scipy.special import expit, logit
from scipy.stats import levene, norm


_REQUIRED_ROW_KEYS = (
    "method",
    "param",
    "estimate",
    "se",
    "p_value",
    "target_kind",
    "status",
    "reason",
    "wall_seconds",
)
_NORMAL_PARAMS = ("a", "c", "variance_derivative")
_MIXTURE_PARAMS = ("gamma", "delta1", "delta2")


class ComparatorInputError(ValueError):
    """Raised when a comparator input cannot define a conditional model."""


class DGLMError(RuntimeError):
    """Raised only for a failed DGLM adapter invocation."""


# A small result object is represented as a dict on purpose: it is convenient
# for callers which want diagnostics, while fit_all only exposes row records.
def _fit_result(
    *,
    status: str,
    reason: str | None,
    estimates: Mapping[str, float] | None = None,
    standard_errors: Mapping[str, float] | None = None,
    p_values: Mapping[str, float] | None = None,
    covariance: np.ndarray | None = None,
    nuisance: Mapping[str, float] | None = None,
    diagnostics: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    return {
        "status": str(status),
        "reason": None if reason in (None, "") else str(reason),
        "estimates": dict(estimates or {}),
        "standard_errors": dict(standard_errors or {}),
        "p_values": dict(p_values or {}),
        "covariance": covariance,
        "nuisance": dict(nuisance or {}),
        "diagnostics": dict(diagnostics or {}),
    }


def _finite_float(value: Any) -> float | None:
    try:
        value = float(value)
    except (TypeError, ValueError):
        return None
    return value if math.isfinite(value) else None


def _json_value(value: Any) -> float | None:
    """Convert a scalar to a JSON-safe Python value."""
    return _finite_float(value)


def _row(
    method: str,
    param: str,
    *,
    estimate: Any = None,
    se: Any = None,
    p_value: Any = None,
    target_kind: str = "scientific",
    status: str = "ok",
    reason: str | None = None,
    wall_seconds: Any = 0.0,
) -> dict[str, Any]:
    # Do not return numpy scalars or NaN: output is intended to be directly
    # serialisable with json.dumps(allow_nan=False).
    return {
        "method": str(method),
        "param": str(param),
        "estimate": _json_value(estimate),
        "se": _json_value(se),
        "p_value": _json_value(p_value),
        "target_kind": str(target_kind),
        "status": str(status),
        "reason": None if reason in (None, "") else str(reason),
        "wall_seconds": _json_value(wall_seconds) if _json_value(wall_seconds) is not None else 0.0,
    }


def _validate_data(g: Sequence[Any], y: Sequence[Any], case: Mapping[str, Any]) -> tuple[np.ndarray, np.ndarray, np.ndarray, float, str]:
    """Validate inputs and return g, y, centred x, MAF and family."""
    try:
        ga = np.asarray(g, dtype=float)
        ya = np.asarray(y, dtype=float)
    except (TypeError, ValueError) as exc:
        raise ComparatorInputError(f"g and y must be numeric: {exc}") from exc
    if ga.ndim != 1 or ya.ndim != 1 or ga.size != ya.size:
        raise ComparatorInputError("g and y must be one-dimensional arrays of equal length")
    if ga.size < 1:
        raise ComparatorInputError("at least one observation is required")
    if not np.all(np.isfinite(ga)) or not np.all(np.isfinite(ya)):
        raise ComparatorInputError("g and y must be finite")
    if not np.all(np.isclose(ga, np.round(ga), rtol=0.0, atol=1e-12)):
        raise ComparatorInputError("g must contain genotype counts 0, 1 or 2")
    ga = np.round(ga)
    if np.any((ga < 0) | (ga > 2)):
        raise ComparatorInputError("g must contain genotype counts 0, 1 or 2")
    try:
        maf = float(case["maf"])
        family = str(case["family"]).lower()
    except (KeyError, TypeError, ValueError) as exc:
        raise ComparatorInputError("case must contain numeric maf and family") from exc
    if not math.isfinite(maf) or not 0.0 < maf < 1.0:
        raise ComparatorInputError("case.maf must lie strictly between zero and one")
    if family not in {"normal", "mixture"}:
        raise ComparatorInputError("case.family must be 'normal' or 'mixture'")
    x = ga - 2.0 * maf
    return ga.astype(float), ya.astype(float), x.astype(float), maf, family


def _normal_nll_grad(theta: np.ndarray, x: np.ndarray, y: np.ndarray) -> tuple[float, np.ndarray]:
    """Normal conditional negative log likelihood and analytic gradient."""
    alpha, a, logvar0, c = np.asarray(theta, dtype=float)
    eta = logvar0 + c * x
    if np.any(~np.isfinite(eta)) or np.any(eta < -745.0) or np.any(eta > 700.0):
        return float("inf"), np.full(4, np.nan)
    mu = alpha + a * x
    residual = y - mu
    invvar = np.exp(-eta)
    nll = 0.5 * float(np.sum(np.log(2.0 * np.pi) + eta + residual * residual * invvar))
    z = np.column_stack((np.ones_like(x), x))
    grad_mean = np.sum((-residual * invvar)[:, None] * z, axis=0)
    grad_var = np.sum((0.5 * (1.0 - residual * residual * invvar))[:, None] * z, axis=0)
    grad = np.r_[grad_mean, grad_var]
    return nll, grad


def _normal_hessian(theta: np.ndarray, x: np.ndarray, y: np.ndarray) -> np.ndarray:
    """Observed Hessian of the normal negative log likelihood."""
    alpha, a, logvar0, c = np.asarray(theta, dtype=float)
    eta = logvar0 + c * x
    mu = alpha + a * x
    residual = y - mu
    invvar = np.exp(-eta)
    z = np.column_stack((np.ones_like(x), x))
    h_mm = z.T @ (invvar[:, None] * z)
    h_mv = z.T @ ((residual * invvar)[:, None] * z)
    h_vv = z.T @ ((0.5 * residual * residual * invvar)[:, None] * z)
    hessian = np.block([[h_mm, h_mv], [h_mv.T, h_vv]])
    return 0.5 * (hessian + hessian.T)


def _normal_target_jacobian(theta: np.ndarray) -> np.ndarray:
    """Jacobian for (a, c, d var(x)/dx at x=0) from likelihood parameters."""
    _, _, logvar0, c = np.asarray(theta, dtype=float)
    baseline_variance = math.exp(float(logvar0))
    # The third row is intentionally not treated as an independent fitted
    # coefficient.  It carries covariance with both logvar0 and c.
    return np.array(
        [
            [0.0, 1.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 1.0],
            [0.0, 0.0, baseline_variance * c, baseline_variance],
        ],
        dtype=float,
    )


def _valid_covariance(hessian: np.ndarray) -> tuple[np.ndarray | None, str | None]:
    """Invert an observed Hessian, rejecting non-identifiable fits."""
    hessian = np.asarray(hessian, dtype=float)
    if hessian.ndim != 2 or hessian.shape[0] != hessian.shape[1] or not np.all(np.isfinite(hessian)):
        return None, "observed Hessian is non-finite or has the wrong shape"
    hessian = 0.5 * (hessian + hessian.T)
    try:
        eigenvalues = np.linalg.eigvalsh(hessian)
    except np.linalg.LinAlgError as exc:
        return None, f"observed Hessian eigendecomposition failed: {exc}"
    if not np.all(np.isfinite(eigenvalues)):
        return None, "observed Hessian eigenvalues are non-finite"
    largest = max(1.0, float(np.max(np.abs(eigenvalues))))
    if float(np.min(eigenvalues)) <= 1e-10 * largest:
        return None, "observed Hessian is not positive definite"
    condition = float(np.max(eigenvalues) / np.min(eigenvalues))
    if not math.isfinite(condition) or condition > 1e14:
        return None, "observed Hessian is numerically ill-conditioned"
    try:
        covariance = np.linalg.inv(hessian)
    except np.linalg.LinAlgError as exc:
        return None, f"observed Hessian inversion failed: {exc}"
    covariance = 0.5 * (covariance + covariance.T)
    if not np.all(np.isfinite(covariance)):
        return None, "inverse observed Hessian is non-finite"
    return covariance, None


def _wald_values(estimates: Mapping[str, float], covariance: np.ndarray, names: Sequence[str]) -> tuple[dict[str, float], dict[str, float]]:
    ses: dict[str, float] = {}
    p_values: dict[str, float] = {}
    diagonal = np.diag(np.asarray(covariance, dtype=float))
    for i, name in enumerate(names):
        variance = float(diagonal[i])
        se = math.sqrt(variance) if math.isfinite(variance) and variance >= 0.0 else float("nan")
        ses[name] = se
        estimate = float(estimates[name])
        p_values[name] = 2.0 * float(norm.sf(abs(estimate / se))) if se > 0.0 and math.isfinite(se) else float("nan")
    return ses, p_values


def _near(value: float, bound: float, *, rtol: float = 2e-6) -> bool:
    return abs(float(value) - float(bound)) <= rtol * max(1.0, abs(float(bound)))


def _normal_starts(x: np.ndarray, y: np.ndarray) -> list[np.ndarray]:
    """Deterministic, data-only starts for the normal likelihood."""
    ymean = float(np.mean(y))
    yvar = max(float(np.var(y)), np.finfo(float).eps)
    xss = float(np.dot(x, x))
    if xss > 0.0:
        a_ols = float(np.dot(x, y - ymean) / xss)
    else:
        a_ols = 0.0
    alpha_ols = float(np.mean(y - a_ols * x))
    resid = y - alpha_ols - a_ols * x
    log_resid2 = np.log(np.maximum(resid * resid, max(yvar * 1e-8, 1e-12)))
    if xss > 0.0:
        c_ols = float(np.dot(x, log_resid2 - np.mean(log_resid2)) / xss)
    else:
        c_ols = 0.0
    logvar = math.log(yvar)
    starts = [
        np.array([ymean, 0.0, logvar, 0.0]),
        np.array([alpha_ols, a_ols, math.log(max(float(np.var(resid)), 1e-12)), c_ols]),
        np.array([ymean, a_ols, logvar, 0.0]),
        np.array([alpha_ols, 0.0, logvar, -0.25]),
    ]
    return [np.asarray(start, dtype=float) for start in starts]


def fit_normal_ml(g: Sequence[Any], y: Sequence[Any], case: Mapping[str, Any]) -> dict[str, Any]:
    """Fit the conditional Normal location/log-variance ML comparator.

    Returned ``covariance`` is for the three reported scientific targets and
    therefore includes the cross-covariance needed by the transformed
    variance derivative.
    """
    started = time.monotonic()
    try:
        _, ya, x, _, family = _validate_data(g, y, case)
        if family != "normal":
            raise ComparatorInputError("normal ML requires case.family == 'normal'")
    except Exception as exc:
        return _fit_result(status="invalid_input", reason=str(exc), diagnostics={"wall_seconds": time.monotonic() - started})

    bounds = [(None, None), (None, None), (-30.0, 30.0), (-30.0, 30.0)]
    best = None
    failures: list[str] = []
    for start in _normal_starts(x, ya):
        try:
            result = minimize(
                lambda theta: _normal_nll_grad(theta, x, ya),
                start,
                jac=True,
                method="L-BFGS-B",
                bounds=bounds,
                options={"maxiter": 2000, "ftol": 1e-12, "gtol": 1e-8, "maxls": 50},
            )
        except Exception as exc:  # keep trying deterministic starts
            failures.append(f"start raised {type(exc).__name__}: {exc}")
            continue
        objective = float(result.fun) if np.isfinite(result.fun) else float("inf")
        if not np.isfinite(objective) or not np.all(np.isfinite(result.x)):
            failures.append(f"start returned non-finite objective ({result.message})")
            continue
        if bool(result.success):
            if best is None or objective < best[0]:
                best = (objective, result)
        else:
            failures.append(f"start did not converge: {result.message}")

    if best is None:
        reason = "; ".join(failures[:3]) or "all deterministic normal-ML starts failed"
        return _fit_result(
            status="fit_failure",
            reason=reason,
            diagnostics={"wall_seconds": time.monotonic() - started, "n_starts": len(_normal_starts(x, ya))},
        )

    theta = np.asarray(best[1].x, dtype=float)
    estimates = {
        "a": float(theta[1]),
        "c": float(theta[3]),
        "variance_derivative": float(math.exp(theta[2]) * theta[3]),
    }
    nuisance = {"mean_intercept": float(theta[0]), "log_variance_intercept": float(theta[2])}
    boundary = any(
        _near(theta[i], bound)
        for i, bound in ((2, -30.0), (2, 30.0), (3, -30.0), (3, 30.0))
    )
    if boundary:
        return _fit_result(
            status="boundary_failure",
            reason="normal ML optimum reached a log-variance/log-slope bound",
            estimates=estimates,
            nuisance=nuisance,
            diagnostics={"wall_seconds": time.monotonic() - started, "theta": theta.tolist(), "objective": best[0]},
        )

    try:
        hessian = _normal_hessian(theta, x, ya)
        covariance_full, hessian_reason = _valid_covariance(hessian)
        if covariance_full is None:
            return _fit_result(
                status="hessian_failure",
                reason=hessian_reason,
                estimates=estimates,
                nuisance=nuisance,
                diagnostics={"wall_seconds": time.monotonic() - started, "theta": theta.tolist(), "objective": best[0]},
            )
        target_jacobian = _normal_target_jacobian(theta)
        covariance = target_jacobian @ covariance_full @ target_jacobian.T
        covariance = 0.5 * (covariance + covariance.T)
        if not np.all(np.isfinite(covariance)):
            raise np.linalg.LinAlgError("transformed covariance is non-finite")
        ses, p_values = _wald_values(estimates, covariance, _NORMAL_PARAMS)
        if not all(math.isfinite(float(value)) for value in ses.values()) or not all(math.isfinite(float(value)) for value in p_values.values()):
            raise np.linalg.LinAlgError("transformed Wald covariance produced non-finite uncertainty")
    except Exception as exc:
        return _fit_result(
            status="hessian_failure",
            reason=f"normal transformed covariance failed: {exc}",
            estimates=estimates,
            nuisance=nuisance,
            diagnostics={"wall_seconds": time.monotonic() - started, "theta": theta.tolist(), "objective": best[0]},
        )
    return _fit_result(
        status="ok",
        reason=None,
        estimates=estimates,
        standard_errors=ses,
        p_values=p_values,
        covariance=covariance,
        nuisance=nuisance,
        diagnostics={
            "wall_seconds": time.monotonic() - started,
            "theta": theta.tolist(),
            "objective": best[0],
            "hessian": hessian.tolist(),
            "full_covariance": covariance_full.tolist(),
            "n_starts": len(_normal_starts(x, ya)),
        },
    )


def _mixture_unpack(u: np.ndarray, span: float, smooth: float = 1e-10) -> tuple[np.ndarray, np.ndarray]:
    """Map unconstrained mixture optimisation coordinates to physical ones.

    Physical order is ``(mu1_0, mu2_0, delta1, delta2, log_sd1, log_sd2,
    logit_pi0, gamma)``.  The positive baseline gap plus the slope-difference
    margin makes component 1 lower than component 2 at every genotype while
    retaining separate delta1 and delta2 parameters.
    """
    u = np.asarray(u, dtype=float)
    if u.size != 8:
        raise ValueError("mixture optimisation vector must have length 8")
    alpha1, delta1, log_gap, delta_difference, log_sd1, log_sd2, lp0, gamma = u
    exp_gap = math.exp(float(log_gap))
    root = math.sqrt(float(delta_difference * delta_difference + smooth * smooth))
    gap = exp_gap + float(span) * root
    q = np.array(
        [alpha1, alpha1 + gap, delta1, delta1 + delta_difference, log_sd1, log_sd2, lp0, gamma],
        dtype=float,
    )
    jacobian = np.zeros((8, 8), dtype=float)
    jacobian[0, 0] = 1.0
    jacobian[1, 0] = 1.0
    jacobian[1, 2] = exp_gap
    jacobian[1, 3] = float(span) * delta_difference / root
    jacobian[2, 1] = 1.0
    jacobian[3, 1] = 1.0
    jacobian[3, 3] = 1.0
    jacobian[4, 4] = 1.0
    jacobian[5, 5] = 1.0
    jacobian[6, 6] = 1.0
    jacobian[7, 7] = 1.0
    return q, jacobian


def _mixture_nll_grad(u: np.ndarray, x: np.ndarray, y: np.ndarray, span: float) -> tuple[float, np.ndarray]:
    q, jacobian = _mixture_unpack(u, span)
    mu1 = q[0] + q[2] * x
    mu2 = q[1] + q[3] * x
    lp = q[6] + q[7] * x
    logp1 = -np.logaddexp(0.0, -lp)
    logp2 = -np.logaddexp(0.0, lp)
    log_sd1, log_sd2 = q[4], q[5]
    var1 = np.exp(2.0 * log_sd1)
    var2 = np.exp(2.0 * log_sd2)
    r1 = y - mu1
    r2 = y - mu2
    logf1 = -0.5 * (np.log(2.0 * np.pi) + 2.0 * log_sd1 + r1 * r1 / var1)
    logf2 = -0.5 * (np.log(2.0 * np.pi) + 2.0 * log_sd2 + r2 * r2 / var2)
    component1 = logp1 + logf1
    component2 = logp2 + logf2
    logmix = np.logaddexp(component1, component2)
    if not np.all(np.isfinite(logmix)):
        return float("inf"), np.full(8, np.nan)
    responsibility1 = np.exp(component1 - logmix)
    responsibility2 = 1.0 - responsibility1
    nll = -float(np.sum(logmix))
    grad_q = np.zeros(8, dtype=float)
    grad_q[0] = float(np.sum(responsibility1 * (mu1 - y) / var1))
    grad_q[1] = float(np.sum(responsibility2 * (mu2 - y) / var2))
    grad_q[2] = float(np.sum(responsibility1 * (mu1 - y) * x / var1))
    grad_q[3] = float(np.sum(responsibility2 * (mu2 - y) * x / var2))
    grad_q[4] = float(np.sum(responsibility1 * (1.0 - r1 * r1 / var1)))
    grad_q[5] = float(np.sum(responsibility2 * (1.0 - r2 * r2 / var2)))
    pi1 = expit(lp)
    grad_q[6] = float(np.sum(pi1 - responsibility1))
    grad_q[7] = float(np.sum((pi1 - responsibility1) * x))
    return nll, jacobian.T @ grad_q


def _numeric_hessian(gradient: Callable[[np.ndarray], np.ndarray], point: np.ndarray) -> np.ndarray:
    """Central-difference Hessian for the smooth mixture coordinates."""
    point = np.asarray(point, dtype=float)
    n = point.size
    hessian = np.empty((n, n), dtype=float)
    for j in range(n):
        step = 2e-5 * max(1.0, abs(float(point[j])))
        plus = point.copy()
        minus = point.copy()
        plus[j] += step
        minus[j] -= step
        gp = np.asarray(gradient(plus), dtype=float)
        gm = np.asarray(gradient(minus), dtype=float)
        if not np.all(np.isfinite(gp)) or not np.all(np.isfinite(gm)):
            raise FloatingPointError("non-finite numerical Hessian difference")
        hessian[:, j] = (gp - gm) / (2.0 * step)
    return 0.5 * (hessian + hessian.T)


def _mixture_start_from_groups(x: np.ndarray, y: np.ndarray, lower: np.ndarray, upper: np.ndarray, span: float, yscale: float) -> np.ndarray:
    def line_fit(mask: np.ndarray) -> tuple[float, float, float]:
        if int(np.sum(mask)) >= 2:
            x_group = x[mask]
            y_group = y[mask]
            x_centered = x_group - np.mean(x_group)
            xss = float(np.dot(x_centered, x_centered))
        else:
            x_group = x[mask]
            y_group = y[mask]
            x_centered = x_group - np.mean(x_group) if np.any(mask) else x_group
            xss = float(np.dot(x_centered, x_centered)) if np.any(mask) else 0.0
        if xss > 1e-12:
            slope = float(np.dot(x_centered, y_group - np.mean(y_group)) / xss)
            intercept = float(np.mean(y_group) - slope * np.mean(x_group))
        else:
            slope = 0.0
            intercept = float(np.mean(y[mask])) if np.any(mask) else float(np.mean(y))
        residual = y[mask] - intercept - slope * x[mask] if np.any(mask) else y - np.mean(y)
        sd = max(float(np.std(residual)), 0.1 * yscale, 1e-6)
        return intercept, slope, sd

    mu1, d1, sd1 = line_fit(lower)
    mu2, d2, sd2 = line_fit(upper)
    dd = d2 - d1
    minimum_gap = max(1e-6 * yscale, 1e-8)
    gap = mu2 - mu1
    if gap <= span * abs(dd) + minimum_gap:
        gap = span * abs(dd) + 2.0 * minimum_gap
    log_gap = math.log(max(gap - span * abs(dd), minimum_gap))
    p1 = float(np.mean(lower)) if lower.size else 0.5
    p1 = min(max(p1, 0.02), 0.98)
    return np.array([mu1, d1, log_gap, dd, math.log(sd1), math.log(sd2), logit(p1), 0.0], dtype=float)


def _mixture_starts(x: np.ndarray, y: np.ndarray, span: float) -> list[np.ndarray]:
    """Deterministic starts based on y and g only, never case effects or labels."""
    ymean = float(np.mean(y))
    yscale = max(float(np.std(y)), 1e-3)
    q25, q50, q75 = [float(v) for v in np.quantile(y, [0.25, 0.50, 0.75])]
    starts: list[np.ndarray] = []
    lower = y <= q50
    upper = ~lower
    if np.any(lower) and np.any(upper):
        starts.append(_mixture_start_from_groups(x, y, lower, upper, span, yscale))
    for spread in (0.6, 1.0, 1.5):
        gap = max(2.0 * spread * yscale, 1e-4)
        starts.append(
            np.array(
                [
                    ymean - spread * yscale,
                    0.0,
                    math.log(gap),
                    0.0,
                    math.log(yscale),
                    math.log(yscale),
                    0.0,
                    0.0,
                ],
                dtype=float,
            )
        )
    # A quantile-separated start and two data-only membership perturbations help
    # avoid selecting a single local optimum without supplying an active axis.
    gap = max(q75 - q25, 1e-4)
    starts.append(np.array([q25, 0.0, math.log(gap), 0.0, math.log(yscale), math.log(yscale), 0.0, 0.35], dtype=float))
    starts.append(np.array([q25, 0.0, math.log(gap), 0.0, math.log(yscale), math.log(yscale), 0.0, -0.35], dtype=float))
    return starts


def fit_mixture_ml(g: Sequence[Any], y: Sequence[Any], case: Mapping[str, Any]) -> dict[str, Any]:
    """Fit a conditional two-Normal mixture with three free genotype slopes.

    The fit estimates the membership intercept and slope, both component
    intercepts and slopes, and both component log-SDs.  Component 1 is the
    lower-mean component at all three genotype values; no memberships are
    accepted from the caller.
    """
    started = time.monotonic()
    try:
        _, ya, x, _, family = _validate_data(g, y, case)
        if family != "mixture":
            raise ComparatorInputError("mixture ML requires case.family == 'mixture'")
    except Exception as exc:
        return _fit_result(status="invalid_input", reason=str(exc), diagnostics={"wall_seconds": time.monotonic() - started})

    genotype_x = np.array([0.0 - 2.0 * float(case["maf"]), 1.0 - 2.0 * float(case["maf"]), 2.0 - 2.0 * float(case["maf"])])
    span = max(1e-12, float(np.max(np.abs(genotype_x))))
    yscale = max(float(np.std(ya)), 1e-3)
    lp_bound = float(logit(0.98))
    logsd_lower = math.log(0.1)
    logsd_upper = math.log(10.0)
    # Match the fixed bounded nuisance parameterisation: a collapsed two-component
    # fit is a reported boundary failure, not an apparently successful ML fit.
    gap_floor = max(0.05, 1e-6 * yscale)
    loggap_lower = math.log(gap_floor)
    loggap_upper = math.log(max(100.0 * yscale, 1.0))
    bounds = [
        (None, None),
        (None, None),
        (loggap_lower, loggap_upper),
        (None, None),
        (logsd_lower, logsd_upper),
        (logsd_lower, logsd_upper),
        (-lp_bound, lp_bound),
        (-20.0, 20.0),
    ]
    starts = _mixture_starts(x, ya, span)
    best = None
    failures: list[str] = []
    for start in starts:
        # Keep starts inside the explicit nuisance bounds without changing any
        # data-derived scientific direction.
        start = np.asarray(start, dtype=float).copy()
        start[2] = min(max(start[2], loggap_lower + 1e-8), loggap_upper - 1e-8)
        start[4:6] = np.clip(start[4:6], logsd_lower + 1e-8, logsd_upper - 1e-8)
        start[6] = np.clip(start[6], -lp_bound + 1e-8, lp_bound - 1e-8)
        start[7] = np.clip(start[7], -19.0, 19.0)
        try:
            result = minimize(
                lambda u: _mixture_nll_grad(u, x, ya, span),
                start,
                jac=True,
                method="L-BFGS-B",
                bounds=bounds,
                options={"maxiter": 3000, "ftol": 1e-12, "gtol": 1e-7, "maxls": 80},
            )
        except Exception as exc:
            failures.append(f"start raised {type(exc).__name__}: {exc}")
            continue
        objective = float(result.fun) if np.isfinite(result.fun) else float("inf")
        if not np.isfinite(objective) or not np.all(np.isfinite(result.x)):
            failures.append(f"start returned non-finite objective ({result.message})")
            continue
        if bool(result.success):
            if best is None or objective < best[0]:
                best = (objective, result)
        else:
            failures.append(f"start did not converge: {result.message}")
    if best is None:
        reason = "; ".join(failures[:3]) or "all deterministic mixture-ML starts failed"
        return _fit_result(
            status="fit_failure",
            reason=reason,
            diagnostics={"wall_seconds": time.monotonic() - started, "n_starts": len(starts)},
        )

    u = np.asarray(best[1].x, dtype=float)
    q, jacobian = _mixture_unpack(u, span)
    estimates = {"gamma": float(q[7]), "delta1": float(q[2]), "delta2": float(q[3])}
    nuisance = {
        "mu1": float(q[0]),
        "mu2": float(q[1]),
        "log_sd1": float(q[4]),
        "log_sd2": float(q[5]),
        "pi0": float(expit(q[6])),
        "logit_pi0": float(q[6]),
    }
    ordered = [float(q[1] + q[3] * value - q[0] - q[2] * value) for value in genotype_x]
    boundary_coordinates = []
    for index, (lower, upper) in enumerate(bounds):
        if lower is not None and _near(u[index], lower):
            boundary_coordinates.append(index)
        if upper is not None and _near(u[index], upper):
            boundary_coordinates.append(index)
    if boundary_coordinates:
        return _fit_result(
            status="boundary_failure",
            reason=f"mixture ML optimum reached bound(s) at coordinates {sorted(set(boundary_coordinates))}",
            estimates=estimates,
            nuisance=nuisance,
            diagnostics={
                "wall_seconds": time.monotonic() - started,
                "u": u.tolist(),
                "objective": best[0],
                "ordered_gaps": ordered,
                "n_starts": len(starts),
            },
        )
    if not np.all(np.asarray(ordered) > 0.0):
        return _fit_result(
            status="hessian_failure",
            reason="ordered component constraint was not satisfied at all genotypes",
            estimates=estimates,
            nuisance=nuisance,
            diagnostics={"wall_seconds": time.monotonic() - started, "ordered_gaps": ordered},
        )

    try:
        hessian = _numeric_hessian(lambda point: _mixture_nll_grad(point, x, ya, span)[1], u)
        covariance_u, hessian_reason = _valid_covariance(hessian)
        if covariance_u is None:
            return _fit_result(
                status="hessian_failure",
                reason=f"mixture Hessian/nonidentifiability failure: {hessian_reason}",
                estimates=estimates,
                nuisance=nuisance,
                diagnostics={"wall_seconds": time.monotonic() - started, "u": u.tolist(), "objective": best[0]},
            )
        target_jacobian = np.array(
            [
                [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0],  # gamma
                [0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],  # delta1
                [0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0],  # delta2
            ],
            dtype=float,
        )
        covariance = target_jacobian @ covariance_u @ target_jacobian.T
        covariance = 0.5 * (covariance + covariance.T)
        if not np.all(np.isfinite(covariance)):
            raise np.linalg.LinAlgError("transformed mixture covariance is non-finite")
        ses, p_values = _wald_values(estimates, covariance, _MIXTURE_PARAMS)
        if not all(math.isfinite(float(value)) for value in ses.values()) or not all(math.isfinite(float(value)) for value in p_values.values()):
            raise np.linalg.LinAlgError("transformed mixture Wald covariance produced non-finite uncertainty")
    except Exception as exc:
        return _fit_result(
            status="hessian_failure",
            reason=f"mixture Hessian/transformed covariance failure: {exc}",
            estimates=estimates,
            nuisance=nuisance,
            diagnostics={"wall_seconds": time.monotonic() - started, "u": u.tolist(), "objective": best[0]},
        )
    return _fit_result(
        status="ok",
        reason=None,
        estimates=estimates,
        standard_errors=ses,
        p_values=p_values,
        covariance=covariance,
        nuisance=nuisance,
        diagnostics={
            "wall_seconds": time.monotonic() - started,
            "u": u.tolist(),
            "objective": best[0],
            "ordered_gaps": ordered,
            "hessian": hessian.tolist(),
            "n_starts": len(starts),
        },
    )


def fit_brown_forsythe(g: Sequence[Any], y: Sequence[Any], case: Mapping[str, Any] | None = None) -> dict[str, Any]:
    """Return the Brown--Forsythe genotype-group omnibus test.

    This is deliberately not represented as a fitted variance slope.  The
    statistic tests equality of group dispersions and has its own omnibus df.
    """
    started = time.monotonic()
    try:
        ga = np.asarray(g, dtype=float)
        ya = np.asarray(y, dtype=float)
        if ga.ndim != 1 or ya.ndim != 1 or ga.size != ya.size:
            raise ComparatorInputError("g and y must be one-dimensional arrays of equal length")
        if not np.all(np.isfinite(ga)) or not np.all(np.isfinite(ya)):
            raise ComparatorInputError("g and y must be finite")
        if not np.all(np.isclose(ga, np.round(ga), rtol=0.0, atol=1e-12)) or np.any((ga < 0.0) | (ga > 2.0)):
            raise ComparatorInputError("g must contain genotype counts 0, 1 or 2")
        ga = np.round(ga)
        groups = [ya[ga == genotype] for genotype in (0.0, 1.0, 2.0)]
        if any(group.size < 2 for group in groups):
            raise ComparatorInputError("Brown-Forsythe requires at least two observations in every genotype group")
        statistic, p_value = levene(*groups, center="median")
        if not math.isfinite(float(statistic)) or not math.isfinite(float(p_value)):
            raise FloatingPointError("Brown-Forsythe returned a non-finite statistic")
    except Exception as exc:
        return _fit_result(
            status="fit_failure",
            reason=f"Brown-Forsythe failed: {exc}",
            diagnostics={"wall_seconds": time.monotonic() - started},
        )
    return _fit_result(
        status="ok",
        reason=None,
        estimates={"group_omnibus": float(statistic)},
        standard_errors={"group_omnibus": float("nan")},
        p_values={"group_omnibus": float(p_value)},
        diagnostics={"wall_seconds": time.monotonic() - started, "test": "Brown-Forsythe", "groups": 3},
    )


def _rows_from_fit(method: str, names: Sequence[str], result: Mapping[str, Any], wall_seconds: float) -> list[dict[str, Any]]:
    status = str(result.get("status", "error"))
    reason = result.get("reason")
    estimates = result.get("estimates", {}) or {}
    ses = result.get("standard_errors", {}) or {}
    p_values = result.get("p_values", {}) or {}
    return [
        _row(
            method,
            name,
            estimate=estimates.get(name),
            se=ses.get(name),
            p_value=p_values.get(name),
            target_kind="scientific",
            status=status,
            reason=reason,
            wall_seconds=wall_seconds,
        )
        for name in names
    ]


def fit_all(g: Sequence[Any], y: Sequence[Any], case: Mapping[str, Any]) -> list[dict[str, Any]]:
    """Fit the family-specific SciPy comparators and Brown--Forsythe.

    Rows are emitted in a stable order and are never dropped when a fit,
    Hessian, or boundary check fails.  DGLM is intentionally not silently
    substituted here; call :func:`fit_dglm` explicitly for the R comparator.
    """
    try:
        _, _, _, _, family = _validate_data(g, y, case)
    except Exception as exc:
        family = str(case.get("family", "unknown")).lower() if isinstance(case, Mapping) else "unknown"
        if family == "normal":
            names = _NORMAL_PARAMS
            method = "normal_ml"
        elif family == "mixture":
            names = _MIXTURE_PARAMS
            method = "mixture_ml"
        else:
            names = ()
            method = "invalid"
        rows = _rows_from_fit(method, names, _fit_result(status="invalid_input", reason=str(exc)), 0.0)
        rows.append(_row("brown_forsythe", "group_omnibus", target_kind="group_omnibus", status="invalid_input", reason=str(exc), wall_seconds=0.0))
        return rows

    rows: list[dict[str, Any]] = []
    if family == "normal":
        started = time.monotonic()
        try:
            result = fit_normal_ml(g, y, case)
        except Exception as exc:
            result = _fit_result(status="fit_failure", reason=f"normal ML raised {type(exc).__name__}: {exc}")
        rows.extend(_rows_from_fit("normal_ml", _NORMAL_PARAMS, result, time.monotonic() - started))
    else:
        started = time.monotonic()
        try:
            result = fit_mixture_ml(g, y, case)
        except Exception as exc:
            result = _fit_result(status="fit_failure", reason=f"mixture ML raised {type(exc).__name__}: {exc}")
        rows.extend(_rows_from_fit("mixture_ml", _MIXTURE_PARAMS, result, time.monotonic() - started))

    started = time.monotonic()
    try:
        bf = fit_brown_forsythe(g, y, case)
    except Exception as exc:
        bf = _fit_result(status="fit_failure", reason=f"Brown-Forsythe raised {type(exc).__name__}: {exc}")
    bf_status = str(bf.get("status", "error"))
    rows.append(
        _row(
            "brown_forsythe",
            "group_omnibus",
            estimate=(bf.get("estimates", {}) or {}).get("group_omnibus"),
            se=None,
            p_value=(bf.get("p_values", {}) or {}).get("group_omnibus"),
            target_kind="group_omnibus",
            status=bf_status,
            reason=bf.get("reason"),
            wall_seconds=time.monotonic() - started,
        )
    )
    return rows


def _dglm_failure_rows(status: str, reason: str, wall_seconds: float) -> list[dict[str, Any]]:
    return [
        _row("dglm", name, target_kind="scientific", status=status, reason=reason, wall_seconds=wall_seconds)
        for name in _NORMAL_PARAMS
    ]


def fit_dglm(
    g: Sequence[Any],
    y: Sequence[Any],
    case: Mapping[str, Any],
    *,
    timeout_seconds: float = 120.0,
    rscript: str = "Rscript",
    adapter: str | os.PathLike[str] | None = None,
) -> list[dict[str, Any]]:
    """Run the explicit R DGLM adapter with a bounded subprocess.

    Timeout is reported as ``timeout``.  Missing R, a non-zero R exit, malformed
    adapter output, and input errors are reported separately as ``error`` or
    ``invalid_input``; none is relabelled as a timeout and no SciPy fallback is
    attempted.
    """
    started = time.monotonic()
    try:
        _, ya, x, _, family = _validate_data(g, y, case)
        if family != "normal":
            raise ComparatorInputError("DGLM adapter requires case.family == 'normal'")
        timeout_seconds = float(timeout_seconds)
        if not math.isfinite(timeout_seconds) or timeout_seconds <= 0.0:
            raise ComparatorInputError("timeout_seconds must be positive and finite")
    except Exception as exc:
        return _dglm_failure_rows("invalid_input", str(exc), time.monotonic() - started)

    if adapter is None:
        adapter_path = Path(__file__).with_name("dglm.R")
    else:
        adapter_path = Path(adapter)
    if not adapter_path.exists():
        return _dglm_failure_rows("error", f"DGLM adapter does not exist: {adapter_path}", time.monotonic() - started)

    try:
        with tempfile.TemporaryDirectory(prefix="gmoq-dglm-") as directory:
            input_path = Path(directory) / "input.tsv"
            output_path = Path(directory) / "output.tsv"
            with input_path.open("w", encoding="utf-8", newline="") as stream:
                writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
                writer.writerow(("g", "y", "x"))
                for gi, yi, xi in zip(np.asarray(g), ya, x):
                    writer.writerow((int(round(float(gi))), repr(float(yi)), repr(float(xi))))
            try:
                completed = subprocess.run(
                    [str(rscript), str(adapter_path), str(input_path), str(output_path)],
                    check=False,
                    capture_output=True,
                    text=True,
                    timeout=timeout_seconds,
                )
            except subprocess.TimeoutExpired as exc:
                elapsed = time.monotonic() - started
                return _dglm_failure_rows("timeout", f"DGLM subprocess exceeded {timeout_seconds:g}s", elapsed)
            except OSError as exc:
                elapsed = time.monotonic() - started
                return _dglm_failure_rows("error", f"DGLM subprocess could not start: {exc}", elapsed)
            if completed.returncode != 0:
                detail = (completed.stderr or completed.stdout or "no adapter diagnostics").strip().replace("\n", " ")
                return _dglm_failure_rows("error", f"DGLM adapter exited {completed.returncode}: {detail[:500]}", time.monotonic() - started)
            if not output_path.exists():
                return _dglm_failure_rows("error", "DGLM adapter exited successfully without output", time.monotonic() - started)
            with output_path.open("r", encoding="utf-8", newline="") as stream:
                records = list(csv.DictReader(stream, delimiter="\t"))
            if len(records) != len(_NORMAL_PARAMS):
                return _dglm_failure_rows("error", f"DGLM adapter returned {len(records)} rows, expected {len(_NORMAL_PARAMS)}", time.monotonic() - started)
            by_param = {record.get("param"): record for record in records}
            if set(by_param) != set(_NORMAL_PARAMS):
                return _dglm_failure_rows("error", "DGLM adapter returned an unexpected parameter schema", time.monotonic() - started)
            elapsed = time.monotonic() - started
            rows = []
            for param in _NORMAL_PARAMS:
                record = by_param[param]
                status = record.get("status") or "error"
                rows.append(
                    _row(
                        "dglm",
                        param,
                        estimate=record.get("estimate"),
                        se=record.get("se"),
                        p_value=record.get("p_value"),
                        target_kind=record.get("target_kind") or "scientific",
                        status=status,
                        reason=record.get("reason") or None,
                        wall_seconds=elapsed,
                    )
                )
            return rows
    except Exception as exc:
        return _dglm_failure_rows("error", f"DGLM adapter output handling failed: {exc}", time.monotonic() - started)


__all__ = [
    "ComparatorInputError",
    "DGLMError",
    "fit_all",
    "fit_brown_forsythe",
    "fit_dglm",
    "fit_mixture_ml",
    "fit_normal_ml",
]
