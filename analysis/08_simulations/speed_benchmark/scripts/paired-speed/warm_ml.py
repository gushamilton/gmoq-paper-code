"""Phenotype-only warm starts for the unchanged direct mixture comparator.

This sensitivity replaces only the per-variant start list. The likelihood,
all eight free parameters, bounds, optimisation, numerical Hessian, standard
errors and failure checks remain those of ``validation.comparators``. Any
unsuccessful warm fit is retried with its canonical deterministic start list.

Use one context in an isolated, single-threaded method process::

    with warm_start(y, {"family": "mixture", "maf": 0.30}) as fit:
        result = fit(g)
        # Optional second argument supplies a different variant's MAF.
        result = fit(g2, {"family": "mixture", "maf": 0.20})

The phenotype-only baseline is fitted once on context entry. Its elapsed time
is available as ``fit.baseline_wall_seconds`` and should be included once in
the method total. ``fit.counts`` reports attempts and fallback outcomes.
``diagnostics['warm_total_wall_seconds']`` includes a failed warm attempt and
its fallback, whereas the original ``wall_seconds`` describes the returned
canonical fit alone. Neither includes the shared baseline fit.

The start-generator patch exists only during each warm call. This adapter is
not intended for concurrent comparator calls in the same process.
"""

from __future__ import annotations

import inspect
import math
import time
from collections.abc import Iterator, Mapping, Sequence
from contextlib import contextmanager
from typing import Any
from unittest.mock import patch

import numpy as np

from validation import comparators, truth


class WarmStartFit:
    """Callable per-phenotype fit with shared baseline and fallback counts."""

    def __init__(self, y: Sequence[float], case: Mapping[str, Any]) -> None:
        self._y = np.array(y, dtype=float, copy=True)
        self._case = dict(case)
        if str(self._case.get("family", "")).lower() != "mixture":
            raise ValueError("warm mixture ML requires case.family == 'mixture'")
        self._active = False
        self.baseline: dict[str, float] | None = None
        self.baseline_error: str | None = None
        started = time.monotonic()
        try:
            # No genotype, case baseline, generating effect or label is used.
            self.baseline = truth.fit_baseline(self._y, "mixture")
        except Exception as exc:
            self.baseline_error = f"{type(exc).__name__}: {exc}"
        self.baseline_wall_seconds = time.monotonic() - started
        self.counts = {
            "variants": 0,
            "warm_attempts": 0,
            "warm_ok": 0,
            "warm_failed": 0,
            "fallbacks": 0,
            "fallback_ok": 0,
            "fallback_failed": 0,
        }

    def _starts(
        self, x: np.ndarray, y: np.ndarray, span: float
    ) -> list[np.ndarray]:
        """Encode the fitted baseline with all genotype slopes initially zero."""
        del x, y
        if self.baseline is None:
            raise RuntimeError("phenotype-only baseline is unavailable")
        b = self.baseline
        # At zero slope difference, the canonical physical mean gap is
        # exp(log_gap) + span * smooth. Subtract that smoothing margin so
        # _mixture_unpack recovers the fitted baseline rather than shifting it.
        smooth = float(
            inspect.signature(comparators._mixture_unpack)
            .parameters["smooth"].default
        )
        gap = float(b["mu2"] - b["mu1"]) - float(span) * abs(smooth)
        if not math.isfinite(gap) or gap <= 0.0:
            raise ValueError("baseline cannot define a positive warm-start gap")
        p1 = float(b["p1"])
        return [
            np.array(
                [
                    b["mu1"],
                    0.0,
                    math.log(gap),
                    0.0,
                    math.log(b["sd1"]),
                    math.log(b["sd2"]),
                    math.log(p1 / (1.0 - p1)),
                    0.0,
                ],
                dtype=float,
            )
        ]

    def __call__(
        self, g: Sequence[float], case: Mapping[str, Any] | None = None
    ) -> dict[str, Any]:
        if not self._active:
            raise RuntimeError("call the fitter inside its warm_start context")
        variant_case = self._case if case is None else dict(case)
        started = time.monotonic()
        self.counts["variants"] += 1
        warm_result: dict[str, Any] | None = None
        warm_wall_seconds = 0.0
        if self.baseline is not None:
            self.counts["warm_attempts"] += 1
            warm_started = time.monotonic()
            try:
                with patch.object(comparators, "_mixture_starts", self._starts):
                    warm_result = comparators.fit_mixture_ml(
                        g, self._y, variant_case
                    )
            except Exception as exc:
                # An adapter/start conversion failure also receives the
                # untouched comparator; its own returned checks decide success.
                warm_result = comparators._fit_result(
                    status="fit_failure",
                    reason=f"warm-start adapter raised {type(exc).__name__}: {exc}",
                )
            warm_wall_seconds = time.monotonic() - warm_started

        fallback = warm_result is None or warm_result["status"] != "ok"
        if fallback:
            if warm_result is not None:
                self.counts["warm_failed"] += 1
            self.counts["fallbacks"] += 1
            # The patch has ended: this uses the original deterministic starts.
            result = comparators.fit_mixture_ml(g, self._y, variant_case)
            outcome = "fallback_ok" if result["status"] == "ok" else "fallback_failed"
            self.counts[outcome] += 1
        else:
            self.counts["warm_ok"] += 1
            result = warm_result

        # Retain every original result and diagnostic field. Additional fields
        # distinguish the sensitivity and account for otherwise hidden retries.
        diagnostics = result["diagnostics"]
        diagnostics["warm_fallback"] = bool(fallback)
        diagnostics["warm_status"] = (
            "baseline_failure" if warm_result is None else warm_result["status"]
        )
        diagnostics["warm_reason"] = (
            self.baseline_error if warm_result is None else warm_result["reason"]
        )
        diagnostics["warm_attempt_wall_seconds"] = warm_wall_seconds
        diagnostics["warm_total_wall_seconds"] = time.monotonic() - started
        if fallback and warm_result is not None:
            diagnostics["warm_failure_diagnostics"] = warm_result["diagnostics"]
        return result


@contextmanager
def warm_start(
    y: Sequence[float], case: Mapping[str, Any]
) -> Iterator[WarmStartFit]:
    """Fit one observed baseline and yield a callable full per-variant refit."""
    fitter = WarmStartFit(y, case)
    fitter._active = True
    try:
        yield fitter
    finally:
        fitter._active = False
