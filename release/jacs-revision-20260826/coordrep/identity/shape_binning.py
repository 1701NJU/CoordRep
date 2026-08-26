"""
shape_binning.py
================

Replace exact CShM decimal values with a discrete, perturbation-tolerant
shape token for use in L1 (ShapeID).

Token grammar::

    <best_shape>/<second_shape>.<s1_bin>.<delta_bin>

Special case – when the top-1 ↔ top-2 gap is below a configurable
boundary threshold the token becomes::

    <best_shape>/<second_shape>_boundary.<s1_bin>.<delta_bin>

Examples::

    Oh/TPr.ideal.D3          # clearly octahedral
    TBP/SPY.good.D1          # decent TBP, moderate gap
    TBP/SPY_boundary.dist.D0 # ambiguous TBP vs SPY, small gap
    Td/SP.ideal.D3           # clear tetrahedral

Configuration keys (DEFAULT_SHAPE_BIN_CONFIG)::

    s1_bins         : list of (upper_bound, label) for best-CShM bucketing
    delta_bins      : list of (upper_bound, label) for Δ(S2−S1) bucketing
    boundary_thresh : Δ below which top-1/top-2 assignment is marked ambiguous
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

import numpy as np


# ──────────────────────────────────────────────────────────────────
# Default binning configuration
# ──────────────────────────────────────────────────────────────────

DEFAULT_SHAPE_BIN_CONFIG: Dict = {
    "s1_bins": [
        (1.0, "ideal"),
        (3.0, "good"),
        (8.0, "dist"),
        (np.inf, "irreg"),
    ],
    "delta_bins": [
        (0.5, "D0"),
        (2.0, "D1"),
        (5.0, "D2"),
        (np.inf, "D3"),
    ],
    "boundary_thresh": 1.0,
}


# ──────────────────────────────────────────────────────────────────
# ShapeBin helpers
# ──────────────────────────────────────────────────────────────────

@dataclass
class ShapeBin:
    """A single bin specification: values < upper_bound get this label."""
    upper_bound: float
    label: str


def _classify(value: float, bins: List[Tuple[float, str]]) -> str:
    for upper, label in bins:
        if value < upper:
            return label
    return bins[-1][1] if bins else "?"


# ──────────────────────────────────────────────────────────────────
# BinnedShape result
# ──────────────────────────────────────────────────────────────────

@dataclass
class BinnedShape:
    """Result of shape binning for one complex."""
    best_shape: str           # e.g. "Oh"
    second_shape: str         # e.g. "TPr"  (or "" if CN has only one ref shape)
    s1_bin: str               # e.g. "ideal"
    delta_bin: str            # e.g. "D3"
    is_boundary: bool         # True if top1-top2 gap < boundary_thresh
    s1_value: float           # raw CShM of best shape
    delta_value: float        # S2 - S1

    @property
    def token(self) -> str:
        """The discrete shape identity token."""
        second = self.second_shape or "?"
        boundary_tag = "_boundary" if self.is_boundary else ""
        return f"{self.best_shape}/{second}{boundary_tag}.{self.s1_bin}.{self.delta_bin}"

    def as_dict(self) -> dict:
        return {
            "best_shape": self.best_shape,
            "second_shape": self.second_shape,
            "s1_bin": self.s1_bin,
            "delta_bin": self.delta_bin,
            "is_boundary": self.is_boundary,
            "s1_value": round(self.s1_value, 4),
            "delta_value": round(self.delta_value, 4),
            "token": self.token,
        }


# ──────────────────────────────────────────────────────────────────
# From ShapeResult (object-level API)
# ──────────────────────────────────────────────────────────────────

def bin_shape(
    ref_shapes: List[str],
    cshm_values: np.ndarray,
    config: Dict | None = None,
) -> BinnedShape:
    """
    Bin a CShM vector into a discrete shape token.

    Parameters
    ----------
    ref_shapes : list[str]
        Reference shape names, same order as *cshm_values*.
    cshm_values : ndarray
        Raw or rounded CShM values.
    config : dict, optional
        Binning thresholds; defaults to DEFAULT_SHAPE_BIN_CONFIG.
    """
    cfg = config or DEFAULT_SHAPE_BIN_CONFIG

    if len(cshm_values) == 0:
        return BinnedShape("?", "", "?", "?", False, 0.0, 0.0)

    order = np.argsort(cshm_values)
    best_idx = order[0]
    best_shape = ref_shapes[best_idx]
    s1 = float(cshm_values[best_idx])

    if len(cshm_values) > 1:
        second_idx = order[1]
        second_shape = ref_shapes[second_idx]
        delta = float(cshm_values[second_idx] - cshm_values[best_idx])
    else:
        second_shape = ""
        delta = float("inf")

    s1_bin = _classify(s1, cfg["s1_bins"])
    delta_bin = _classify(delta, cfg["delta_bins"])
    is_boundary = delta < cfg["boundary_thresh"]

    return BinnedShape(
        best_shape=best_shape,
        second_shape=second_shape,
        s1_bin=s1_bin,
        delta_bin=delta_bin,
        is_boundary=is_boundary,
        s1_value=s1,
        delta_value=delta,
    )


def bin_shape_from_result(shape_result, config: Dict | None = None) -> BinnedShape:
    """
    Convenience wrapper for ``coordrep.geometry.shape.ShapeResult``.
    """
    vals = shape_result.values_rounded
    if vals is None:
        vals = shape_result.values
    return bin_shape(shape_result.ref_shapes, vals, config)


# ──────────────────────────────────────────────────────────────────
# From CoordRep string (string-level API)
# ──────────────────────────────────────────────────────────────────

_RE_SHAPE_BLOCK = re.compile(
    r'<ShapeBest:(?P<best>\w[\w-]*)'
    r'\|Class:(?P<cls>\w+)'
    r'\|Delta:(?P<dbin>\d+)'
    r'\|V:(?P<vals>[^>]+)>'
)

# Shape name order per CN, matching IDEAL_GEOMETRIES in shape.py
_SHAPES_BY_CN = {
    2: ["L"],
    3: ["TP"],
    4: ["Td", "SP"],
    5: ["TBP", "SPY"],
    6: ["Oh", "TPr"],
}

_ALL_SHAPES = {s for shapes in _SHAPES_BY_CN.values() for s in shapes}


def _infer_ref_shapes_from_vals(best_shape: str, n_vals: int) -> List[str]:
    """
    Reconstruct the reference-shape name list from the best-shape label
    and number of CShM values.
    """
    for cn, shapes in _SHAPES_BY_CN.items():
        if best_shape in shapes and len(shapes) == n_vals:
            return shapes
    # Fallback: best_shape first, then unknowns
    return [best_shape] + [f"S{i}" for i in range(1, n_vals)]


def bin_shape_from_string(
    coordrep_str: str,
    config: Dict | None = None,
) -> BinnedShape:
    """
    Extract and bin shape information directly from a CoordRep string.
    """
    cfg = config or DEFAULT_SHAPE_BIN_CONFIG
    m = _RE_SHAPE_BLOCK.search(coordrep_str)
    if not m:
        return BinnedShape("?", "", "?", "?", False, 0.0, 0.0)

    best = m.group("best")
    vals_str = m.group("vals")
    try:
        vals = np.array([float(v) for v in vals_str.split(",")])
    except ValueError:
        return BinnedShape(best, "", "?", "?", False, 0.0, 0.0)

    ref_shapes = _infer_ref_shapes_from_vals(best, len(vals))
    return bin_shape(ref_shapes, vals, cfg)
