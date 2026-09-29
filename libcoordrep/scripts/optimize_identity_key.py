#!/usr/bin/env python3
"""
optimize_identity_key.py
========================

Grid-search over S1-bin boundaries and delta-threshold to find the
``CoordRep-ID-v1`` configuration that **minimises the false non-duplicate
rate at L1** while **preserving fac/mer and cis/trans discrimination**.

Search space
------------
- S1 bin boundary-1 (ideal ↔ good):  0.5, 1.0, 1.5, 2.0
- S1 bin boundary-2 (good ↔ dist):   2.0, 3.0, 4.0, 5.0
- Delta boundary-thresh (boundary flag): 0.3, 0.5, 1.0, 1.5, 2.0
- Delta bin boundary-1 (D0 ↔ D1):    0.3, 0.5, 1.0
- Delta bin boundary-2 (D1 ↔ D2):    1.5, 2.0, 3.0

Evaluation
----------
For each config point, encode N molecules, perturb them at σ=0.01 Å
(20 trials each), and measure:

  1.  L1 retention rate (higher is better)
  2.  Stereo preservation: do molecules with *different* constraint
      signatures still get *different* L1 keys?  (must stay high)
  3.  Combined score = L1_retention − λ × stereo_loss

Outputs (→ output_dir/):
  grid_search_results.csv   – all grid points with scores
  recommended_config.json   – the winning CoordRep-ID-v1 configuration
"""

from __future__ import annotations

import argparse
import itertools
import json
import sys
import warnings
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")

sys.path.insert(0, str(Path(__file__).parent.parent))

from coordrep import encode_molecule, CoordRepConfig
from coordrep.io.tmqm_reader import TMQMReader, Atom, RawMolecule
from coordrep.identity.identity_keys import extract_identity_keys
from coordrep.identity.shape_binning import DEFAULT_SHAPE_BIN_CONFIG


# ──────────────────────────────────────────────────────────────────
# Noise helper
# ──────────────────────────────────────────────────────────────────

def _add_noise(mol: RawMolecule, sigma: float, seed: int) -> RawMolecule:
    rng = np.random.RandomState(seed)
    new_atoms = [
        Atom(a.index, a.element, a.x + dx, a.y + dy, a.z + dz)
        for a, (dx, dy, dz) in zip(mol.atoms, rng.normal(0, sigma, (len(mol.atoms), 3)))
    ]
    return RawMolecule(mol.mol_id, new_atoms, mol.bond_orders, mol.properties)


# ──────────────────────────────────────────────────────────────────
# Build config dict from grid parameters
# ──────────────────────────────────────────────────────────────────

def _make_config(s1_b1, s1_b2, d_thresh, d_b1, d_b2) -> dict:
    return {
        "s1_bins": [
            (s1_b1, "ideal"),
            (s1_b2, "good"),
            (8.0, "dist"),
            (np.inf, "irreg"),
        ],
        "delta_bins": [
            (d_b1, "D0"),
            (d_b2, "D1"),
            (5.0, "D2"),
            (np.inf, "D3"),
        ],
        "boundary_thresh": d_thresh,
    }


# ──────────────────────────────────────────────────────────────────
# Evaluation at one grid point
# ──────────────────────────────────────────────────────────────────

def precompute_noisy_strings(
    refs: List[dict],
    mol_map: Dict[str, RawMolecule],
    sigma: float,
    n_trials: int,
    enc_config=None,
) -> List[Tuple[str, str]]:
    """
    Pre-compute all (ref_string, noisy_string) pairs once.
    Returns list of (ref_coordrep, noisy_coordrep) tuples.
    """
    enc_config = enc_config or CoordRepConfig.default()
    pairs = []
    for ref in refs:
        mol = mol_map[ref["mol_id"]]
        for t in range(n_trials):
            noisy = _add_noise(mol, sigma, seed=t)
            try:
                cc = encode_molecule(noisy, enc_config)
                s = cc.canonicalize().to_string()
            except Exception:
                continue
            if not s:
                continue
            pairs.append((ref["string"], s))
    return pairs


def evaluate_config(
    precomputed_pairs: List[Tuple[str, str]],
    bin_cfg: dict,
) -> dict:
    """
    Return L1 retention and stereo-preservation for one bin_cfg.
    Uses pre-computed (ref_string, noisy_string) pairs to avoid
    re-running encode_molecule for every grid point.
    """
    l1_matches = 0
    l1_total = 0
    stereo_distinct_ok = 0
    stereo_distinct_total = 0

    for ref_str, noisy_str in precomputed_pairs:
        ref_keys = extract_identity_keys(ref_str, bin_cfg)
        nk = extract_identity_keys(noisy_str, bin_cfg)
        ref_csig = ref_keys.constraint_signature

        l1_total += 1
        if nk.L1_ShapeID == ref_keys.L1_ShapeID:
            l1_matches += 1

        if ref_csig != "none":
            stereo_distinct_total += 1
            if nk.constraint_signature == ref_csig:
                stereo_distinct_ok += 1

    l1_ret = l1_matches / l1_total if l1_total else 0
    stereo_pres = stereo_distinct_ok / stereo_distinct_total if stereo_distinct_total else 1.0

    return {
        "l1_retention": l1_ret,
        "stereo_preservation": stereo_pres,
        "l1_total": l1_total,
        "stereo_total": stereo_distinct_total,
    }


# ──────────────────────────────────────────────────────────────────
# Main
# ──────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--tmqm-dir", default="/data/CoordRep/CoordSMILES/tmQM-master/tmQM")
    parser.add_argument("--n-molecules", type=int, default=100)
    parser.add_argument("--n-trials", type=int, default=20)
    parser.add_argument("--sigma", type=float, default=0.01)
    parser.add_argument("--lambda-stereo", type=float, default=2.0,
                        help="Penalty weight for stereo-discrimination loss")
    parser.add_argument("--output-dir",
                        default="/data/CoordRep/coordrep-release/libcoordrep/outputs/identity_optimize")
    args = parser.parse_args()

    out = Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)
    enc_config = CoordRepConfig.default()

    # ── Load & encode references ──────────────────────────────────
    print("Loading tmQM …")
    reader = TMQMReader(args.tmqm_dir)
    reader.load()
    all_ids = reader.get_mol_ids()

    rng = np.random.RandomState(42)
    sample_ids = rng.choice(all_ids, min(args.n_molecules, len(all_ids)), replace=False)

    mol_map: Dict[str, RawMolecule] = {}
    refs: List[dict] = []
    for mid in sample_ids:
        mol = reader.get_molecule(mid)
        if mol is None:
            continue
        try:
            cc = encode_molecule(mol, enc_config)
            s = cc.canonicalize().to_string()
        except Exception:
            continue
        if not s:
            continue
        keys = extract_identity_keys(s)
        refs.append({"mol_id": mid, "string": s, "keys": keys})
        mol_map[mid] = mol

    print(f"  Encoded {len(refs)} molecules")

    # ── Pre-compute noisy strings (encode_molecule is the bottleneck) ──
    print(f"\nPre-computing noisy encodings: {len(refs)} mols × {args.n_trials} trials at σ={args.sigma} …")
    pairs = precompute_noisy_strings(refs, mol_map, args.sigma, args.n_trials, enc_config)
    print(f"  {len(pairs)} valid (ref, noisy) pairs")

    # ── Grid search ───────────────────────────────────────────────
    s1_b1_vals = [0.5, 1.0, 1.5, 2.0]
    s1_b2_vals = [2.0, 3.0, 4.0, 5.0]
    d_thresh_vals = [0.3, 0.5, 1.0, 1.5, 2.0]
    d_b1_vals = [0.3, 0.5, 1.0]
    d_b2_vals = [1.5, 2.0, 3.0]

    grid = list(itertools.product(s1_b1_vals, s1_b2_vals, d_thresh_vals, d_b1_vals, d_b2_vals))
    # Filter invalid combos: s1_b1 < s1_b2, d_b1 < d_b2
    grid = [(a, b, c, d, e) for a, b, c, d, e in grid if a < b and d < e]

    print(f"\nGrid search: {len(grid)} configurations (string-level only, fast) …")

    results = []
    for i, (s1_b1, s1_b2, d_thresh, d_b1, d_b2) in enumerate(grid):
        cfg = _make_config(s1_b1, s1_b2, d_thresh, d_b1, d_b2)
        ev = evaluate_config(pairs, cfg)

        stereo_loss = 1.0 - ev["stereo_preservation"]
        score = ev["l1_retention"] - args.lambda_stereo * stereo_loss

        row = {
            "s1_b1": s1_b1, "s1_b2": s1_b2,
            "d_thresh": d_thresh, "d_b1": d_b1, "d_b2": d_b2,
            "l1_retention": round(ev["l1_retention"], 4),
            "stereo_preservation": round(ev["stereo_preservation"], 4),
            "stereo_loss": round(stereo_loss, 4),
            "score": round(score, 4),
            "false_nondup_l1": round(1 - ev["l1_retention"], 4),
        }
        results.append(row)

        if (i + 1) % 50 == 0:
            print(f"  {i+1}/{len(grid)} configs evaluated …")

    res_df = pd.DataFrame(results).sort_values("score", ascending=False)
    res_df.to_csv(out / "grid_search_results.csv", index=False)

    # ── Best config ───────────────────────────────────────────────
    best = res_df.iloc[0]
    print(f"\n{'='*60}")
    print("BEST CONFIGURATION  (CoordRep-ID-v1)")
    print(f"{'='*60}")
    print(f"  S1 bin boundaries:    {best['s1_b1']} / {best['s1_b2']} / 8.0")
    print(f"  Delta bin boundaries: {best['d_b1']} / {best['d_b2']} / 5.0")
    print(f"  Boundary threshold:   {best['d_thresh']}")
    print(f"  L1 retention:         {best['l1_retention']:.1%}")
    print(f"  Stereo preservation:  {best['stereo_preservation']:.1%}")
    print(f"  False non-dup (L1):   {best['false_nondup_l1']:.1%}")
    print(f"  Score:                {best['score']:.4f}")

    recommended = {
        "version": "CoordRep-ID-v1",
        "s1_bins": [
            [float(best["s1_b1"]), "ideal"],
            [float(best["s1_b2"]), "good"],
            [8.0, "dist"],
            [float("inf"), "irreg"],
        ],
        "delta_bins": [
            [float(best["d_b1"]), "D0"],
            [float(best["d_b2"]), "D1"],
            [5.0, "D2"],
            [float("inf"), "D3"],
        ],
        "boundary_thresh": float(best["d_thresh"]),
        "evaluation": {
            "sigma": args.sigma,
            "n_trials": args.n_trials,
            "n_molecules": len(refs),
            "l1_retention": float(best["l1_retention"]),
            "stereo_preservation": float(best["stereo_preservation"]),
            "false_nondup_l1": float(best["false_nondup_l1"]),
            "score": float(best["score"]),
        },
    }
    with open(out / "recommended_config.json", "w") as f:
        json.dump(recommended, f, indent=2)

    # ── Top-5 configs ─────────────────────────────────────────────
    print(f"\nTop-5 configurations:")
    print(res_df[["s1_b1", "s1_b2", "d_thresh", "d_b1", "d_b2",
                   "l1_retention", "stereo_preservation", "false_nondup_l1", "score"]].head().to_string(index=False))

    print(f"\nAll outputs saved to {out}")


if __name__ == "__main__":
    main()
