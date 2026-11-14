#!/usr/bin/env python3
"""
benchmark_identity_layers.py
============================

Benchmark L0/L1/L2/L3 identity-key retention under:
  A.  Gaussian noise perturbations  (σ = 0.001 … 0.20 Å)
  B.  Systematic M-L bond scaling   (±1 … ±5 %)

For each perturbation, measure:
  - retention rate per identity level
  - false non-duplicate rate  (1 − retention)

Stratify results by:
  - coordination number (CN)
  - best-shape label
  - delta = S2 − S1 quartile
  - maximum ligand denticity

Outputs (→ output_dir/):
  noise_retention.csv        – per-σ retention
  scaling_retention.csv      – per-Δ% retention
  noise_stratified.csv       – per-σ × stratum
  scaling_stratified.csv     – per-Δ% × stratum
  false_nondup_summary.csv   – headline false-non-dup rates
  consensus_results.csv      – geometry-boundary flags
"""

from __future__ import annotations

import argparse
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
from coordrep.io.tmqm_reader import TMQMReader, RawMolecule, Atom, TRANSITION_METALS
from coordrep.identity import extract_identity_keys
from coordrep.identity.shape_binning import bin_shape_from_string, DEFAULT_SHAPE_BIN_CONFIG
from coordrep.identity.stability_consensus import compute_consensus_batch


# ──────────────────────────────────────────────────────────────────
# Noise & scaling helpers  (shared with test_geometric_robustness)
# ──────────────────────────────────────────────────────────────────

def _add_noise(mol: RawMolecule, sigma: float, seed: int) -> RawMolecule:
    rng = np.random.RandomState(seed)
    new_atoms = []
    for a in mol.atoms:
        dx, dy, dz = rng.normal(0, sigma, 3)
        new_atoms.append(Atom(index=a.index, element=a.element,
                              x=a.x + dx, y=a.y + dy, z=a.z + dz))
    return RawMolecule(mol_id=mol.mol_id, atoms=new_atoms,
                       bond_orders=mol.bond_orders, properties=mol.properties)


def _scale_bonds(mol: RawMolecule, sf: float) -> RawMolecule:
    from coordrep.io.tmqm_reader import find_metal_index
    midx = find_metal_index(mol.atoms)
    if midx is None:
        return mol
    mc = mol.atoms[midx].coords
    new_atoms = []
    for a in mol.atoms:
        if a.element in TRANSITION_METALS:
            new_atoms.append(a); continue
        vec = a.coords - mc
        d = np.linalg.norm(vec)
        if 0.5 < d < 3.0:
            disp = vec * (sf - 1.0)
            new_atoms.append(Atom(index=a.index, element=a.element,
                                  x=a.x + disp[0], y=a.y + disp[1], z=a.z + disp[2]))
        else:
            new_atoms.append(a)
    return RawMolecule(mol_id=mol.mol_id, atoms=new_atoms,
                       bond_orders=mol.bond_orders, properties=mol.properties)


# ──────────────────────────────────────────────────────────────────
# Encode reference data for each molecule
# ──────────────────────────────────────────────────────────────────

def _encode_ref(mol, config) -> dict | None:
    """Encode the undisturbed molecule and extract metadata + keys."""
    try:
        cc = encode_molecule(mol, config)
        s = cc.canonicalize().to_string()
    except Exception:
        return None
    if not s:
        return None

    keys = extract_identity_keys(s)

    # Denticity info
    max_dent = 1
    for lig in cc.ligands:
        if lig.dent > max_dent:
            max_dent = lig.dent

    return {
        "mol_id": mol.mol_id,
        "string": s,
        "keys": keys,
        "cn": keys.cn,
        "best_shape": keys.best_shape,
        "delta_value": keys.binned_shape.delta_value if keys.binned_shape else 0.0,
        "max_dent": max_dent,
    }


# ──────────────────────────────────────────────────────────────────
# Retention computation
# ──────────────────────────────────────────────────────────────────

def _retention_row(ref_keys, noisy_keys, extra: dict) -> dict:
    row = dict(extra)
    row["L0_match"] = int(noisy_keys.L0_StateKey == ref_keys.L0_StateKey)
    row["L1_match"] = int(noisy_keys.L1_ShapeID == ref_keys.L1_ShapeID)
    row["L2_match"] = int(noisy_keys.L2_TopoID == ref_keys.L2_TopoID)
    row["L3_match"] = int(noisy_keys.L3_ConnID == ref_keys.L3_ConnID)
    return row


def run_noise_benchmark(
    refs: List[dict],
    molecules: Dict[str, RawMolecule],
    sigmas: List[float],
    n_trials: int,
    config,
) -> pd.DataFrame:
    rows = []
    for ref in refs:
        mol = molecules[ref["mol_id"]]
        for sigma in sigmas:
            for trial in range(n_trials):
                noisy = _add_noise(mol, sigma, seed=trial * 1000 + int(sigma * 10000))
                try:
                    cc = encode_molecule(noisy, config)
                    s = cc.canonicalize().to_string()
                except Exception:
                    continue
                if not s:
                    continue
                nk = extract_identity_keys(s)
                rows.append(_retention_row(ref["keys"], nk, {
                    "mol_id": ref["mol_id"],
                    "sigma": sigma,
                    "cn": ref["cn"],
                    "best_shape": ref["best_shape"],
                    "delta_value": ref["delta_value"],
                    "max_dent": ref["max_dent"],
                }))
    return pd.DataFrame(rows)


def run_scaling_benchmark(
    refs: List[dict],
    molecules: Dict[str, RawMolecule],
    scale_factors: List[float],
    config,
) -> pd.DataFrame:
    rows = []
    for ref in refs:
        mol = molecules[ref["mol_id"]]
        for sf in scale_factors:
            scaled = _scale_bonds(mol, sf)
            try:
                cc = encode_molecule(scaled, config)
                s = cc.canonicalize().to_string()
            except Exception:
                continue
            if not s:
                continue
            nk = extract_identity_keys(s)
            rows.append(_retention_row(ref["keys"], nk, {
                "mol_id": ref["mol_id"],
                "pct_change": round((sf - 1.0) * 100, 1),
                "cn": ref["cn"],
                "best_shape": ref["best_shape"],
                "delta_value": ref["delta_value"],
                "max_dent": ref["max_dent"],
            }))
    return pd.DataFrame(rows)


# ──────────────────────────────────────────────────────────────────
# Stratification helpers
# ──────────────────────────────────────────────────────────────────

def _delta_quartile(v: float) -> str:
    if v < 1.0:
        return "Q1(<1)"
    if v < 3.0:
        return "Q2(1-3)"
    if v < 6.0:
        return "Q3(3-6)"
    return "Q4(≥6)"


def _dent_group(d: int) -> str:
    if d == 1:
        return "mono"
    if d == 2:
        return "bi"
    return "multi"


def stratified_summary(df: pd.DataFrame, group_col: str) -> pd.DataFrame:
    """Aggregate retention by *group_col* × perturbation level."""
    df = df.copy()
    df["delta_q"] = df["delta_value"].apply(_delta_quartile)
    df["dent_grp"] = df["max_dent"].apply(_dent_group)

    strata = ["cn", "best_shape", "delta_q", "dent_grp"]
    agg_cols = ["L0_match", "L1_match", "L2_match", "L3_match"]

    all_rows = []
    for stratum in strata:
        g = df.groupby([group_col, stratum])[agg_cols].mean().reset_index()
        g["stratum_type"] = stratum
        g = g.rename(columns={stratum: "stratum_value"})
        all_rows.append(g)

    return pd.concat(all_rows, ignore_index=True)


# ──────────────────────────────────────────────────────────────────
# Main
# ──────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--tmqm-dir", default="/data/CoordRep/CoordSMILES/tmQM-master/tmQM")
    parser.add_argument("--n-molecules", type=int, default=200)
    parser.add_argument("--n-trials", type=int, default=30)
    parser.add_argument("--output-dir", default="/data/CoordRep/coordrep-release/libcoordrep/outputs/identity_benchmark")
    args = parser.parse_args()

    out = Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)
    config = CoordRepConfig.default()

    # ── Load ──────────────────────────────────────────────────────
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
        r = _encode_ref(mol, config)
        if r is None:
            continue
        mol_map[mid] = mol
        refs.append(r)

    print(f"  Encoded {len(refs)} / {len(sample_ids)} molecules")

    # ── A. Noise benchmark ────────────────────────────────────────
    sigmas = [0.001, 0.005, 0.01, 0.02, 0.05, 0.10]
    print(f"\nNoise benchmark: {len(refs)} mols × {len(sigmas)} σ × {args.n_trials} trials …")
    noise_df = run_noise_benchmark(refs, mol_map, sigmas, args.n_trials, config)
    noise_df.to_csv(out / "noise_detail.csv", index=False)

    noise_ret = noise_df.groupby("sigma")[["L0_match", "L1_match", "L2_match", "L3_match"]].mean().reset_index()
    noise_ret.to_csv(out / "noise_retention.csv", index=False)

    print("\n  Noise retention (mean):")
    print(f"  {'σ':>8}  {'L0':>7}  {'L1':>7}  {'L2':>7}  {'L3':>7}")
    print("  " + "-" * 42)
    for _, r in noise_ret.iterrows():
        print(f"  {r['sigma']:>8.3f}  {r['L0_match']:>6.1%}  {r['L1_match']:>6.1%}  "
              f"{r['L2_match']:>6.1%}  {r['L3_match']:>6.1%}")

    noise_strat = stratified_summary(noise_df, "sigma")
    noise_strat.to_csv(out / "noise_stratified.csv", index=False)

    # ── B. Scaling benchmark ──────────────────────────────────────
    sfs = [0.95, 0.97, 0.98, 0.99, 1.01, 1.02, 1.03, 1.05]
    print(f"\nBond-scaling benchmark: {len(refs)} mols × {len(sfs)} factors …")
    scale_df = run_scaling_benchmark(refs, mol_map, sfs, config)
    scale_df.to_csv(out / "scaling_detail.csv", index=False)

    scale_ret = scale_df.groupby("pct_change")[["L0_match", "L1_match", "L2_match", "L3_match"]].mean().reset_index()
    scale_ret.to_csv(out / "scaling_retention.csv", index=False)

    print("\n  Scaling retention (mean):")
    print(f"  {'Δ%':>8}  {'L0':>7}  {'L1':>7}  {'L2':>7}  {'L3':>7}")
    print("  " + "-" * 42)
    for _, r in scale_ret.iterrows():
        print(f"  {r['pct_change']:>+7.1f}%  {r['L0_match']:>6.1%}  {r['L1_match']:>6.1%}  "
              f"{r['L2_match']:>6.1%}  {r['L3_match']:>6.1%}")

    scale_strat = stratified_summary(scale_df, "pct_change")
    scale_strat.to_csv(out / "scaling_stratified.csv", index=False)

    # ── C. Consensus analysis ─────────────────────────────────────
    print(f"\nConsensus analysis: {min(100, len(refs))} mols × σ=[0.005, 0.01] × 50 trials …")
    consensus_mols = [mol_map[r["mol_id"]] for r in refs[:100]]
    cons_results = compute_consensus_batch(
        consensus_mols, sigmas=[0.005, 0.01], n_trials=50, config=config,
    )
    cons_df = pd.DataFrame([cr.as_dict() for cr in cons_results])
    cons_df.to_csv(out / "consensus_results.csv", index=False)

    n_boundary = cons_df["is_geometry_boundary"].sum()
    n_total = len(cons_df)
    print(f"  Geometry-boundary structures: {n_boundary}/{n_total} "
          f"({n_boundary/n_total:.1%})")

    # ── D. False non-duplicate summary ────────────────────────────
    fndr = []
    for _, r in noise_ret.iterrows():
        fndr.append({
            "perturbation": f"noise_σ={r['sigma']:.3f}",
            "L0_false_nondup": round(1 - r["L0_match"], 4),
            "L1_false_nondup": round(1 - r["L1_match"], 4),
            "L2_false_nondup": round(1 - r["L2_match"], 4),
            "L3_false_nondup": round(1 - r["L3_match"], 4),
        })
    for _, r in scale_ret.iterrows():
        fndr.append({
            "perturbation": f"scale_{r['pct_change']:+.1f}%",
            "L0_false_nondup": round(1 - r["L0_match"], 4),
            "L1_false_nondup": round(1 - r["L1_match"], 4),
            "L2_false_nondup": round(1 - r["L2_match"], 4),
            "L3_false_nondup": round(1 - r["L3_match"], 4),
        })
    fndr_df = pd.DataFrame(fndr)
    fndr_df.to_csv(out / "false_nondup_summary.csv", index=False)

    print(f"\nAll outputs saved to {out}")


if __name__ == "__main__":
    main()
