#!/usr/bin/env python3
"""
tool_a_ablation_eval.py
=======================

Evaluate Tool A donor-marker prediction under ablation modes,
with stratified metrics by denticity, CN, metal row, donor element,
and ligand frequency.

Public entry point: ``evaluate_tool_a_ablation()``.
"""

from __future__ import annotations

import math
import re
from collections import Counter, defaultdict
from typing import Dict, List, Optional, Tuple

import numpy as np

from .ablation_masking import make_ablation_input, ABLATION_MODES
from .tool_a_baselines import _extract_ligand_smiles

# Donor atoms list
DONOR_ATOMS = ['C', 'N', 'O', 'S', 'P', 'F', 'Cl', 'Br', 'I',
               'Se', 'Te', 'As', 'Si', 'B', 'H']

_METAL_RE = re.compile(r'\[Metal:(\w+)\|.*?CN:(\d+)\]')

_3D_METALS = {"Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn"}
_4D_METALS = {"Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd"}
_5D_METALS = {"La", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg"}


def _metal_row(metal: str) -> str:
    if metal in _3D_METALS:
        return "3d"
    if metal in _4D_METALS:
        return "4d"
    if metal in _5D_METALS:
        return "5d"
    return "other"


def _cn_group(cn) -> str:
    cn = int(cn) if cn else 0
    if cn in (4, 5, 6):
        return str(cn)
    return "other"


def _donor_element_group(donor: str) -> str:
    if donor in ("N", "O", "S", "P", "C"):
        return donor
    if donor in ("F", "Cl", "Br", "I"):
        return "halide"
    return "other"


def _get_ligand_id_for_pos(tokens: list, mask_pos: int) -> str:
    """Find the ligand ID (L1, L2, ...) containing *mask_pos*.

    Handles single-token (``L1``) and split (``L``, ``2``) variants.
    """
    if not isinstance(mask_pos, int):
        return "?"
    for i in range(mask_pos, -1, -1):
        if re.match(r'^L\d+$', tokens[i]):
            return tokens[i]
        # split variant: "L" followed by digit
        if (tokens[i] == "L" and i + 1 < len(tokens)
                and tokens[i + 1].isdigit()):
            return f"L{tokens[i + 1]}"
    return "?"


def _get_complex_id(sid: str) -> str:
    m = re.match(r'^(tmqm_\d+)_donor_\d+$', sid)
    return m.group(1) if m else sid


def build_denticity_map(samples: List[dict]) -> Dict[str, str]:
    """Pre-compute denticity by grouping samples by complex + ligand."""
    from collections import defaultdict
    complex_ligs = defaultdict(lambda: defaultdict(list))
    for s in samples:
        sid = s.get("id", "")
        cid = _get_complex_id(sid)
        tokens = s.get("tokens", [])
        mp = s.get("mask_pos")
        lid = _get_ligand_id_for_pos(tokens, mp) if isinstance(mp, int) else "?"
        complex_ligs[cid][lid].append(sid)

    dent_map = {}
    for cid, ligs in complex_ligs.items():
        for lid, sids in ligs.items():
            d = len(sids)
            label = f"dent={d}" if d <= 2 else "dent>=3"
            for sid in sids:
                dent_map[sid] = label
    return dent_map


# ──────────────────────────────────────────────────────────────────
# Main evaluation function
# ──────────────────────────────────────────────────────────────────

def evaluate_tool_a_ablation(
    samples: List[dict],
    predict_fn,
    ablation_modes: List[str],
    top_k: List[int] = (1, 5),
    confidence_thresholds: List[float] = (0.5, 0.7, 0.85, 0.95),
    ligand_freq_counter: Optional[Counter] = None,
    rare_ligand_threshold: int = 5,
) -> Dict:
    """
    Parameters
    ----------
    samples : list[dict]
        Each has ``tokens``, ``mask_pos``, ``y_true``, ``group``.
    predict_fn : callable
        ``predict_fn(ablated_tokens, mask_positions) -> list[list[(token, prob)]]``
        Returns top-k predictions per mask position.
    ablation_modes : list[str]
        Subset of ABLATION_MODES to evaluate.
    top_k : list[int]
    confidence_thresholds : list[float]
    ligand_freq_counter : Counter, optional
        Maps ligand SMILES → count in training set.
    rare_ligand_threshold : int

    Returns
    -------
    dict with keys:
        per_mode : dict[mode -> aggregate metrics]
        per_prediction : list[dict]  (one row per sample × mode)
        stratified : dict[stratum_name -> DataFrame-ready rows]
        ligand_level : list[dict]  per-ligand all-correct
    """
    per_prediction = []
    per_mode_accum = defaultdict(lambda: defaultdict(list))

    # Pre-compute ligand frequency bucket
    _lig_freq = ligand_freq_counter or Counter()

    # Pre-compute denticity map
    _dent_map = build_denticity_map(samples)

    for s in samples:
        tokens = s["tokens"]
        mask_pos = s.get("mask_pos")
        if isinstance(mask_pos, int):
            donor_positions = [mask_pos]
        elif isinstance(mask_pos, list):
            donor_positions = mask_pos
        else:
            continue

        y_true = s.get("y_true", {}).get("donor_atom")
        if not y_true:
            continue

        g = s.get("group", {})
        metal = g.get("metal_element", "?")
        cn = g.get("cn", 0)
        dent = _dent_map.get(s.get("id", ""), "?")
        mrow = _metal_row(metal)
        delem = _donor_element_group(y_true)
        cn_grp = _cn_group(cn)

        # Ligand frequency bucket
        lig_smi = _extract_ligand_smiles(tokens, donor_positions[0]) if donor_positions else "?"
        lig_freq_val = _lig_freq.get(lig_smi, 0)
        lig_bucket = "common" if lig_freq_val >= rare_ligand_threshold else "rare"

        for mode in ablation_modes:
            try:
                abl_tokens, abl_pos = make_ablation_input(
                    tokens, donor_positions, mode, shuffle_seed=42,
                )
            except Exception:
                continue

            if not abl_pos:
                continue

            try:
                preds_per_pos = predict_fn(abl_tokens, abl_pos)
            except Exception:
                continue

            if not preds_per_pos:
                continue

            for pos_idx, topk_preds in enumerate(preds_per_pos):
                if not topk_preds:
                    continue

                pred_tokens = [t for t, p in topk_preds]
                pred_probs = [p for t, p in topk_preds]
                top1_token = pred_tokens[0] if pred_tokens else ""
                top1_prob = pred_probs[0] if pred_probs else 0.0

                is_top1 = top1_token == y_true
                is_top5 = y_true in pred_tokens[:5]
                rank = (pred_tokens.index(y_true) + 1) if y_true in pred_tokens else 0
                rr = 1.0 / rank if rank > 0 else 0.0

                # NLL
                true_prob = 0.0
                for t, p in topk_preds:
                    if t == y_true:
                        true_prob = p
                        break
                nll = -math.log(true_prob + 1e-12)

                row = {
                    "id": s.get("id", ""),
                    "ablation_mode": mode,
                    "ligand_id": "",
                    "donor_position": abl_pos[pos_idx] if pos_idx < len(abl_pos) else -1,
                    "true_token": y_true,
                    "top1_pred": top1_token,
                    "top5_pred": pred_tokens[:5],
                    "probability": top1_prob,
                    "true_prob": true_prob,
                    "nll": nll,
                    "rank": rank,
                    "rr": rr,
                    "correct_top1": is_top1,
                    "correct_top5": is_top5,
                    "denticity": dent,
                    "CN": cn,
                    "cn_group": cn_grp,
                    "metal": metal,
                    "metal_row": mrow,
                    "donor_element": delem,
                    "ligand_smiles": lig_smi,
                    "ligand_frequency_bucket": lig_bucket,
                    "ligand_smiles_available": mode not in (
                        "no_ligand_smiles_keep_length",
                        "no_ligand_smiles_collapsed",
                        "metal_cn_only",
                    ),
                    "geometry_available": mode not in (
                        "no_geometry_stereo", "metal_cn_only",
                    ),
                    "stereo_available": mode not in (
                        "no_geometry_stereo", "metal_cn_only",
                    ),
                }

                per_prediction.append(row)

                # Accum for aggregate
                per_mode_accum[mode]["top1"].append(int(is_top1))
                per_mode_accum[mode]["top5"].append(int(is_top5))
                per_mode_accum[mode]["rr"].append(rr)
                per_mode_accum[mode]["nll"].append(nll)

                # Confidence buckets
                for thr in confidence_thresholds:
                    key = f"conf_{thr}"
                    if top1_prob >= thr:
                        per_mode_accum[mode][f"{key}_covered"].append(1)
                        per_mode_accum[mode][f"{key}_correct"].append(int(is_top1))
                    else:
                        per_mode_accum[mode][f"{key}_covered"].append(0)

    # ── Aggregate per mode ────────────────────────────────────────
    per_mode = {}
    for mode, acc in per_mode_accum.items():
        n = len(acc["top1"])
        result = {
            "n": n,
            "top1_acc": float(np.mean(acc["top1"])) if n else 0,
            "top5_acc": float(np.mean(acc["top5"])) if n else 0,
            "mrr": float(np.mean(acc["rr"])) if n else 0,
            "mean_nll": float(np.mean(acc["nll"])) if n else 0,
        }
        for thr in confidence_thresholds:
            key = f"conf_{thr}"
            covered = acc.get(f"{key}_covered", [])
            correct = acc.get(f"{key}_correct", [])
            cov_rate = float(np.mean(covered)) if covered else 0
            acc_at_cov = float(np.mean(correct)) if correct else 0
            result[f"coverage_{thr}"] = cov_rate
            result[f"accuracy_at_{thr}"] = acc_at_cov
        per_mode[mode] = result

    # ── Stratified summaries ──────────────────────────────────────
    strata_keys = {
        "denticity": "denticity",
        "CN": "cn_group",
        "metal_row": "metal_row",
        "donor_element": "donor_element",
        "ligand_frequency": "ligand_frequency_bucket",
    }
    stratified = {}
    for stratum_name, field in strata_keys.items():
        rows = []
        groups = defaultdict(lambda: defaultdict(list))
        for r in per_prediction:
            key = (r["ablation_mode"], r[field])
            groups[key]["top1"].append(int(r["correct_top1"]))
            groups[key]["top5"].append(int(r["correct_top5"]))
            groups[key]["rr"].append(r["rr"])
        for (mode, val), acc in groups.items():
            n = len(acc["top1"])
            rows.append({
                "ablation_mode": mode,
                stratum_name: val,
                "n": n,
                "top1_acc": round(float(np.mean(acc["top1"])), 4),
                "top5_acc": round(float(np.mean(acc["top5"])), 4),
                "mrr": round(float(np.mean(acc["rr"])), 4),
            })
        stratified[stratum_name] = rows

    # ── Per-ligand all-correct ────────────────────────────────────
    # Group predictions by (id, mode, ligand_smiles) and check if all
    # donors in that ligand are correct
    lig_groups = defaultdict(list)
    for r in per_prediction:
        key = (r["id"], r["ablation_mode"], r["ligand_smiles"])
        lig_groups[key].append(r)

    ligand_level = []
    for (cid, mode, smi), preds in lig_groups.items():
        all_correct = all(p["correct_top1"] for p in preds)
        n_donors = len(preds)
        dent = preds[0]["denticity"] if preds else "?"
        ligand_level.append({
            "id": cid,
            "ablation_mode": mode,
            "ligand_smiles": smi,
            "denticity": dent,
            "n_donors": n_donors,
            "all_correct": all_correct,
        })

    return {
        "per_mode": per_mode,
        "per_prediction": per_prediction,
        "stratified": stratified,
        "ligand_level": ligand_level,
    }
