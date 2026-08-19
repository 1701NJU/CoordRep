#!/usr/bin/env python
"""Run a PaiNN E(3)-equivariant baseline on the frozen CoordStatePairs cohort.

This entry point reuses the audited SchNet OOF data/split/statistics runner and
changes only the model factory and model label. Consequently, cohort IDs,
graph-grouped outer folds, deterministic validation groups, target scaling,
seeds, stopping rule, pair manifest, metrics, and output schemas remain
identical to the SchNet comparison.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import torch
from torch_geometric.data import Data

import run_schnet_3d_coordstatepairs_oof as runner
from painn_equivariant_model import MultiTargetPaiNN


MODEL_LABEL = "painn_3d"
DEFAULT_OUT = Path("revision_experiments/results/coordstatepairs_painn_3d_oof")
_LAST_ARGS = None


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def patched_parse_args():
    global _LAST_ARGS
    args = runner._ORIGINAL_PARSE_ARGS()
    schnet_default = Path("revision_experiments/results/coordstatepairs_schnet_3d_oof")
    if args.out == schnet_default:
        args.out = DEFAULT_OUT
    _LAST_ARGS = args
    return args


def random_rotation(seed: int = 20260803) -> torch.Tensor:
    generator = torch.Generator().manual_seed(seed)
    matrix = torch.randn(3, 3, generator=generator)
    q, _ = torch.linalg.qr(matrix)
    if torch.det(q) < 0:
        q[:, 0] *= -1
    return q


@torch.no_grad()
def equivariance_preflight(out: Path, device: torch.device) -> dict:
    torch.manual_seed(20260803)
    z = torch.tensor([26, 7, 7, 8, 6, 6, 1, 1], dtype=torch.long)
    pos = torch.tensor(
        [
            [0.00, 0.00, 0.00],
            [1.85, 0.05, 0.10],
            [-1.82, -0.08, 0.04],
            [0.12, 1.91, -0.07],
            [-0.10, -1.88, 0.12],
            [0.20, 0.15, 2.02],
            [0.55, 0.42, 2.89],
            [-0.50, -0.35, -1.10],
        ],
        dtype=torch.float32,
    )
    model = MultiTargetPaiNN(
        hidden=64,
        interactions=3,
        num_gaussians=24,
        cutoff=5.0,
        max_neighbors=32,
        n_targets=2,
        target_mean=np.asarray([2.9, 5.7], dtype=np.float32),
        target_scale=np.asarray([0.93, 4.14], dtype=np.float32),
    ).to(device)
    model.eval()

    def evaluate(z_value: torch.Tensor, pos_value: torch.Tensor) -> np.ndarray:
        data = Data(z=z_value, pos=pos_value)
        data.batch = torch.zeros(len(z_value), dtype=torch.long)
        return model(data.to(device)).detach().cpu().numpy().reshape(-1)

    base = evaluate(z, pos)
    rotation = random_rotation()
    rotated = evaluate(z, pos @ rotation.T)
    reflected_pos = pos.clone()
    reflected_pos[:, 0] *= -1
    reflected = evaluate(z, reflected_pos)
    translated = evaluate(z, pos + torch.tensor([3.1, -2.7, 5.4]))
    permutation = torch.tensor([5, 2, 7, 0, 4, 1, 6, 3], dtype=torch.long)
    permuted = evaluate(z[permutation], pos[permutation])
    comparisons = {
        "proper_rotation": rotated,
        "reflection": reflected,
        "translation": translated,
        "atom_permutation": permuted,
    }
    differences = {
        name: float(np.max(np.abs(value - base))) for name, value in comparisons.items()
    }
    audit = {
        "model": MODEL_LABEL,
        "device": str(device),
        "torch_version": torch.__version__,
        "base_prediction": base.tolist(),
        "transformed_predictions": {name: value.tolist() for name, value in comparisons.items()},
        "maximum_absolute_differences": differences,
        "tolerance": 1.0e-4,
        "passed": bool(max(differences.values()) <= 1.0e-4),
        "scope": (
            "Numerical invariance of the scalar outputs under one proper rotation, one "
            "reflection, one translation, and one atom permutation; this is an implementation "
            "regression test, not an empirical property benchmark."
        ),
    }
    out.mkdir(parents=True, exist_ok=True)
    (out / "EQUIVARIANCE_PREFLIGHT.json").write_text(
        json.dumps(audit, indent=2), encoding="utf-8"
    )
    if not audit["passed"]:
        raise RuntimeError(f"PaiNN invariance preflight failed: {differences}")
    return audit


def finalize_metadata(out: Path, preflight: dict) -> None:
    model_path = Path(__file__).with_name("painn_equivariant_model.py").resolve()
    runner_path = Path(runner.__file__).resolve()
    entry_path = Path(__file__).resolve()
    card = {
        "model_label": MODEL_LABEL,
        "architecture": "PaiNN-style scalar/vector E(3)-equivariant message passing",
        "reference": (
            "Schuett, Unke, Gastegger. Equivariant message passing for the prediction of "
            "tensorial properties and molecular spectra. ICML 2021."
        ),
        "inputs": ["atomic number", "Cartesian coordinates in angstrom"],
        "readouts": {
            "hl_gap_ev": "invariant mean-pooled scalar head",
            "dipole_moment": "norm of an equivariant COM-centered latent-charge vector",
        },
        "comparison_role": (
            "Established E(3)-equivariant common-cohort comparator; not an equal-information "
            "ablation against 2D graph or CoordRep-field models."
        ),
        "equivariance_preflight": preflight,
        "source_sha256": {
            model_path.name: sha256(model_path),
            runner_path.name: sha256(runner_path),
            entry_path.name: sha256(entry_path),
        },
    }
    (out / "MODEL_CARD.json").write_text(json.dumps(card, indent=2), encoding="utf-8")

    metadata_path = out / "run_metadata.json"
    if metadata_path.exists():
        metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
        metadata["model"] = MODEL_LABEL
        metadata["architecture"] = card["architecture"]
        metadata["model_source_sha256"] = card["source_sha256"]
        metadata_path.write_text(json.dumps(metadata, indent=2), encoding="utf-8")


def main() -> None:
    runner._ORIGINAL_PARSE_ARGS = runner.parse_args
    runner.parse_args = patched_parse_args
    runner.MultiTargetSchNet = MultiTargetPaiNN
    runner.MODEL_LABEL = MODEL_LABEL
    runner.MODEL_GUARDRAIL = (
        "PaiNN receives optimized Cartesian coordinates as extra information relative to 2D "
        "graph/record inputs. It is an established E(3)-equivariant common-cohort comparator, "
        "not an equal-information representation ablation or a uniformly retuned leaderboard."
    )

    preview_args = runner._ORIGINAL_PARSE_ARGS()
    out = DEFAULT_OUT if preview_args.out == Path(
        "revision_experiments/results/coordstatepairs_schnet_3d_oof"
    ) else preview_args.out
    # Restore argv parsing for the actual runner: patched_parse_args will parse once more.
    preflight = equivariance_preflight(out.resolve(), torch.device(preview_args.device))
    runner.main()
    if _LAST_ARGS is None:
        raise RuntimeError("Runner arguments were not captured")
    finalize_metadata(Path(_LAST_ARGS.out).resolve(), preflight)


if __name__ == "__main__":
    main()
