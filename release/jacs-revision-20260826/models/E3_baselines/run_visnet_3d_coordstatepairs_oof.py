#!/usr/bin/env python
"""Run an official-PyG-trunk ViSNet CoordStatePairs baseline."""

from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import run_schnet_3d_coordstatepairs_oof as runner
from visnet_equivariant_model import MultiTargetViSNet


MODEL_LABEL = "visnet_3d"
MODEL_ARCHITECTURE = (
    "official PyG ViSNet vector-scalar representation trunk "
    "(SE(3) terminology in the source paper) with custom invariant readouts"
)
DEFAULT_OUT = Path("revision_experiments/results/coordstatepairs_visnet_3d_oof")
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
    if args.out == Path("revision_experiments/results/coordstatepairs_schnet_3d_oof"):
        args.out = DEFAULT_OUT
    # The frozen production configuration is parameter matched to SchNet:
    # 282,882 ViSNet parameters versus 301,827 SchNet parameters.
    if not option_was_supplied("--hidden"):
        args.hidden = 64
    if not option_was_supplied("--interactions"):
        args.interactions = 4
    if not option_was_supplied("--precision"):
        args.precision = "bfloat16"
    _LAST_ARGS = args
    return args


def option_was_supplied(option: str) -> bool:
    return any(token == option or token.startswith(f"{option}=") for token in sys.argv[1:])


def finalize_metadata(out: Path, preflight: dict, args: object) -> None:
    model_path = Path(__file__).with_name("visnet_equivariant_model.py").resolve()
    runner_path = Path(runner.__file__).resolve()
    entry_path = Path(__file__).resolve()
    sources = {
        model_path.name: sha256(model_path),
        runner_path.name: sha256(runner_path),
        entry_path.name: sha256(entry_path),
    }
    card = {
        "model_label": MODEL_LABEL,
        "architecture": MODEL_ARCHITECTURE,
        "reference": {
            "citation": (
                "Wang et al. Enhancing geometric representations for molecules with "
                "equivariant vector-scalar interactive message passing. Nature "
                "Communications 2024, 15, 313."
            ),
            "doi": "10.1038/s41467-023-43720-2",
        },
        "production_configuration": {
            "lmax": 1,
            "layers": int(args.interactions),
            "attention_heads": 8,
            "hidden_channels": int(args.hidden),
            "radial_basis_functions": int(args.num_gaussians),
            "cutoff_angstrom": float(args.cutoff),
            "max_neighbors": int(args.max_neighbors),
            "requested_training_precision": str(args.precision),
            "effective_training_precision": (
                str(args.precision)
                if str(preflight.get("device", "")).startswith("cuda")
                else "float32"
            ),
        },
        "inputs": ["atomic number", "optimized Cartesian coordinates in angstrom"],
        "readouts": {
            "hl_gap_ev": "invariant atomwise-mean scalar head",
            "dipole_moment": (
                "norm of a COM-centered latent-scalar-weight vector; latent weights are "
                "not interpreted as physical charges"
            ),
        },
        "comparison_role": (
            "Modern equivariant, approximately parameter-matched common-cohort "
            "comparator using z+coordinates only; not an equal-information ablation "
            "against the chemically attributed 2D graph or CoordRep-field models."
        ),
        "equivariance_preflight": preflight,
        "source_sha256": sources,
    }
    (out / "MODEL_CARD.json").write_text(
        json.dumps(card, indent=2, allow_nan=False), encoding="utf-8"
    )
    metadata_path = out / "run_metadata.json"
    if metadata_path.exists():
        metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
        metadata["model"] = MODEL_LABEL
        metadata["architecture"] = card["architecture"]
        metadata["model_source_sha256"] = sources
        metadata_path.write_text(
            json.dumps(metadata, indent=2, allow_nan=False), encoding="utf-8"
        )


def main() -> None:
    runner._ORIGINAL_PARSE_ARGS = runner.parse_args
    runner.parse_args = patched_parse_args
    runner.MultiTargetSchNet = MultiTargetViSNet
    runner.MODEL_LABEL = MODEL_LABEL
    runner.MODEL_ARCHITECTURE = MODEL_ARCHITECTURE
    runner.MODEL_GUARDRAIL = (
        "ViSNet receives atomic numbers and optimized Cartesian coordinates but omits "
        "bond order and formal/total charge attributes available to the chemically "
        "attributed 2D graph or CoordRep fields. It is a modern equivariant, non-nested "
        "common-cohort comparator, not an equal-information representation ablation or a "
        "uniformly retuned leaderboard."
    )
    runner.main()
    if _LAST_ARGS is None:
        raise RuntimeError("Runner arguments were not captured")
    out = Path(_LAST_ARGS.out).resolve()
    preflight_path = out / "INVARIANCE_PREFLIGHT.json"
    if not preflight_path.is_file():
        raise RuntimeError("Shared dense invariance preflight was not written")
    preflight = json.loads(preflight_path.read_text(encoding="utf-8"))
    finalize_metadata(out, preflight, _LAST_ARGS)


if __name__ == "__main__":
    main()
