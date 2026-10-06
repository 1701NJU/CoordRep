#!/usr/bin/env python3
"""Inference-only access to the three archived CoordRep diagnostic models.

Inputs are explicit, pretokenized token lists, including [CLS] and [SEP].
This entry point does not tokenize records, add vocabulary entries, truncate
inputs, train models, or reproduce a full experimental cohort.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass, fields
import json
from pathlib import Path
import sys
from typing import Any, Sequence
import warnings

import torch

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from brain.model import CoordRepForMLM, CoordRepModelConfig
from ranker import CoordRepRanker


MODEL_KEYS = ("production_mlm", "factorized_mlm", "semantic_ranker")
MAX_TOKENS = 768
SPECIAL_TOKENS = ("[CLS]", "[SEP]", "[PAD]", "[MASK]", "[UNK]")


@dataclass
class LoadedModel:
    """A loaded model and its fixed vocabulary, for CLI or Python callers."""

    model_key: str
    model: torch.nn.Module
    config: CoordRepModelConfig
    token2id: dict[str, int]
    id2token: dict[int, str]
    device: torch.device


def _read_json(path: Path) -> Any:
    with path.open(encoding="utf-8") as handle:
        return json.load(handle)


def load_model(
    bundle: str | Path, model_key: str, device: str = "cpu"
) -> LoadedModel:
    """Safely load a tensor-only state_dict with an explicit architecture.

    Each bundle/<model_key>/ directory must contain model.pt, config.json and
    tokenizer.json. Architecture fields may be at the config root or in an
    ``architecture`` object; all dataclass fields must be present. For the
    ranker, ``pooling`` is read from the config root, never inferred.
    """
    if model_key not in MODEL_KEYS:
        raise ValueError(f"Unknown model key: {model_key!r}")
    model_dir = Path(bundle).resolve() / model_key
    payload = _read_json(model_dir / "config.json")
    if not isinstance(payload, dict):
        raise ValueError("config.json must contain an object")
    architecture = payload.get("architecture", payload)
    if not isinstance(architecture, dict):
        raise ValueError("The architecture configuration must be an object")
    field_names = {field.name for field in fields(CoordRepModelConfig)}
    missing = sorted(field_names - architecture.keys())
    if missing:
        raise ValueError(f"Explicit architecture fields are missing: {missing}")
    config = CoordRepModelConfig(
        **{name: architecture[name] for name in field_names}
    )
    if config.max_position_embeddings != MAX_TOKENS:
        raise ValueError("These archived models require 768 position embeddings")
    if config.vocab_size != 10000:
        raise ValueError("These archived models require vocab_size=10000")

    tokenizer = _read_json(model_dir / "tokenizer.json")
    token2id = tokenizer.get("token2id") if isinstance(tokenizer, dict) else None
    if not isinstance(token2id, dict) or not token2id:
        raise ValueError("tokenizer.json must contain a nonempty token2id object")
    for token, token_id in token2id.items():
        if not isinstance(token, str) or type(token_id) is not int:
            raise ValueError("Vocabulary entries must map strings to integer IDs")
        if not 0 <= token_id < config.vocab_size:
            raise ValueError(f"Vocabulary ID out of range for token {token!r}")
    if len(set(token2id.values())) != len(token2id):
        raise ValueError("Vocabulary IDs must be unique")
    missing_special = [token for token in SPECIAL_TOKENS if token not in token2id]
    if missing_special:
        raise ValueError(f"Special tokens are missing: {missing_special}")
    if token2id["[PAD]"] != config.pad_token_id:
        raise ValueError("Tokenizer PAD ID does not match the architecture")

    # Restricted tensor loading only: do not fall back to unsafe pickle loading.
    state = torch.load(
        model_dir / "model.pt", map_location="cpu", weights_only=True
    )
    if not isinstance(state, dict) or not state:
        raise ValueError("model.pt must contain a nonempty, pure state_dict")
    if not all(isinstance(key, str) and torch.is_tensor(value)
               for key, value in state.items()):
        raise ValueError("model.pt must contain only string-to-tensor entries")
    if model_key == "semantic_ranker":
        pooling = payload.get("pooling")
        if pooling not in ("mean", "cls"):
            raise ValueError("Ranker config must explicitly specify pooling")
        model = CoordRepRanker(config, pooling=pooling)
    else:
        model = CoordRepForMLM(config)
    model.load_state_dict(state, strict=True)
    target_device = torch.device(device)
    model.eval().to(target_device)
    return LoadedModel(
        model_key, model, config, dict(token2id),
        {token_id: token for token, token_id in token2id.items()}, target_device
    )


def read_tokens(path: str | Path) -> list[str]:
    """Read a JSON token list, or an object containing a ``tokens`` list."""
    data = _read_json(Path(path))
    tokens = data.get("tokens") if isinstance(data, dict) else data
    if not isinstance(tokens, list):
        raise ValueError("Token JSON must be a list or an object with a tokens list")
    if not all(isinstance(token, str) and token for token in tokens):
        raise ValueError("Each input token must be a nonempty string")
    return tokens


def token_ids(loaded: LoadedModel, tokens: Sequence[str]) -> list[int]:
    """Map against the frozen vocabulary, warning on UNK substitution."""
    if isinstance(tokens, (str, bytes)) or not tokens:
        raise ValueError("Input must be a nonempty sequence of token strings")
    if not all(isinstance(token, str) and token for token in tokens):
        raise ValueError("Each input token must be a nonempty string")
    if len(tokens) > MAX_TOKENS:
        raise ValueError(f"Input has {len(tokens)} tokens; maximum is {MAX_TOKENS}")
    if tokens[0] != "[CLS]" or tokens[-1] != "[SEP]":
        raise ValueError("Input must explicitly start with [CLS] and end with [SEP]")
    unknown = sorted({token for token in tokens if token not in loaded.token2id})
    if unknown:
        warnings.warn(
            f"Unknown tokens mapped to [UNK], without vocabulary insertion: {unknown}",
            RuntimeWarning, stacklevel=2
        )
    unk_id = loaded.token2id["[UNK]"]
    return [loaded.token2id.get(token, unk_id) for token in tokens]


@torch.inference_mode()
def run_inference(loaded: LoadedModel, tokens: Sequence[str]) -> torch.Tensor:
    """Return raw CPU logits [1,L,V] or a raw ranker score tensor [1].

    The raw tensor is useful for checking parity against original checkpoints.
    Ranker values are compatibility scores, not calibrated probabilities.
    """
    ids = torch.tensor([token_ids(loaded, tokens)], dtype=torch.long,
                       device=loaded.device)
    attention_mask = (ids != loaded.token2id["[PAD]"]).long()
    result = loaded.model(ids, attention_mask)
    tensor = result if loaded.model_key == "semantic_ranker" else result[0]
    if not torch.isfinite(tensor).all():
        raise RuntimeError("Inference returned a non-finite value")
    return tensor.detach().cpu()


def prediction_json(
    loaded: LoadedModel, tokens: Sequence[str], output: torch.Tensor, top_k: int = 5
) -> dict[str, Any]:
    """Present ranker scores or full-vocabulary MLM top-k at [MASK] positions."""
    result: dict[str, Any] = {
        "model": loaded.model_key,
        "n_tokens": len(tokens),
        "unknown_tokens": sorted(set(tokens) - loaded.token2id.keys()),
    }
    if loaded.model_key == "semantic_ranker":
        result["compatibility_score"] = float(output[0])
        result["score_is_probability"] = False
        return result
    if not 1 <= top_k <= loaded.config.vocab_size:
        raise ValueError(f"top_k must be between 1 and {loaded.config.vocab_size}")
    mask_positions = [index for index, token in enumerate(tokens) if token == "[MASK]"]
    if not mask_positions:
        raise ValueError("MLM input must contain at least one explicit [MASK] token")
    predictions = []
    for position in mask_positions:
        values, indices = torch.softmax(output[0, position], dim=-1).topk(top_k)
        predictions.append({
            "position_0based": position,
            "top_k": [
                {"token_id": int(token_id),
                 "token": loaded.id2token.get(int(token_id)),
                 "probability": float(probability)}
                for probability, token_id in zip(values, indices)
            ],
        })
    result["masked_predictions"] = predictions
    result["unassigned_vocabulary_token_is_null"] = True
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bundle", type=Path, default=HERE / "bundle",
                        help="Directory containing the three model subdirectories")
    parser.add_argument("--model", choices=MODEL_KEYS)
    parser.add_argument("--tokens-json", type=Path,
                        help="Pretokenized input JSON; no record tokenization is performed")
    parser.add_argument("--top-k", type=int, default=5)
    parser.add_argument("--device", default="cpu", choices=("cpu", "cuda"))
    parser.add_argument("--smoke-test", action="store_true",
                        help="Run the supplied synthetic example for each selected model")
    args = parser.parse_args()
    if args.top_k < 1:
        parser.error("--top-k must be positive")
    if args.smoke_test and args.tokens_json:
        parser.error("--smoke-test cannot be combined with --tokens-json")
    if not args.smoke_test and (args.model is None or args.tokens_json is None):
        parser.error("Supply --model and --tokens-json, or use --smoke-test")
    keys = MODEL_KEYS if args.smoke_test and args.model is None else (args.model,)
    results = []
    try:
        for key in keys:
            path = HERE / "examples" / f"{key}.json" if args.smoke_test else args.tokens_json
            tokens = read_tokens(path)
            loaded = load_model(args.bundle, key, args.device)
            output = run_inference(loaded, tokens)
            results.append(prediction_json(loaded, tokens, output, args.top_k))
            del loaded, output
    except (OSError, ValueError, RuntimeError) as exc:
        parser.error(str(exc))
    print(json.dumps(results if args.smoke_test else results[0], indent=2,
                     ensure_ascii=False, allow_nan=False))


if __name__ == "__main__":
    main()
