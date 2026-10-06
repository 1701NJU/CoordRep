# Released CoordRep pretrained weights

This additive release publishes three archived CoordRep inference checkpoints,
their exact vocabularies, explicit architectures, a safe inference entry point,
and synthetic examples. No model was retrained and no reported result was
changed. The original checkpoint SHA-256 digests match the frozen field-order
analysis manifest. Every exported tensor was checked for exact equality with
its archived source. Optimizer states and training metadata are omitted.

| Model | Published vocabulary | Verified analysis mapping |
|---|---:|---|
| `production_mlm` | 657 tokens | Supplementary Table S5 production-MLM frozen inference; Table S14B conditional donor-field inference |
| `factorized_mlm` | 376 tokens | Table S5 archived factorized-MLM frozen inference; **not** the separately retrained three-seed Table S4 models |
| `semantic_ranker` | 657 tokens | Table S5 semantic-decoy ranking frozen inference |

The tensor-only checkpoints are downloadable assets on the
[pretrained-model release](https://github.com/1701NJU/CoordRep/releases/tag/coordrep-pretrained-20261006).
The code, configurations, matching tokenizer files and manifests are in this
directory on the existing `jacs-revision` branch. Checkpoint assets are not
stored inside the frozen 2026-08-26 or 2026-09-29 evidence packages.

## Download and inference

Use Python 3.10 or later and PyTorch with `weights_only=True` support. CPU
verification used Python 3.12.14 and PyTorch 2.5.1. From this directory:

```bash
python download_weights.py --output downloaded
python inference.py --bundle downloaded/models --smoke-test
python inference.py --bundle downloaded/models --model production_mlm --tokens-json examples/production_mlm.json --top-k 5
python inference.py --bundle downloaded/models --model factorized_mlm --tokens-json examples/factorized_mlm.json --top-k 5
python inference.py --bundle downloaded/models --model semantic_ranker --tokens-json examples/semantic_ranker.json
```

The downloader checks the archive SHA-256 before extraction and all model,
configuration and tokenizer digests afterwards. Existing output directories
are not overwritten. `manifest.json` records both source and exported
checkpoint digests; `download_manifest.json` records the archive digest.

`inference.py` loads only tensor state dictionaries and requires all
architecture fields, including 768 position embeddings. It takes explicit
pretokenized lists with `[CLS]` and `[SEP]`. It does not insert vocabulary
entries or silently truncate inputs. Unknown tokens map to `[UNK]` with a
warning. MLM outputs are masked-token predictions; the ranker output is an
uncalibrated compatibility score, not a molecular property prediction.

## Frozen sequence interface

These models use their archived sequence interface, not retroactively
re-encoded rc3 inputs. The production and factorized vocabularies are not
interchangeable. `brain/tokenizer.py` is the original production tokenizer:
it retains composite metal-header tokens. The later tokenizer in the frozen
2026-08-26 code package factorizes such headers and must not replace this
production interface. Explicit token lists are used in the new inference
entry point to avoid changing the frozen input definition.

`brain/model.py` and `ranker.py` preserve the archived architectures. The
historical training/loading helpers in those modules are not used by the
new inference entry point. Use `inference.py` for restricted checkpoint
loading, rather than historical helpers.

## Scope and provenance

The archived production training configuration identifies tmQM as its input
source (`cod_dir` is null). This release distributes learned parameters and
code under the accompanying MIT License, not a license to redistribute the
underlying datasets. No CIF, coordinate file, licensed CSD record, row-level
property input, fold/pair manifest, or row-level prediction is included.

The separately retrained three-seed tokenizer models in Table S4 and the
property-prediction models underlying Figure 4 and Tables S8/S8A did not save
checkpoint files. They are **not** included in this release. Their existing
configurations, protocols and permitted result summaries remain available in
the earlier evidence packages. Released pretrained weights permit inference
with the published models; they do not by themselves reproduce every cohort
analysis or all property predictions from aggregate CSV files.

See `verification.json` for strict model loading, archived/exported tensor
equality and synthetic CPU inference-parity checks. Those checks establish
faithful packaging, not a new benchmark result.
