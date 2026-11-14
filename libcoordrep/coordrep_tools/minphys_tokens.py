"""
minphys_tokens.py
=================
Generate auxiliary MinPhys tokens and append them to CoordRep token lists.

Per-ligand token:
  <Aux:L1|HA:small|Charge:neutral|Donor:borderline>

Global token:
  <AuxGlobal|SizeLoad:medium|ChargeSum:neutral|HardSoft:mixed>

These are NOT part of the CoordRep core grammar.  They are only used
for the MinPhys diagnostic ablation in the Ranker experiments.
"""

from __future__ import annotations

from typing import List

from .minphys_features import MinPhysBinner, extract_ligand_features


def make_minphys_tokens(tokens: List[str], binner: MinPhysBinner) -> List[str]:
    """
    Given a CoordRep token list and a fitted MinPhysBinner, return
    the list of auxiliary MinPhys tokens (per-ligand + global).

    Important: for hard decoys, call this on the **decoy** token list
    so that Aux tokens reflect the *replaced* ligands.
    """
    feat = binner.featurise_complex(tokens)

    aux_tokens = []
    for pl in feat["per_ligand"]:
        tok = (f"<Aux:{pl['lid']}"
               f"|HA:{pl['ha_bin']}"
               f"|Charge:{pl['charge_bin']}"
               f"|Donor:{pl['donor_hs']}>")
        aux_tokens.append(tok)

    g = feat["global"]
    aux_tokens.append(
        f"<AuxGlobal"
        f"|SizeLoad:{g['size_load']}"
        f"|ChargeSum:{g['charge_sum']}"
        f"|HardSoft:{g['hardsoft']}>"
    )

    return aux_tokens


def append_minphys(tokens: List[str], binner: MinPhysBinner) -> List[str]:
    """Return a new token list with MinPhys aux tokens appended."""
    aux = make_minphys_tokens(tokens, binner)
    return list(tokens) + aux


def register_minphys_tokens_in_tokenizer(tokenizer, train_token_lists: List[List[str]],
                                          binner: MinPhysBinner):
    """
    Pre-generate all possible MinPhys tokens from training data and
    add them to the tokenizer vocabulary.  This ensures consistent
    token IDs for the ranker.
    """
    seen = set()
    for toks in train_token_lists:
        aux = make_minphys_tokens(toks, binner)
        seen.update(aux)

    # Also add canonical aux tokens that may not appear in training data
    for ha in ["small", "medium", "large", "unk"]:
        for ch in ["anionic", "neutral", "cationic", "mixed", "unk"]:
            for dn in ["hard", "borderline", "soft", "mixed", "unk"]:
                for lid in ["L1", "L2", "L3", "L4", "L5", "L6", "L7", "L8"]:
                    seen.add(f"<Aux:{lid}|HA:{ha}|Charge:{ch}|Donor:{dn}>")
                seen.add(f"<AuxGlobal|SizeLoad:{ha}|ChargeSum:{ch}|HardSoft:{dn}>")

    for tok in sorted(seen):
        if tok not in tokenizer.token2id:
            tid = tokenizer.next_id
            if tid < tokenizer.config.vocab_size:
                tokenizer.token2id[tok] = tid
                tokenizer.id2token[tid] = tok
                tokenizer.next_id = tid + 1

    return len(seen)
