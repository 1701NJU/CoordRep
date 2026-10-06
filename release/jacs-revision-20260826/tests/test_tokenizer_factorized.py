#!/usr/bin/env python3
"""
test_tokenizer_factorized.py
============================
Tests that the tokenizer always emits factorized metal tokens,
regardless of input format (composite pipe or factorized semicolon).
"""

import sys
import os
import importlib

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

# Import tokenizer module directly to avoid brain/__init__.py pulling in torch
_spec = importlib.util.spec_from_file_location(
    "brain.tokenizer",
    os.path.join(os.path.dirname(__file__), "..", "brain", "tokenizer.py"),
)
_mod = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_mod)
CoordRepTokenizer = _mod.CoordRepTokenizer
TokenizerConfig = _mod.TokenizerConfig


def test_composite_input_produces_factorized_tokens():
    """Composite [Metal:Fe|ox:+2|d:d6|CN:6] → [Fe] ;ox=+2 ;d=6 ;CN=6"""
    tok = CoordRepTokenizer()
    block = "[Metal:Fe|ox:+2|d:d6|CN:6]"
    result = tok._tokenize_metal(block)
    assert result[0] == "[Fe]", f"Expected [Fe], got {result[0]}"
    assert ";ox=+2" in result, f"Missing ;ox=+2 in {result}"
    assert ";d=6" in result, f"Missing ;d=6 in {result}"
    assert ";CN=6" in result, f"Missing ;CN=6 in {result}"
    # Must NOT contain the composite token
    assert "[Metal:Fe|ox:+2|d:d6|CN:6]" not in result, "Composite token should not appear"
    print(f"  PASS: composite → {result}")


def test_semicolon_input_produces_factorized_tokens():
    """Semicolon [Fe;ox=+2;d=6;CN=6] → [Fe] ;ox=+2 ;d=6 ;CN=6"""
    tok = CoordRepTokenizer()
    block = "[Fe;ox=+2;d=6;CN=6]"
    result = tok._tokenize_metal(block)
    assert result[0] == "[Fe]", f"Expected [Fe], got {result[0]}"
    assert ";ox=+2" in result, f"Missing ;ox=+2 in {result}"
    assert ";d=6" in result, f"Missing ;d=6 in {result}"
    assert ";CN=6" in result, f"Missing ;CN=6 in {result}"
    print(f"  PASS: semicolon → {result}")


def test_composite_minimal():
    """Composite with CN only: [Metal:Ir|CN:6] → [Ir] ;CN=6"""
    tok = CoordRepTokenizer()
    block = "[Metal:Ir|CN:6]"
    result = tok._tokenize_metal(block)
    assert result[0] == "[Ir]", f"Expected [Ir], got {result[0]}"
    assert ";CN=6" in result, f"Missing ;CN=6 in {result}"
    assert len(result) == 2, f"Expected 2 tokens, got {len(result)}: {result}"
    print(f"  PASS: minimal composite → {result}")


def test_factorized_tokens_in_vocab():
    """All factorized tokens should be in SPECIAL_TOKENS / vocab."""
    tok = CoordRepTokenizer()
    assert "[Fe]" in tok.token2id, "[Fe] not in vocab"
    assert ";ox=+2" in tok.token2id, ";ox=+2 not in vocab"
    assert ";CN=6" in tok.token2id, ";CN=6 not in vocab"
    assert ";d=6" in tok.token2id, ";d=6 not in vocab"
    print("  PASS: factorized tokens exist in vocab")


def test_composite_tokens_not_generated():
    """Tokenizing a composite metal block should NOT produce composite token."""
    tok = CoordRepTokenizer()
    block = "[Metal:Pd|ox:+2|d:d8|CN:4]"
    result = tok._tokenize_metal(block)
    for r in result:
        assert not r.startswith("[Metal:"), f"Composite token found: {r}"
    print(f"  PASS: no composite tokens in output: {result}")


def test_full_coordrep_string():
    """Full CoordRep string should have factorized metal tokens."""
    tok = CoordRepTokenizer()
    s = "[Metal:Fe|ox:+2|d:d6|CN:6]<Shape:Oh=1.23>{trans:L1:N:1--L2:Cl:1}|L1=NCC|L2=Cl|"
    tokens = tok.tokenize(s)
    assert "[Fe]" in tokens, f"[Fe] not in tokens: {tokens[:10]}"
    assert ";CN=6" in tokens, f";CN=6 not in tokens: {tokens[:10]}"
    assert not any(t.startswith("[Metal:") for t in tokens), \
        f"Composite token found in full tokenization: {tokens[:10]}"
    print(f"  PASS: full string → first 8 tokens: {tokens[:8]}")


def test_metal_and_cn_are_independent_tokens():
    """Metal and CN must be distinct, independently maskable tokens."""
    tok = CoordRepTokenizer()
    block = "[Metal:Fe|ox:+2|d:d6|CN:6]"
    result = tok._tokenize_metal(block)
    metal_tokens = [t for t in result if t.startswith('[') and t.endswith(']')]
    cn_tokens = [t for t in result if t.startswith(';CN=')]
    assert len(metal_tokens) == 1, f"Expected 1 metal token, got {metal_tokens}"
    assert len(cn_tokens) == 1, f"Expected 1 CN token, got {cn_tokens}"
    assert metal_tokens[0] != cn_tokens[0], "Metal and CN should be different tokens"
    print(f"  PASS: metal={metal_tokens[0]}, CN={cn_tokens[0]} (independent)")


if __name__ == "__main__":
    print("Running factorized tokenizer tests …\n")
    test_composite_input_produces_factorized_tokens()
    test_semicolon_input_produces_factorized_tokens()
    test_composite_minimal()
    test_factorized_tokens_in_vocab()
    test_composite_tokens_not_generated()
    test_full_coordrep_string()
    test_metal_and_cn_are_independent_tokens()
    print("\nAll tests passed ✓")
