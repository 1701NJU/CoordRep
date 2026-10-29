#!/usr/bin/env python3
"""
test_basic_encoding.py
======================
Tests for basic CoordRep encoding roundtrip and structure.
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import re
import pytest


def test_coordrep_structure_regex():
    """A valid CoordRep string should match the expected structural pattern."""
    example = (
        "[Metal:Fe|ox:+2|d:d6|CN:6]"
        "<Shape:Oh=1.20,TPr=8.50>"
        "{trans:L1:N:1--L4:N:1}"
        "{trans:L2:N:1--L5:N:1}"
        "{trans:L3:N:1--L6:N:1}"
        "|L1=[H]N([H])[H]|L2=[H]N([H])[H]|L3=[H]N([H])[H]|"
        "L4=[H]N([H])[H]|L5=[H]N([H])[H]|L6=[H]N([H])[H]|"
    )

    # Metal block
    assert re.search(r'\[Metal:[A-Z][a-z]?\|', example), "Missing metal block"

    # Shape block
    assert re.search(r'<Shape:[A-Za-z]+=[\d.]+', example), "Missing shape block"

    # At least one constraint
    assert re.search(r'\{(trans|cis):', example), "Missing constraint block"

    # At least one ligand
    assert re.search(r'\|L\d+=', example), "Missing ligand block"

    print("  PASS: structural regex matches")


def test_metal_extraction():
    """Metal can be extracted from the standard format."""
    example = "[Metal:Ir|ox:+3|d:d6|CN:6]<Shape:Oh=0.80>|L1=Cl|"
    m = re.search(r'\[Metal:([A-Z][a-z]?)\|', example)
    assert m is not None, "Metal pattern not found"
    assert m.group(1) == "Ir", f"Expected Ir, got {m.group(1)}"
    print(f"  PASS: metal = {m.group(1)}")


def test_cn_extraction():
    """CN can be extracted from the standard format."""
    example = "[Metal:Pt|ox:+2|d:d8|CN:4]<Shape:SP=0.50>|L1=Cl|"
    m = re.search(r'CN:(\d+)', example)
    assert m is not None, "CN pattern not found"
    assert m.group(1) == "4", f"Expected 4, got {m.group(1)}"
    print(f"  PASS: CN = {m.group(1)}")


def test_ligand_count_matches_cn():
    """Number of ligand blocks should be >= CN (may equal for monodentate)."""
    cn = 4
    example = (
        f"[Metal:Pt|ox:+2|d:d8|CN:{cn}]"
        "<Shape:SP=0.50,Td=12.30>"
        "{trans:L1:Cl:1--L2:Cl:1}"
        "|L1=Cl|L2=Cl|L3=[H]N([H])[H]|L4=[H]N([H])[H]|"
    )
    ligands = re.findall(r'\|L(\d+)=', example)
    assert len(ligands) >= cn, f"Expected >= {cn} ligands, got {len(ligands)}"
    print(f"  PASS: {len(ligands)} ligands for CN={cn}")


def test_shape_values_are_numeric():
    """CShM values in shape block should be parseable as floats."""
    example = "<Shape:Oh=1.20,TPr=8.50,Td=15.30>"
    values = re.findall(r'=(\d+\.\d+)', example)
    assert len(values) >= 2, f"Expected >= 2 shape values, got {len(values)}"
    for v in values:
        assert float(v) >= 0, f"Negative CShM value: {v}"
    print(f"  PASS: shape values = {values}")


def test_constraint_format():
    """Constraint blocks should follow {type:Lx:Elem:idx--Ly:Elem:idx}."""
    example = "{trans:L1:N:1--L2:Cl:1}"
    m = re.match(
        r'\{(trans|cis|fac|mer):(L\d+:[A-Z][a-z]?:\d+)--(L\d+:[A-Z][a-z]?:\d+)\}',
        example
    )
    assert m is not None, f"Constraint format mismatch: {example}"
    print(f"  PASS: constraint = {m.group(0)}")


def test_tokenizer_roundtrip():
    """Tokenize → encode → decode should preserve structure."""
    import importlib
    _spec = importlib.util.spec_from_file_location(
        "brain.tokenizer",
        os.path.join(os.path.dirname(__file__), "..", "brain", "tokenizer.py"),
    )
    _mod = importlib.util.module_from_spec(_spec)
    _spec.loader.exec_module(_mod)
    tok = _mod.CoordRepTokenizer()
    s = "[Metal:Fe|ox:+2|d:d6|CN:6]<Shape:Oh=V_006>{trans:L1:N:1--L2:N:1}|L1=NCC|"
    ids = tok.encode(s)
    decoded = tok.decode(ids)
    # Decoded should contain key elements
    assert "[Fe]" in decoded, f"Missing [Fe] in decoded: {decoded[:50]}"
    assert ";CN=6" in decoded, f"Missing ;CN=6 in decoded: {decoded[:50]}"
    print(f"  PASS: roundtrip preserved key tokens")


if __name__ == "__main__":
    print("Running basic encoding tests …\n")
    test_coordrep_structure_regex()
    test_metal_extraction()
    test_cn_extraction()
    test_ligand_count_matches_cn()
    test_shape_values_are_numeric()
    test_constraint_format()
    test_tokenizer_roundtrip()
    print("\nAll tests passed ✓")
