#!/usr/bin/env python3
"""
test_identity_keys.py
=====================
Tests that extract_identity_keys produces L0/L1/L2/L3 correctly.
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import pytest


def _get_extract_fn():
    """Import identity key extraction; skip if dependencies missing."""
    try:
        from coordrep.identity.identity_keys import extract_identity_keys
        return extract_identity_keys
    except ImportError as e:
        pytest.skip(f"identity_keys not importable: {e}")


EXAMPLE_COORDREP = (
    "[Metal:Pt|ox:+2|d:d8|CN:4]"
    "<ShapeBest:SP|Class:dist|Delta:3|V:0.50,12.30>"
    "{trans:L1:Cl:1--L2:Cl:1}"
    "|L1=Cl|L2=Cl|L3=[H]N([H])[H]|L4=[H]N([H])[H]|"
)


def test_extract_produces_all_layers():
    """extract_identity_keys should return L0_StateKey, L1_ShapeID, L2_TopoID, L3_ConnID."""
    extract = _get_extract_fn()
    keys = extract(EXAMPLE_COORDREP)
    assert hasattr(keys, 'L0_StateKey'), "Missing L0_StateKey"
    assert hasattr(keys, 'L1_ShapeID'), "Missing L1_ShapeID"
    assert hasattr(keys, 'L2_TopoID'), "Missing L2_TopoID"
    assert hasattr(keys, 'L3_ConnID'), "Missing L3_ConnID"
    print(f"  L0: {keys.L0_StateKey[:60]}...")
    print(f"  L1: {keys.L1_ShapeID[:60]}...")
    print(f"  L2: {keys.L2_TopoID[:60]}...")
    print(f"  L3: {keys.L3_ConnID[:60]}...")


def test_L0_is_full_string():
    """L0 should be the full CoordRep string (geometry state)."""
    extract = _get_extract_fn()
    keys = extract(EXAMPLE_COORDREP)
    assert keys.L0_StateKey == EXAMPLE_COORDREP, "L0 should equal the input string"


def test_L3_is_coarsest():
    """L3 (ConnID) should be the coarsest: metal + CN + ligands only."""
    extract = _get_extract_fn()
    keys = extract(EXAMPLE_COORDREP)
    assert "Pt" in keys.L3_ConnID, "L3 should contain metal"
    assert "CN4" in keys.L3_ConnID or "CN=4" in keys.L3_ConnID, "L3 should contain CN"


def test_L2_contains_shape():
    """L2 (TopoID) should include metal + CN + best shape + ligands."""
    extract = _get_extract_fn()
    keys = extract(EXAMPLE_COORDREP)
    assert "Pt" in keys.L2_TopoID
    assert "SP" in keys.L2_TopoID or "Td" in keys.L2_TopoID, "L2 should contain best shape"


def test_different_cshm_gives_different_L0_same_L3():
    """Two strings differing only in CShM values should have different L0 but same L3."""
    extract = _get_extract_fn()
    s1 = "[Metal:Pt|ox:+2|d:d8|CN:4]<ShapeBest:SP|Class:dist|Delta:3|V:0.50,12.30>{trans:L1:Cl:1--L2:Cl:1}|L1=Cl|L2=Cl|L3=[H]N([H])[H]|L4=[H]N([H])[H]|"
    s2 = "[Metal:Pt|ox:+2|d:d8|CN:4]<ShapeBest:SP|Class:dist|Delta:5|V:1.20,11.00>{trans:L1:Cl:1--L2:Cl:1}|L1=Cl|L2=Cl|L3=[H]N([H])[H]|L4=[H]N([H])[H]|"
    k1 = extract(s1)
    k2 = extract(s2)
    assert k1.L0_StateKey != k2.L0_StateKey, "L0 should differ when CShM differs"
    assert k1.L3_ConnID == k2.L3_ConnID, "L3 should be same when only CShM differs"


def test_metal_field():
    """Extracted metal should be correct."""
    extract = _get_extract_fn()
    keys = extract(EXAMPLE_COORDREP)
    assert keys.metal == "Pt", f"Expected Pt, got {keys.metal}"


if __name__ == "__main__":
    print("Running identity key tests …\n")
    test_extract_produces_all_layers()
    test_L0_is_full_string()
    test_L3_is_coarsest()
    test_L2_contains_shape()
    test_different_cshm_gives_different_L0_same_L3()
    test_metal_field()
    print("\nAll tests passed ✓")
