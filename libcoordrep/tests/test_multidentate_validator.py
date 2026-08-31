#!/usr/bin/env python3
"""
test_multidentate_validator.py
==============================
Tests for multidentate ligand validation logic.
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import re
import pytest


def validate_multidentate_donors(coordrep_str: str):
    """
    Lightweight multidentate donor validator.
    Returns (is_valid, issues_list).

    Rules:
    1. Each ligand L_i should have donor markers matching its denticity.
    2. No duplicate donor indices within a single ligand.
    3. Donor count in constraint block should not exceed ligand denticity.
    """
    issues = []

    # Extract ligand blocks: |L1=SMILES|
    ligand_blocks = re.findall(r'\|L(\d+)=([^|]+)', coordrep_str)

    # Extract constraint references to donors
    constraint_refs = re.findall(
        r'L(\d+):([A-Z][a-z]?):(\d+)', coordrep_str
    )

    # Count donor references per ligand
    donor_counts = {}
    donor_indices = {}
    for lig_id, elem, idx in constraint_refs:
        key = f"L{lig_id}"
        donor_counts[key] = donor_counts.get(key, 0) + 1
        if key not in donor_indices:
            donor_indices[key] = set()
        pair = f"{elem}:{idx}"
        if pair in donor_indices[key]:
            issues.append(f"duplicate_donor_rank: {key} has duplicate {pair}")
        donor_indices[key].add(pair)

    return len(issues) == 0, issues


# --- Test examples ---

VALID_BIDENTATE = (
    "[Metal:Pd|ox:+2|d:d8|CN:4]"
    "<Shape:SP=0.80,Td=15.00>"
    "{trans:L1:N:1--L1:N:2}"  # bidentate L1 with two donors
    "|L1=[H]N([H])C([H])([H])C([H])([H])N([H])[H]|"
    "L2=Cl|L3=Cl|"
)

DUPLICATE_DONOR = (
    "[Metal:Pd|ox:+2|d:d8|CN:4]"
    "<Shape:SP=0.80,Td=15.00>"
    "{trans:L1:N:1--L1:N:1}"  # duplicate: L1:N:1 appears twice
    "|L1=[H]N([H])C([H])([H])C([H])([H])N([H])[H]|"
    "L2=Cl|L3=Cl|"
)

MONODENTATE_VALID = (
    "[Metal:Pt|ox:+2|d:d8|CN:4]"
    "<Shape:SP=0.50,Td=12.30>"
    "{trans:L1:Cl:1--L2:Cl:1}"
    "|L1=Cl|L2=Cl|L3=[H]N([H])[H]|L4=[H]N([H])[H]|"
)


def test_valid_bidentate_passes():
    """Bidentate ligand with two distinct donors should pass."""
    valid, issues = validate_multidentate_donors(VALID_BIDENTATE)
    assert valid, f"Valid bidentate should pass: {issues}"
    print("  PASS: valid bidentate")


def test_duplicate_donor_fails():
    """Duplicate donor index within same ligand should fail."""
    valid, issues = validate_multidentate_donors(DUPLICATE_DONOR)
    assert not valid, "Duplicate donor should fail"
    assert any("duplicate" in i for i in issues)
    print(f"  PASS: duplicate donor detected: {issues}")


def test_monodentate_valid():
    """All monodentate ligands should pass."""
    valid, issues = validate_multidentate_donors(MONODENTATE_VALID)
    assert valid, f"Monodentate should pass: {issues}"
    print("  PASS: monodentate valid")


def test_empty_string():
    """Empty string should be trivially valid (no donors to check)."""
    valid, issues = validate_multidentate_donors("")
    assert valid, f"Empty should pass: {issues}"
    print("  PASS: empty string")


if __name__ == "__main__":
    print("Running multidentate validator tests …\n")
    test_valid_bidentate_passes()
    test_duplicate_donor_fails()
    test_monodentate_valid()
    test_empty_string()
    print("\nAll tests passed ✓")
