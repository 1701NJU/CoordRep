#!/usr/bin/env python3
"""
test_coordrep_validator.py
==========================
Tests for CoordRep grammar validation and rule-based repair.
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import pytest


def _get_validator():
    try:
        from coordrep_tools.validate import is_valid_coordrep
        return is_valid_coordrep
    except ImportError as e:
        pytest.skip(f"validate not importable: {e}")


def _get_repair():
    try:
        from scripts.run_toolb_baselines import repair_rule_only
        return repair_rule_only
    except ImportError as e:
        pytest.skip(f"repair not importable: {e}")


VALID_EXAMPLE = (
    "[Metal:Co|ox:+3|d:d6|CN:6]"
    "<Shape:Oh=1.20,TPr=8.50>"
    "{trans:L1:N:1--L2:N:1}"
    "|L1=[H]N([H])[H]|L2=[H]N([H])[H]|L3=Cl|L4=Cl|L5=Cl|L6=Cl|"
)

MISSING_BRACKET = (
    "[Metal:Co|ox:+3|d:d6|CN:6]"
    "<Shape:Oh=1.20,TPr=8.50>"
    "{trans:L1:N:1--L2:N:1}"
    "|L1=[H]N([H][H]|L2=[H]N([H])[H]|L3=Cl|L4=Cl|L5=Cl|L6=Cl|"
)

MISSING_SHAPE_CLOSE = (
    "[Metal:Co|ox:+3|d:d6|CN:6]"
    "<Shape:Oh=1.20,TPr=8.50"
    "{trans:L1:N:1--L2:N:1}"
    "|L1=[H]N([H])[H]|L2=[H]N([H])[H]|L3=Cl|L4=Cl|L5=Cl|L6=Cl|"
)


def test_valid_example_passes():
    """A correctly formed CoordRep string should pass validation."""
    is_valid = _get_validator()
    assert is_valid(VALID_EXAMPLE), "Valid example should pass"
    print("  PASS: valid example passes")


def test_missing_bracket_fails():
    """A string with a missing parenthesis should fail validation."""
    is_valid = _get_validator()
    assert not is_valid(MISSING_BRACKET), "Missing bracket should fail"
    print("  PASS: missing bracket fails")


def test_missing_shape_close_fails():
    """A string with unclosed <Shape: should fail."""
    is_valid = _get_validator()
    assert not is_valid(MISSING_SHAPE_CLOSE), "Missing > should fail"
    print("  PASS: missing shape close fails")


def test_empty_string_fails():
    """Empty string should not be valid."""
    is_valid = _get_validator()
    assert not is_valid(""), "Empty string should fail"
    print("  PASS: empty string fails")


def test_rule_repair_restores_validity():
    """Rule-based repair should restore a bracket-corrupted string to valid."""
    is_valid = _get_validator()
    repair = _get_repair()
    repaired = repair(MISSING_BRACKET)
    assert is_valid(repaired), f"Repaired string should be valid: {repaired[:80]}..."
    print(f"  PASS: rule repair restores validity (len={len(repaired)})")


def test_rule_repair_on_valid_is_noop():
    """Rule repair on an already valid string should not corrupt it."""
    is_valid = _get_validator()
    repair = _get_repair()
    repaired = repair(VALID_EXAMPLE)
    assert is_valid(repaired), "Repairing valid string should keep it valid"
    print("  PASS: repair on valid is safe")


if __name__ == "__main__":
    print("Running CoordRep validator tests …\n")
    test_valid_example_passes()
    test_missing_bracket_fails()
    test_missing_shape_close_fails()
    test_empty_string_fails()
    test_rule_repair_restores_validity()
    test_rule_repair_on_valid_is_noop()
    print("\nAll tests passed ✓")
