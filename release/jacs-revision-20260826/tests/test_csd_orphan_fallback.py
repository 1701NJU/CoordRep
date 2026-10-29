from __future__ import annotations

import pytest

pytest.importorskip("ccdc", reason="licensed CCDC API is optional")

from scripts.audit_csd_orphan_geometry_fallback import _other_bond_endpoint


class _AtomWrapper:
    def __init__(self, index: int):
        self.index = index


class _FreshWrapperBond:
    @property
    def atoms(self):
        # Mimic the CCDC API: accessing bond.atoms returns wrappers that are
        # equal by molecule-local identity but are not the caller's object.
        return (_AtomWrapper(7), _AtomWrapper(11))


def test_other_endpoint_does_not_depend_on_python_object_identity():
    source_wrapper = _AtomWrapper(7)
    other = _other_bond_endpoint(_FreshWrapperBond(), source_wrapper)
    assert other is not None
    assert other.index == 11


class _MalformedBond:
    @property
    def atoms(self):
        return (_AtomWrapper(7), _AtomWrapper(7))


def test_other_endpoint_rejects_nonunique_endpoint_mapping():
    assert _other_bond_endpoint(_MalformedBond(), _AtomWrapper(7)) is None
