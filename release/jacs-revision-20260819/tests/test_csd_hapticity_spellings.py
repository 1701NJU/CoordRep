from types import SimpleNamespace

import pytest

from coordrep_tools.csd_adapter import _has_hapticity


@pytest.mark.parametrize(
    "bond_type",
    ["Pi", "pi", "Delocalized", "delocalized", "Delocalised", "delocalised"],
)
def test_hapticity_accepts_ccdc_spelling_variants(bond_type):
    metal = SimpleNamespace(
        bonds=[SimpleNamespace(bond_type=bond_type)]
    )
    assert _has_hapticity(None, metal)


@pytest.mark.parametrize("bond_type", ["Single", "Double", "Triple", None])
def test_hapticity_does_not_reject_atom_resolved_bonds(bond_type):
    metal = SimpleNamespace(
        bonds=[SimpleNamespace(bond_type=bond_type)]
    )
    assert not _has_hapticity(None, metal)
