from coordrep.audit.metal_policy import (
    ALL_METALS,
    ALL_METAL_SET,
    D_BLOCK_METALS,
    F_BLOCK_METALS,
    METAL_POLICY_ID,
    METAL_POLICY_SHA256,
    P_BLOCK_METALS,
    S_BLOCK_METALS,
    is_target_metal_symbol,
    metal_block,
    metal_blocks,
)


def test_frozen_ccdc_360_policy_is_closed_and_unique():
    assert METAL_POLICY_ID == "ccdc-3.6.0-is_metal-frozen-20260825-v1"
    assert len(ALL_METALS) == 96
    assert len(ALL_METAL_SET) == 96
    assert len(METAL_POLICY_SHA256) == 64
    assert S_BLOCK_METALS | P_BLOCK_METALS | D_BLOCK_METALS | F_BLOCK_METALS == ALL_METAL_SET
    assert not S_BLOCK_METALS & P_BLOCK_METALS
    assert not S_BLOCK_METALS & D_BLOCK_METALS
    assert not S_BLOCK_METALS & F_BLOCK_METALS
    assert not P_BLOCK_METALS & D_BLOCK_METALS
    assert not P_BLOCK_METALS & F_BLOCK_METALS
    assert not D_BLOCK_METALS & F_BLOCK_METALS


def test_policy_promotes_s_p_d_and_f_block_metals():
    assert is_target_metal_symbol("Mg")
    assert is_target_metal_symbol("Al")
    assert is_target_metal_symbol("Fe")
    assert is_target_metal_symbol("Ce")
    assert not is_target_metal_symbol("Si")
    assert not is_target_metal_symbol("C")
    assert metal_block("Mg") == "s"
    assert metal_block("Al") == "p"
    assert metal_block("Fe") == "d"
    assert metal_block("Ce") == "f"
    assert metal_blocks(["Fe", "Al", "Fe"]) == ("p", "d")
