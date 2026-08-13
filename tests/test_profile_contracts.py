"""Scientific profile contracts shared by DockLens and the PyMOL plug-in."""

from __future__ import annotations

import pytest

from docklens import interaction_core as core


def test_supported_scientific_profiles_have_immutable_snapshots():
    assert tuple(core.HBOND_PRESETS) == ("plip", "luna", "dsv", "luna_dsv")

    for profile in core.HBOND_PRESETS:
        snapshot = core.cutoffs_for_preset(profile)
        with pytest.raises(TypeError):
            snapshot["hbond_dist"] = 999


def test_dsv_profile_uses_runtime_discovery_studio_2024_defaults():
    dsv = core.cutoffs_for_preset("dsv")

    assert dsv["hbond_dist"] == pytest.approx(3.4)
    assert dsv["hbond_angle"] == pytest.approx(90)
    assert dsv["carbon_hbond_dist"] == pytest.approx(3.8)
    assert dsv["saltbridge_dist"] == pytest.approx(4.0)
    assert dsv["charge_dist"] == pytest.approx(5.6)
    assert dsv["pication_dist"] == pytest.approx(5.0)
    assert dsv["pication_angle_max"] == pytest.approx(40)
    assert dsv["pipi_dist"] == pytest.approx(6.0)
    assert dsv["pipi_closest_atom_dist"] == pytest.approx(4.5)
    assert dsv["pipi_stacked_theta_max"] == pytest.approx(50)
    assert dsv["pipi_stacked_gamma_max"] == pytest.approx(35)
    assert dsv["pipi_t_theta_deviation_max"] == pytest.approx(30)
    assert dsv["pipi_t_gamma_min"] == pytest.approx(55)
    assert dsv["halogen_f_dist"] == pytest.approx(3.7)
    assert dsv["halogen_vdw_fraction"] == pytest.approx(1.0)
    assert dsv["halogen_donor_angle_min"] == pytest.approx(120)
    assert dsv["halogen_acceptor_angle_min"] == pytest.approx(75)
    assert dsv["alkyl_dist"] == pytest.approx(5.5)
    assert dsv["pialkyl_dist"] == pytest.approx(5.5)
    assert dsv["pi_lone_pair_dist"] == pytest.approx(3.0)
    assert dsv["pi_lone_pair_angle"] == pytest.approx(45)


def test_luna_profile_uses_native_luna_defaults():
    luna = core.cutoffs_for_preset("luna")

    assert luna["hbond_dist"] == pytest.approx(3.9)
    assert luna["hbond_h_a_dist"] == pytest.approx(2.8)
    assert luna["hbond_angle"] == pytest.approx(90)
    assert luna["carbon_hbond_dist"] == pytest.approx(4.0)
    assert luna["carbon_hbond_h_a_dist"] == pytest.approx(3.0)
    assert luna["carbon_hbond_angle"] == pytest.approx(110)
    assert luna["saltbridge_dist"] == pytest.approx(6.0)
    assert luna["pipi_dist"] == pytest.approx(6.0)
    assert luna["pipi_stacked_theta_max"] == pytest.approx(30)
    assert luna["pipi_t_theta_deviation_max"] == pytest.approx(30)
    assert luna["alkyl_dist"] == pytest.approx(4.5)
    assert luna["pication_dist"] == pytest.approx(6.0)
    assert luna["halogen_dist"] == pytest.approx(4.0)
    assert luna["halogen_donor_angle_min"] == pytest.approx(120)
    assert luna["halogen_acceptor_angle_min"] == pytest.approx(80)
    assert luna["chalcogen_dist"] == pytest.approx(4.0)
    assert luna["chalcogen_donor_angle_min"] == pytest.approx(120)
    assert luna["chalcogen_acceptor_angle_min"] == pytest.approx(80)
    assert luna["metal_dist"] == pytest.approx(2.8)


def test_luna_dsv_profile_uses_conservative_overlap_and_union_of_families():
    hybrid = core.cutoffs_for_preset("luna_dsv")

    assert hybrid["hbond_dist"] == pytest.approx(3.4)
    assert hybrid["hbond_h_a_dist"] == pytest.approx(2.8)
    assert hybrid["carbon_hbond_dist"] == pytest.approx(3.8)
    assert hybrid["carbon_hbond_angle"] == pytest.approx(110)
    assert hybrid["saltbridge_dist"] == pytest.approx(4.0)
    assert hybrid["pication_dist"] == pytest.approx(5.0)
    assert hybrid["alkyl_dist"] == pytest.approx(4.5)
    assert hybrid["halogen_acceptor_angle_min"] == pytest.approx(80)
    assert hybrid["metal_dist"] == pytest.approx(2.8)

    families = core.default_types_for_profile("luna_dsv")
    assert "chalcogen" in families
    assert set(core.default_types_for_profile("dsv")) <= set(families)
    assert set(core.default_types_for_profile("luna")) <= set(families)
    assert "proximal" not in families
    assert "vdw" not in families
    assert "vdw_clash" not in families


def test_palette_contract_is_complete_and_stable_for_every_type():
    assert "chalcogen" in core.VALID_TYPES
    assert "attractive_charge" in core.VALID_TYPES
    assert "charge_repulsion" in core.VALID_TYPES
    assert set(core.VALID_TYPES) == set(core.INTERACTION_COLORS)
    for interaction_type in core.VALID_TYPES:
        value = core.color_hex(interaction_type)
        assert value.startswith("#")
        assert len(value) == 7
        int(value[1:], 16)


def _atom(index, element, name, coord, *, side, charge=0):
    value = core.Atom(
        index,
        element,
        name,
        "LIG" if side == "ligand" else "REC",
        "1",
        coord=coord,
        fcharge=charge,
    )
    value.side = side
    return value


def _bond(first, second):
    first.neighbors.append(second)
    second.neighbors.append(first)
    first.bond_orders[second.idx] = "1"
    second.bond_orders[first.idx] = "1"


def test_luna_requires_explicit_donor_hydrogen_and_applies_both_distances():
    donor = _atom(1, "N", "N1", (0, 0, 0), side="receptor")
    hydrogen = _atom(2, "H", "H1", (1, 0, 0), side="receptor")
    acceptor = _atom(3, "O", "O1", (3.7, 0, 0), side="ligand")
    base = _atom(4, "C", "C1", (4.7, 0, 0), side="ligand")
    _bond(donor, hydrogen)
    _bond(acceptor, base)

    assert core.compute_interactions(
        [donor, hydrogen], [acceptor, base], types=["hbond"], chemistry_profile="luna"
    )

    donor.neighbors.clear()
    hydrogen.neighbors.clear()
    assert not core.compute_interactions(
        [donor], [acceptor, base], types=["hbond"], chemistry_profile="luna"
    )


def test_hybrid_halogen_and_chalcogen_require_both_native_distance_policies():
    carbon = _atom(1, "C", "C1", (-1, 0, 0), side="receptor")
    chlorine = _atom(2, "Cl", "CL", (0, 0, 0), side="receptor")
    sulfur = _atom(3, "S", "S1", (0, 0, 0), side="receptor")
    acceptor = _atom(4, "O", "O1", (3.5, 0, 0), side="ligand")
    base = _atom(5, "C", "C2", (4.5, 0, 0), side="ligand")
    _bond(carbon, chlorine)
    _bond(carbon, sulfur)
    _bond(acceptor, base)

    assert core.compute_interactions(
        [carbon, chlorine, sulfur],
        [acceptor, base],
        types=["halogen"],
        chemistry_profile="luna",
    )
    assert not core.compute_interactions(
        [carbon, chlorine, sulfur],
        [acceptor, base],
        types=["halogen"],
        chemistry_profile="luna_dsv",
    )
    assert core.compute_interactions(
        [carbon, chlorine, sulfur],
        [acceptor, base],
        types=["chalcogen"],
        chemistry_profile="luna",
    )
    assert not core.compute_interactions(
        [carbon, chlorine, sulfur],
        [acceptor, base],
        types=["chalcogen"],
        chemistry_profile="luna_dsv",
    )


def test_charge_families_expand_coverage_without_duplicating_salt_bridges():
    cation = _atom(1, "N", "N1", (0, 0, 0), side="receptor", charge=1)
    anion = _atom(2, "O", "O1", (3.0, 0, 0), side="ligand", charge=-1)

    close = core.compute_interactions(
        [cation],
        [anion],
        types=["saltbridge", "attractive_charge"],
        chemistry_profile="luna_dsv",
    )
    assert [record["type"] for record in close] == ["saltbridge"]

    anion.coord = core._v((5.0, 0, 0))
    extended = core.compute_interactions(
        [cation],
        [anion],
        types=["saltbridge", "attractive_charge"],
        chemistry_profile="luna_dsv",
    )
    assert [record["type"] for record in extended] == ["attractive_charge"]

    second_cation = _atom(
        3, "N", "N2", (5.0, 0, 0), side="ligand", charge=1
    )
    repulsive = core.compute_interactions(
        [cation],
        [second_cation],
        types=["charge_repulsion"],
        chemistry_profile="luna_dsv",
    )
    assert len(repulsive) == 1


def test_luna_explicit_charge_families_do_not_double_count_one_ionic_pair():
    cation = _atom(1, "N", "N1", (0, 0, 0), side="receptor", charge=1)
    anion = _atom(2, "O", "O1", (3.0, 0, 0), side="ligand", charge=-1)

    records = core.compute_interactions(
        [cation],
        [anion],
        types=["saltbridge", "attractive_charge"],
        chemistry_profile="luna",
    )

    assert [record["type"] for record in records] == ["saltbridge"]


def test_dsv_hydrogen_fallback_is_side_local_and_auditable():
    donor = _atom(1, "N", "N1", (0, 0, 0), side="ligand")
    acceptor = _atom(2, "O", "O1", (3.0, 0, 0), side="receptor")
    acceptor_base = _atom(3, "C", "C1", (4.0, 0, 0), side="receptor")
    _bond(acceptor, acceptor_base)

    fallback = core.compute_interactions(
        [acceptor, acceptor_base],
        [donor],
        types=["hbond"],
        chemistry_profile="dsv",
    )
    assert [record["chemistry_basis"] for record in fallback] == [
        "inferred_hydrogen"
    ]

    explicit_hydrogen = _atom(4, "H", "H1", (10.0, 0, 0), side="receptor")
    spectator = _atom(5, "C", "C2", (11.0, 0, 0), side="receptor")
    _bond(spectator, explicit_hydrogen)

    mixed_protonation = core.compute_interactions(
        [acceptor, acceptor_base, spectator, explicit_hydrogen],
        [donor],
        types=["hbond"],
        chemistry_profile="dsv",
    )
    assert [record["chemistry_basis"] for record in mixed_protonation] == [
        "inferred_hydrogen"
    ]
