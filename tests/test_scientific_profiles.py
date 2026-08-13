"""Scientific profile contracts derived from the local DSV/LUNA defaults."""

from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor
import re

import pytest

from docklens import interaction_core as core


VALID_SCIENTIFIC_PROFILES = ("plip", "luna", "dsv", "luna_dsv")


# Prefer the established DockLens names while allowing a more explicit public
# spelling for criteria whose geometry was not represented before this feature.
CRITERION_KEYS = {
    "strong_da": ("hbond_dist", "hbond_da_dist"),
    "strong_ha": ("hbond_h_a_dist", "hbond_ha_dist"),
    "strong_dha_min": ("hbond_angle", "hbond_dha_angle_min"),
    "strong_har_min": ("hbond_acceptor_angle", "hbond_har_angle_min"),
    "strong_dar_min": ("hbond_dar_angle", "hbond_dar_angle_min"),
    "weak_da": ("carbon_hbond_dist", "weak_hbond_da_dist"),
    "weak_ha": ("carbon_hbond_h_a_dist", "weak_hbond_ha_dist"),
    "weak_dha_min": ("carbon_hbond_angle", "weak_hbond_dha_angle_min"),
    "ionic": ("saltbridge_dist", "ionic_dist"),
    "charge": ("charge_dist", "charge_charge_dist"),
    "cation_pi_dist": ("pication_dist", "cation_pi_dist"),
    "cation_pi_angle_max": (
        "pication_angle_max",
        "pication_angle",
        "cation_pi_angle_max",
    ),
    "pipi_center": ("pipi_dist", "pipi_center_dist"),
    "pipi_closest_atom": ("pipi_closest_atom_dist",),
    "pipi_stack_theta_max": ("pipi_stacked_theta_max",),
    "pipi_stack_gamma_max": ("pipi_stacked_gamma_max",),
    "pipi_t_theta_max": ("pipi_t_theta_max", "pipi_t_theta_deviation_max"),
    "pipi_t_gamma_min": ("pipi_t_gamma_min",),
    "pipi_slope_min": ("pipi_slope_angle_min",),
    "pipi_slope_max": ("pipi_slope_angle_max",),
    "pipi_offset_min": ("pipi_offset_angle_min",),
    "pipi_offset_max": ("pipi_offset_angle_max",),
    "halogen_f": ("halogen_f_dist", "halogen_dist"),
    "halogen_clbri_vdw": (
        "halogen_clbri_vdw_fraction",
        "halogen_vdw_fraction",
    ),
    "halogen_centroid": ("halogen_centroid_dist", "halogen_ring_dist"),
    "halogen_cxa_min": (
        "halogen_cxa_angle_min",
        "halogen_donor_angle_min",
        "halogen_angle",
    ),
    "halogen_xay_min": (
        "halogen_xay_angle_min",
        "halogen_xar_angle_min",
        "halogen_acceptor_angle_min",
    ),
    "halogen_displacement_max": ("halogen_displacement_angle_max",),
    "alkyl_center": ("alkyl_dist", "alkyl_center_dist"),
    "sulfur_pi_edge_dist": ("pi_sulfur_edge_dist", "pi_sulfur_dist"),
    "sulfur_pi_edge_angle_min": ("pi_sulfur_edge_angle_min",),
    "sulfur_pi_face_dist": ("pi_sulfur_face_dist",),
    "sulfur_pi_face_angle_max": ("pi_sulfur_face_angle_max",),
    "lone_pair_pi_dist": ("pi_lone_pair_dist",),
    "lone_pair_pi_angle_max": ("pi_lone_pair_angle",),
    "chalcogen_acceptor": ("chalcogen_dist", "chalcogen_acceptor_dist"),
    "chalcogen_centroid": ("chalcogen_centroid_dist", "chalcogen_ring_dist"),
    "chalcogen_rya_min": (
        "chalcogen_rya_angle_min",
        "chalcogen_donor_angle_min",
    ),
    "chalcogen_yan_min": (
        "chalcogen_yan_angle_min",
        "chalcogen_acceptor_angle_min",
    ),
    "chalcogen_displacement_max": ("chalcogen_displacement_angle_max",),
    "metal": ("metal_dist",),
}


def _project_criteria(snapshot, expected):
    projected = {}
    for semantic_name in expected:
        aliases = CRITERION_KEYS[semantic_name]
        projected[semantic_name] = next(
            (snapshot[name] for name in aliases if name in snapshot),
            "<missing>",
        )
    return projected


@pytest.mark.parametrize("profile", VALID_SCIENTIFIC_PROFILES)
def test_scientific_profiles_are_accepted_by_detection(profile):
    assert core.compute_interactions(
        [], [], types=["hbond"], chemistry_profile=profile
    ) == []


def test_unknown_scientific_profile_is_rejected_at_snapshot_boundary():
    with pytest.raises(ValueError, match="[Uu]nknown.*profile|[Uu]nsupported.*profile"):
        core.cutoffs_for_preset("not-a-scientific-profile")


def test_dsv_profile_matches_discovery_studio_2024_defaults():
    expected = {
        "strong_da": 3.4,
        "weak_da": 3.8,
        "ionic": 4.0,
        "cation_pi_dist": 5.0,
        "cation_pi_angle_max": 40.0,
        "pipi_center": 6.0,
        "pipi_closest_atom": 4.5,
        "pipi_stack_theta_max": 50.0,
        "pipi_stack_gamma_max": 35.0,
        "pipi_t_theta_max": 30.0,
        "pipi_t_gamma_min": 55.0,
        "halogen_f": 3.7,
        "halogen_clbri_vdw": 1.0,
        "halogen_cxa_min": 120.0,
        "halogen_xay_min": 75.0,
        "alkyl_center": 5.5,
        "sulfur_pi_edge_dist": 6.0,
        "sulfur_pi_edge_angle_min": 70.0,
        "sulfur_pi_face_dist": 4.5,
        "sulfur_pi_face_angle_max": 25.0,
        "lone_pair_pi_dist": 3.0,
        "lone_pair_pi_angle_max": 45.0,
    }
    snapshot = core.cutoffs_for_preset("dsv")

    assert _project_criteria(snapshot, expected) == expected
    if any(name in snapshot for name in CRITERION_KEYS["charge"]):
        assert _project_criteria(snapshot, {"charge": 5.6}) == {"charge": 5.6}


def test_luna_profile_matches_local_luna_defaults_tsv():
    expected = {
        "strong_da": 3.9,
        "strong_ha": 2.8,
        "strong_dha_min": 90.0,
        "strong_har_min": 90.0,
        "strong_dar_min": 90.0,
        "weak_da": 4.0,
        "weak_ha": 3.0,
        "weak_dha_min": 110.0,
        "ionic": 6.0,
        "pipi_center": 6.0,
        "pipi_slope_min": 30.0,
        "pipi_slope_max": 60.0,
        "pipi_offset_min": 30.0,
        "pipi_offset_max": 60.0,
        "alkyl_center": 4.5,
        "cation_pi_dist": 6.0,
        "halogen_f": 4.0,
        "halogen_centroid": 4.5,
        "halogen_cxa_min": 120.0,
        "halogen_xay_min": 80.0,
        "halogen_displacement_max": 60.0,
        "chalcogen_acceptor": 4.0,
        "chalcogen_centroid": 4.5,
        "chalcogen_rya_min": 120.0,
        "chalcogen_yan_min": 80.0,
        "chalcogen_displacement_max": 60.0,
        "metal": 2.8,
    }
    snapshot = core.cutoffs_for_preset("luna")

    assert _project_criteria(snapshot, expected) == expected


def test_luna_dsv_is_conservative_and_enables_only_represented_families():
    conservative = core.cutoffs_for_preset("luna_dsv")
    expected_overlap = {
        # Minimum of native maximum-distance criteria.
        "strong_da": 3.4,
        "weak_da": 3.8,
        "ionic": 4.0,
        "pipi_center": 6.0,
        "cation_pi_dist": 5.0,
        "alkyl_center": 4.5,
        "halogen_f": 3.7,
        # Maximum of native minimum-angle criteria.
        "strong_dha_min": 90.0,
        "weak_dha_min": 110.0,
        "halogen_cxa_min": 120.0,
        "halogen_xay_min": 80.0,
    }

    assert _project_criteria(conservative, expected_overlap) == expected_overlap

    represented_union = {
        "hbond",
        "carbon_hbond",
        "saltbridge",
        "attractive_charge",
        "charge_repulsion",
        "pipi",
        "pication",
        "alkyl",
        "halogen",
        "metal",
        "pi_sulfur",
        "pi_lone_pair",
        "chalcogen",
    }
    assert represented_union <= set(core.VALID_TYPES)
    assert {"proximal", "vdw", "van_der_waals", "vdw_clash"}.isdisjoint(
        core.VALID_TYPES
    )


def test_profile_cutoff_snapshots_are_immutable_and_independent():
    dsv = core.cutoffs_for_preset("dsv")
    original = dsv["hbond_dist"]

    with pytest.raises(TypeError):
        dsv["hbond_dist"] = 999.0

    assert core.cutoffs_for_preset("dsv")["hbond_dist"] == original
    assert core.cutoffs_for_preset("luna") is not dsv


def _charged_pair_has_salt_bridge(profile):
    cation = core.Atom(
        1, "N", "N1", "LIG", "1", coord=(0.0, 0.0, 0.0), fcharge=1
    )
    anion = core.Atom(
        2, "O", "O1", "ASP", "2", coord=(5.0, 0.0, 0.0), fcharge=-1
    )
    cation.side = "receptor"
    anion.side = "ligand"
    return bool(
        core.compute_interactions(
            [cation],
            [anion],
            types=["saltbridge"],
            chemistry_profile=profile,
        )
    )


def test_concurrent_profile_contexts_do_not_leak_cutoffs():
    profiles = ("dsv", "luna") * 16

    with ThreadPoolExecutor(max_workers=8) as executor:
        observed = tuple(executor.map(_charged_pair_has_salt_bridge, profiles))

    assert observed == tuple(profile == "luna" for profile in profiles)


def test_canonical_color_registry_covers_every_represented_type_stably():
    first = {kind: core.color_hex(kind) for kind in core.VALID_TYPES}
    second = {kind: core.color_hex(kind) for kind in reversed(core.VALID_TYPES)}

    assert first == second
    assert all(re.fullmatch(r"#[0-9A-F]{6}", color) for color in first.values())
    assert "chalcogen" in first
