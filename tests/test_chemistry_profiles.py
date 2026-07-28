"""Chemistry-aware and initial Discovery Studio calibration tests."""

from __future__ import annotations

import pytest

from docklens import batch_runner, export
from docklens import interaction_core as core
from docklens.parser_mol2 import parse_mol2


def _atom(
    index,
    element,
    name,
    coord,
    *,
    sybyl_type="",
    resname="LIG",
    side="receptor",
):
    atom = core.Atom(
        index,
        element,
        name,
        resname,
        "1",
        "A",
        coord=coord,
        serial=index + 1,
        sybyl_type=sybyl_type,
    )
    atom.side = side
    return atom


def _bond(left, right, order="1"):
    left.neighbors.append(right)
    right.neighbors.append(left)
    left.bond_orders[right.idx] = order
    right.bond_orders[left.idx] = order


def test_mol2_preserves_sybyl_partial_charge_and_bond_order(tmp_path):
    source = tmp_path / "chemistry.mol2"
    source.write_text(
        """@<TRIPOS>MOLECULE
chemistry
2 1 1 0 0
SMALL
USER_CHARGES

@<TRIPOS>ATOM
1 N1 0.0 0.0 0.0 N.am 1 LIG1 -0.321
2 C1 1.3 0.0 0.0 C.2  1 LIG1  0.321
@<TRIPOS>BOND
1 1 2 am
@<TRIPOS>SUBSTRUCTURE
1 LIG1 1 GROUP
""",
        encoding="utf-8",
    )

    pose = parse_mol2(source)[0]
    nitrogen, carbon = pose.atoms

    assert nitrogen.sybyl_type == "N.am"
    assert nitrogen.partial_charge == -0.321
    assert nitrogen.bond_orders[carbon.idx] == "am"
    assert carbon.bond_orders[nitrogen.idx] == "am"


def test_dsv_excludes_amide_nitrogen_acceptor_without_changing_plip():
    nitrogen = _atom(0, "N", "N1", (0, 0, 0), sybyl_type="N.am")

    legacy = core.classify([nitrogen], [], has_h=False)
    chemistry_aware = core.classify(
        [nitrogen], [], has_h=False, chemistry_profile="dsv"
    )

    assert legacy["acceptors"] == [nitrogen]
    assert chemistry_aware["acceptors"] == []


def test_dsv_uses_sybyl_formal_charge_types_for_ligand_ionic_centres():
    quaternary_n = _atom(0, "N", "NQ", (0, 0, 0), sybyl_type="N.4")
    carbon = _atom(1, "C", "C1", (5, 0, 0), sybyl_type="C.2")
    oxygen_1 = _atom(2, "O", "O1", (4, 1, 0), sybyl_type="O.co2")
    oxygen_2 = _atom(3, "O", "O2", (4, -1, 0), sybyl_type="O.co2")
    _bond(carbon, oxygen_1)
    _bond(carbon, oxygen_2)

    legacy = core.classify(
        [quaternary_n, carbon, oxygen_1, oxygen_2],
        [],
        has_h=False,
    )
    chemistry_aware = core.classify(
        [quaternary_n, carbon, oxygen_1, oxygen_2],
        [],
        has_h=False,
        chemistry_profile="dsv",
    )

    assert legacy["cations"] == []
    assert legacy["anions"] == []
    assert len(chemistry_aware["cations"]) == 1
    assert chemistry_aware["cations"][0][2] is quaternary_n
    assert len(chemistry_aware["anions"]) == 1
    assert chemistry_aware["anions"][0][2] in {oxygen_1, oxygen_2}


def test_dsv_never_treats_carbonyl_or_carboxylate_oxygen_as_donor():
    carbonyl = _atom(0, "O", "O1", (0, 0, 0), sybyl_type="O.2")
    carboxylate = _atom(1, "O", "O2", (2, 0, 0), sybyl_type="O.co2")
    explicit_h = _atom(2, "H", "H1", (0, 1, 0))
    _bond(carbonyl, explicit_h)

    legacy = core.classify([carbonyl, carboxylate, explicit_h], [], has_h=True)
    chemistry_aware = core.classify(
        [carbonyl, carboxylate, explicit_h],
        [],
        has_h=True,
        chemistry_profile="dsv",
    )

    assert [donor for donor, _hydrogens in legacy["donors"]] == [carbonyl]
    assert chemistry_aware["donors"] == []


def test_dsv_excludes_protonated_hydroxyl_from_acceptors():
    hydroxyl = _atom(0, "O", "OG", (0, 0, 0), sybyl_type="O.3")
    hydrogen = _atom(1, "H", "HG", (1, 0, 0))
    _bond(hydroxyl, hydrogen)

    features = core.classify(
        [hydroxyl, hydrogen],
        [],
        has_h=True,
        chemistry_profile="dsv",
    )

    assert hydroxyl not in features["acceptors"]


def test_dsv_pdb_fallback_only_infers_known_protein_hydroxyl_donor():
    asparagine_oxygen = _atom(0, "O", "OD1", (0, 0, 0), resname="ASN")
    asparagine_carbon = _atom(1, "C", "CG", (1.2, 0, 0), resname="ASN")
    serine_oxygen = _atom(2, "O", "OG", (4, 0, 0), resname="SER")
    serine_carbon = _atom(3, "C", "CB", (5.2, 0, 0), resname="SER")
    ligand_oxygen = _atom(4, "O", "O1", (8, 0, 0))
    ligand_carbon = _atom(5, "C", "C1", (9.2, 0, 0))
    _bond(asparagine_oxygen, asparagine_carbon)
    _bond(serine_oxygen, serine_carbon)
    _bond(ligand_oxygen, ligand_carbon)

    chemistry_aware = core.classify(
        [
            asparagine_oxygen,
            asparagine_carbon,
            serine_oxygen,
            serine_carbon,
            ligand_oxygen,
            ligand_carbon,
        ],
        [],
        has_h=False,
        chemistry_profile="dsv",
    )

    assert chemistry_aware["donors"] == [(serine_oxygen, [])]


def test_dsv_uses_bond_order_before_inferring_a_missing_hydrogen():
    nitrile_nitrogen = _atom(0, "N", "NX", (0, 0, 0))
    nitrile_carbon = _atom(1, "C", "CX", (1.2, 0, 0))
    _bond(nitrile_nitrogen, nitrile_carbon, "3")

    chemistry_aware = core.classify(
        [nitrile_nitrogen, nitrile_carbon],
        [],
        has_h=False,
        chemistry_profile="dsv",
    )

    assert chemistry_aware["donors"] == []


def test_dsv_infers_donor_per_atom_despite_unrelated_explicit_hydrogen():
    donor = _atom(0, "O", "OH", (0, 0, 0), sybyl_type="O.3")
    donor_carbon = _atom(1, "C", "C1", (-1.3, 0, 0), sybyl_type="C.3")
    _bond(donor, donor_carbon)

    acceptor = _atom(2, "O", "O1", (3.0, 0, 0), sybyl_type="O.2", side="ligand")
    acceptor_carbon = _atom(3, "C", "C2", (4.2, 0, 0), sybyl_type="C.2", side="ligand")
    unrelated_h = _atom(4, "H", "HX", (8, 8, 8), side="ligand")
    _bond(acceptor, acceptor_carbon, "2")

    records = core.compute_interactions(
        [donor, donor_carbon],
        [acceptor, acceptor_carbon, unrelated_h],
        types=["hbond"],
        cutoffs=core.cutoffs_for_preset("dsv"),
        chemistry_profile="dsv",
    )

    assert len(records) == 1
    assert records[0]["a_obj"] is donor
    assert records[0]["chemistry_basis"] == "inferred_hydrogen"
    assert records[0]["confidence"] == "medium"


def test_dsv_explicit_hydrogen_has_high_confidence_basis():
    donor = _atom(0, "N", "N1", (0, 0, 0), sybyl_type="N.am")
    hydrogen = _atom(1, "H", "H1", (1, 0, 0))
    acceptor = _atom(2, "O", "O1", (3.0, 0, 0), sybyl_type="O.2", side="ligand")
    acceptor_base = _atom(3, "C", "C1", (3.0, 1.0, 0), sybyl_type="C.2", side="ligand")
    _bond(donor, hydrogen)
    _bond(acceptor, acceptor_base, "2")

    records = core.compute_interactions(
        [donor, hydrogen],
        [acceptor, acceptor_base],
        types=["hbond"],
        cutoffs=core.cutoffs_for_preset("dsv"),
        chemistry_profile="dsv",
    )

    assert len(records) == 1
    assert records[0]["chemistry_basis"] == "explicit_hydrogen"
    assert records[0]["confidence"] == "high"

    detail = batch_runner._detail_from_interaction(
        records[0],
        "PEP",
        "pose.mol2",
        frozenset(),
    )
    frame = export.detail_dataframe(batch_runner.make_result(details=(detail,)))

    assert detail.chemistry_basis == "explicit_hydrogen"
    assert detail.chemistry_confidence == "high"
    assert detail.hydrogen_atom == "H1"
    assert detail.hydrogen_atom_serial == 2
    assert detail.hydrogen_acceptor_distance_A == pytest.approx(2.0)
    assert detail.donor_hydrogen_acceptor_angle_deg == pytest.approx(180.0)
    assert frame.loc[0, "chemistry_basis"] == "explicit_hydrogen"
    assert frame.loc[0, "chemistry_confidence"] == "high"
    assert frame.loc[0, "hydrogen_atom"] == "H1"
    assert frame.loc[0, "hydrogen_acceptor_distance_A"] == pytest.approx(2.0)


def test_dsv_rejects_explicit_hbond_without_acceptor_base_geometry():
    donor = _atom(0, "N", "N1", (0, 0, 0), sybyl_type="N.am")
    hydrogen = _atom(1, "H", "H1", (1, 0, 0))
    isolated_acceptor = _atom(
        2, "O", "O1", (3.0, 0, 0), sybyl_type="O.2", side="ligand"
    )
    _bond(donor, hydrogen)

    records = core.compute_interactions(
        [donor, hydrogen],
        [isolated_acceptor],
        types=["hbond"],
        chemistry_profile="dsv",
    )

    assert records == []


def test_dsv_uses_explicit_hydrogen_geometry_instead_of_heavy_atom_proxy():
    donor = _atom(0, "N", "N", (0, 0, 0), sybyl_type="N.am")
    hydrogen = _atom(1, "H", "HN", (1, 0, 0))
    acceptor = _atom(2, "O", "O", (3.9, 0, 0), sybyl_type="O.2", side="ligand")
    acceptor_base = _atom(3, "C", "C", (3.9, 1, 0), sybyl_type="C.2", side="ligand")
    _bond(donor, hydrogen)
    _bond(acceptor, acceptor_base, "2")

    records = core.compute_interactions(
        [donor, hydrogen],
        [acceptor, acceptor_base],
        types=["hbond"],
        cutoffs=core.cutoffs_for_preset("dsv"),
        chemistry_profile="dsv",
    )

    assert len(records) == 1
    assert records[0]["hydrogen_obj"] is hydrogen
    assert records[0]["hydrogen_acceptor_distance"] == pytest.approx(2.9)
    assert records[0]["donor_hydrogen_acceptor_angle"] == pytest.approx(180.0)
    assert records[0]["hydrogen_acceptor_base_angle"] == pytest.approx(90.0)


def test_dsv_rejects_explicit_hydrogen_with_invalid_acceptor_geometry():
    donor = _atom(0, "N", "N", (0, 0, 0), sybyl_type="N.am")
    hydrogen = _atom(1, "H", "HN", (1, 0, 0))
    acceptor = _atom(2, "O", "O", (3.0, 0, 0), sybyl_type="O.2", side="ligand")
    acceptor_base = _atom(3, "C", "C", (2.0, 0, 0), sybyl_type="C.2", side="ligand")
    _bond(donor, hydrogen)
    _bond(acceptor, acceptor_base, "2")

    records = core.compute_interactions(
        [donor, hydrogen],
        [acceptor, acceptor_base],
        types=["hbond"],
        cutoffs=core.cutoffs_for_preset("dsv"),
        chemistry_profile="dsv",
    )

    assert records == []


def test_dsv_emits_one_auditable_record_per_qualifying_hydrogen():
    donor = _atom(0, "N", "ND2", (0, 0, 0), sybyl_type="N.am")
    hydrogen_1 = _atom(1, "H", "HD21", (1, 0.1, 0))
    hydrogen_2 = _atom(2, "H", "HD22", (1, -0.1, 0))
    acceptor = _atom(3, "O", "O", (3.5, 0, 0), sybyl_type="O.2", side="ligand")
    acceptor_base = _atom(4, "C", "C", (3.5, 0, 1), sybyl_type="C.2", side="ligand")
    _bond(donor, hydrogen_1)
    _bond(donor, hydrogen_2)
    _bond(acceptor, acceptor_base, "2")

    records = core.compute_interactions(
        [donor, hydrogen_1, hydrogen_2],
        [acceptor, acceptor_base],
        types=["hbond"],
        cutoffs=core.cutoffs_for_preset("dsv"),
        chemistry_profile="dsv",
    )

    assert {record["hydrogen_obj"].name for record in records} == {
        "HD21",
        "HD22",
    }


def test_dsv_carbon_hbond_requires_a_polarized_carbon_donor():
    unpolarized = _atom(0, "C", "CB", (0, 0, 0), sybyl_type="C.3")
    carbon_neighbor = _atom(1, "C", "CA", (-1, 0, 0), sybyl_type="C.3")
    unpolarized_h = _atom(2, "H", "HB", (1, 0, 0))
    polarized = _atom(3, "C", "CA", (0, 4, 0), sybyl_type="C.3")
    nitrogen_neighbor = _atom(4, "N", "N", (-1, 4, 0), sybyl_type="N.am")
    polarized_h = _atom(5, "H", "HA", (1, 4, 0))
    _bond(unpolarized, carbon_neighbor)
    _bond(unpolarized, unpolarized_h)
    _bond(polarized, nitrogen_neighbor)
    _bond(polarized, polarized_h)

    features = core.classify(
        [
            unpolarized,
            carbon_neighbor,
            unpolarized_h,
            polarized,
            nitrogen_neighbor,
            polarized_h,
        ],
        [],
        has_h=True,
        chemistry_profile="dsv",
    )

    assert [atom for atom, _hydrogens in features["carbon_donors"]] == [polarized]


def test_dsv_carbon_hbond_accepts_corpus_calibrated_hydrogen_distance():
    donor = _atom(0, "C", "CA", (0, 0, 0), sybyl_type="C.3")
    nitrogen_neighbor = _atom(1, "N", "N", (-1, 0, 0), sybyl_type="N.am")
    hydrogen = _atom(2, "H", "HA", (1, 0, 0))
    acceptor = _atom(3, "O", "O", (3.6, 0, 0), sybyl_type="O.2", side="ligand")
    acceptor_base = _atom(4, "C", "C", (3.6, 1, 0), sybyl_type="C.2", side="ligand")
    _bond(donor, nitrogen_neighbor)
    _bond(donor, hydrogen)
    _bond(acceptor, acceptor_base, "2")

    records = core.compute_interactions(
        [donor, nitrogen_neighbor, hydrogen],
        [acceptor, acceptor_base],
        types=["carbon_hbond"],
        cutoffs=core.cutoffs_for_preset("dsv"),
        chemistry_profile="dsv",
    )

    assert len(records) == 1


def test_dsv_does_not_infer_missing_hydrogen_on_an_explicitly_protonated_side():
    donor = _atom(0, "O", "OG", (0, 0, 0), sybyl_type="O.3")
    donor_carbon = _atom(1, "C", "CB", (-1.2, 0, 0), sybyl_type="C.3")
    unrelated_hydrogen = _atom(2, "H", "H", (8, 8, 8))
    _bond(donor, donor_carbon)

    features = core.classify(
        [donor, donor_carbon, unrelated_hydrogen],
        [],
        has_h=True,
        chemistry_profile="dsv",
    )

    assert features["donors"] == []


def test_dsv_alkyl_feature_uses_hydrophobic_protein_side_chains():
    histidine_cb = _atom(0, "C", "CB", (0, 0, 0), sybyl_type="C.3", resname="HIS")
    histidine_ca = _atom(1, "C", "CA", (1.5, 0, 0), sybyl_type="C.3", resname="HIS")
    cysteine_cb = _atom(2, "C", "CB", (4, 0, 0), sybyl_type="C.3", resname="CYS")
    cysteine_sg = _atom(3, "S", "SG", (5.5, 0, 0), sybyl_type="S.3", resname="CYS")
    _bond(histidine_cb, histidine_ca)
    _bond(cysteine_cb, cysteine_sg)

    features = core.classify(
        [histidine_cb, histidine_ca, cysteine_cb, cysteine_sg],
        [],
        has_h=False,
        chemistry_profile="dsv",
    )

    assert features["alkyl"] == [cysteine_cb]


def test_pi_lone_pair_detects_axial_acceptor_above_aromatic_ring():
    ring_atoms = [
        _atom(
            index,
            "C",
            f"C{index}",
            (
                float(__import__("math").cos(index * __import__("math").pi / 3)),
                float(__import__("math").sin(index * __import__("math").pi / 3)),
                0.0,
            ),
            sybyl_type="C.ar",
            side="ligand",
        )
        for index in range(6)
    ]
    for index, atom in enumerate(ring_atoms):
        _bond(atom, ring_atoms[(index + 1) % len(ring_atoms)], "ar")
    acceptor = _atom(10, "O", "O", (0, 0, 3.0), sybyl_type="O.2", side="receptor")

    legacy_records = core.compute_interactions(
        [acceptor],
        ring_atoms,
        types=["pi_lone_pair"],
        chemistry_profile="plip",
    )
    records = core.compute_interactions(
        [acceptor],
        ring_atoms,
        types=["pi_lone_pair"],
        chemistry_profile="dsv",
    )

    assert legacy_records == []
    assert len(records) == 1
    assert records[0]["type"] == "pi_lone_pair"
    assert records[0]["dist"] == pytest.approx(3.0)
    assert records[0]["theta"] == pytest.approx(0.0)

    detail = batch_runner._detail_from_interaction(
        records[0],
        "PEP",
        "pose.mol2",
        frozenset(),
    )
    frame = export.detail_dataframe(batch_runner.make_result(details=(detail,)))

    assert detail.theta_deg == pytest.approx(0.0)
    assert frame.loc[0, "theta_deg"] == pytest.approx(0.0)


def test_pi_lone_pair_rejects_planar_aliphatic_ring():
    ring_atoms = [
        _atom(
            index,
            "C",
            f"C{index}",
            (
                float(__import__("math").cos(index * __import__("math").pi / 3)),
                float(__import__("math").sin(index * __import__("math").pi / 3)),
                0.0,
            ),
            sybyl_type="C.3",
            side="ligand",
        )
        for index in range(6)
    ]
    for index, atom in enumerate(ring_atoms):
        _bond(atom, ring_atoms[(index + 1) % len(ring_atoms)], "1")
    acceptor = _atom(10, "O", "O", (0, 0, 3.0), sybyl_type="O.2", side="receptor")

    records = core.compute_interactions(
        [acceptor],
        ring_atoms,
        types=["pi_lone_pair"],
        chemistry_profile="dsv",
    )

    assert records == []


def test_dsv_does_not_treat_neutral_histidine_as_cation():
    neutral_nd1 = _atom(0, "N", "ND1", (0, 0, 0), sybyl_type="N.ar", resname="HIS")
    neutral_ne2 = _atom(1, "N", "NE2", (1, 0, 0), sybyl_type="N.ar", resname="HIS")
    one_hydrogen = _atom(2, "H", "HD1", (-0.5, 0, 0), resname="HIS")
    _bond(neutral_nd1, neutral_ne2, "ar")
    _bond(neutral_nd1, one_hydrogen)

    protonated_nd1 = _atom(3, "N", "ND1", (0, 3, 0), sybyl_type="N.ar", resname="HIP")
    protonated_ne2 = _atom(4, "N", "NE2", (1, 3, 0), sybyl_type="N.ar", resname="HIP")
    _bond(protonated_nd1, protonated_ne2, "ar")

    neutral = core.classify(
        [neutral_nd1, neutral_ne2, one_hydrogen],
        [],
        has_h=True,
        chemistry_profile="dsv",
    )
    protonated = core.classify(
        [protonated_nd1, protonated_ne2],
        [],
        has_h=False,
        chemistry_profile="dsv",
    )

    assert neutral["cations"] == []
    assert protonated["cations"]


def test_cutoff_values_do_not_implicitly_switch_chemistry_profile():
    donor = _atom(0, "O", "OH", (0, 0, 0), sybyl_type="O.3")
    amide_acceptor = _atom(1, "N", "N1", (3.0, 0, 0), sybyl_type="N.am", side="ligand")
    amide_neighbors = [
        _atom(
            index,
            "C",
            "C%d" % index,
            (8.0 + index, 0, 0),
            sybyl_type="C.3",
            side="ligand",
        )
        for index in range(2, 5)
    ]
    for neighbor in amide_neighbors:
        _bond(amide_acceptor, neighbor)
    dsv_geometry = dict(core.cutoffs_for_preset("dsv"))

    legacy_records = core.compute_interactions(
        [donor],
        [amide_acceptor, *amide_neighbors],
        types=["hbond"],
        cutoffs=dsv_geometry,
    )
    chemistry_records = core.compute_interactions(
        [donor],
        [amide_acceptor, *amide_neighbors],
        types=["hbond"],
        cutoffs={**dsv_geometry, "hbond_dist": 3.6},
        chemistry_profile="dsv",
    )

    assert any(record["b_obj"] is amide_acceptor for record in legacy_records)
    assert all("chemistry_basis" not in record for record in legacy_records)
    assert chemistry_records == []


def test_batch_flows_forward_the_selected_chemistry_profile(monkeypatch, fixture_path):
    profiles = []

    def capture_profile(*args, **kwargs):
        profiles.append(kwargs.get("chemistry_profile"))
        return []

    monkeypatch.setattr(batch_runner, "compute_interactions", capture_profile)

    batch_runner.run(
        [fixture_path("minimal_complex.pdb")],
        hbond_preset="dsv",
    )
    batch_runner.run_paired(
        fixture_path("minimal_complex.pdb"),
        fixture_path("two_poses_sol3.pdbqt"),
        hbond_preset="dsv",
    )

    assert profiles
    assert set(profiles) == {"dsv"}


def test_mol2_batch_export_preserves_explicit_hydrogen_geometry(tmp_path):
    source = tmp_path / "reference.mol2"
    source.write_text(
        """@<TRIPOS>MOLECULE
reference
4 2 2 0 0
BIOPOLYMER
USER_CHARGES

@<TRIPOS>ATOM
1 N  0.0 0.0 0.0 N.am 1 ASN1  0.0
2 HN 1.0 0.0 0.0 H    1 ASN1  0.0
3 O  3.9 0.0 0.0 O.2  2 LIG1 -0.4
4 C  3.9 1.0 0.0 C.2  2 LIG1  0.4
@<TRIPOS>BOND
1 1 2 1
2 3 4 2
@<TRIPOS>SUBSTRUCTURE
1 ASN1 1 RESIDUE 0 A ASN 0
2 LIG1 3 GROUP 0 B LIG 0
""",
        encoding="utf-8",
    )

    result = batch_runner.run(
        [source],
        types=["hbond"],
        hbond_preset="dsv",
    )
    frame = export.detail_dataframe(result)

    assert len(result.details) == 1
    assert frame.loc[0, "hydrogen_atom"] == "HN"
    assert frame.loc[0, "hydrogen_acceptor_distance_A"] == pytest.approx(2.9)
    assert frame.loc[0, "donor_hydrogen_acceptor_angle_deg"] == pytest.approx(180.0)
    assert frame.loc[0, "hydrogen_acceptor_base_angle_deg"] == pytest.approx(90.0)


def _aromatic_ring(side="receptor"):
    import math

    atoms = [
        _atom(
            index,
            "C",
            f"C{index}",
            (math.cos(index * math.pi / 3), math.sin(index * math.pi / 3), 0.0),
            sybyl_type="C.ar",
            resname="PHE",
            side=side,
        )
        for index in range(6)
    ]
    for index, atom in enumerate(atoms):
        _bond(atom, atoms[(index + 1) % len(atoms)], "ar")
    return atoms


def test_dsv_detects_pi_sigma_with_explicit_hydrogen_geometry():
    ring_atoms = _aromatic_ring()
    carbon = _atom(10, "C", "C1", (0, 0, 3.5), sybyl_type="C.3", side="ligand")
    hydrogen = _atom(11, "H", "H1", (0, 0, 2.5), side="ligand")
    anchor = _atom(12, "C", "C2", (1.2, 0, 3.5), sybyl_type="C.3", side="ligand")
    _bond(carbon, hydrogen)
    _bond(carbon, anchor)

    records = core.compute_interactions(
        ring_atoms,
        [carbon, hydrogen, anchor],
        types=["pi_sigma"],
        chemistry_profile="dsv",
    )

    assert len(records) == 1
    assert records[0]["type"] == "pi_sigma"
    assert records[0]["hydrogen_obj"] is hydrogen
    assert records[0]["theta"] == pytest.approx(0.0)


def test_dsv_detects_pi_donor_hydrogen_bond():
    ring_atoms = _aromatic_ring()
    donor = _atom(10, "N", "N1", (0, 0, 4.0), sybyl_type="N.am", side="ligand")
    hydrogen = _atom(11, "H", "H1", (0, 0, 3.0), side="ligand")
    _bond(donor, hydrogen)

    records = core.compute_interactions(
        ring_atoms,
        [donor, hydrogen],
        types=["pi_donor_hbond"],
        chemistry_profile="dsv",
    )

    assert len(records) == 1
    assert records[0]["type"] == "pi_donor_hbond"
    assert records[0]["hydrogen_obj"] is hydrogen
    assert records[0]["theta"] == pytest.approx(0.0)


def test_dsv_pi_alkyl_counts_one_semantic_ring_group_pair():
    ring_atoms = _aromatic_ring()
    carbon_1 = _atom(10, "C", "C1", (0, 0, 4.5), sybyl_type="C.3", side="ligand")
    carbon_2 = _atom(11, "C", "C2", (1.0, 0, 4.5), sybyl_type="C.3", side="ligand")
    _bond(carbon_1, carbon_2)

    records = core.compute_interactions(
        ring_atoms,
        [carbon_1, carbon_2],
        types=["pialkyl"],
        chemistry_profile="dsv",
    )

    assert len(records) == 1


def test_dsv_profile_uses_the_discovery_studio_geometry_observed_in_2m5d():
    cutoffs = core.cutoffs_for_preset("dsv")

    assert cutoffs["pialkyl_dist"] == pytest.approx(4.9)
    assert cutoffs["alkyl_dist"] == pytest.approx(4.2)
    assert cutoffs["metal_dist"] == pytest.approx(3.0)
    assert cutoffs["pi_sulfur_dist"] == pytest.approx(5.3)
    assert cutoffs["pi_sigma_h_centroid_dist"] == pytest.approx(4.3)
    assert cutoffs["pi_sigma_axis_angle"] == pytest.approx(40.0)
    assert cutoffs["pi_sigma_dha_angle"] == pytest.approx(160.0)
    assert cutoffs["pi_donor_dist"] == pytest.approx(5.2)
    assert cutoffs["pi_donor_h_centroid_dist"] == pytest.approx(4.1)
    assert cutoffs["pi_donor_axis_angle"] == pytest.approx(45.0)
    assert cutoffs["pi_donor_dha_angle"] == pytest.approx(145.0)


def test_dsv_pi_sigma_accepts_observed_hydrogen_centroid_distance():
    ring_atoms = _aromatic_ring()
    carbon = _atom(10, "C", "C1", (0, 0, 3.8), sybyl_type="C.3", side="ligand")
    hydrogen = _atom(11, "H", "H1", (0, 0, 2.8), side="ligand")
    anchor = _atom(12, "C", "C2", (1.2, 0, 3.8), sybyl_type="C.3", side="ligand")
    _bond(carbon, hydrogen)
    _bond(carbon, anchor)

    records = core.compute_interactions(
        ring_atoms,
        [carbon, hydrogen, anchor],
        types=["pi_sigma"],
        chemistry_profile="dsv",
    )

    assert len(records) == 1


def test_dsv_pi_donor_accepts_observed_theta_limit():
    import math

    ring_atoms = _aromatic_ring()
    theta = math.radians(44.0)
    direction = (math.sin(theta), 0.0, math.cos(theta))
    hydrogen = _atom(
        11,
        "H",
        "H1",
        tuple(3.9 * component for component in direction),
        side="ligand",
    )
    donor = _atom(
        10,
        "N",
        "N1",
        tuple(4.9 * component for component in direction),
        sybyl_type="N.am",
        side="ligand",
    )
    _bond(donor, hydrogen)

    records = core.compute_interactions(
        ring_atoms,
        [donor, hydrogen],
        types=["pi_donor_hbond"],
        chemistry_profile="dsv",
    )

    assert len(records) == 1


def test_new_dsv_interaction_types_are_exportable():
    assert "pi_sigma" in core.VALID_TYPES
    assert "pi_donor_hbond" in core.VALID_TYPES
