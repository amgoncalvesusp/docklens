"""Focused geometry/feature tests for the dependency-free interaction core."""

from __future__ import annotations

from types import SimpleNamespace

import numpy as np
import pytest

from docklens import interaction_core as core


def atom(index, element, name, coord, charge=0, resname="LIG"):
    value = core.Atom(
        index,
        element,
        name,
        resname,
        "1",
        "A",
        coord=coord,
        fcharge=charge,
        serial=index + 1,
    )
    value.side = "receptor"
    return value


def feature(**updates):
    values = {
        "donors": [],
        "carbon_donors": [],
        "acceptors": [],
        "cations": [],
        "anions": [],
        "halogens": [],
        "alkyl": [],
        "metals": [],
        "sulfurs": [],
        "rings": [],
    }
    values.update(updates)
    return values


def ring(center, normal=(0.0, 0.0, 1.0), tag="PHE1A_ring"):
    representative = atom(90, "C", "CG", center, resname="PHE")
    return SimpleNamespace(
        centroid=np.asarray(center, dtype=float),
        normal=np.asarray(normal, dtype=float),
        tag=tag,
        atoms=[representative],
    )


def test_vector_helpers_and_cutoff_snapshots_are_isolated():
    assert core._dist((0, 0, 0), (3, 4, 0)) == 5.0
    assert core._angle_at((0, 0, 0), (0, 0, 0), (1, 0, 0)) == 0.0
    assert core._angle_at((0, 0, 0), (1, 0, 0), (0, 1, 0)) == 90.0
    assert core._plane_angle((0, 0, 1), (0, 1, 0)) == 90.0
    assert core._proj_offset((1, 0, 2), (0, 0, 0), (0, 0, 1)) == 1.0
    assert core._planar_deviation([(0, 0, 0), (1, 0, 0), (0, 1, 0)]) == pytest.approx(0)
    assert core.cutoffs_for_preset("dsv")["hbond_dist"] == 3.5
    assert core.cutoffs_for_preset("unknown")["hbond_dist"] == 4.1


def test_ring_perception_and_labels():
    atoms = []
    for index in range(6):
        angle = index * np.pi / 3
        atoms.append(atom(index, "C", "C%d" % index, (np.cos(angle), np.sin(angle), 0)))
    for index, value in enumerate(atoms):
        value.neighbors.extend((atoms[(index - 1) % 6], atoms[(index + 1) % 6]))

    rings = core._build_rings(atoms)

    assert len(rings) == 1
    assert rings[0].tag.endswith("_ring")
    assert core.endpoint_name(rings[0]) == "ring"
    assert core.endpoint_resid(rings[0]) == "LIG1A"


def test_classify_recognizes_charges_donors_halogens_metals_and_sulfur():
    positive = atom(0, "N", "NZ", (0, 0, 0), 1, "LYS")
    negative = atom(1, "O", "OD1", (1, 0, 0), -1, "ASP")
    carbon = atom(2, "C", "C1", (2, 0, 0))
    hydrogen = atom(3, "H", "H1", (2, 1, 0))
    neighbor = atom(4, "C", "C2", (3, 0, 0))
    chlorine = atom(5, "Cl", "CL", (4, 0, 0))
    metal = atom(6, "Zn", "ZN", (5, 0, 0))
    sulfur = atom(7, "S", "S1", (6, 0, 0))
    carbon.neighbors.extend((hydrogen, neighbor))
    hydrogen.neighbors.append(carbon)
    neighbor.neighbors.extend((carbon, chlorine))
    chlorine.neighbors.append(neighbor)

    result = core.classify(
        [positive, negative, carbon, hydrogen, neighbor, chlorine, metal, sulfur],
        [],
        has_h=True,
    )

    assert result["cations"] and result["anions"]
    assert result["carbon_donors"]
    assert result["halogens"]
    assert result["metals"] == [metal]
    assert result["sulfurs"] == [sulfur]


def test_mixed_hydrogen_input_uses_one_conservative_hbond_policy():
    """One protonated side must not enable distance-only donors on the other."""
    receptor_donor = atom(20, "N", "N", (0, 0, 0))
    receptor_acceptor = atom(21, "O", "O", (0, 3, 0))
    ligand_acceptor = atom(22, "O", "O", (0, 0, 3))
    ligand_donor = atom(23, "N", "N", (0, 3, 3))
    ligand_hydrogen = atom(24, "H", "H", (0, 3, 2))
    ligand_donor.neighbors.append(ligand_hydrogen)
    ligand_hydrogen.neighbors.append(ligand_donor)
    for value in (ligand_acceptor, ligand_donor, ligand_hydrogen):
        value.side = "ligand"

    records = core.compute_interactions(
        [receptor_donor, receptor_acceptor],
        [ligand_acceptor, ligand_donor, ligand_hydrogen],
        types=["hbond"],
    )

    assert all(record["a_obj"] is not receptor_donor for record in records)


def test_pairwise_detectors_emit_roles_and_subtypes():
    cation = atom(0, "N", "N+", (0, 0, 0), 1)
    anion = atom(1, "O", "O-", (3, 0, 0), -1)
    alkyl_a = atom(2, "C", "C1", (0, 0, 0))
    alkyl_b = atom(3, "C", "C2", (3, 0, 0))
    metal = atom(4, "Zn", "ZN", (0, 0, 0))
    acceptor = atom(5, "O", "O1", (2, 0, 0))
    sulfur = atom(6, "S", "S1", (2, 0, 0))
    halogen = atom(7, "Cl", "CL", (0, 0, 0))
    bonded = atom(8, "C", "CX", (-1, 0, 0))
    r1, r2 = ring((0, 0, 0)), ring((0, 0, 4), tag="TYR2A_ring")
    pi_anion = atom(9, "O", "OPI", (0, 0, 3), -1)
    fa = feature(
        cations=[(cation.coord, cation.label(), cation)],
        alkyl=[alkyl_a],
        metals=[metal],
        halogens=[(halogen, bonded)],
        rings=[r1],
    )
    fb = feature(
        anions=[(anion.coord, anion.label(), anion)],
        alkyl=[alkyl_b],
        acceptors=[acceptor],
        sulfurs=[sulfur],
        rings=[r2],
    )

    assert core.detect_saltbridge(fa, fb)[0]["a_role"] == "cation"
    assert core.detect_alkyl(fa, fb)[0]["type"] == "alkyl"
    assert core.detect_metal(fa, fb)[0]["type"] == "metal"
    assert core.detect_halogen(fa, fb)[0]["type"] == "halogen"
    assert core.detect_pipi(fa, fb)[0]["subtype"] == "sandwich"
    assert core.detect_pication(feature(rings=[r1]), feature(cations=fa["cations"]))
    assert core.detect_pialkyl(feature(rings=[r1]), feature(alkyl=[alkyl_b]))
    assert core.detect_pi_sulfur(feature(rings=[r1]), feature(sulfurs=[sulfur]))
    assert core.detect_pi_anion(
        feature(rings=[r1]),
        feature(anions=[(pi_anion.coord, pi_anion.label(), pi_anion)]),
    )
