"""Atom numbering must not change physical features or duplicate ionic centres."""

import numpy as np
import pytest

from docklens import interaction_core as core


def two_rings(offset):
    atoms = []
    for ring in range(2):
        for i in range(6):
            angle = i * np.pi / 3
            atoms.append(
                core.Atom(
                    offset + ring * 6 + i,
                    "C",
                    "C",
                    "LIG",
                    "1",
                    coord=(ring * 5 + np.cos(angle), np.sin(angle), 0),
                    serial=ring * 6 + i + 1,
                    sybyl_type="C.ar",
                )
            )
        members = atoms[-6:]
        for i, atom in enumerate(members):
            atom.neighbors = [members[(i - 1) % 6], members[(i + 1) % 6]]
    return atoms


def test_ring_names_and_member_order_survive_pymol_index_offset():
    def signature(offset):
        return [
            (r.tag, tuple(a.serial for a in r.atoms))
            for r in core._build_rings(two_rings(offset))
        ]

    expected = signature(0)
    for offset in (1, 27, 3876, 10000):
        assert signature(offset) == expected


def test_ring_identity_survives_reserialization_and_reversed_input():
    def signature(atoms):
        return [
            (r.tag, tuple(tuple(a.coord) for a in r.atoms))
            for r in core._build_rings(atoms)
        ]

    original = two_rings(0)
    renumbered = two_rings(0)
    for atom in renumbered:
        atom.serial = 13 - atom.serial
        atom.neighbors = list(reversed(atom.neighbors))
    assert signature(list(reversed(renumbered))) == signature(original)


@pytest.mark.parametrize("profile", ["plip", "luna", "dsv", "luna_dsv"])
@pytest.mark.parametrize(
    "resn,names,charge,feature",
    [
        ("ARG", ("NH1", "NH2", "NE"), 1, "cations"),
        ("GLU", ("OE1", "OE2"), -1, "anions"),
    ],
)
def test_formal_charge_does_not_duplicate_protein_group(
    profile, resn, names, charge, feature
):
    atoms = [
        core.Atom(
            i,
            "N" if charge > 0 else "O",
            name,
            resn,
            "1",
            coord=(i, 0, 0),
            fcharge=charge if i == 0 else 0,
        )
        for i, name in enumerate(names)
    ]
    result = core.classify(atoms, [], False, chemistry_profile=profile)
    assert len(result[feature]) == 1
    assert result[feature][0][0] == pytest.approx(
        np.mean([a.coord for a in atoms], axis=0)
    )


@pytest.mark.parametrize("profile", ["luna", "dsv", "luna_dsv"])
def test_untyped_pdb_tyrosine_keeps_pi_sulfur(profile):
    atoms = two_rings(0)[:6]
    for atom, name in zip(atoms, ("CG", "CD1", "CE1", "CZ", "CE2", "CD2")):
        atom.resn, atom.name, atom.sybyl_type = "TYR", name, ""
        atom.side = "receptor"
    sulfur = core.Atom(100, "S", "S1", "LIG", "2", coord=(0, 0, 4))
    sulfur.side = "ligand"
    assert (
        len(
            core.compute_interactions(
                atoms, [sulfur], types=["pi_sulfur"], chemistry_profile=profile
            )
        )
        == 1
    )
    atoms[0].resn = "LIG"
    assert not core.compute_interactions(
        atoms, [sulfur], types=["pi_sulfur"], chemistry_profile=profile
    )


def test_chalcogen_checks_both_bond_extensions():
    sulfur = core.Atom(0, "S", "S", "LIG", "1", coord=(0, 0, 0))
    wrong = core.Atom(1, "C", "C1", "LIG", "1", coord=(0, 1, 0))
    aligned = core.Atom(2, "C", "C2", "LIG", "1", coord=(-1, 0, 0))
    acceptor = core.Atom(3, "O", "O", "ASP", "2", coord=(3, 0, 0))
    base = core.Atom(4, "C", "CG", "ASP", "2", coord=(4, 0, 0))
    sulfur.neighbors = [wrong, aligned]
    wrong.neighbors = aligned.neighbors = [sulfur]
    acceptor.neighbors, base.neighbors = [base], [acceptor]
    result = core.compute_interactions(
        [sulfur, wrong, aligned],
        [acceptor, base],
        types=["chalcogen"],
        chemistry_profile="luna",
    )
    assert len(result) == 1
    assert result[0]["donor_angle"] == pytest.approx(180)
