"""Integration contracts for shared-receptor and combined multi-input runs."""

from __future__ import annotations

from pathlib import Path

from docklens import batch_runner as br
from docklens.entity_resolver import resolve
from docklens.input_plan import InputJob, InputPlan
from docklens.parser_mol2 import parse_mol2


def _pdb_atom(record, serial, name, resn, chain, resi, x, y, z, element):
    return (
        f"{record:<6}{serial:>5} {name:^4} {resn:>3} {chain}{resi:>4}    "
        f"{x:>8.3f}{y:>8.3f}{z:>8.3f}  1.00  0.00          {element:>2}\n"
    )


def _write_receptor(path: Path, x=0.0):
    path.write_text(
        _pdb_atom("ATOM", 1, "O", "SER", "A", 70, x, 0, 0, "O") + "END\n",
        encoding="utf-8",
    )


def _write_pdbqt(path: Path, offsets):
    blocks = []
    for index, offset in enumerate(offsets, 1):
        blocks.append(
            "MODEL %8d\n" % index
            + _pdb_atom("HETATM", 1, "N1", "LIG", "B", 1, offset, 0, 0, "N")
            .rstrip("\n")
            + "    -0.350 N\n"
            + _pdb_atom("HETATM", 2, "C1", "LIG", "B", 1, offset + 1.2, 0, 0, "C")
            .rstrip("\n")
            + "     0.100 C\n"
            + "ENDMDL\n"
        )
    path.write_text("".join(blocks) + "END\n", encoding="utf-8")


def _write_grouped_mol2(path: Path):
    path.write_text(
        """@<TRIPOS>MOLECULE
multi_group
7 4 3 0 0
SMALL
NO_CHARGES

@<TRIPOS>ATOM
1 CA 0.000 0.000 0.000 C.3 1 ALA1 0.000
2 C1 3.000 0.000 0.000 C.3 2 LIGA 0.000
3 C2 4.400 0.000 0.000 C.3 2 LIGA 0.000
4 N1 0.000 3.000 0.000 N.3 3 LIGB 0.000
5 C3 1.400 3.000 0.000 C.3 3 LIGB 0.000
6 O1 0.000 0.000 3.000 O.2 1 ALA1 0.000
7 H1 0.000 0.960 3.000 H 1 ALA1 0.000
@<TRIPOS>BOND
1 1 6 1
2 2 3 1
3 4 5 1
4 6 7 1
@<TRIPOS>SUBSTRUCTURE
1 ALA1 1 RESIDUE 0 A ROOT
2 LIGA 2 GROUP 0 ligand|A|LIGA ****
3 LIGB 4 GROUP 0 ligand|B|LIGB ****
""",
        encoding="utf-8",
    )


def test_shared_receptor_two_multipose_files_produces_distinct_sources(tmp_path):
    receptor = tmp_path / "protein.pdb"
    dock_a = tmp_path / "dock_A.pdbqt"
    dock_b = tmp_path / "dock_B.pdbqt"
    _write_receptor(receptor)
    _write_pdbqt(dock_a, (2.9, 8.0, 9.0))
    _write_pdbqt(dock_b, (2.9, 8.0))

    result = br.run_shared_receptor(receptor, [dock_a, dock_b])

    assert len(result.summaries) == 5
    assert {summary.source_file for summary in result.summaries} == {
        "dock_A.pdbqt",
        "dock_B.pdbqt",
    }
    assert {summary.source_id for summary in result.summaries} == {
        "S000001",
        "S000002",
    }
    assert all(record.input_mode == "paired" for record in result.input_qc)
    assert {
        record.receptor_source_file for record in result.input_qc
    } == {"protein.pdb"}


def test_shared_receptor_parses_receptor_once(tmp_path, monkeypatch):
    receptor = tmp_path / "protein.pdb"
    dock_a = tmp_path / "dock_A.pdbqt"
    dock_b = tmp_path / "dock_B.pdbqt"
    _write_receptor(receptor)
    _write_pdbqt(dock_a, (2.9,))
    _write_pdbqt(dock_b, (2.9,))
    original = br.parse_file
    calls = []

    def counted(path):
        if Path(path).resolve() == receptor.resolve():
            calls.append(path)
        return original(path)

    monkeypatch.setattr(br, "parse_file", counted)

    br.run_shared_receptor(receptor, [dock_a, dock_b])

    assert len(calls) == 1


def test_shared_receptor_accepts_pdbqt_and_mol2_in_one_plan(tmp_path, fixture_path):
    receptor = tmp_path / "protein.pdb"
    dock = tmp_path / "dock.pdbqt"
    _write_receptor(receptor)
    _write_pdbqt(dock, (2.9,))
    mol2 = Path(fixture_path("minimal_ccdc_sol16.mol2"))

    result = br.run_shared_receptor(receptor, [dock, mol2])

    assert len(result.summaries) == 2
    assert {item.source_file for item in result.summaries} == {
        "dock.pdbqt",
        "minimal_ccdc_sol16.mol2",
    }
    assert not [item for item in result.input_qc if item.status == "error"]


def test_two_receptor_groups_keep_their_provenance_separate(tmp_path):
    receptor_a = tmp_path / "protein_A.pdb"
    receptor_b = tmp_path / "protein_B.pdb"
    dock_a = tmp_path / "dock_A.pdbqt"
    dock_b = tmp_path / "dock_B.pdbqt"
    _write_receptor(receptor_a, 0.0)
    _write_receptor(receptor_b, 10.0)
    _write_pdbqt(dock_a, (2.9,))
    _write_pdbqt(dock_b, (12.9,))

    plan = InputPlan(
        jobs=(
            InputJob("paired", str(dock_a), str(receptor_a), "group-a"),
            InputJob("paired", str(dock_b), str(receptor_b), "group-b"),
        )
    )
    result = br.run_plan(plan)

    assert {
        (record.source_file, record.receptor_source_file)
        for record in result.input_qc
    } == {("dock_A.pdbqt", "protein_A.pdb"), ("dock_B.pdbqt", "protein_B.pdb")}


def test_paired_plan_rejects_multi_model_receptor_without_choosing_one(tmp_path):
    receptor = tmp_path / "multi_model.pdb"
    ligand = tmp_path / "ligand.pdbqt"
    receptor.write_text(
        "MODEL        1\n"
        + _pdb_atom("ATOM", 1, "O", "SER", "A", 70, 0, 0, 0, "O")
        + "ENDMDL\n"
        + "MODEL        2\n"
        + _pdb_atom("ATOM", 1, "O", "SER", "A", 70, 1, 0, 0, "O")
        + "ENDMDL\n",
        encoding="utf-8",
    )
    _write_pdbqt(ligand, (2.9,))
    result = br.run_shared_receptor(receptor, [ligand])

    assert result.summaries == ()
    assert [item.code for item in result.input_qc] == [
        "receptor_multiple_models"
    ]


def test_combined_pdb_two_disconnected_ligands_yields_two_resolutions(tmp_path):
    source = tmp_path / "two_ligands.pdb"
    source.write_text(
        "".join(
            (
                _pdb_atom("ATOM", 1, "CA", "ALA", "A", 1, 0, 0, 0, "C"),
                _pdb_atom("HETATM", 2, "C1", "LGA", "B", 1, 3, 0, 0, "C"),
                _pdb_atom("HETATM", 3, "C2", "LGA", "B", 1, 4.4, 0, 0, "C"),
                _pdb_atom("HETATM", 4, "N1", "LGB", "C", 1, 0, 3, 0, "N"),
                _pdb_atom("HETATM", 5, "C3", "LGB", "C", 1, 1.4, 3, 0, "C"),
                "END\n",
            )
        ),
        encoding="utf-8",
    )

    result = br.run([source])

    assert len(result.summaries) == 2
    assert {summary.ligand_id for summary in result.summaries} == {"LGA1B", "LGB1C"}
    assert all(summary.resolution_method == "hetatm-components" for summary in result.summaries)


def test_combined_mol2_groups_never_enter_each_other_receptor(tmp_path):
    source = tmp_path / "grouped.mol2"
    _write_grouped_mol2(source)

    poses = parse_mol2(str(source))
    resolutions = resolve(poses[0])

    assert len(resolutions) == 2
    assert all(len(item.receptor_atoms) == 3 for item in resolutions)
    assert all(
        {atom.resn for atom in item.receptor_atoms} == {"ALA"}
        for item in resolutions
    )


def test_invalid_combined_file_does_not_cancel_valid_sources(tmp_path):
    valid = tmp_path / "valid.pdb"
    invalid = tmp_path / "invalid.pdb"
    _write_receptor(valid)
    invalid.write_text(
        "ATOM      1  CA  ALA A   1      bad     0.000   0.000  1.00  0.00           C\n",
        encoding="utf-8",
    )

    result = br.run([valid, invalid])

    assert result.summaries == ()
    assert any(item.code == "parse_error" for item in result.input_qc)
    assert any(item.source_file == "valid.pdb" for item in result.input_qc)


def test_pdbqt_parser_preserves_atom_record_classification(tmp_path):
    source = tmp_path / "combined.pdbqt"
    source.write_text(
        "MODEL 1\n"
        + _pdb_atom("ATOM", 1, "CA", "ALA", "A", 1, 0, 0, 0, "C")
        .rstrip("\n")
        + "     0.000 C\n"
        + _pdb_atom("HETATM", 2, "C1", "LIG", "B", 1, 3, 0, 0, "C")
        .rstrip("\n")
        + "     0.000 C\nENDMDL\n",
        encoding="utf-8",
    )

    pose = br.parse_file(str(source))[0]

    assert pose.is_hetatm == {1: False, 2: True}


def test_combined_pdbqt_uses_atom_record_classification(tmp_path):
    source = tmp_path / "combined_classified.pdbqt"
    source.write_text(
        _pdb_atom("ATOM", 1, "OE1", "GLU", "A", 166, 0, 0, 0, "O")
        .rstrip("\n")
        + "     0.000 OA\n"
        + _pdb_atom("HETATM", 2, "N1", "LIG", "B", 1, 2.9, 0, 0, "N")
        .rstrip("\n")
        + "     0.000 N\nEND\n",
        encoding="utf-8",
    )

    result = br.run([source])

    assert len(result.summaries) == 1
    assert result.summaries[0].resolution_method == "hetatm"
    assert any(
        item.input_mode == "combined" and item.ligand_atoms == 1
        for item in result.input_qc
    )
