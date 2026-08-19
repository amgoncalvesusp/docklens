"""Qt contracts for the multi-group input-plan dialog."""

from __future__ import annotations

from docklens.input_dialog import InputPlanDialog


def test_paired_dialog_appends_files_and_keeps_receptor_groups(qtbot, tmp_path):
    dialog = InputPlanDialog(paired_only=True)
    qtbot.addWidget(dialog)
    receptor_a = tmp_path / "protein_a.pdb"
    receptor_b = tmp_path / "protein_b.pdb"
    ligand_a = tmp_path / "dock_a.pdbqt"
    ligand_b = tmp_path / "dock_b.pdbqt"

    first = dialog._groups[0]
    first.receptor_edit.setText(str(receptor_a))
    first.add_paths([ligand_a, ligand_a])
    dialog.add_receptor_group()
    second = dialog._groups[1]
    second.receptor_edit.setText(str(receptor_b))
    second.add_paths([ligand_b])

    plan = dialog.plan()

    assert [(job.group_id, job.receptor_source_file) for job in plan.jobs] == [
        ("receptor_1", "protein_a.pdb"),
        ("receptor_2", "protein_b.pdb"),
    ]
    assert [job.source_path for job in plan.jobs] == [
        str(ligand_a.resolve()),
        str(ligand_b.resolve()),
    ]
