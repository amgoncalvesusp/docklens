"""Contracts for the version 1.2 input-plan abstraction."""

from __future__ import annotations

import pytest

from docklens.input_plan import (
    InputJob,
    InputPlan,
    plan_from_paths,
    run_shared_receptor_plan,
)


def test_shared_receptor_plan_creates_one_paired_job_per_unique_ligand(tmp_path):
    receptor = tmp_path / "protein.pdb"
    ligand_a = tmp_path / "dock_A.pdbqt"
    ligand_b = tmp_path / "dock_B.pdbqt"

    plan = run_shared_receptor_plan(receptor, [ligand_a, ligand_a, ligand_b])

    assert plan.jobs == (
        InputJob(
            kind="paired",
            receptor_path=str(receptor),
            source_path=str(ligand_a),
            group_id="receptor_1",
        ),
        InputJob(
            kind="paired",
            receptor_path=str(receptor),
            source_path=str(ligand_b),
            group_id="receptor_1",
        ),
    )


def test_plan_from_paths_keeps_combined_sources_separate(tmp_path):
    first = tmp_path / "complex1.pdb"
    second = tmp_path / "complex2.mol2"

    plan = plan_from_paths([first, second])

    assert [job.kind for job in plan.jobs] == ["combined", "combined"]
    assert [job.source_path for job in plan.jobs] == [
        str(first.resolve()),
        str(second.resolve()),
    ]
    assert [job.group_id for job in plan.jobs] == ["combined_1", "combined_2"]


def test_input_plan_rejects_invalid_job_shape(tmp_path):
    with pytest.raises(ValueError, match="kind"):
        InputJob(kind="unknown", source_path=str(tmp_path / "x.pdb"))
    with pytest.raises(ValueError, match="receptor_path"):
        InputJob(kind="paired", source_path=str(tmp_path / "x.pdb"))
    with pytest.raises(ValueError, match="combined"):
        InputJob(
            kind="combined",
            receptor_path=str(tmp_path / "protein.pdb"),
            source_path=str(tmp_path / "x.pdb"),
        )


def test_input_plan_deduplicates_only_within_the_same_group(tmp_path):
    receptor_a = tmp_path / "protein_a.pdb"
    receptor_b = tmp_path / "protein_b.pdb"
    ligand = tmp_path / "dock.pdbqt"

    plan = InputPlan(
        jobs=(
            InputJob("paired", str(ligand), str(receptor_a), "group-a"),
            InputJob("paired", str(ligand), str(receptor_a), "group-a"),
            InputJob("paired", str(ligand), str(receptor_b), "group-b"),
        )
    )

    assert len(plan.jobs) == 2
    assert {job.group_id for job in plan.jobs} == {"group-a", "group-b"}
