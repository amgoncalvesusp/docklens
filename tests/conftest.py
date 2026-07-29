from __future__ import annotations

import os
from pathlib import Path

import pytest


# GUI tests must not initialize the native Windows platform integration inside
# an automated, non-interactive session. Set this before pytest-qt creates its
# QApplication so teardown follows the same headless path as CI/self-check.
os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

FIXTURE_DIR = Path(__file__).parent / "fixtures"


@pytest.fixture
def fixture_path():
    """Return absolute paths rooted in the repository-owned fixture folder."""

    def resolve(name: str) -> str:
        path = (FIXTURE_DIR / name).resolve()
        assert path.is_file(), f"Missing test fixture: {path}"
        return str(path)

    return resolve


@pytest.fixture
def multi_source_result():
    """Three observations across two uploaded ligand files."""
    from docklens.results import Detail, Endpoint, InputQC, Summary, make_result

    ligand = Endpoint("ligand", "atom", "C1", (1,), "LIG", "1")
    receptor = Endpoint("receptor", "atom", "OE1", (2,), "GLU", "166")
    summaries = (
        Summary(
            ligand_id="LIG-A",
            source_file="ligand_a.mol2",
            sol=1,
            pose=1,
            docking_score=-8.0,
            n_total_interactions=1,
            n_key_residue_interactions=1,
            counts={"hbond": 1},
            source_id="S000001",
            pose_id="S000001:P0001:R001",
            source_path="C:/data/ligand_a.mol2",
        ),
        Summary(
            ligand_id="LIG-A",
            source_file="ligand_a.mol2",
            sol=2,
            pose=2,
            docking_score=-7.5,
            n_total_interactions=0,
            n_key_residue_interactions=0,
            counts={},
            source_id="S000001",
            pose_id="S000001:P0002:R001",
            source_path="C:/data/ligand_a.mol2",
        ),
        Summary(
            ligand_id="LIG-B",
            source_file="ligand_b.mol2",
            sol=1,
            pose=1,
            docking_score=-9.0,
            n_total_interactions=1,
            n_key_residue_interactions=1,
            counts={"hbond": 1},
            source_id="S000002",
            pose_id="S000002:P0001:R001",
            source_path="C:/data/ligand_b.mol2",
        ),
    )
    details = tuple(
        Detail(
            ligand_id=ligand_id,
            source_file=source_file,
            interaction_type="hbond",
            subtype="Conventional Hydrogen Bond",
            ligand=ligand,
            receptor=receptor,
            distance_A=distance,
            source_id=source_id,
            pose_id=pose_id,
            interaction_id=interaction_id,
            pose=1,
            source_path=source_path,
            is_key_residue=True,
        )
        for (
            ligand_id,
            source_file,
            source_id,
            pose_id,
            interaction_id,
            source_path,
            distance,
        ) in (
            (
                "LIG-A",
                "ligand_a.mol2",
                "S000001",
                "S000001:P0001:R001",
                "I-A",
                "C:/data/ligand_a.mol2",
                2.8,
            ),
            (
                "LIG-B",
                "ligand_b.mol2",
                "S000002",
                "S000002:P0001:R001",
                "I-B",
                "C:/data/ligand_b.mol2",
                3.0,
            ),
        )
    )
    qc = tuple(
        InputQC(
            source_id=source_id,
            source_file=source_file,
            source_path=source_path,
            poses_found=pose_count,
            poses_processed=pose_count,
        )
        for source_id, source_file, source_path, pose_count in (
            ("S000001", "ligand_a.mol2", "C:/data/ligand_a.mol2", 2),
            ("S000002", "ligand_b.mol2", "C:/data/ligand_b.mol2", 1),
        )
    )
    return make_result(
        summaries=summaries,
        details=details,
        key_residues=("GLU166",),
        receptor_residues=("GLU166",),
        input_qc=qc,
    )
