from __future__ import annotations

import pytest

from docklens import export
from docklens.results import AnalysisParameters, Detail, Endpoint, Summary, make_result


def _endpoint(side, resname, resseq, chain="", atom_name="X"):
    return Endpoint(
        side=side,
        kind="atom",
        atom_name=atom_name,
        atom_serials=(1,),
        resname=resname,
        resseq=resseq,
        chain=chain,
    )


def _detail(interaction_id, residue, interaction_type):
    resname, resseq, chain = residue
    return Detail(
        ligand_id="PEP",
        source_file="poses.mol2",
        interaction_type=interaction_type,
        subtype="",
        ligand=_endpoint("ligand", "LIG", "1"),
        receptor=_endpoint("receptor", resname, resseq, chain),
        distance_A=3.0,
        source_id="source",
        pose_id="pose-1",
        interaction_id=interaction_id,
        pose=1,
    )


def _summary(pose_id, pose):
    return Summary(
        ligand_id="PEP",
        source_file="poses.mol2",
        sol=None,
        pose=pose,
        docking_score=-7.0,
        n_total_interactions=0,
        n_key_residue_interactions=0,
        counts={},
        pose_id=pose_id,
    )


def test_key_residue_coverage_dataframe_separates_pairs_from_coverage():
    result = make_result(
        details=(
            _detail("i1", ("SER", "70", "A"), "hbond"),
            _detail("i2", ("SER", "70", "A"), "saltbridge"),
            _detail("i3", ("LYS", "73", "A"), "carbon_hbond"),
        ),
        summaries=(_summary("pose-1", 1), _summary("pose-2", 2)),
        key_residues=("SER70", "LYS73", "GLU166"),
    )

    frame = export.key_residue_coverage_dataframe(result)

    assert len(frame) == 2
    first = frame.iloc[0]
    assert first["raw_key_pair_count"] == 3
    assert first["distinct_key_residue_count"] == 2
    assert first["configured_key_count"] == 3
    assert first["key_residue_coverage"] == pytest.approx(2 / 3)
    assert first["conventional_hbond_residue_count"] == 1
    assert first["distinct_key_residues"] == "LYS73; SER70"
    assert first["conventional_hbond_residues"] == "SER70"
    assert frame.iloc[1]["raw_key_pair_count"] == 0


def test_parameters_report_key_matching_audit():
    result = make_result(
        details=(),
        summaries=(),
        receptor_residues=("SER70A", "SER70B", "LYS73A"),
        parameters=AnalysisParameters(key_residues=("SER70", "GLU166")),
        key_residues=("SER70", "GLU166"),
    )

    frame = export.parameters_dataframe(result)
    values = dict(frame.itertuples(index=False, name=None))

    assert values["key_residues_matched"] == "SER70"
    assert values["key_residues_unmatched"] == "GLU166"
    assert values["key_residues_ambiguous"] == "SER70"
