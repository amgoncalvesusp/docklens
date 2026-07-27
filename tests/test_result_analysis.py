from dataclasses import FrozenInstanceError

import pytest

from docklens.result_analysis import (
    analyze_key_residues,
    build_key_residue_matrix,
    detail_key_residues,
)
from docklens.results import AnalysisParameters, Detail, Endpoint, Summary, make_result


def endpoint(
    *,
    side: str,
    resname: str,
    resseq: str,
    chain: str = "",
    atom_name: str = "X",
) -> Endpoint:
    return Endpoint(
        side=side,
        kind="atom",
        atom_name=atom_name,
        atom_serials=(1,),
        resname=resname,
        resseq=resseq,
        chain=chain,
    )


def detail(
    *,
    pose_id: str,
    residue: tuple[str, str, str],
    interaction_type: str,
    interaction_id: str,
) -> Detail:
    resname, resseq, chain = residue
    return Detail(
        ligand_id="PEP",
        source_file="poses.mol2",
        interaction_type=interaction_type,
        subtype="",
        ligand=endpoint(side="ligand", resname="LIG", resseq="1"),
        receptor=endpoint(
            side="receptor",
            resname=resname,
            resseq=resseq,
            chain=chain,
        ),
        distance_A=3.0,
        source_id="source",
        pose_id=pose_id,
        interaction_id=interaction_id,
        pose=int(pose_id[-1]),
    )


def summary(pose_id: str) -> Summary:
    return Summary(
        ligand_id="PEP",
        source_file="poses.mol2",
        sol=None,
        pose=int(pose_id[-1]),
        docking_score=None,
        n_total_interactions=0,
        n_key_residue_interactions=0,
        counts={},
        pose_id=pose_id,
    )


def test_detail_key_residues_matches_chain_specific_and_chainless_keys():
    item = detail(
        pose_id="pose-1",
        residue=("SER", "70", "A"),
        interaction_type="hbond",
        interaction_id="i1",
    )

    assert detail_key_residues(item, ("LYS73", "SER70A")) == ("SER70A",)
    assert detail_key_residues(item, ("SER70",)) == ("SER70",)
    assert detail_key_residues(item, ("SER70B",)) == ()


def test_analysis_reports_auditable_metrics_for_every_pose():
    details = (
        detail(
            pose_id="pose-1",
            residue=("SER", "70", "A"),
            interaction_type="hbond",
            interaction_id="i1",
        ),
        detail(
            pose_id="pose-1",
            residue=("SER", "70", "A"),
            interaction_type="saltbridge",
            interaction_id="i2",
        ),
        detail(
            pose_id="pose-1",
            residue=("LYS", "73", "A"),
            interaction_type="carbon_hbond",
            interaction_id="i3",
        ),
        detail(
            pose_id="pose-1",
            residue=("ASN", "132", "A"),
            interaction_type="hbond",
            interaction_id="i4",
        ),
        detail(
            pose_id="pose-1",
            residue=("GLY", "71", "A"),
            interaction_type="hbond",
            interaction_id="i5",
        ),
    )
    result = make_result(
        details=details,
        summaries=(summary("pose-1"), summary("pose-2")),
        key_residues=("SER70", "LYS73A", "GLU166"),
    )

    pose_1, pose_2 = analyze_key_residues(result)

    assert pose_1.pose_id == "pose-1"
    assert pose_1.raw_key_pair_count == 3
    assert pose_1.distinct_key_residues == ("LYS73A", "SER70")
    assert pose_1.distinct_key_residue_count == 2
    assert pose_1.configured_key_count == 3
    assert pose_1.coverage == pytest.approx(2 / 3)
    assert pose_1.conventional_hbond_residues == ("SER70",)
    assert pose_1.conventional_hbond_residue_count == 1
    assert pose_1.interaction_types == ("carbon_hbond", "hbond", "saltbridge")
    assert pose_1.interaction_type_diversity == 3

    assert pose_2.pose_id == "pose-2"
    assert pose_2.raw_key_pair_count == 0
    assert pose_2.distinct_key_residues == ()
    assert pose_2.configured_key_count == 3
    assert pose_2.coverage == 0.0
    assert pose_2.conventional_hbond_residues == ()
    assert pose_2.interaction_types == ()


def test_analysis_uses_parameter_keys_when_result_key_set_is_empty():
    item = detail(
        pose_id="pose-1",
        residue=("GLU", "166", "A"),
        interaction_type="hbond",
        interaction_id="i1",
    )
    result = make_result(
        details=(item,),
        summaries=(summary("pose-1"),),
        parameters=AnalysisParameters(key_residues=("GLU166",)),
    )

    metrics = analyze_key_residues(result)

    assert metrics[0].distinct_key_residues == ("GLU166",)
    assert metrics[0].coverage == 1.0


def test_analysis_normalizes_legacy_parameter_key_separators():
    item = detail(
        pose_id="pose-1",
        residue=("SER", "70", "A"),
        interaction_type="hbond",
        interaction_id="i1",
    )
    result = make_result(
        details=(item,),
        summaries=(summary("pose-1"),),
        parameters=AnalysisParameters(
            key_residues=("SER70;", "LYS 73", "GLU166,\nASN132")
        ),
    )

    metrics = analyze_key_residues(result)[0]
    matrix = build_key_residue_matrix(result)

    assert metrics.configured_key_count == 4
    assert metrics.distinct_key_residues == ("SER70",)
    assert matrix.residues == ("ASN132", "GLU166", "LYS73", "SER70")


def test_chain_specific_key_takes_precedence_over_overlapping_chainless_key():
    item = detail(
        pose_id="pose-1",
        residue=("SER", "70", "A"),
        interaction_type="hbond",
        interaction_id="i1",
    )
    result = make_result(
        details=(item,),
        summaries=(summary("pose-1"),),
        key_residues=("SER70", "SER70A"),
    )

    metrics = analyze_key_residues(result)[0]
    matrix = build_key_residue_matrix(result)

    assert detail_key_residues(item, ("SER70", "SER70A")) == ("SER70A",)
    assert metrics.distinct_key_residues == ("SER70A",)
    assert metrics.coverage == 0.5
    assert dict(matrix.rows[0].counts) == {"SER70": 0, "SER70A": 1}


def test_zero_configured_keys_has_defined_zero_coverage():
    item = detail(
        pose_id="pose-1",
        residue=("SER", "70", "A"),
        interaction_type="hbond",
        interaction_id="i1",
    )
    result = make_result(details=(item,), summaries=(summary("pose-1"),))

    metrics = analyze_key_residues(result)
    matrix = build_key_residue_matrix(result)

    assert metrics[0].raw_key_pair_count == 0
    assert metrics[0].configured_key_count == 0
    assert metrics[0].coverage == 0.0
    assert matrix.residues == ()
    assert matrix.rows[0].counts == {}


def test_key_residue_matrix_counts_pairs_and_keeps_empty_poses():
    result = make_result(
        details=(
            detail(
                pose_id="pose-1",
                residue=("SER", "70", "A"),
                interaction_type="hbond",
                interaction_id="i1",
            ),
            detail(
                pose_id="pose-1",
                residue=("SER", "70", "A"),
                interaction_type="saltbridge",
                interaction_id="i2",
            ),
            detail(
                pose_id="pose-1",
                residue=("LYS", "73", "A"),
                interaction_type="hbond",
                interaction_id="i3",
            ),
        ),
        summaries=(summary("pose-1"), summary("pose-2")),
        key_residues=("SER70", "LYS73A", "GLU166"),
    )

    matrix = build_key_residue_matrix(result)

    assert matrix.residues == ("GLU166", "LYS73A", "SER70")
    assert matrix.rows[0].pose_id == "pose-1"
    assert dict(matrix.rows[0].counts) == {
        "GLU166": 0,
        "LYS73A": 1,
        "SER70": 2,
    }
    assert matrix.rows[1].pose_id == "pose-2"
    assert dict(matrix.rows[1].counts) == {
        "GLU166": 0,
        "LYS73A": 0,
        "SER70": 0,
    }


def test_analysis_contracts_are_immutable():
    result = make_result(
        summaries=(summary("pose-1"),),
        key_residues=("SER70",),
    )
    metrics = analyze_key_residues(result)[0]
    matrix = build_key_residue_matrix(result)

    with pytest.raises(FrozenInstanceError):
        metrics.coverage = 1.0
    with pytest.raises(TypeError):
        matrix.rows[0].counts["SER70"] = 2
