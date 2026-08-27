from __future__ import annotations

from docklens.source_identity import parse_source_pose_identity


def test_gold_name_preserves_verbatim_name_and_extracts_identity():
    identity = parse_source_pose_identity(
        "EOS100588|Actives_100588|mol2|1|dock48"
    )

    assert identity.molecule_name == (
        "EOS100588|Actives_100588|mol2|1|dock48"
    )
    assert identity.ligand_id == "EOS100588"
    assert identity.source_pose_label == "dock48"


def test_source_pose_label_accepts_clear_vendor_patterns():
    assert parse_source_pose_identity("EOS|pose_12").source_pose_label == "pose_12"
    assert parse_source_pose_identity("EOS|sol3").source_pose_label == "sol3"
    assert (
        parse_source_pose_identity("EOS|solution3").source_pose_label
        == "solution3"
    )


def test_unknown_pipe_name_keeps_first_nonempty_field_as_ligand():
    identity = parse_source_pose_identity("  compound-7 | vendor | batch-A  ")

    assert identity.molecule_name == "  compound-7 | vendor | batch-A  "
    assert identity.ligand_id == "compound-7"
    assert identity.source_pose_label == ""


def test_empty_source_name_has_empty_identity():
    assert parse_source_pose_identity("").ligand_id == ""
