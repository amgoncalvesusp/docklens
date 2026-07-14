"""Security and resilience contracts for batch input processing."""

from __future__ import annotations

import pytest

from docklens import batch_runner as br


def test_missing_and_unsupported_direct_inputs_are_manifested(fixture_path, tmp_path):
    missing = tmp_path / "missing.pdb"
    unsupported = tmp_path / "notes.txt"
    unsupported.write_text("not a structure", encoding="utf-8")

    result = br.run(
        [missing, unsupported, fixture_path("minimal_complex.pdb")]
    )

    assert result.summaries
    by_code = {record.code: record for record in result.input_qc if record.code}
    assert by_code["missing_input"].source_file == "missing.pdb"
    assert by_code["unsupported_input"].source_file == "notes.txt"
    assert {record.source_id for record in result.input_qc} == {
        "S000001",
        "S000002",
        "S000003",
    }


def test_oversized_input_is_rejected_before_parser_runs(tmp_path, monkeypatch):
    source = tmp_path / "large.pdb"
    source.write_bytes(b"123456789")
    called = False

    def unexpected_parser(_path):
        nonlocal called
        called = True
        return []

    monkeypatch.setitem(br._PARSERS, ".pdb", unexpected_parser)

    result = br.run([source], max_file_size_bytes=8)

    assert not called
    assert not result.summaries
    assert result.input_qc[0].code == "file_too_large"
    assert result.input_qc[0].status == "error"


def test_input_stat_failure_is_audited(tmp_path, monkeypatch):
    source = tmp_path / "unreadable.pdb"
    source.write_text("HEADER", encoding="utf-8")

    def denied(_path):
        raise PermissionError("private operating-system detail")

    monkeypatch.setattr(br.os.path, "getsize", denied)

    result = br.run(source)

    record = result.input_qc[0]
    assert record.code == "file_stat_error"
    assert record.message == "Could not inspect input file (PermissionError)."
    assert "private" not in record.message


def test_pose_failure_is_audited_without_aborting_later_poses(
    fixture_path, monkeypatch
):
    original_resolve = br.resolve

    def fail_first_pose(pose, import_stem):
        if pose.pose_index == 0:
            raise RuntimeError("sensitive\n=attacker-controlled")
        return original_resolve(pose, import_stem=import_stem)

    monkeypatch.setattr(br, "resolve", fail_first_pose)

    result = br.run([fixture_path("two_pose_complex_sol7.pdb")])

    assert {summary.pose for summary in result.summaries} == {2}
    failure = next(record for record in result.input_qc if record.code == "pose_error")
    assert failure.pose_id.endswith(":P0001:R000")
    assert "sensitive" not in failure.message
    assert "\n" not in failure.message


def test_types_and_preset_are_normalized_and_deduplicated(fixture_path):
    result = br.run(
        [fixture_path("minimal_complex.pdb")],
        types=[" HBOND ", "hbond", "ALKYL"],
        hbond_preset=" DSV ",
    )

    assert result.parameters.interaction_types == ("hbond", "alkyl")
    assert result.parameters.hbond_preset == "dsv"


@pytest.mark.parametrize(
    ("kwargs", "message"),
    [
        ({"types": ["hbond", "made_up"]}, "Unknown interaction type"),
        ({"types": []}, "At least one interaction type"),
        ({"hbond_preset": "unknown"}, "Unknown H-bond preset"),
        ({"max_file_size_bytes": 0}, "max_file_size_bytes"),
    ],
)
def test_invalid_run_options_fail_fast(fixture_path, kwargs, message):
    with pytest.raises(ValueError, match=message):
        br.run([fixture_path("minimal_complex.pdb")], **kwargs)


def test_parse_exception_text_is_not_exposed_in_qc(tmp_path, monkeypatch):
    source = tmp_path / "broken.pdb"
    source.write_text("HEADER", encoding="utf-8")

    def hostile_parser(_path):
        raise ValueError("secret-path=C:\\private\\sample\n=FORMULA\x00")

    monkeypatch.setitem(br._PARSERS, ".pdb", hostile_parser)

    result = br.run([source])

    record = result.input_qc[0]
    assert record.code == "parse_error"
    assert record.message == "Could not parse input file (ValueError)."
    assert "private" not in record.message


def test_empty_directory_is_present_in_input_manifest(tmp_path):
    source_dir = tmp_path / "empty"
    source_dir.mkdir()

    result = br.run([source_dir])

    assert result.input_qc[0].code == "no_supported_files"
    assert result.input_qc[0].source_file == "empty"


def test_empty_input_collection_fails_fast():
    with pytest.raises(ValueError, match="at least one input"):
        br.run([])
