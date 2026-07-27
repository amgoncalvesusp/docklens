from __future__ import annotations

from dataclasses import replace

import pytest

from docklens import batch_runner
from docklens.export_views import build_export_view
from docklens.residue_keys import match_key_residues, parse_key_residues
from docklens.results import (
    AnalysisParameters,
    Detail,
    Endpoint,
    ExportFilter,
    Summary,
    make_result,
    with_key_residues,
)


def _endpoint(side, resname, resseq, chain, atom):
    return Endpoint(
        side=side,
        kind="atom",
        atom_name=atom,
        atom_serials=(1,),
        resname=resname,
        resseq=resseq,
        chain=chain,
        element=atom[:1],
    )


def _result():
    detail = Detail(
        ligand_id="pep",
        source_file="pose.mol2",
        interaction_type="hbond",
        subtype="",
        ligand=_endpoint("ligand", "LIG", "1", "", "N"),
        receptor=_endpoint("receptor", "ASN", "170", "A", "OD1"),
        distance_A=2.9,
        source_id="S1",
        pose_id="P1",
        interaction_id="I1",
        pose=1,
    )
    summary = Summary(
        ligand_id="pep",
        source_file="pose.mol2",
        sol=None,
        pose=1,
        docking_score=None,
        n_total_interactions=1,
        n_key_residue_interactions=0,
        counts={"hbond": 1},
        source_id="S1",
        pose_id="P1",
    )
    return make_result(
        details=(detail,),
        summaries=(summary,),
        receptor_residues=("ASN170A", "SER70A", "SER70B"),
        parameters=AnalysisParameters(),
    )


def test_key_residues_accept_semicolon_newline_and_spaced_number():
    raw = "SER70; LYS73,\nASN132; SER130; ASN 170; GLU166;"

    parsed = parse_key_residues(raw)

    assert parsed.keys == (
        "ASN132",
        "ASN170",
        "GLU166",
        "LYS73",
        "SER130",
        "SER70",
    )
    assert parsed.invalid == ()
    assert batch_runner.normalize_key_residues(raw) == set(parsed.keys)


def test_key_residue_parser_reports_invalid_fragments():
    parsed = parse_key_residues("SER70; residue?; ASN")

    assert parsed.keys == ("SER70",)
    assert parsed.invalid == ("ASN", "RESIDUE?")


def test_matching_reports_chain_ambiguity_and_unmatched_keys():
    match = match_key_residues(
        ("SER70", "ASN170A", "GLU166"),
        ("SER70A", "SER70B", "ASN170A"),
    )

    assert match.matched_keys == ("ASN170A", "SER70")
    assert match.matched_residues == ("ASN170A", "SER70A", "SER70B")
    assert match.ambiguous_keys == ("SER70",)
    assert match.unmatched_keys == ("GLU166",)


def test_recompute_uses_the_same_canonical_key_parser():
    updated = with_key_residues(_result(), "ASN 170;")

    assert updated.key_residues == frozenset({"ASN170"})
    assert updated.details[0].is_key_residue is True
    assert updated.summaries[0].n_key_residue_interactions == 1


def test_key_only_export_rejects_when_no_configured_key_exists_in_receptor():
    result = replace(
        _result(),
        key_residues=frozenset({"GLU166"}),
        parameters=replace(_result().parameters, key_residues=("GLU166",)),
    )

    with pytest.raises(ValueError, match="No configured key residue"):
        build_export_view(
            result,
            ExportFilter(scope="filtered", key_only=True),
        )


def test_key_only_export_allows_existing_residue_without_interactions():
    result = replace(
        _result(),
        details=(),
        key_residues=frozenset({"SER70"}),
        parameters=replace(_result().parameters, key_residues=("SER70",)),
    )

    view = build_export_view(
        result,
        ExportFilter(scope="filtered", key_only=True),
    )

    assert view.details == ()
    assert view.summaries[0].n_key_residue_interactions == 0
