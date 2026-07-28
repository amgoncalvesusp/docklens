from __future__ import annotations

from docklens.analysis_profiles import build_analysis_view, detail_matches_profile
from docklens import export
from docklens.export_views import build_export_view
from docklens.results import Detail, Endpoint, ExportFilter, Summary, make_result


def _endpoint(side, atom_name="X", resname="LIG", resseq="1", chain=""):
    return Endpoint(
        side=side,
        kind="atom",
        atom_name=atom_name,
        atom_serials=(1,),
        resname=resname,
        resseq=resseq,
        chain=chain,
    )


def _detail(interaction_id, interaction_type, distance):
    return Detail(
        ligand_id="PEP",
        source_file="pose.mol2",
        interaction_type=interaction_type,
        subtype="",
        ligand=_endpoint("ligand", atom_name="C1"),
        receptor=_endpoint(
            "receptor", atom_name="NZ", resname="LYS", resseq="76", chain="A"
        ),
        distance_A=distance,
        source_id="source-1",
        pose_id="pose-1",
        interaction_id=interaction_id,
        pose=1,
    )


def _summary(counts):
    return Summary(
        ligand_id="PEP",
        source_file="pose.mol2",
        sol=None,
        pose=1,
        docking_score=-8.0,
        n_total_interactions=sum(counts.values()),
        n_key_residue_interactions=0,
        counts=counts,
        source_id="source-1",
        pose_id="pose-1",
    )


def test_ds_like_profile_keeps_dsv_hydrophobics_and_excludes_long_saltbridge():
    result = make_result(
        details=(
            _detail("i1", "hbond", 3.0),
            _detail("i2", "alkyl", 3.7),
            _detail("i3", "pialkyl", 4.6),
            _detail("i4", "saltbridge", 4.3),
            _detail("i5", "pi_lone_pair", 3.0),
        ),
        summaries=(
            _summary(
                {
                    "hbond": 1,
                    "alkyl": 1,
                    "pialkyl": 1,
                    "saltbridge": 1,
                    "pi_lone_pair": 1,
                }
            ),
        ),
    )

    view = build_analysis_view(result, "ds_like")

    assert [detail.interaction_type for detail in view.details] == [
        "hbond",
        "alkyl",
        "pialkyl",
        "pi_lone_pair",
    ]
    assert view.summaries[0].n_total_interactions == 4
    assert view.summaries[0].counts["hbond"] == 1
    assert view.summaries[0].counts["pi_lone_pair"] == 1
    assert view.summaries[0].counts["alkyl"] == 1
    assert view.summaries[0].counts["pialkyl"] == 1
    assert view.summaries[0].counts["saltbridge"] == 0


def test_ds_like_profile_keeps_short_saltbridge():
    assert detail_matches_profile(_detail("i1", "saltbridge", 3.8), "ds_like")


def test_export_parameters_record_analysis_profile():
    result = make_result()

    frame = export.parameters_dataframe(
        result,
        ExportFilter(scope="filtered", analysis_profile="ds_like"),
    )
    values = dict(frame.itertuples(index=False, name=None))

    assert values["analysis_profile"] == "ds_like"


def test_all_scope_exports_complete_result_even_when_view_is_conservative():
    result = make_result(
        details=(_detail("i1", "hbond", 3.0), _detail("i2", "alkyl", 3.7)),
        summaries=(_summary({"hbond": 1, "alkyl": 1}),),
    )

    view = build_export_view(
        result,
        ExportFilter(scope="all", analysis_profile="ds_like"),
    )

    assert view is result
