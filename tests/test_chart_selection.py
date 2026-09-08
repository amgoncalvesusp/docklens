from dataclasses import replace

from docklens.results import Summary, make_result
from tests.test_analytics import _result


def summary(ligand, pose, count, source="library"):
    return Summary(ligand, "library.sdf", None, pose, None, count, 0, {},
                   source_id=source, pose_id=f"{source}:{ligand}:{pose}",
                   source_molecule_name=f"{ligand} dock{pose}")


def test_top_100_counts_ligands_not_poses_or_files_and_uses_mean():
    from docklens.chart_selection import select_chart_ligands
    records = [summary(f"L{i:03}", 1, i) for i in range(102)]
    records += [summary("L000", i, 0) for i in range(2, 150)]
    result = make_result(summaries=records)
    selected = select_chart_ligands(result)
    assert selected.metadata["ligands_before_limit"] == 102
    assert selected.metadata["ligands_displayed"] == 100
    assert {s.ligand_id for s in selected.result.summaries} == {f"L{i:03}" for i in range(2, 102)}
    assert len(result.summaries) == 250


def test_identity_preserves_distinct_sources_and_groups_raw_pose_names():
    from docklens.chart_selection import select_chart_ligands
    result = make_result(summaries=[summary("A", 1, 2), summary("A", 2, 4),
                                   summary("A", 1, 5, "other")])
    selected = select_chart_ligands(result)
    assert selected.metadata["ligands_displayed"] == 2
    assert len(selected.result.summaries) == 3


def test_canonical_compound_names_are_not_guessed_from_pose_suffixes():
    from docklens.chart_selection import select_chart_ligands
    result = make_result(summaries=[summary("A_pose1", 1, 1), summary("A_pose2", 1, 1)])
    assert select_chart_ligands(result).metadata["ligands_displayed"] == 2


def test_fingerprint_ranking_deduplicates_atom_pairs_and_includes_empty_poses():
    from docklens.chart_selection import select_chart_ligands
    result = _result(["p1", "p2", "p3"], [
        ("p1", "GLU166", "hbond", "H1", "OE1", 2.8),
        ("p1", "GLU166", "hbond", "H2", "OE2", 2.9),
        ("p3", "SER70", "saltbridge", "N1", "OG", 3.8),
    ])
    records = tuple(replace(s, source_id="file", ligand_id="A" if s.pose_id != "p3" else "B")
                    for s in result.summaries)
    selected = select_chart_ligands(replace(result, summaries=records), "fingerprint_features")
    assert [r["ligand_id"] for r in selected.metadata["ligand_ranking"]] == ["B", "A"]
    assert [r["mean_per_observation"] for r in selected.metadata["ligand_ranking"]] == [1.0, .5]


def test_observation_display_includes_every_selected_ligand_and_preserves_order():
    from docklens.chart_selection import select_chart_ligands
    records = [summary("A", i, i) for i in range(1, 250)]
    records += [summary(f"L{i:03}", 1, 1) for i in range(99)]
    selected = select_chart_ligands(make_result(summaries=records))
    shown = selected.observation_result()
    assert len(shown.summaries) == 100
    assert len({s.ligand_id for s in shown.summaries}) == 100
    assert shown.summaries[0].pose == 249
    assert len(selected.result.summaries) == 348


def test_ranking_ties_are_independent_of_input_order():
    from docklens.chart_selection import select_chart_ligands
    records = [summary(f"L{i:03}", 1, 1) for i in range(101)]
    first = select_chart_ligands(make_result(summaries=records))
    second = select_chart_ligands(make_result(summaries=reversed(records)))
    assert first.metadata["ligand_ranking"] == second.metadata["ligand_ranking"]


def test_plot_builders_limit_ligands_before_fingerprint_materialization(monkeypatch):
    from docklens import plotting
    from docklens.analytics import fingerprint_matrix
    records = [summary(f"L{i:03}", 1, i) for i in range(120)]
    result = make_result(summaries=records)
    calls = []
    def matrix(selected):
        calls.append(len(selected.summaries))
        return fingerprint_matrix(selected)
    monkeypatch.setattr(plotting, "fingerprint_matrix", matrix)
    artifact = plotting.build_fingerprint_chart(result)
    assert calls == [100]
    assert artifact.metadata["chart_selection"]["ligands_before_limit"] == 120
    assert any("Top 100 of 120" in text.get_text() for text in artifact.figure.texts)


def test_keyword_chart_calls_preserve_public_api():
    from docklens.plotting import build_comparison_chart, build_retention_chart
    result = make_result(summaries=[summary("A", 1, 1)])
    difference = build_comparison_chart(system_a=result, system_b=result)
    retention = build_retention_chart(docking=result, md=result)
    assert difference.metadata["system_a_observations"] == 1
    assert retention.metadata["md_observations"] == 1


def test_barcode_feature_cap_does_not_reduce_similarity_fingerprint():
    import pandas as pd
    from docklens.plotting import build_fingerprint_chart, build_similarity_chart
    result = make_result(summaries=[summary("A", 1, 1), summary("B", 1, 1)])
    ids = [s.pose_id for s in result.summaries]
    matrix = pd.DataFrame([[True]*45+[False]*5, [True]*50], index=ids,
                          columns=pd.MultiIndex.from_tuples([(f"R{i}", "hbond") for i in range(50)]))
    barcode = build_fingerprint_chart(result, matrix=matrix)
    similarity = build_similarity_chart(result, matrix=matrix)
    assert barcode.metadata["features_before_limit"] == 50
    assert barcode.metadata["features_displayed"] == 40
    assert len(barcode.data) == 80
    pair = similarity.data.query("observation_a != observation_b")
    assert pair["tanimoto_similarity"].tolist() == [.9, .9]
    assert all(0 <= t.get_position()[1] <= 1 for t in barcode.figure.texts)
    assert all("\n" not in tick.get_text() for tick in barcode.figure.axes[0].get_xticklabels())


def test_heatmap_labels_every_one_of_the_100_selected_ligands():
    from docklens.plotting import build_interaction_heatmap_chart
    ids = [f"p{i}" for i in range(100)]
    result = _result(ids, [(pose, "GLU166", "hbond", "H1", "OE1", 2.8) for pose in ids])
    artifact = build_interaction_heatmap_chart(result)
    assert artifact.metadata["rows_displayed"] == 100
    assert len(artifact.figure.axes[0].get_yticklabels()) == 100
