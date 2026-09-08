from unittest.mock import patch

from docklens.analytics_widgets import AnalyticsWorkspace
from docklens import analytics_widgets as widgets
from docklens.results import make_result
from tests.test_chart_selection import summary


def test_updates_are_batched_and_only_active_page_builds_figures(qtbot):
    workspace = AnalyticsWorkspace()
    qtbot.addWidget(workspace.residue_page)
    try:
        result = make_result(summaries=[summary("A", 1, 1)])
        with patch.object(widgets, "build_residue_chart", wraps=widgets.build_residue_chart) as residue, \
             patch.object(widgets, "build_fingerprint_chart", wraps=widgets.build_fingerprint_chart) as barcode:
            workspace.set_result(result, refresh=False)
            workspace.set_comparison(None, refresh=False)
            workspace.set_ligand_group(None, refresh=False)
            assert residue.call_count == barcode.call_count == 0
            workspace.refresh()
            assert residue.call_count == 1
            assert barcode.call_count == 0
            workspace.activate(1)
            assert barcode.call_count == 1
            workspace.activate(0)
            workspace.activate(1)
            assert residue.call_count == barcode.call_count == 1
    finally:
        workspace.dispose()


def test_ui_exports_disclose_ranked_ligands_and_do_not_trim_source_tables(qtbot):
    workspace = AnalyticsWorkspace()
    qtbot.addWidget(workspace.residue_page)
    try:
        original = make_result(summaries=[summary(f"L{i:03}", 1, i) for i in range(105)])
        workspace.set_result(original)
        assert len(workspace._source_result.summaries) == 105
        assert len(workspace._result.summaries) == 100
        assert "Top 100 of 105" in workspace._selection_notices[0].text()
        artifact = workspace.current_artifact(1)
        assert artifact.metadata["chart_selection"]["ligands_before_limit"] == 105
        assert workspace.fingerprint_panel.minimumHeight() >= 1400
    finally:
        workspace.dispose()


def test_uncertainty_callback_artifacts_disclose_chart_ranking(qtbot):
    from docklens.analytics import AnalysisContext
    from docklens.uncertainty import block_bootstrap_occupancy, block_bootstrap_difference
    from tests.test_dynamic_plotting import _matrix
    workspace = AnalyticsWorkspace()
    qtbot.addWidget(workspace.fingerprint_page)
    try:
        original = make_result(summaries=[summary(f"L{i:03}", 1, i) for i in range(105)])
        workspace.set_result(original)
        workspace.set_comparison(original)
        workspace.set_mode("md")
        context = AnalysisContext(mode="md")
        occupancy = block_bootstrap_occupancy(_matrix(), context, iterations=100, block_size=2)
        difference = block_bootstrap_difference(_matrix(), _matrix(), context,
                                               iterations=100, block_size_a=2, block_size_b=2)
        workspace._uncertainty_ready(occupancy)
        workspace._comparison_uncertainty_ready(difference)
        for panel in (workspace.uncertainty_panel, workspace.compare_uncertainty_panel):
            assert panel.artifact.metadata["chart_selection"]["ligands_displayed"] == 100
            assert "Top 100 of 105" in panel.artifact.metadata["selection_notice"]
        assert workspace.compare_uncertainty_panel.artifact.metadata["comparison_chart_selection"]["ligands_displayed"] == 100
    finally:
        workspace.dispose()


def test_evidence_is_current_when_residue_chart_is_hidden(qtbot):
    from tests.test_analytics import _result
    workspace = AnalyticsWorkspace()
    qtbot.addWidget(workspace.fingerprint_page)
    try:
        workspace.activate(3)
        workspace.residue_panel.artifact = None
        result = _result(["p1"], [("p1", "GLU166", "hbond", "H1", "OE1", 2.8)])
        workspace.set_result(result)
        assert "100.0%" in workspace.evidence_text("GLU166")
        assert workspace.residue_panel.artifact is None
    finally:
        workspace.dispose()
