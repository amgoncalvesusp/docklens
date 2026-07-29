"""Offscreen contracts for the DockLens 1.0 analytical workspace."""

from __future__ import annotations

from dataclasses import replace
import time

from PyQt5 import QtCore, QtWidgets

from docklens import batch_runner as br
from docklens.analytics_widgets import AnalyticsWorkspace, ChartPanel
from docklens.main_window import MainWindow
from docklens.plotting import build_residue_chart
from docklens.project_session import load_project, save_project
from docklens.observation_series import ObservationPoint, ObservationSeries
from docklens.results import make_result


def test_v1_workspace_navigation_and_ds_like_profile_are_always_available(qtbot):
    window = MainWindow()
    qtbot.addWidget(window)

    labels = [button.text() for button in window.workspace_buttons]

    assert labels == ["Residues", "Fingerprint", "Compare", "Tables"]
    assert window.workspace_stack.count() == 4
    assert window.analysis_combo.findData("complete") >= 0
    assert window.analysis_combo.findData("ds_like") >= 0
    assert "Discovery Studio" in window.analysis_combo.itemText(
        window.analysis_combo.findData("ds_like")
    )


def test_navigation_switches_peer_workspaces(qtbot):
    window = MainWindow()
    qtbot.addWidget(window)

    window.workspace_buttons[1].click()

    assert window.workspace_stack.currentIndex() == 1
    assert window.workspace_buttons[1].isChecked()


def test_analytics_refreshes_from_the_same_profiled_result_as_tables(
    qtbot, fixture_path
):
    window = MainWindow()
    qtbot.addWidget(window)
    window._result = br.run([fixture_path("minimal_complex.pdb")])

    window._refresh_tables()

    artifact = window.analytics_workspace.residue_panel.artifact
    assert artifact is not None
    assert artifact.metadata["total_observations"] == len(window._result.summaries)
    assert len(artifact.data) > 0


def test_mode_changes_scientific_labels_without_reinterpreting_raw_rows(
    qtbot, fixture_path
):
    window = MainWindow()
    qtbot.addWidget(window)
    window._result = br.run([fixture_path("minimal_complex.pdb")])
    window._refresh_tables()
    detail_count = window.detail_model.rowCount()

    window.mode_combo.setCurrentIndex(window.mode_combo.findData("md"))

    artifact = window.analytics_workspace.residue_panel.artifact
    assert artifact.metadata["mode"] == "md"
    assert artifact.figure.axes[0].get_xlabel() == (
        "MD occupancy (% of saved frames)"
    )
    assert window.detail_model.rowCount() == detail_count


def test_docklens_lens_selection_updates_inspector(qtbot, fixture_path):
    window = MainWindow()
    qtbot.addWidget(window)
    window._result = br.run([fixture_path("minimal_complex.pdb")])
    window._refresh_tables()
    residue = window._result.details[0].receptor_residue

    window.analytics_workspace.select_residue(residue)

    assert window.lens_residue.text() == residue
    assert residue in window.lens_evidence.text()


def test_major_workspaces_are_scrollable_on_compact_screens(qtbot):
    window = MainWindow()
    qtbot.addWidget(window)
    window.resize(900, 650)
    window.show()
    QtWidgets.QApplication.processEvents()

    for area in window.workspace_scroll_areas:
        assert area.horizontalScrollBarPolicy() == QtCore.Qt.ScrollBarAsNeeded
        assert area.verticalScrollBarPolicy() == QtCore.Qt.ScrollBarAsNeeded


def test_retention_requires_explicit_docking_a_and_md_b_roles(
    qtbot, fixture_path
):
    window = MainWindow()
    qtbot.addWidget(window)
    result = br.run([fixture_path("minimal_complex.pdb")])
    window._result = result
    window._refresh_tables()
    workspace = window.analytics_workspace
    workspace.set_comparison(result)
    retention_item = workspace.compare_analysis_combo.model().item(1)

    assert not retention_item.isEnabled()

    workspace.system_b_role.setCurrentIndex(
        workspace.system_b_role.findData("md")
    )
    assert retention_item.isEnabled()

    workspace.compare_analysis_combo.setCurrentIndex(1)
    assert not workspace.retention_panel.isHidden()

    workspace.system_a_role.setCurrentIndex(
        workspace.system_a_role.findData("md")
    )
    assert workspace.compare_analysis_combo.currentData() == "difference"


def test_chart_panel_disposes_qt_canvas_on_close(qtbot):
    panel = ChartPanel()
    qtbot.addWidget(panel)
    panel.set_artifact(build_residue_chart(make_result()))

    panel.close()
    QtCore.QCoreApplication.sendPostedEvents(
        None, QtCore.QEvent.DeferredDelete
    )

    assert panel._canvas is None


def test_main_window_disposes_all_chart_canvases_before_shutdown(qtbot):
    window = MainWindow()
    qtbot.addWidget(window)
    panels = (
        window.analytics_workspace.residue_panel,
        window.analytics_workspace.residue_barcode_panel,
        window.analytics_workspace.fingerprint_panel,
        window.analytics_workspace.similarity_panel,
        window.analytics_workspace.compare_panel,
        window.analytics_workspace.retention_panel,
        window.analytics_workspace.state_population_panel,
        window.analytics_workspace.state_timeline_panel,
        window.analytics_workspace.transition_panel,
        window.analytics_workspace.uncertainty_panel,
    )

    window.close()

    assert all(panel._canvas is None for panel in panels)


def test_fingerprint_workspace_contains_pose_family_and_dynamic_state_tools(
    qtbot, fixture_path
):
    window = MainWindow()
    qtbot.addWidget(window)
    window._result = br.run([fixture_path("minimal_complex.pdb")])
    window._refresh_tables()
    workspace = window.analytics_workspace

    assert [
        workspace.dynamic_tabs.tabText(index)
        for index in range(workspace.dynamic_tabs.count())
    ] == ["Population", "Pose order", "Transitions", "Confidence"]
    assert workspace.state_analysis.group_term == "pose family"
    assert workspace.dynamic_group_box.title() == "Pose families"
    assert not workspace.dynamic_tabs.isTabEnabled(2)
    assert not workspace.dynamic_tabs.isTabEnabled(3)
    assert workspace.state_table.rowCount() >= 1
    assert workspace.current_artifact(1) is workspace.fingerprint_panel.artifact
    workspace.fingerprint_export_combo.setCurrentIndex(
        workspace.fingerprint_export_combo.findData("similarity")
    )
    assert workspace.current_artifact(1) is workspace.similarity_panel.artifact
    for data, panel in (
        ("population", workspace.state_population_panel),
        ("timeline", workspace.state_timeline_panel),
    ):
        workspace.fingerprint_export_combo.setCurrentIndex(
            workspace.fingerprint_export_combo.findData(data)
        )
        assert workspace.current_artifact(1) is panel.artifact

    window.mode_combo.setCurrentIndex(window.mode_combo.findData("md"))

    assert workspace.state_analysis.group_term == "interaction state"
    assert workspace.dynamic_tabs.tabText(1) == "Timeline"
    assert workspace.dynamic_group_box.title() == "Dynamic interaction states"
    assert workspace.dynamic_tabs.isTabEnabled(2)
    assert workspace.dynamic_tabs.isTabEnabled(3)
    assert workspace.compute_uncertainty_button.isEnabled()
    workspace.fingerprint_export_combo.setCurrentIndex(
        workspace.fingerprint_export_combo.findData("transitions")
    )
    assert workspace.current_artifact(1) is workspace.transition_panel.artifact
    workspace.fingerprint_export_combo.setCurrentIndex(
        workspace.fingerprint_export_combo.findData("confidence")
    )
    assert workspace.current_artifact(1) is workspace.uncertainty_panel.artifact


def test_project_controls_are_visible_and_save_requires_a_result(qtbot, fixture_path):
    window = MainWindow()
    qtbot.addWidget(window)

    assert window.open_project_button.isEnabled()
    assert not window.save_project_button.isEnabled()

    window._result = br.run([fixture_path("minimal_complex.pdb")])
    window._files = [fixture_path("minimal_complex.pdb")]
    window._refresh_tables()

    assert window.save_project_button.isEnabled()


def test_project_round_trip_restores_profile_mode_workspace_and_cached_result(
    qtbot, fixture_path, tmp_path
):
    source = fixture_path("minimal_complex.pdb")
    original = MainWindow()
    qtbot.addWidget(original)
    original._files = [source]
    original._result = br.run([source])
    original.analysis_combo.setCurrentIndex(
        original.analysis_combo.findData("ds_like")
    )
    original.mode_combo.setCurrentIndex(original.mode_combo.findData("md"))
    original.analytics_workspace.state_threshold_spin.setValue(0.8)
    original._refresh_tables()
    observation_id = original._result.summaries[0].pose_id
    mapped_series = ObservationSeries(
        mode="md",
        points=(
            ObservationPoint(
                observation_id,
                0,
                frame_index=10,
                time_ns=2.0,
                replica_id="replica-B",
            ),
        ),
        time_step_ns=0.2,
    )
    original.analytics_workspace.set_observation_series(mapped_series)
    original.workspace_stack.setCurrentIndex(1)
    project_path = tmp_path / "roundtrip.docklens"
    save_project(original._project_state(), project_path)

    restored = MainWindow()
    qtbot.addWidget(restored)
    restored._restore_project(load_project(project_path))

    assert restored.analysis_combo.currentData() == "ds_like"
    assert restored.mode_combo.currentData() == "md"
    assert restored.workspace_stack.currentIndex() == 1
    assert restored.analytics_workspace.state_threshold_spin.value() == 0.8
    assert restored._result == original._result
    assert restored.analytics_workspace.state_analysis.series == mapped_series
    assert restored.save_project_button.isEnabled()


def test_explicit_series_reaches_states_bootstrap_and_evidence(qtbot, fixture_path):
    result = br.run([fixture_path("minimal_complex.pdb")])
    workspace = AnalyticsWorkspace()
    qtbot.addWidget(workspace.fingerprint_page)
    workspace.set_mode("md")
    workspace.set_result(result)
    observation_id = result.summaries[0].pose_id
    series = ObservationSeries(
        mode="md",
        points=(
            ObservationPoint(
                observation_id,
                0,
                frame_index=3,
                time_ns=0.3,
                replica_id="run-A",
            ),
        ),
        time_step_ns=0.1,
    )

    workspace.set_observation_series(series)

    assert workspace.state_analysis.series == series
    assert "explicit trajectory map" in workspace.state_status.text().lower()
    assert "1 replica" in workspace.state_status.text()


def test_md_comparison_exposes_independent_block_bootstrap_intervals(
    qtbot, fixture_path
):
    result = br.run([fixture_path("minimal_complex.pdb")])
    window = MainWindow()
    qtbot.addWidget(window)
    window._result = result
    window._refresh_tables()
    window.mode_combo.setCurrentIndex(window.mode_combo.findData("md"))
    window.analytics_workspace.set_comparison(result)
    workspace = window.analytics_workspace

    assert workspace.compare_confidence_button.isEnabled()
    workspace.compare_confidence_button.click()
    qtbot.waitUntil(
        lambda: workspace.compare_uncertainty_panel.artifact is not None
        and workspace.compare_uncertainty_panel.artifact.kind
        == "md-differential-occupancy-confidence-intervals",
        timeout=10_000,
    )

    assert workspace.compare_uncertainty_panel.artifact.kind == (
        "md-differential-occupancy-confidence-intervals"
    )
    assert workspace.compare_uncertainty_panel.artifact.metadata["design"] == (
        "independent B - A"
    )


def test_representative_path_must_be_a_declared_supported_structure(
    qtbot, fixture_path, tmp_path
):
    source = fixture_path("minimal_complex.pdb")
    window = MainWindow()
    qtbot.addWidget(window)
    window._files = [source]
    window._result = br.run([source])
    observation_id = window._result.summaries[0].pose_id

    assert window._trusted_representative_path(observation_id) is not None

    undeclared = tmp_path / "undeclared.pdb"
    undeclared.write_text("END\n", encoding="utf-8")
    summary = window._result.summaries[0]
    window._result = make_result(
        summaries=(replace(summary, source_path=str(undeclared)),),
        details=window._result.details,
    )

    assert window._trusted_representative_path(observation_id) is None

    window._result = make_result(
        summaries=(
            replace(
                summary,
                source_path=r"\\server\share\malicious.pdb",
            ),
        ),
        details=window._result.details,
    )
    assert window._trusted_representative_path(observation_id) is None


def test_project_actions_remain_scrollable_in_a_short_window(qtbot):
    window = MainWindow()
    qtbot.addWidget(window)
    window.resize(900, 300)
    window.show()
    QtWidgets.QApplication.processEvents()

    assert window.nav_scroll_area.verticalScrollBarPolicy() == (
        QtCore.Qt.ScrollBarAsNeeded
    )
    assert window.nav_scroll_area.verticalScrollBar().maximum() > 0


def test_stale_comparison_task_cannot_reopen_interval_panel(
    qtbot, fixture_path, monkeypatch
):
    import docklens.analytics_widgets as widgets

    original = widgets.block_bootstrap_difference

    def delayed(*args, **kwargs):
        time.sleep(0.15)
        return original(*args, **kwargs)

    monkeypatch.setattr(widgets, "block_bootstrap_difference", delayed)
    result = br.run([fixture_path("minimal_complex.pdb")])
    window = MainWindow()
    qtbot.addWidget(window)
    window._result = result
    window._refresh_tables()
    window.mode_combo.setCurrentIndex(window.mode_combo.findData("md"))
    workspace = window.analytics_workspace
    workspace.set_comparison(result)
    workspace.compare_confidence_button.click()
    workspace.system_a_role.setCurrentIndex(
        workspace.system_a_role.findData("docking")
    )
    workspace.compare_analysis_combo.setCurrentIndex(
        workspace.compare_analysis_combo.findData("retention")
    )
    qtbot.wait(400)

    assert workspace.compare_analysis_combo.currentData() == "retention"
    assert workspace.compare_uncertainty_panel.isHidden()
    assert workspace.compare_confidence_button.text() == (
        "Compute Δ confidence intervals"
    )
    assert workspace.current_artifact(2) is workspace.retention_panel.artifact
