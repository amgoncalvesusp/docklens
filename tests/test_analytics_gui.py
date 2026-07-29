"""Offscreen contracts for the DockLens 1.0 analytical workspace."""

from __future__ import annotations

from PyQt5 import QtCore, QtWidgets

from docklens import batch_runner as br
from docklens.analytics_widgets import ChartPanel
from docklens.main_window import MainWindow
from docklens.plotting import build_residue_chart
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
    )

    window.close()

    assert all(panel._canvas is None for panel in panels)
