"""Offscreen UI integration tests for the v0.3 workflow."""

from __future__ import annotations

from pathlib import Path

from PyQt5 import QtWidgets

from docklens import batch_runner as br
from docklens.main_window import MainWindow


def test_window_runs_filters_and_recomputes_keys(qtbot, fixture_path):
    window = MainWindow()
    qtbot.addWidget(window)
    window._files = [fixture_path("minimal_complex.pdb")]

    window._run()

    assert window._result.summaries
    residue = window._result.details[0].receptor_residue
    window.key_edit.setText(residue)
    window._key_text_changed()
    assert window._result.key_residues == frozenset({residue})
    assert window.detail_proxy.rowCount() > 0


def test_window_exports_selected_scope(qtbot, fixture_path, tmp_path, monkeypatch):
    window = MainWindow()
    qtbot.addWidget(window)
    window._result = br.run([fixture_path("minimal_complex.pdb")])
    output = tmp_path / "ui-export.xlsx"
    monkeypatch.setattr(
        QtWidgets.QFileDialog,
        "getSaveFileName",
        lambda *args, **kwargs: (str(output), "Excel (*.xlsx)"),
    )
    monkeypatch.setattr(
        window,
        "_choose_export_filter",
        lambda include_matrix_mode: br.ExportFilter(
            scope="filtered", interaction_types=frozenset({"hbond"})
        ),
    )
    monkeypatch.setattr(QtWidgets.QMessageBox, "information", lambda *args: None)

    window._export_xlsx()

    assert Path(output).is_file()
