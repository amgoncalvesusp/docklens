"""Key selections remain drafts until one explicit, atomic application."""

from dataclasses import replace
from threading import Event, Timer
from unittest.mock import Mock

from PyQt5 import QtCore, QtWidgets
import pytest

from docklens import batch_runner as br
from docklens.main_window import MainWindow
from docklens.analysis_tasks import AnalysisTaskRunner
from docklens.project_session import load_project, save_project


@pytest.fixture
def window(qtbot, multi_source_result, monkeypatch):
    widget = MainWindow()
    qtbot.addWidget(widget)
    widget._result = replace(
        multi_source_result,
        receptor_residues=frozenset({"GLU166", "SER70A", "SER70B", "ASN170"}),
    )
    widget.key_edit.setText("GLU166")
    widget._populate_residue_list()
    monkeypatch.setattr(widget, "_refresh_tables", Mock())
    monkeypatch.setattr(widget, "_write_dockinghub_result", Mock(return_value=True))
    return widget


def _item(window, text):
    return next(window.res_list.item(i) for i in range(window.res_list.count())
                if window.res_list.item(i).text() == text)


def test_text_and_checkbox_drafts_do_not_compute_or_write(window, monkeypatch):
    original = window._result
    compute = Mock(wraps=br.recompute_key)
    monkeypatch.setattr(br, "recompute_key", compute)
    window.key_edit.setText("SER70")
    window._key_text_changed()
    _item(window, "ASN170").setCheckState(QtCore.Qt.Checked)
    assert window._result is original
    compute.assert_not_called()
    window._refresh_tables.assert_not_called()
    window._write_dockinghub_result.assert_not_called()
    assert window.recalculate_keys_button.isEnabled()
    assert "pending" in window.key_status.text().lower()


def test_one_apply_updates_both_systems_and_writes_once(window, monkeypatch):
    window._comparison_result = window._result
    compute = Mock(wraps=br.recompute_key)
    monkeypatch.setattr(br, "recompute_key", compute)
    window.key_edit.setText("ASN170; SER70")
    window._key_text_changed()
    window.recalculate_keys_button.click()
    assert window._result.key_residues == frozenset({"ASN170", "SER70"})
    assert window._comparison_result.key_residues == window._result.key_residues
    assert compute.call_count == 2
    window._refresh_tables.assert_called_once()
    window._write_dockinghub_result.assert_called_once()
    assert not window.recalculate_keys_button.isEnabled()
    window._recompute_key()
    assert compute.call_count == 2


def test_comparison_failure_preserves_both_applied_results(window, monkeypatch):
    original = window._result
    window._comparison_result = original
    window.key_edit.setText("SER70")
    error = Mock(side_effect=[br.recompute_key(original, "SER70"), ValueError("bad B")])
    monkeypatch.setattr(br, "recompute_key", error)
    monkeypatch.setattr(QtWidgets.QMessageBox, "critical", Mock())
    window._recompute_key()
    assert window._result is original
    assert window._comparison_result is original
    assert window.recalculate_keys_button.isEnabled()
    window._write_dockinghub_result.assert_not_called()


def test_checkboxes_preserve_chainless_and_unmatched_draft_keys(window):
    window.key_edit.setText("SER70; GLU999")
    window._key_text_changed()
    assert _item(window, "SER70A").checkState() == QtCore.Qt.Checked
    assert _item(window, "SER70B").checkState() == QtCore.Qt.Checked
    _item(window, "ASN170").setCheckState(QtCore.Qt.Checked)
    assert br.normalize_key_residues(window.key_edit.text()) == {"SER70", "GLU999", "ASN170"}
    _item(window, "SER70A").setCheckState(QtCore.Qt.Unchecked)
    assert br.normalize_key_residues(window.key_edit.text()) == {"SER70B", "GLU999", "ASN170"}


def test_discard_restores_applied_keys_without_calculation(window):
    original = window._result
    window.key_edit.setText("SER70")
    window._key_text_changed()
    window.discard_keys_button.click()
    assert window.key_edit.text() == "GLU166"
    assert window._result is original
    assert not window.recalculate_keys_button.isEnabled()
    window._refresh_tables.assert_not_called()
    window._write_dockinghub_result.assert_not_called()


@pytest.mark.parametrize("action", ["_export_csv", "_export_xlsx", "_export_figure", "_export_all_tiff", "_save_project"])
def test_pending_keys_block_exports_and_project_dialogs(window, monkeypatch, action):
    window.key_edit.setText("SER70")
    dialog = Mock(side_effect=AssertionError("No destination dialog while pending"))
    monkeypatch.setattr(QtWidgets.QFileDialog, "getSaveFileName", dialog)
    monkeypatch.setattr(QtWidgets.QFileDialog, "getExistingDirectory", dialog)
    warning = Mock()
    monkeypatch.setattr(QtWidgets.QMessageBox, "warning", warning)
    getattr(window, action)()
    warning.assert_called_once()
    dialog.assert_not_called()


def test_project_snapshot_cannot_claim_pending_keys(window):
    window.key_edit.setText("SER70")
    with pytest.raises(ValueError, match="pending"):
        window._project_state()


def test_large_recalculation_keeps_gui_responsive_and_applies_once(window, qtbot, monkeypatch):
    assert hasattr(window, "_key_tasks")
    window._key_background_threshold = 0
    started, release = Event(), Event()
    real_compute = br.recompute_key

    def compute(result, keys):
        started.set()
        assert release.wait(5)
        return real_compute(result, keys)

    monkeypatch.setattr(br, "recompute_key", compute)
    window.key_edit.setText("SER70")
    window._recompute_key()
    try:
        qtbot.waitUntil(started.is_set)
        assert window._key_recalculation_running
        assert "recalculating" in window.key_status.text().lower()
        release.set()
        qtbot.waitUntil(lambda: not window._key_recalculation_running)
        assert window._result.key_residues == frozenset({"SER70"})
        window._write_dockinghub_result.assert_called_once()
    finally:
        release.set()


@pytest.mark.parametrize("cancel_action", ["reset", "discard", "edit", "restore", "close"])
def test_context_change_invalidates_an_inflight_key_recalculation(window, qtbot, monkeypatch, cancel_action):
    assert hasattr(window, "_key_tasks")
    window._key_background_threshold = 0
    started, release = Event(), Event()
    real_compute = br.recompute_key
    original = window._result
    project = window._project_state()

    def compute(result, keys):
        started.set()
        assert release.wait(5)
        return real_compute(result, keys)

    monkeypatch.setattr(br, "recompute_key", compute)
    window.key_edit.setText("SER70")
    window._recompute_key()
    try:
        qtbot.waitUntil(started.is_set)
        if cancel_action == "reset":
            window._reset()
        elif cancel_action == "discard":
            window._discard_key_changes()
        elif cancel_action == "edit":
            window.key_edit.setText("ASN170")
        elif cancel_action == "restore":
            window._restore_project(project)
        else:
            timer = Timer(0.05, release.set)
            timer.start()
            window.close()
            timer.join()
        release.set()
        window._key_tasks.wait_for_done()
        QtWidgets.QApplication.processEvents()
        if cancel_action == "reset":
            assert window._result is None
        else:
            assert window._result is original
            assert window._result.key_residues == frozenset({"GLU166"})
        window._write_dockinghub_result.assert_not_called()
    finally:
        release.set()


def test_applied_keys_round_trip_in_schema5(window, tmp_path, qtbot):
    window.analytics_workspace.set_ranking_criterion("fingerprint_features")
    window.key_edit.setText("SER70")
    window._recompute_key()
    path = tmp_path / "applied.docklens"
    save_project(window._project_state(), path)
    restored = MainWindow()
    qtbot.addWidget(restored)
    restored._restore_project(load_project(path))
    assert restored._result.key_residues == frozenset({"SER70"})
    assert restored.key_edit.text() == "SER70"
    assert not restored.recalculate_keys_button.isEnabled()
    assert restored.analytics_workspace.ranking_criterion == "fingerprint_features"


def test_cancelled_queue_does_not_delay_latest_key_job(qtbot):
    runner = AnalysisTaskRunner(max_threads=1, clear_pending_on_invalidate=True)
    started, release, obsolete_ran, latest_ran = Event(), Event(), Event(), Event()
    completed = []

    def first():
        started.set()
        assert release.wait(5)

    try:
        runner.start(first, on_success=completed.append)
        qtbot.waitUntil(started.is_set)
        runner.start(obsolete_ran.set, on_success=completed.append)
        runner.invalidate()
        runner.start(latest_ran.set, on_success=lambda value: completed.append("latest"))
        release.set()
        qtbot.waitUntil(lambda: completed == ["latest"])
        assert latest_ran.is_set()
        assert not obsolete_ran.is_set()
    finally:
        release.set()
        runner.wait_for_done()
