"""
main_window.py — PyQt5 desktop UI for DockLens.

Load a file / list / folder, resolve ligand vs. receptor, detect interactions and
show two sortable/filterable tables (Summary, Detail) with distinct per-type
colouring. Key residues are editable as free text AND pickable from a checkbox
list of the detected protein residues; counts recompute without re-detection.
An H-bond criteria preset switches between PLIP (default) and a chemistry-aware
strict profile. Export to CSV / XLSX. Reset starts fresh.
"""

from __future__ import annotations

import logging
import os
import sys

from PyQt5 import QtCore, QtGui, QtWidgets

from . import __version__, batch_runner as br
from . import export
from .analysis_profiles import build_analysis_view
from .figure_export import export_figure_bundle
from .integration_result import write_integration_result
from .main_window_ui import build_main_window_ui
from .residue_keys import match_key_residues, parse_key_residues

INVENTOR = "Adriano Marques Gonçalves — Universidade de Araraquara (UNIARA)"
LOGGER = logging.getLogger(__name__)

# ---- Reciprocal-space atlas palette -----------------------------------------
BLUE = "#071A2E"
ACCENT = "#0072B2"
GRAY = "#5B6976"
OFFWHITE = "#F7F9FA"
CREAM = "#EEF2F4"
STRUCT = "#D8E0E5"

_STYLE = f"""
* {{ font-family: 'Aptos', 'Segoe UI Variable', 'Segoe UI', 'Noto Sans', sans-serif; }}
QMainWindow, QWidget {{ background: {OFFWHITE}; color: #142330; font-size: 12px; }}
QWidget#topBar, QWidget#navigationRail {{ background: {BLUE}; }}
QLabel#title {{ color: white; font-size: 20px; font-weight: 700; }}
QLabel#subtitle {{ color: #B7C9D8; font-size: 11px; font-weight: 400; }}
QLabel#credit {{ color: {GRAY}; font-size: 11px; }}
QLabel#workspaceTitle {{ color: #101C26; font-size: 20px; font-weight: 700; }}
QLabel#workspaceDescription {{ color: {GRAY}; font-size: 12px; }}
QLabel#sectionTitle {{ color: #101C26; font-size: 14px; font-weight: 700; }}
QLabel#metricStrip {{
    background: #EDF3F6; border: 1px solid #CAD7DE; border-radius: 3px;
    color: #294151; padding: 7px 10px;
}}
QGroupBox {{
    background: white; border: 1px solid {STRUCT}; border-radius: 4px;
    margin-top: 14px; padding: 9px; font-weight: 700; color: {BLUE};
}}
QGroupBox::title {{ subcontrol-origin: margin; left: 10px; padding: 0 4px; }}
QPushButton {{
    background: white; border: 1px solid #AEBBC4; border-radius: 4px;
    padding: 7px 11px; color: #172733;
}}
QPushButton:hover {{ background: {STRUCT}; }}
QPushButton:disabled {{ background: #EDF0F2; color: #87939B; border-color: #D4DBDF; }}
QPushButton#primary {{ background: {BLUE}; color: white; border: 1px solid {BLUE};
    font-weight: 700; }}
QPushButton#primary:hover {{ background: #123654; }}
QPushButton#navButton {{
    background: transparent; border: none; border-radius: 3px; color: #DCE7EE;
    padding: 10px 12px; text-align: left; font-weight: 600;
}}
QPushButton#navButton:hover {{ background: #102E49; color: white; }}
QPushButton#navButton:checked {{ background: {ACCENT}; color: white; }}
QLineEdit, QComboBox {{
    background: white; border: 1px solid #AEBBC4; border-radius: 4px;
    padding: 5px 7px; min-height: 20px;
}}
QPushButton:focus, QLineEdit:focus, QComboBox:focus, QListWidget:focus,
QTableView:focus {{ border: 2px solid {ACCENT}; }}
QTableView {{
    background: white; alternate-background-color: #F4F7F8;
    gridline-color: {STRUCT}; selection-background-color: #CDE8F7;
    selection-color: #142330; border: 1px solid {STRUCT};
}}
QHeaderView::section {{
    background: #E9EFF2; color: #203644; padding: 6px; border: none;
    border-right: 1px solid #D5DEE3; border-bottom: 1px solid #CBD6DC;
    font-weight: 700;
}}
QTabWidget::pane {{ border: 1px solid {STRUCT}; background: white; }}
QTabBar::tab {{ background: #E9EFF2; padding: 8px 14px;
    border: 1px solid {STRUCT}; border-bottom: none; }}
QTabBar::tab:selected {{ background: white; color: {BLUE}; font-weight: 700; }}
QStatusBar {{ background: {BLUE}; color: white; }}
QListWidget {{ background: white; border: 1px solid {STRUCT}; border-radius: 4px; }}
QWidget#analysisField {{ background: white; border: 1px solid {STRUCT};
    border-radius: 3px; }}
QWidget#lensPanel {{ background: white; border-left: 1px solid {STRUCT}; }}
QScrollArea#workspaceScroll {{ background: {OFFWHITE}; border: none; }}
"""


def resource_path(name):
    """Path to a bundled asset, working both from source and PyInstaller."""
    base = getattr(sys, "_MEIPASS", None)
    if base:
        return os.path.join(base, "docklens", "assets", name)
    return os.path.join(os.path.dirname(__file__), "assets", name)


class MainWindow(QtWidgets.QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("DockLens")
        self.setWindowIcon(QtGui.QIcon(resource_path("docklens_icon.png")))
        self.resize(1200, 800)
        self._files = []
        self._result = None
        self._comparison_result = None
        self._launch_manifest = None
        self._syncing = False
        self._key_invalid_tokens = ()
        self.setStyleSheet(_STYLE)
        self._build_ui()

    # ------------------------------------------------------------------ UI
    def _build_ui(self):
        build_main_window_ui(self, resource_path)

    def _workspace_changed(self, index):
        if 0 <= index < len(self.workspace_buttons):
            self.workspace_buttons[index].setChecked(True)

    def _update_lens(self, residue):
        self.lens_residue.setText(residue or "No selection")
        self.lens_evidence.setText(
            self.analytics_workspace.evidence_text(residue)
        )

    # -------------------------------------------------------------- actions
    def _open_files(self):
        files, _ = QtWidgets.QFileDialog.getOpenFileNames(
            self,
            "Select docking files",
            "",
            "Structures (*.mol2 *.pdb *.pdbqt);;All files (*)",
        )
        if files:
            self._files = list(files)
            self.status.showMessage("%d file(s) selected." % len(files))

    def _open_folder(self):
        folder = QtWidgets.QFileDialog.getExistingDirectory(self, "Select folder")
        if folder:
            self._files = [folder]
            self.status.showMessage("Folder selected: %s" % folder)

    def _selected_types(self):
        return [t for t, cb in self.type_boxes.items() if cb.isChecked()]

    def _hbond_preset(self):
        return self.preset_combo.currentData()

    def _analysis_profile(self):
        return self.analysis_combo.currentData()

    def _mode_changed(self):
        mode = self.mode_combo.currentData() or "docking"
        self.analytics_workspace.set_mode(mode)
        self._update_lens(self.analytics_workspace.selected_residue)

    def _load_comparison(self):
        files, _ = QtWidgets.QFileDialog.getOpenFileNames(
            self,
            "Select System B structures",
            "",
            "Structures (*.mol2 *.pdb *.pdbqt);;All files (*)",
        )
        if not files:
            return
        try:
            comparison = br.run(
                files,
                key_residues=self.key_edit.text(),
                hbond_preset=self._hbond_preset(),
            )
        except Exception:  # noqa: BLE001 - contain parser failures at GUI boundary
            LOGGER.exception("Comparison dataset analysis failed")
            QtWidgets.QMessageBox.critical(
                self,
                "System B could not be analyzed",
                "Verify the selected structures and analysis criteria, then "
                "try again.",
            )
            return
        self._comparison_result = comparison
        self._refresh_tables()
        self.workspace_stack.setCurrentIndex(2)
        self.status.showMessage(
            "System B loaded: %d observation(s), %d interaction row(s)."
            % (len(comparison.summaries), len(comparison.details))
        )

    def _export_figure(self):
        if not self._require_result():
            return
        artifact = self.analytics_workspace.current_artifact(
            self.workspace_stack.currentIndex()
        )
        if artifact is None:
            QtWidgets.QMessageBox.information(
                self,
                "Choose an analytical workspace",
                "Open Residues, Fingerprint or Compare before exporting a figure.",
            )
            return
        path, _ = QtWidgets.QFileDialog.getSaveFileName(
            self,
            "Export publication figure bundle",
            artifact.kind,
            "PNG image (*.png);;SVG vector (*.svg);;PDF vector (*.pdf)",
        )
        if not path:
            return
        suffix = os.path.splitext(path)[1].lower().lstrip(".")
        file_format = suffix if suffix in {"png", "svg", "pdf"} else "png"
        try:
            outputs = export_figure_bundle(
                artifact,
                path,
                formats=(file_format,),
                dpi=300,
                extra_metadata={
                    "analysis_profile": self._analysis_profile(),
                    "hbond_preset": self._hbond_preset(),
                },
            )
        except Exception:  # noqa: BLE001 - contain library errors at UI boundary
            LOGGER.exception("Publication figure export failed")
            QtWidgets.QMessageBox.critical(
                self,
                "Figure export failed",
                "DockLens could not write the figure bundle. Verify the "
                "destination and available disk space, then try again.",
            )
            return
        QtWidgets.QMessageBox.information(
            self,
            "Figure bundle exported",
            "Figure, source rows and reproducibility manifest were written:\n"
            + "\n".join(outputs),
        )

    def _run(self):
        if not self._files:
            QtWidgets.QMessageBox.warning(
                self, "No input", "Open a file or folder first."
            )
            return
        if self._launch_manifest is not None:
            result = br.run_paired(
                self._launch_manifest.receptor_path,
                self._launch_manifest.poses_path,
                key_residues=self.key_edit.text(),
                hbond_preset=self._hbond_preset(),
            )
        else:
            result = br.run(
                self._files,
                key_residues=self.key_edit.text(),
                hbond_preset=self._hbond_preset(),
            )
        if result.pending and self._launch_manifest is None:
            if not self._confirm_pending(result):
                self.status.showMessage(
                    "Run cancelled at ligand/receptor confirmation."
                )
                return
            result = br.run(
                self._files,
                key_residues=self.key_edit.text(),
                confirm_fallback=True,
                hbond_preset=self._hbond_preset(),
            )
        self._result = result
        self.dataset_context.setText(
            f"{len(result.summaries)} observation(s) · "
            f"{len(result.details)} raw interaction row(s)"
        )
        self._write_dockinghub_result(show_warning=True)
        self._populate_residue_list()
        self._refresh_tables()
        self.dataset_context.setText(
            f"DockingHub · {len(self._result.summaries)} pose(s)"
        )
        errors = sum(record.status == "error" for record in result.input_qc)
        self.status.showMessage(
            "%d ligand/pose row(s), %d interaction(s), %d error(s) [%s]."
            % (
                len(result.summaries),
                len(result.details),
                errors,
                self.preset_combo.currentText(),
            )
        )
        if errors:
            QtWidgets.QMessageBox.warning(
                self,
                "Completed with errors",
                "%d input file(s) could not be processed. Details will be included "
                "in the Input QC export sheet." % errors,
            )

    def load_manifest(self, manifest):
        """Load and immediately analyze an explicit DockingHub receptor/poses pair."""
        preset_index = self.preset_combo.findData(manifest.hbond_preset)
        if preset_index >= 0:
            self.preset_combo.setCurrentIndex(preset_index)
        self._key_invalid_tokens = ()
        self.key_edit.setText(" ".join(manifest.key_residues))
        self._files = [str(manifest.receptor_path), str(manifest.poses_path)]
        self._launch_manifest = manifest
        try:
            self._result = br.run_paired(
                manifest.receptor_path,
                manifest.poses_path,
                key_residues=manifest.key_residues,
                hbond_preset=manifest.hbond_preset,
            )
        except Exception:  # noqa: BLE001 - GUI boundary must contain parser failures
            LOGGER.exception("Unexpected DockingHub paired-analysis failure")
            self._result = None
            self.status.showMessage("DockingHub pair could not be analyzed.")
            QtWidgets.QMessageBox.critical(
                self,
                "DockingHub integration error",
                "The receptor/poses pair could not be analyzed. Verify both files in DockingHub and try again.",
            )
            return
        self._populate_residue_list()
        self._refresh_tables()
        roundtrip_failed = not self._write_dockinghub_result(show_warning=True)
        errors = sum(record.status == "error" for record in self._result.input_qc)
        self.status.showMessage(
            "DockingHub pair: %d pose(s), %d interaction(s), %d error(s)%s."
            % (
                len(self._result.summaries),
                len(self._result.details),
                errors,
                "; result export failed" if roundtrip_failed else "; result ready",
            )
        )

    def _confirm_pending(self, result):
        previews = "\n\n".join(
            "%s:\n%s" % (os.path.basename(p.source_file), p.preview)
            for p in result.pending
        )
        box = QtWidgets.QMessageBox(self)
        box.setWindowTitle("Confirm ligand/receptor split")
        box.setIcon(QtWidgets.QMessageBox.Question)
        box.setText(
            "Some files could not be split unambiguously and used a heuristic "
            "fallback. Confirm using the suggested split?"
        )
        box.setDetailedText(previews)
        confirm = box.addButton("Use suggested", QtWidgets.QMessageBox.AcceptRole)
        box.addButton("Cancel", QtWidgets.QMessageBox.RejectRole)
        box.exec_()
        if box.clickedButton() is not confirm:
            return False
        return True

    def _refresh_tables(self):
        if self._result is None:
            return
        view = build_analysis_view(self._result, self._analysis_profile())
        self.summary_model.set_dataframe(export.summary_dataframe(view))
        self.coverage_model.set_dataframe(export.key_residue_coverage_dataframe(view))
        self.detail_model.set_dataframe(export.detail_dataframe(view))
        self.analytics_workspace.set_result(view)
        if self._comparison_result is not None:
            comparison_view = build_analysis_view(
                self._comparison_result, self._analysis_profile()
            )
            self.analytics_workspace.set_comparison(comparison_view)
        else:
            self.analytics_workspace.set_comparison(None)
        self.lens_profile.setText(
            "Profile: %s\nCounting unit: one presence per observation, "
            "residue and interaction type"
            % self.analysis_combo.currentText()
        )
        self._update_lens(self.analytics_workspace.selected_residue)
        self._apply_filters()
        if len(view.details) < 500:
            self.summary_view.resizeColumnsToContents()
            self.coverage_view.resizeColumnsToContents()
            self.detail_view.resizeColumnsToContents()

    # ---- key residues: text field <-> checkbox list stay in sync ----
    def _populate_residue_list(self):
        self._syncing = True
        self.res_list.clear()
        key_set = br.normalize_key_residues(self.key_edit.text())
        for res in sorted(self._result.receptor_residues):
            item = QtWidgets.QListWidgetItem(res)
            item.setFlags(item.flags() | QtCore.Qt.ItemIsUserCheckable)
            checked = res.upper() in key_set or res.rstrip("_").upper() in key_set
            item.setCheckState(QtCore.Qt.Checked if checked else QtCore.Qt.Unchecked)
            self.res_list.addItem(item)
        self._syncing = False
        self._update_key_status()

    def _residue_checks_changed(self, _item):
        if self._syncing:
            return
        checked = [
            self.res_list.item(i).text()
            for i in range(self.res_list.count())
            if self.res_list.item(i).checkState() == QtCore.Qt.Checked
        ]
        self._syncing = True
        self.key_edit.setText(" ".join(checked))
        self._syncing = False
        self._key_invalid_tokens = ()
        self._recompute_key()

    def _key_text_changed(self):
        if self._syncing:
            return
        parsed = parse_key_residues(self.key_edit.text())
        key_set = frozenset(parsed.keys)
        self._key_invalid_tokens = parsed.invalid
        self._syncing = True
        self.key_edit.setText(" ".join(parsed.keys))
        for i in range(self.res_list.count()):
            it = self.res_list.item(i)
            on = (
                it.text().upper() in key_set or it.text().rstrip("_").upper() in key_set
            )
            it.setCheckState(QtCore.Qt.Checked if on else QtCore.Qt.Unchecked)
        self._syncing = False
        self._recompute_key()

    def _update_key_status(self):
        keys = tuple(sorted(br.normalize_key_residues(self.key_edit.text())))
        if not keys:
            message = "No key residues configured."
        elif self._result is None:
            message = (
                f"{len(keys)} key residue(s) configured; matching will be "
                "checked after detection."
            )
        else:
            match = match_key_residues(keys, self._result.receptor_residues)
            message = (
                f"{len(match.matched_keys)} of {len(keys)} key residue "
                "identifier(s) matched"
            )
            if match.matched_residues:
                message += (
                    f" ({len(match.matched_residues)} concrete receptor residue(s))"
                )
            if match.unmatched_keys:
                message += "; unmatched: " + ", ".join(match.unmatched_keys)
            if match.ambiguous_keys:
                message += "; chain-ambiguous: " + ", ".join(match.ambiguous_keys)
            message += "."
        if self._key_invalid_tokens:
            message += " Invalid: " + ", ".join(self._key_invalid_tokens) + "."
        self.key_status.setText(message)

    def _recompute_key(self):
        if self._result is None:
            self._update_key_status()
            return
        self._result = br.recompute_key(self._result, self.key_edit.text())
        self._refresh_tables()
        self._update_key_status()
        self._write_dockinghub_result(show_warning=False)

    def _write_dockinghub_result(self, *, show_warning):
        if (
            self._result is None
            or self._launch_manifest is None
            or getattr(self._launch_manifest, "result_path", None) is None
        ):
            return True
        try:
            write_integration_result(self._launch_manifest, self._result)
            return True
        except (OSError, ValueError):
            self.status.showMessage("DockingHub result export failed.")
            if show_warning:
                QtWidgets.QMessageBox.warning(
                    self,
                    "DockingHub result warning",
                    "The analysis is available in DockLens, but its DockingHub result file could not be written.",
                )
            return False

    def _filter_residue_list(self, text):
        needle = text.lower()
        for i in range(self.res_list.count()):
            it = self.res_list.item(i)
            it.setHidden(needle not in it.text().lower())

    def _apply_filters(self):
        types = set(self._selected_types())
        self.detail_proxy.type_filter = types
        self.detail_proxy.key_only = self.key_only_cb.isChecked()
        self.detail_proxy.text_filter = self.search_edit.text()
        self.summary_proxy.text_filter = self.search_edit.text()
        self.coverage_proxy.text_filter = self.search_edit.text()
        self.detail_proxy.invalidateFilter()
        self.summary_proxy.invalidateFilter()
        self.coverage_proxy.invalidateFilter()

    def _reset(self):
        """Clear everything for a fresh analysis."""
        self._files = []
        self._result = None
        self._comparison_result = None
        self._launch_manifest = None
        self._syncing = True
        self._key_invalid_tokens = ()
        self.key_edit.clear()
        self.search_edit.clear()
        self.res_filter.clear()
        self.res_list.clear()
        self.key_only_cb.setChecked(False)
        for cb in self.type_boxes.values():
            cb.setChecked(True)
        self.preset_combo.setCurrentIndex(0)
        self.analysis_combo.setCurrentIndex(0)
        self.mode_combo.setCurrentIndex(0)
        self._syncing = False
        import pandas as pd

        self.summary_model.set_dataframe(pd.DataFrame())
        self.coverage_model.set_dataframe(pd.DataFrame())
        self.detail_model.set_dataframe(pd.DataFrame())
        self.analytics_workspace.set_comparison(None)
        self.analytics_workspace.set_result(None)
        self.dataset_context.setText("No dataset loaded")
        self._update_lens("")
        self._update_key_status()
        self.status.showMessage(
            "Reset. Open a file, list or folder, then Run detection."
        )

    def _export_csv(self):
        if not self._require_result():
            return
        prefix, _ = QtWidgets.QFileDialog.getSaveFileName(
            self, "Export CSV (prefix)", "interactions", "CSV (*.csv)"
        )
        if not prefix:
            return
        export_filter = self._choose_export_filter(include_matrix_mode=False)
        if export_filter is None:
            return
        try:
            paths = export.export_csv(self._result, prefix, export_filter)
        except Exception as exc:  # noqa: BLE001 - present a safe UI error
            QtWidgets.QMessageBox.critical(self, "Export failed", str(exc))
            return
        QtWidgets.QMessageBox.information(
            self, "Exported", "Wrote:\n" + "\n".join(paths)
        )

    def _export_xlsx(self):
        if not self._require_result():
            return
        path, _ = QtWidgets.QFileDialog.getSaveFileName(
            self, "Export XLSX", "interactions.xlsx", "Excel (*.xlsx)"
        )
        if not path:
            return
        export_filter = self._choose_export_filter(include_matrix_mode=True)
        if export_filter is None:
            return
        try:
            out = export.export_xlsx(self._result, path, export_filter)
        except Exception as exc:  # noqa: BLE001 - present a safe UI error
            QtWidgets.QMessageBox.critical(self, "Export failed", str(exc))
            return
        QtWidgets.QMessageBox.information(self, "Exported", "Wrote:\n" + out)

    def _choose_export_filter(self, include_matrix_mode):
        box = QtWidgets.QMessageBox(self)
        box.setWindowTitle("Export scope")
        box.setText("Choose which results to export.")
        all_button = box.addButton("All results", QtWidgets.QMessageBox.AcceptRole)
        filtered_button = box.addButton(
            "Filtered interactions (all poses)", QtWidgets.QMessageBox.ActionRole
        )
        box.addButton(QtWidgets.QMessageBox.Cancel)
        box.exec_()
        clicked = box.clickedButton()
        if clicked not in (all_button, filtered_button):
            return None
        scope = "filtered" if clicked is filtered_button else "all"
        matrix_mode = "count"
        if include_matrix_mode:
            label, accepted = QtWidgets.QInputDialog.getItem(
                self,
                "Residue Matrix",
                "Cell values:",
                ["Count", "Presence (0/1)"],
                0,
                False,
            )
            if not accepted:
                return None
            matrix_mode = "presence" if label.startswith("Presence") else "count"
        return br.ExportFilter(
            scope=scope,
            interaction_types=(
                frozenset(self._selected_types()) if scope == "filtered" else None
            ),
            text=self.search_edit.text() if scope == "filtered" else "",
            key_only=(self.key_only_cb.isChecked() if scope == "filtered" else False),
            matrix_mode=matrix_mode,
            analysis_profile=self._analysis_profile(),
        )

    def _require_result(self):
        if self._result is None or (
            not self._result.summaries and not self._result.input_qc
        ):
            QtWidgets.QMessageBox.warning(
                self, "Nothing to export", "Run detection first."
            )
            return False
        return True

    def resizeEvent(self, event):
        super().resizeEvent(event)
        if hasattr(self, "lens_panel"):
            self.lens_panel.setVisible(self.width() >= 1080)
        if hasattr(self, "nav_rail"):
            compact = self.width() < 820
            self.nav_rail.setMinimumWidth(104 if compact else 132)
            self.nav_rail.setMaximumWidth(122 if compact else 172)

    def closeEvent(self, event):
        if hasattr(self, "analytics_workspace"):
            self.analytics_workspace.dispose()
        super().closeEvent(event)

    def _about(self):
        dlg = QtWidgets.QDialog(self)
        dlg.setWindowTitle("About DockLens")
        lay = QtWidgets.QVBoxLayout(dlg)
        logo = QtWidgets.QLabel()
        pix = QtGui.QPixmap(resource_path("docklens_logo.png"))
        if not pix.isNull():
            logo.setPixmap(pix.scaledToWidth(360, QtCore.Qt.SmoothTransformation))
            logo.setAlignment(QtCore.Qt.AlignCenter)
        lay.addWidget(logo)
        for text, obj in [
            ("DockLens %s" % __version__, "title"),
            ("Visual Intermolecular Interaction Analytics", "subtitle"),
            ("Inventor: %s" % INVENTOR, "credit"),
        ]:
            lb = QtWidgets.QLabel(text)
            lb.setObjectName(obj)
            lb.setAlignment(QtCore.Qt.AlignCenter)
            lb.setWordWrap(True)
            lay.addWidget(lb)
        btn = QtWidgets.QPushButton("Close")
        btn.clicked.connect(dlg.accept)
        lay.addWidget(btn)
        dlg.exec_()


def launch(launch_manifest=None):
    app = QtWidgets.QApplication(sys.argv)
    app.setApplicationName("DockLens")
    app.setStyleSheet(_STYLE)
    app.setWindowIcon(QtGui.QIcon(resource_path("docklens_icon.png")))

    splash_pix = QtGui.QPixmap(resource_path("docklens_logo.png"))
    splash = None
    if not splash_pix.isNull():
        splash = QtWidgets.QSplashScreen(
            splash_pix.scaledToWidth(560, QtCore.Qt.SmoothTransformation)
        )
        splash.showMessage(
            "  " + INVENTOR,
            QtCore.Qt.AlignBottom | QtCore.Qt.AlignHCenter,
            QtGui.QColor(GRAY),
        )
        splash.show()
        app.processEvents()

    win = MainWindow()
    win.show()
    if launch_manifest is not None:
        win.load_manifest(launch_manifest)
    if splash is not None:
        splash.finish(win)
    return app.exec_()
