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
from .integration_result import write_integration_result
from .interaction_core import (
    INTERACTION_COLORS,
    VALID_TYPES,
    color_hex,
)
from .residue_keys import match_key_residues, parse_key_residues

INVENTOR = "Adriano Marques Gonçalves — Universidade de Araraquara (UNIARA)"
LOGGER = logging.getLogger(__name__)

# ---- BioLens / DockLens palette ----------------------------------------------
BLUE = "#003B5C"  # primary: titles, primary buttons, status
GRAY = "#60707A"  # secondary text, borders
OFFWHITE = "#FBFBF9"  # main background
CREAM = "#F0EFEA"  # secondary panels
STRUCT = "#E0E4E6"  # structure base / subtle fills

_STYLE = f"""
* {{ font-family: 'Inter', 'Open Sans', 'Segoe UI', Arial, sans-serif; }}
QMainWindow, QWidget {{ background: {OFFWHITE}; color: #1c2b33; }}
QLabel#title {{ color: {BLUE}; font-size: 22px; font-weight: 700; }}
QLabel#subtitle {{ color: {GRAY}; font-size: 12px; font-weight: 300; }}
QLabel#credit {{ color: {GRAY}; font-size: 11px; }}
QGroupBox {{
    background: {OFFWHITE}; border: 1px solid {GRAY}; border-radius: 8px;
    margin-top: 14px; padding: 8px; font-weight: 700; color: {BLUE};
}}
QGroupBox::title {{ subcontrol-origin: margin; left: 10px; padding: 0 4px; }}
QPushButton {{
    background: {CREAM}; border: 1px solid {GRAY}; border-radius: 6px;
    padding: 6px 12px; color: #1c2b33;
}}
QPushButton:hover {{ background: {STRUCT}; }}
QPushButton#primary {{ background: {BLUE}; color: white; border: 1px solid {BLUE};
    font-weight: 700; }}
QPushButton#primary:hover {{ background: #00527d; }}
QLineEdit, QComboBox {{
    background: white; border: 1px solid {GRAY}; border-radius: 6px; padding: 4px;
}}
QTableView {{
    background: {OFFWHITE}; alternate-background-color: {CREAM};
    gridline-color: {STRUCT}; selection-background-color: #cfe3ef;
    selection-color: #1c2b33;
}}
QHeaderView::section {{
    background: {BLUE}; color: white; padding: 4px; border: none; font-weight: 700;
}}
QTabBar::tab {{ background: {CREAM}; padding: 6px 14px; border: 1px solid {GRAY};
    border-bottom: none; border-top-left-radius: 6px; border-top-right-radius: 6px; }}
QTabBar::tab:selected {{ background: {BLUE}; color: white; }}
QStatusBar {{ background: {BLUE}; color: white; }}
QListWidget {{ background: white; border: 1px solid {GRAY}; border-radius: 6px; }}
"""


def resource_path(name):
    """Path to a bundled asset, working both from source and PyInstaller."""
    base = getattr(sys, "_MEIPASS", None)
    if base:
        return os.path.join(base, "docklens", "assets", name)
    return os.path.join(os.path.dirname(__file__), "assets", name)


class DataFrameModel(QtCore.QAbstractTableModel):
    """Table model over a pandas DataFrame with type-aware sort + colour."""

    def __init__(self, df, colour_type_col=None):
        super().__init__()
        self._df = df.reset_index(drop=True)
        self._colour_type_col = colour_type_col

    def set_dataframe(self, df):
        self.beginResetModel()
        self._df = df.reset_index(drop=True)
        self.endResetModel()

    def rowCount(self, parent=QtCore.QModelIndex()):
        return 0 if parent.isValid() else len(self._df)

    def columnCount(self, parent=QtCore.QModelIndex()):
        return 0 if parent.isValid() else self._df.shape[1]

    def headerData(self, section, orientation, role=QtCore.Qt.DisplayRole):
        if role != QtCore.Qt.DisplayRole:
            return None
        if orientation == QtCore.Qt.Horizontal:
            return str(self._df.columns[section])
        return str(section + 1)

    def data(self, index, role=QtCore.Qt.DisplayRole):
        if not index.isValid():
            return None
        col = self._df.columns[index.column()]
        val = self._df.iat[index.row(), index.column()]
        if role in (QtCore.Qt.DisplayRole, QtCore.Qt.EditRole):
            if val is None or (isinstance(val, float) and val != val):
                return ""
            return str(val)
        if role == QtCore.Qt.UserRole:
            return val
        if role == QtCore.Qt.BackgroundRole and self._colour_type_col == col:
            if val in INTERACTION_COLORS:
                return QtGui.QBrush(QtGui.QColor(color_hex(val)))
        if role == QtCore.Qt.ForegroundRole and self._colour_type_col == col:
            if val in INTERACTION_COLORS:
                return QtGui.QBrush(QtGui.QColor("white"))
        return None


class MultiFilterProxy(QtCore.QSortFilterProxyModel):
    """Sort numerically where possible; filter by type set / text / key-only."""

    def __init__(self):
        super().__init__()
        self.type_filter = None
        self.text_filter = ""
        self.key_only = False
        self.setSortRole(QtCore.Qt.UserRole)

    def _col(self, name):
        model = self.sourceModel()
        for c in range(model.columnCount()):
            if model.headerData(c, QtCore.Qt.Horizontal) == name:
                return c
        return -1

    def lessThan(self, left, right):
        lv = self.sourceModel().data(left, QtCore.Qt.UserRole)
        rv = self.sourceModel().data(right, QtCore.Qt.UserRole)
        try:
            return float(lv) < float(rv)
        except (TypeError, ValueError):
            return str(lv) < str(rv)

    def filterAcceptsRow(self, row, parent):
        model = self.sourceModel()

        def cell(name):
            c = self._col(name)
            return "" if c < 0 else (model.data(model.index(row, c)) or "")

        if self.type_filter is not None and self._col("interaction_type") >= 0:
            if cell("interaction_type") not in self.type_filter:
                return False
        if self.key_only and self._col("is_key_residue") >= 0:
            if cell("is_key_residue") not in ("True", "true", "1"):
                return False
        if self.text_filter:
            needle = self.text_filter.lower()
            hay = " ".join(
                cell(n) for n in ("receptor_residue", "source_file", "ligand_id")
            ).lower()
            if needle not in hay:
                return False
        return True


class MainWindow(QtWidgets.QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("DockLens")
        self.setWindowIcon(QtGui.QIcon(resource_path("docklens_icon.png")))
        self.resize(1200, 800)
        self._files = []
        self._result = None
        self._launch_manifest = None
        self._syncing = False
        self._key_invalid_tokens = ()
        self._build_ui()

    # ------------------------------------------------------------------ UI
    def _build_ui(self):
        central = QtWidgets.QWidget()
        self.setCentralWidget(central)
        root = QtWidgets.QVBoxLayout(central)

        root.addLayout(self._header())
        root.addLayout(self._toolbar())
        root.addWidget(self._filters_card())
        root.addWidget(self._type_card())
        root.addWidget(self._tables(), 1)

        self.status = self.statusBar()
        self.status.showMessage("Open a file, list or folder, then Run detection.")
        self._init_models()

    def _header(self):
        bar = QtWidgets.QHBoxLayout()
        logo = QtWidgets.QLabel()
        pix = QtGui.QPixmap(resource_path("docklens_icon.png"))
        if not pix.isNull():
            logo.setPixmap(pix.scaledToHeight(48, QtCore.Qt.SmoothTransformation))
        bar.addWidget(logo)
        titles = QtWidgets.QVBoxLayout()
        t = QtWidgets.QLabel("DockLens")
        t.setObjectName("title")
        s = QtWidgets.QLabel("Visual Intermolecular Interaction Analytics")
        s.setObjectName("subtitle")
        titles.addWidget(t)
        titles.addWidget(s)
        bar.addLayout(titles)
        bar.addStretch(1)
        about = QtWidgets.QPushButton("About")
        about.clicked.connect(self._about)
        bar.addWidget(about)
        return bar

    def _toolbar(self):
        top = QtWidgets.QHBoxLayout()
        b_files = QtWidgets.QPushButton("Open file(s)...")
        b_folder = QtWidgets.QPushButton("Open folder...")
        b_run = QtWidgets.QPushButton("Run detection")
        b_run.setObjectName("primary")
        b_reset = QtWidgets.QPushButton("Reset")
        b_files.clicked.connect(self._open_files)
        b_folder.clicked.connect(self._open_folder)
        b_run.clicked.connect(self._run)
        b_reset.clicked.connect(self._reset)
        for b in (b_files, b_folder, b_run, b_reset):
            top.addWidget(b)
        top.addSpacing(16)
        top.addWidget(QtWidgets.QLabel("H-bond criteria:"))
        self.preset_combo = QtWidgets.QComboBox()
        self.preset_combo.addItem("PLIP (default)", "plip")
        self.preset_combo.addItem(
            "DS-calibrated beta (explicit-H geometry)", "dsv"
        )
        top.addWidget(self.preset_combo)
        top.addSpacing(12)
        top.addWidget(QtWidgets.QLabel("Analysis view:"))
        self.analysis_combo = QtWidgets.QComboBox()
        self.analysis_combo.addItem("Complete", "complete")
        self.analysis_combo.addItem("Conservative polar/specific", "ds_like")
        self.analysis_combo.currentIndexChanged.connect(self._refresh_tables)
        top.addWidget(self.analysis_combo)
        top.addStretch(1)
        b_csv = QtWidgets.QPushButton("Export CSV")
        b_xlsx = QtWidgets.QPushButton("Export XLSX")
        b_csv.clicked.connect(self._export_csv)
        b_xlsx.clicked.connect(self._export_xlsx)
        top.addWidget(b_csv)
        top.addWidget(b_xlsx)
        return top

    def _filters_card(self):
        card = QtWidgets.QGroupBox("Filters & key residues")
        lay = QtWidgets.QGridLayout(card)

        lay.addWidget(QtWidgets.QLabel("Key residues:"), 0, 0)
        self.key_edit = QtWidgets.QLineEdit()
        self.key_edit.setPlaceholderText(
            "e.g. SER70; LYS73; GLU166 (space, comma, semicolon or line break)"
        )
        self.key_edit.editingFinished.connect(self._key_text_changed)
        lay.addWidget(self.key_edit, 0, 1, 1, 2)

        lay.addWidget(QtWidgets.QLabel("Search:"), 0, 3)
        self.search_edit = QtWidgets.QLineEdit()
        self.search_edit.setPlaceholderText("residue / file / ligand")
        self.search_edit.textChanged.connect(self._apply_filters)
        lay.addWidget(self.search_edit, 0, 4)
        self.key_only_cb = QtWidgets.QCheckBox("Key residues only")
        self.key_only_cb.stateChanged.connect(self._apply_filters)
        lay.addWidget(self.key_only_cb, 0, 5)

        lay.addWidget(
            QtWidgets.QLabel("Pick key residues from detected protein:"), 1, 0, 1, 2
        )
        self.res_filter = QtWidgets.QLineEdit()
        self.res_filter.setPlaceholderText("filter list...")
        self.res_filter.textChanged.connect(self._filter_residue_list)
        lay.addWidget(self.res_filter, 1, 3, 1, 3)

        self.res_list = QtWidgets.QListWidget()
        self.res_list.setMaximumHeight(140)
        self.res_list.itemChanged.connect(self._residue_checks_changed)
        lay.addWidget(self.res_list, 2, 0, 1, 6)
        self.key_status = QtWidgets.QLabel("No key residues configured.")
        self.key_status.setWordWrap(True)
        self.key_status.setObjectName("subtitle")
        lay.addWidget(self.key_status, 3, 0, 1, 6)
        return card

    def _type_card(self):
        card = QtWidgets.QGroupBox("Interaction types (filter Detail table)")
        grid = QtWidgets.QGridLayout(card)
        self.type_boxes = {}
        for i, t in enumerate(VALID_TYPES):
            cb = QtWidgets.QCheckBox(t)
            cb.setChecked(True)
            cb.stateChanged.connect(self._apply_filters)
            self.type_boxes[t] = cb
            grid.addWidget(cb, i // 6, i % 6)
        return card

    def _tables(self):
        self.tabs = QtWidgets.QTabWidget()
        self.summary_view = QtWidgets.QTableView()
        self.coverage_view = QtWidgets.QTableView()
        self.detail_view = QtWidgets.QTableView()
        for v in (self.summary_view, self.coverage_view, self.detail_view):
            v.setSortingEnabled(True)
            v.setAlternatingRowColors(True)
            v.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectRows)
            v.horizontalHeader().setStretchLastSection(True)
        self.tabs.addTab(self.summary_view, "Summary")
        self.tabs.addTab(self.coverage_view, "Key Residue Coverage")
        self.tabs.addTab(self.detail_view, "Detail")
        return self.tabs

    def _init_models(self):
        import pandas as pd

        self.summary_model = DataFrameModel(pd.DataFrame())
        self.coverage_model = DataFrameModel(pd.DataFrame())
        self.detail_model = DataFrameModel(
            pd.DataFrame(), colour_type_col="interaction_type"
        )
        self.summary_proxy = MultiFilterProxy()
        self.summary_proxy.setSourceModel(self.summary_model)
        self.coverage_proxy = MultiFilterProxy()
        self.coverage_proxy.setSourceModel(self.coverage_model)
        self.detail_proxy = MultiFilterProxy()
        self.detail_proxy.setSourceModel(self.detail_model)
        self.summary_view.setModel(self.summary_proxy)
        self.coverage_view.setModel(self.coverage_proxy)
        self.detail_view.setModel(self.detail_proxy)

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
                self.status.showMessage("Run cancelled at ligand/receptor confirmation.")
                return
            result = br.run(
                self._files,
                key_residues=self.key_edit.text(),
                confirm_fallback=True,
                hbond_preset=self._hbond_preset(),
            )
        self._result = result
        self._write_dockinghub_result(show_warning=True)
        self._populate_residue_list()
        self._refresh_tables()
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
        self.coverage_model.set_dataframe(
            export.key_residue_coverage_dataframe(view)
        )
        self.detail_model.set_dataframe(export.detail_dataframe(view))
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
                    f" ({len(match.matched_residues)} concrete receptor "
                    "residue(s))"
                )
            if match.unmatched_keys:
                message += "; unmatched: " + ", ".join(match.unmatched_keys)
            if match.ambiguous_keys:
                message += (
                    "; chain-ambiguous: " + ", ".join(match.ambiguous_keys)
                )
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
        self._syncing = False
        import pandas as pd

        self.summary_model.set_dataframe(pd.DataFrame())
        self.coverage_model.set_dataframe(pd.DataFrame())
        self.detail_model.set_dataframe(pd.DataFrame())
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
