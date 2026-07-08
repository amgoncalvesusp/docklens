"""
main_window.py — PyQt5 desktop UI for the interaction explorer.

Load a file / list / folder, resolve ligand vs. receptor, detect interactions and
show two sortable/filterable tables (Summary, Detail) with Okabe-Ito colouring.
Key residues are editable before running and re-applied afterwards without
re-running detection. Export to CSV / XLSX. Fallback ligand/receptor splits ask
for confirmation before use.
"""

from __future__ import annotations

import os

from PyQt5 import QtCore, QtGui, QtWidgets

from . import batch_runner as br
from . import export
from .entity_resolver import set_sides
from .interaction_core import (
    INTERACTION_COLORS,
    VALID_TYPES,
    color_hex,
    compute_interactions,
)

_NUMERIC_COLS = {
    "sol",
    "pose",
    "docking_score",
    "distance_A",
    "n_total_interactions",
    "n_key_residue_interactions",
    *VALID_TYPES,
}


class DataFrameModel(QtCore.QAbstractTableModel):
    """Table model over a pandas DataFrame with type-aware sort + colour."""

    def __init__(self, df, colour_type_col=None):
        super().__init__()
        self._df = df.reset_index(drop=True)
        self._colour_type_col = colour_type_col  # column name to shade by type

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
        if role == QtCore.Qt.UserRole:  # raw value for numeric sorting
            return val
        if role == QtCore.Qt.BackgroundRole and self._colour_type_col == col:
            if val in INTERACTION_COLORS:
                return QtGui.QBrush(QtGui.QColor(color_hex(val)))
        return None


class MultiFilterProxy(QtCore.QSortFilterProxyModel):
    """Sort numerically where possible; filter by type set / text / key-only."""

    def __init__(self):
        super().__init__()
        self.type_filter = None  # set of allowed types, or None
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
            if c < 0:
                return ""
            return model.data(model.index(row, c), QtCore.Qt.DisplayRole) or ""

        if self.type_filter is not None:
            if self._col("interaction_type") >= 0 and (
                cell("interaction_type") not in self.type_filter
            ):
                return False
        if self.key_only:
            if self._col("is_key_residue") >= 0 and (
                cell("is_key_residue") not in ("True", "true", "1")
            ):
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
        self.resize(1150, 720)
        self._files = []
        self._result = None
        self._build_ui()

    # ------------------------------------------------------------------ UI
    def _build_ui(self):
        central = QtWidgets.QWidget()
        self.setCentralWidget(central)
        root = QtWidgets.QVBoxLayout(central)

        top = QtWidgets.QHBoxLayout()
        b_files = QtWidgets.QPushButton("Open file(s)...")
        b_folder = QtWidgets.QPushButton("Open folder...")
        b_run = QtWidgets.QPushButton("Run detection")
        b_run.setStyleSheet("font-weight:bold;")
        b_files.clicked.connect(self._open_files)
        b_folder.clicked.connect(self._open_folder)
        b_run.clicked.connect(self._run)
        top.addWidget(b_files)
        top.addWidget(b_folder)
        top.addWidget(b_run)
        top.addStretch(1)
        b_csv = QtWidgets.QPushButton("Export CSV")
        b_xlsx = QtWidgets.QPushButton("Export XLSX")
        b_csv.clicked.connect(self._export_csv)
        b_xlsx.clicked.connect(self._export_xlsx)
        top.addWidget(b_csv)
        top.addWidget(b_xlsx)
        root.addLayout(top)

        row2 = QtWidgets.QHBoxLayout()
        row2.addWidget(QtWidgets.QLabel("Key residues:"))
        self.key_edit = QtWidgets.QLineEdit()
        self.key_edit.setPlaceholderText("e.g. ASP32 LYS36 SER70 (space/comma)")
        self.key_edit.editingFinished.connect(self._apply_key_residues)
        row2.addWidget(self.key_edit, 2)
        row2.addWidget(QtWidgets.QLabel("Search:"))
        self.search_edit = QtWidgets.QLineEdit()
        self.search_edit.setPlaceholderText("residue / file / ligand")
        self.search_edit.textChanged.connect(self._apply_filters)
        row2.addWidget(self.search_edit, 1)
        self.key_only_cb = QtWidgets.QCheckBox("Key residues only")
        self.key_only_cb.stateChanged.connect(self._apply_filters)
        row2.addWidget(self.key_only_cb)
        root.addLayout(row2)

        self.type_boxes = {}
        tbox = QtWidgets.QGroupBox("Interaction types (filter Detail table)")
        tgrid = QtWidgets.QGridLayout(tbox)
        for i, t in enumerate(VALID_TYPES):
            cb = QtWidgets.QCheckBox(t)
            cb.setChecked(True)
            cb.stateChanged.connect(self._apply_filters)
            self.type_boxes[t] = cb
            tgrid.addWidget(cb, i // 6, i % 6)
        root.addWidget(tbox)

        self.tabs = QtWidgets.QTabWidget()
        self.summary_view = QtWidgets.QTableView()
        self.detail_view = QtWidgets.QTableView()
        for v in (self.summary_view, self.detail_view):
            v.setSortingEnabled(True)
            v.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectRows)
            v.horizontalHeader().setStretchLastSection(True)
        self.tabs.addTab(self.summary_view, "Summary")
        self.tabs.addTab(self.detail_view, "Detail")
        root.addWidget(self.tabs, 1)

        self.status = self.statusBar()
        self.status.showMessage("Open a file, list or folder, then Run detection.")

        import pandas as pd

        self.summary_model = DataFrameModel(pd.DataFrame())
        self.detail_model = DataFrameModel(
            pd.DataFrame(), colour_type_col="interaction_type"
        )
        self.summary_proxy = MultiFilterProxy()
        self.summary_proxy.setSourceModel(self.summary_model)
        self.detail_proxy = MultiFilterProxy()
        self.detail_proxy.setSourceModel(self.detail_model)
        self.summary_view.setModel(self.summary_proxy)
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

    def _run(self):
        if not self._files:
            QtWidgets.QMessageBox.warning(
                self, "No input", "Open a file or folder first."
            )
            return
        result = br.run(self._files, key_residues=self.key_edit.text())
        if result.pending and not self._confirm_pending(result):
            self.status.showMessage("Run cancelled at ligand/receptor confirmation.")
            return
        self._result = result
        self._refresh_tables()
        self.status.showMessage(
            "%d ligand/pose row(s), %d interaction(s)."
            % (len(result.summaries), len(result.details))
        )

    def _confirm_pending(self, result):
        """Ask the user to confirm fallback ligand/receptor splits."""
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
        key_set = br.normalize_key_residues(self.key_edit.text())
        for pend in result.pending:
            res = pend.resolution
            if not res.ligand_atoms:
                continue
            set_sides(res)
            inters = compute_interactions(
                res.receptor_atoms,
                res.ligand_atoms,
                res.waters,
                self._selected_types(),
            )
            details = [
                br._detail_from_interaction(
                    it, res.ligand_id, pend.source_file, key_set
                )
                for it in inters
            ]
            summ = br._summarize(
                res.ligand_id,
                pend.source_file,
                pend.pose.sol,
                pend.pose.pose_index + 1,
                pend.pose.score,
                details,
                key_set,
            )
            result.details.extend(details)
            result.summaries.append(summ)
        result.pending = []
        return True

    def _refresh_tables(self):
        if self._result is None:
            return
        self.summary_model.set_dataframe(export.summary_dataframe(self._result))
        self.detail_model.set_dataframe(export.detail_dataframe(self._result))
        self._apply_filters()
        self.summary_view.resizeColumnsToContents()
        self.detail_view.resizeColumnsToContents()

    def _apply_key_residues(self):
        if self._result is None:
            return
        br.recompute_key(self._result, self.key_edit.text())
        self._refresh_tables()

    def _apply_filters(self):
        types = set(self._selected_types())
        self.detail_proxy.type_filter = types
        self.detail_proxy.key_only = self.key_only_cb.isChecked()
        self.detail_proxy.text_filter = self.search_edit.text()
        self.summary_proxy.text_filter = self.search_edit.text()
        self.detail_proxy.invalidateFilter()
        self.summary_proxy.invalidateFilter()

    def _export_csv(self):
        if not self._require_result():
            return
        prefix, _ = QtWidgets.QFileDialog.getSaveFileName(
            self, "Export CSV (prefix)", "interactions", "CSV (*.csv)"
        )
        if not prefix:
            return
        paths = export.export_csv(self._result, prefix)
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
        out = export.export_xlsx(self._result, path)
        QtWidgets.QMessageBox.information(self, "Exported", "Wrote:\n" + out)

    def _require_result(self):
        if self._result is None or not self._result.summaries:
            QtWidgets.QMessageBox.warning(
                self, "Nothing to export", "Run detection first."
            )
            return False
        return True


def launch():
    import sys

    app = QtWidgets.QApplication(sys.argv)
    win = MainWindow()
    win.show()
    return app.exec_()
