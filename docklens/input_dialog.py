"""Qt dialog for building paired or combined DockLens input plans."""

from __future__ import annotations

import os
from pathlib import Path

from PyQt5 import QtCore, QtWidgets

from .input_plan import InputJob, InputPlan, SUPPORTED_SUFFIXES


def _structure_filter():
    return "Structures (*.mol2 *.pdb *.pdbqt);;All files (*)"


def _scan_folder(folder):
    root = Path(folder)
    return [
        str(path)
        for path in sorted(root.rglob("*"), key=lambda item: str(item).casefold())
        if path.is_file() and path.suffix.lower() in SUPPORTED_SUFFIXES
    ]


class _ReceptorGroupEditor(QtWidgets.QGroupBox):
    """One receptor and an append-only ligand source list."""

    changed = QtCore.pyqtSignal()

    def __init__(self, parent=None, index=1):
        super().__init__(f"Receptor group {index}", parent)
        self.group_id = f"receptor_{index}"
        layout = QtWidgets.QGridLayout(self)
        layout.addWidget(QtWidgets.QLabel("Receptor / protein:"), 0, 0)
        self.receptor_edit = QtWidgets.QLineEdit()
        self.receptor_edit.setPlaceholderText("Choose one PDB/PDBQT/MOL2 receptor")
        layout.addWidget(self.receptor_edit, 0, 1)
        choose = QtWidgets.QPushButton("Choose protein")
        choose.clicked.connect(self._choose_receptor)
        layout.addWidget(choose, 0, 2)
        layout.addWidget(QtWidgets.QLabel("Ligand / docking files:"), 1, 0)
        self.ligand_list = QtWidgets.QListWidget()
        self.ligand_list.setSelectionMode(QtWidgets.QAbstractItemView.ExtendedSelection)
        self.ligand_list.setMinimumHeight(82)
        layout.addWidget(self.ligand_list, 1, 1, 3, 1)
        add_files = QtWidgets.QPushButton("Add files")
        add_files.clicked.connect(self._add_files)
        layout.addWidget(add_files, 1, 2)
        add_folder = QtWidgets.QPushButton("Add folder")
        add_folder.clicked.connect(self._add_folder)
        layout.addWidget(add_folder, 2, 2)
        remove = QtWidgets.QPushButton("Remove selected")
        remove.clicked.connect(self._remove_selected)
        layout.addWidget(remove, 3, 2)
        self.receptor_edit.textChanged.connect(self.changed)
        self.ligand_list.itemChanged.connect(lambda _item: self.changed.emit())

    def _choose_receptor(self):
        path, _ = QtWidgets.QFileDialog.getOpenFileName(
            self, "Choose receptor / protein", "", _structure_filter()
        )
        if path:
            self.receptor_edit.setText(path)

    def _add_files(self):
        paths, _ = QtWidgets.QFileDialog.getOpenFileNames(
            self, "Add ligand / docking files", "", _structure_filter()
        )
        self.add_paths(paths)

    def _add_folder(self):
        folder = QtWidgets.QFileDialog.getExistingDirectory(
            self, "Add ligand / docking folder"
        )
        if folder:
            self.add_paths(_scan_folder(folder))

    def add_paths(self, paths):
        existing = {
            os.path.normcase(self.ligand_list.item(i).text())
            for i in range(self.ligand_list.count())
        }
        for path in paths:
            normalized = os.path.abspath(os.fspath(path))
            if os.path.normcase(normalized) in existing:
                continue
            self.ligand_list.addItem(normalized)
            existing.add(os.path.normcase(normalized))
        self.changed.emit()

    def _remove_selected(self):
        for item in self.ligand_list.selectedItems():
            self.ligand_list.takeItem(self.ligand_list.row(item))
        self.changed.emit()

    def paths(self):
        return tuple(self.ligand_list.item(i).text() for i in range(self.ligand_list.count()))

    def jobs(self):
        receptor = self.receptor_edit.text().strip()
        if not receptor:
            raise ValueError(f"{self.title()} requires a receptor file")
        paths = self.paths()
        if not paths:
            raise ValueError(f"{self.title()} requires at least one ligand file")
        return tuple(
            InputJob("paired", source, receptor, self.group_id)
            for source in paths
        )


class InputPlanDialog(QtWidgets.QDialog):
    """Build a multi-group paired plan or a multi-file combined plan."""

    def __init__(self, parent=None, *, title="Select DockLens input plan", paired_only=False):
        super().__init__(parent)
        self.setWindowTitle(title)
        self.resize(760, 560)
        self._groups = []
        layout = QtWidgets.QVBoxLayout(self)
        intro = QtWidgets.QLabel(
            "Choose either external receptor groups or complex files. "
            "Adding files again appends them to the current list."
        )
        intro.setWordWrap(True)
        layout.addWidget(intro)
        self.tabs = QtWidgets.QTabWidget()
        self.paired_tab = QtWidgets.QWidget()
        self._build_paired_tab()
        self.tabs.addTab(self.paired_tab, "Protein + ligands")
        if not paired_only:
            self.combined_tab = QtWidgets.QWidget()
            self._build_combined_tab()
            self.tabs.addTab(self.combined_tab, "Complex files")
        layout.addWidget(self.tabs, 1)
        buttons = QtWidgets.QDialogButtonBox(
            QtWidgets.QDialogButtonBox.Ok | QtWidgets.QDialogButtonBox.Cancel
        )
        buttons.accepted.connect(self._accept_checked)
        buttons.rejected.connect(self.reject)
        layout.addWidget(buttons)
        self.add_receptor_group()

    def _build_paired_tab(self):
        layout = QtWidgets.QVBoxLayout(self.paired_tab)
        add_group = QtWidgets.QPushButton("Add another receptor group")
        add_group.clicked.connect(self.add_receptor_group)
        layout.addWidget(add_group)
        scroll = QtWidgets.QScrollArea()
        scroll.setWidgetResizable(True)
        host = QtWidgets.QWidget()
        self._group_layout = QtWidgets.QVBoxLayout(host)
        self._group_layout.addStretch(1)
        scroll.setWidget(host)
        layout.addWidget(scroll, 1)

    def _build_combined_tab(self):
        layout = QtWidgets.QVBoxLayout(self.combined_tab)
        self.combined_list = QtWidgets.QListWidget()
        self.combined_list.setSelectionMode(QtWidgets.QAbstractItemView.ExtendedSelection)
        layout.addWidget(self.combined_list, 1)
        buttons = QtWidgets.QHBoxLayout()
        add_files = QtWidgets.QPushButton("Add files")
        add_files.clicked.connect(self._add_combined_files)
        add_folder = QtWidgets.QPushButton("Add folder")
        add_folder.clicked.connect(self._add_combined_folder)
        remove = QtWidgets.QPushButton("Remove selected")
        remove.clicked.connect(self._remove_combined)
        for button in (add_files, add_folder, remove):
            buttons.addWidget(button)
        buttons.addStretch(1)
        layout.addLayout(buttons)

    def add_receptor_group(self):
        editor = _ReceptorGroupEditor(index=len(self._groups) + 1)
        self._groups.append(editor)
        self._group_layout.insertWidget(self._group_layout.count() - 1, editor)

    def _add_combined_files(self):
        paths, _ = QtWidgets.QFileDialog.getOpenFileNames(
            self, "Add complex files", "", _structure_filter()
        )
        self._add_combined_paths(paths)

    def _add_combined_folder(self):
        folder = QtWidgets.QFileDialog.getExistingDirectory(self, "Add complex folder")
        if folder:
            self._add_combined_paths(_scan_folder(folder))

    def _add_combined_paths(self, paths):
        existing = {
            os.path.normcase(self.combined_list.item(i).text())
            for i in range(self.combined_list.count())
        }
        for path in paths:
            normalized = os.path.abspath(os.fspath(path))
            if os.path.normcase(normalized) in existing:
                continue
            self.combined_list.addItem(normalized)
            existing.add(os.path.normcase(normalized))

    def _remove_combined(self):
        for item in self.combined_list.selectedItems():
            self.combined_list.takeItem(self.combined_list.row(item))

    def _accept_checked(self):
        try:
            self.plan()
        except ValueError as exc:
            QtWidgets.QMessageBox.warning(self, "Incomplete input plan", str(exc))
            return
        self.accept()

    def plan(self) -> InputPlan:
        if self.tabs.currentWidget() is self.paired_tab:
            jobs = tuple(job for group in self._groups for job in group.jobs())
            return InputPlan(jobs)
        paths = tuple(
            self.combined_list.item(i).text()
            for i in range(self.combined_list.count())
        )
        if not paths:
            raise ValueError("Select at least one complex file")
        return InputPlan(
            InputJob("combined", path, group_id="combined_%d" % index)
            for index, path in enumerate(paths, 1)
        )


__all__ = ["InputPlanDialog"]
