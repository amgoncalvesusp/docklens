"""Desktop project persistence and representative-opening behavior."""

from __future__ import annotations

import logging
from pathlib import Path

from PyQt5 import QtCore, QtGui, QtWidgets

from . import __version__, batch_runner as br
from .project_session import (
    ProjectDataset,
    ProjectState,
    build_project_input,
    is_local_filesystem_path,
    load_project,
    save_project,
    validate_project_inputs,
)


LOGGER = logging.getLogger(__name__)
_WORKSPACES = ("residues", "fingerprint", "compare", "tables")


class ProjectControllerMixin:
    """Add transactional project save/load behavior to ``MainWindow``."""

    def _project_source_files(self, paths):
        suffixes = {".mol2", ".pdb", ".pdbqt"}
        resolved = []
        for value in paths:
            path = Path(value).expanduser()
            if path.is_dir():
                resolved.extend(
                    candidate
                    for candidate in path.rglob("*")
                    if candidate.is_file()
                    and candidate.suffix.lower() in suffixes
                )
            elif path.is_file():
                resolved.append(path)
        return tuple(
            build_project_input(path)
            for path in sorted(
                {path.resolve() for path in resolved},
                key=lambda item: str(item).casefold(),
            )
        )

    def _project_state(self):
        active_index = self.workspace_stack.currentIndex()
        active = (
            _WORKSPACES[active_index]
            if 0 <= active_index < len(_WORKSPACES)
            else "residues"
        )
        workspace = self.analytics_workspace
        primary = ProjectDataset(
            label="System A",
            mode=self.mode_combo.currentData() or "docking",
            time_step_ns=workspace.time_step_spin.value(),
            inputs=self._project_source_files(self._files),
            result=self._result,
            observation_series=workspace._series,
        )
        comparison = None
        if self._comparison_result is not None:
            comparison = ProjectDataset(
                label="System B",
                mode=(
                    workspace.system_b_role.currentData()
                    or self.mode_combo.currentData()
                    or "docking"
                ),
                time_step_ns=workspace.time_step_spin.value(),
                inputs=self._project_source_files(self._comparison_files),
                result=self._comparison_result,
                observation_series=workspace._comparison_series,
            )
        return ProjectState(
            app_version=__version__,
            analysis_profile=self._analysis_profile(),
            hbond_preset=self._hbond_preset(),
            key_residues=tuple(
                sorted(br.normalize_key_residues(self.key_edit.text()))
            ),
            selected_types=tuple(self._selected_types()),
            active_workspace=active,
            selected_residue=workspace.selected_residue,
            state_threshold=workspace.state_threshold_spin.value(),
            bootstrap_iterations=int(
                workspace.bootstrap_iterations_combo.currentData()
            ),
            primary=primary,
            comparison=comparison,
            bootstrap_block_size=(
                workspace.bootstrap_block_size_spin.value() or None
            ),
            bootstrap_seed=workspace.bootstrap_seed,
            confidence_level=workspace.confidence_level,
        )

    def _save_project(self):
        if not self._require_result():
            return
        path, _ = QtWidgets.QFileDialog.getSaveFileName(
            self,
            "Save reproducible DockLens project",
            "analysis.docklens",
            "DockLens project (*.docklens)",
        )
        if not path:
            return
        if not path.lower().endswith(".docklens"):
            path += ".docklens"
        share_warning = QtWidgets.QMessageBox.question(
            self,
            "Save project provenance",
            "The project records absolute local source paths and SHA-256 "
            "digests so inputs can be checked when reopened. These paths may "
            "contain local user or folder names. Continue?",
            QtWidgets.QMessageBox.Yes | QtWidgets.QMessageBox.No,
            QtWidgets.QMessageBox.Yes,
        )
        if share_warning != QtWidgets.QMessageBox.Yes:
            return
        try:
            outputs = save_project(self._project_state(), path)
        except Exception:  # noqa: BLE001 - contain filesystem/codec errors
            LOGGER.exception("DockLens project save failed")
            QtWidgets.QMessageBox.critical(
                self,
                "Project could not be saved",
                "DockLens could not write a verified project. Check that the "
                "source files and destination remain accessible.",
            )
            return
        self.status.showMessage("Project saved: %s" % outputs[0])
        QtWidgets.QMessageBox.information(
            self,
            "Project saved",
            "The analysis and its methods record were saved:\n"
            + "\n".join(str(output) for output in outputs),
        )

    def _open_project(self):
        path, _ = QtWidgets.QFileDialog.getOpenFileName(
            self,
            "Open DockLens project",
            "",
            "DockLens project (*.docklens)",
        )
        if not path:
            return
        try:
            project = load_project(path)
            if project.primary.result is None:
                raise ValueError("project does not contain a cached result")
            stale_messages = validate_project_inputs(project)
        except (OSError, ValueError):
            LOGGER.exception("DockLens project load failed")
            QtWidgets.QMessageBox.critical(
                self,
                "Project could not be opened",
                "The project is invalid, incomplete or failed its integrity "
                "checks. The current analysis was not changed.",
            )
            return
        self._restore_project(project, stale_messages)

    @staticmethod
    def _set_combo_data(combo, value):
        index = combo.findData(value)
        if index >= 0:
            combo.setCurrentIndex(index)

    def _restore_project(self, project, stale_messages=()):
        """Apply a fully decoded project after all integrity checks succeed."""
        controls = (
            self.preset_combo,
            self.analysis_combo,
            self.mode_combo,
        )
        for control in controls:
            control.blockSignals(True)
        try:
            self._set_combo_data(self.preset_combo, project.hbond_preset)
            self._set_combo_data(
                self.analysis_combo, project.analysis_profile
            )
            self._set_combo_data(self.mode_combo, project.primary.mode)
        finally:
            for control in controls:
                control.blockSignals(False)
        self._files = [item.path for item in project.primary.inputs]
        self._result = project.primary.result
        self._comparison_result = (
            project.comparison.result
            if project.comparison is not None
            else None
        )
        self._comparison_files = (
            [item.path for item in project.comparison.inputs]
            if project.comparison is not None
            else []
        )
        self.key_edit.setText(" ".join(project.key_residues))
        selected_types = set(project.selected_types)
        for kind, checkbox in self.type_boxes.items():
            checkbox.setChecked(kind in selected_types)
        workspace = self.analytics_workspace
        workspace.state_threshold_spin.setValue(project.state_threshold)
        if project.primary.time_step_ns is not None:
            workspace.time_step_spin.setValue(project.primary.time_step_ns)
        index = workspace.bootstrap_iterations_combo.findData(
            project.bootstrap_iterations
        )
        if index >= 0:
            workspace.bootstrap_iterations_combo.setCurrentIndex(index)
        workspace.bootstrap_block_size_spin.setValue(
            project.bootstrap_block_size or 0
        )
        workspace.bootstrap_seed = project.bootstrap_seed
        workspace.confidence_level = project.confidence_level
        workspace.set_mode(project.primary.mode)
        workspace.set_observation_series(
            project.primary.observation_series,
            refresh=False,
            result=project.primary.result,
        )
        workspace.set_observation_series(
            (
                project.comparison.observation_series
                if project.comparison is not None
                else None
            ),
            comparison=True,
            refresh=False,
            result=(
                project.comparison.result
                if project.comparison is not None
                else None
            ),
        )
        workspace.system_a_role.setCurrentIndex(
            workspace.system_a_role.findData(project.primary.mode)
        )
        if project.comparison is not None:
            workspace.system_b_role.setCurrentIndex(
                workspace.system_b_role.findData(project.comparison.mode)
            )
        self._populate_residue_list()
        self._refresh_tables()
        workspace.select_residue(project.selected_residue)
        workspace_map = {
            name: index for index, name in enumerate(_WORKSPACES)
        }
        self.workspace_stack.setCurrentIndex(
            workspace_map.get(project.active_workspace, 0)
        )
        self._project_stale = bool(stale_messages)
        self.run_detection_button.setEnabled(not self._project_stale)
        self.dataset_context.setText(
            f"{len(self._result.summaries)} cached observation(s) · "
            f"{len(self._result.details)} raw interaction row(s)"
        )
        if stale_messages:
            self.status.showMessage(
                "Project opened from verified cached results; source files "
                "changed or are unavailable, so re-run is disabled."
            )
            QtWidgets.QMessageBox.warning(
                self,
                "Project sources changed",
                "The cached analysis passed internal integrity checks and was "
                "opened. Re-running is disabled because:\n"
                + "\n".join(stale_messages),
            )
        else:
            self.status.showMessage(
                "Verified DockLens project opened from cached results."
            )

    def _open_representative(self, observation_id):
        source_path = self._trusted_representative_path(observation_id)
        if source_path is None:
            QtWidgets.QMessageBox.information(
                self,
                "Representative source unavailable",
                "The representative remains identified in the project, but "
                "its declared external structure file is not available.",
            )
            return
        opened = QtGui.QDesktopServices.openUrl(
            QtCore.QUrl.fromLocalFile(str(source_path))
        )
        if opened:
            self.status.showMessage(
                "Opened representative %s with the system structure viewer. "
                "This is a medoid, not a 'best pose'." % observation_id
            )
        else:
            QtWidgets.QMessageBox.warning(
                self,
                "Could not open representative",
                "Associate the structure format with PyMOL or another "
                "molecular viewer and try again.",
            )

    def _trusted_representative_path(self, observation_id):
        """Return a declared local structure path, without opening arbitrary data."""
        if self._result is None:
            return None
        source_path = next(
            (
                summary.source_path
                for summary in self._result.summaries
                if summary.pose_id == observation_id and summary.source_path
            ),
            "",
        )
        if not source_path:
            source_path = next(
                (
                    detail.source_path
                    for detail in self._result.details
                    if detail.pose_id == observation_id and detail.source_path
                ),
                "",
            )
        if not source_path:
            return None
        if not is_local_filesystem_path(source_path):
            return None
        candidate = Path(source_path).expanduser()
        if candidate.suffix.lower() not in {".mol2", ".pdb", ".pdbqt"}:
            return None
        candidate = candidate.resolve()
        declared = tuple(self._files) + tuple(self._comparison_files)
        trusted = False
        for value in declared:
            root = Path(value).expanduser()
            if root.is_dir():
                try:
                    candidate.relative_to(root.resolve())
                    trusted = True
                except ValueError:
                    pass
            elif root.resolve() == candidate:
                trusted = True
            if trusted:
                break
        return candidate if trusted and candidate.is_file() else None


__all__ = ["ProjectControllerMixin"]
