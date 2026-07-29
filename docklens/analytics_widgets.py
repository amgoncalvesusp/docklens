"""Qt widgets for the DockLens 1.0 analytical workspaces."""

from __future__ import annotations

import pandas as pd
from PyQt5 import QtCore, QtWidgets
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg

from .analysis_tasks import AnalysisTaskRunner
from .analytics import (
    AnalysisContext,
    episode_statistics,
    fingerprint_clusters,
    fingerprint_matrix,
    fingerprint_similarity,
    residue_type_prevalence,
)
from .dynamic_plotting import (
    build_difference_uncertainty_chart,
    build_state_population_chart,
    build_state_timeline_chart,
    build_transition_chart,
    build_uncertainty_chart,
)
from .dynamic_states import (
    interaction_state_analysis,
    state_summary_frame,
)
from .dynamic_ui import build_dynamic_analysis
from .ligand_selection import (
    default_md_series_for_result,
    ligand_groups,
    subset_observation_series,
    subset_run_result,
)
from .observation_series import (
    ObservationSeries,
    observation_series_from_dataframe,
)
from .plotting import (
    ChartArtifact,
    build_comparison_chart,
    build_fingerprint_chart,
    build_residue_chart,
    build_retention_chart,
    build_similarity_chart,
)
from .results import RunResult, make_result
from .uncertainty import (
    block_bootstrap_difference,
    block_bootstrap_occupancy,
)


_MAX_SIMILARITY_OBSERVATIONS = 300

class ChartPanel(QtWidgets.QWidget):
    """Replaceable Matplotlib canvas with a safe empty state."""

    def __init__(self, minimum_height=360, parent=None):
        super().__init__(parent)
        self.artifact: ChartArtifact | None = None
        self._canvas = None
        self._layout = QtWidgets.QVBoxLayout(self)
        self._layout.setContentsMargins(0, 0, 0, 0)
        self.setMinimumHeight(minimum_height)

    def set_artifact(self, artifact: ChartArtifact):
        self._dispose_canvas()
        self.artifact = artifact
        self._canvas = FigureCanvasQTAgg(artifact.figure)
        self._canvas.setFocusPolicy(QtCore.Qt.StrongFocus)
        self._canvas.setMinimumHeight(self.minimumHeight())
        self._layout.addWidget(self._canvas)
        self._canvas.draw()

    def _dispose_canvas(self):
        canvas = self._canvas
        if canvas is None:
            return
        self._canvas = None
        self._layout.removeWidget(canvas)
        canvas.close()
        canvas.setParent(None)
        canvas.figure.set_canvas(None)
        canvas.deleteLater()

    def closeEvent(self, event):
        self._dispose_canvas()
        super().closeEvent(event)

def _heading(title, description):
    widget = QtWidgets.QWidget()
    layout = QtWidgets.QVBoxLayout(widget)
    layout.setContentsMargins(0, 0, 0, 4)
    title_label = QtWidgets.QLabel(title)
    title_label.setObjectName("workspaceTitle")
    description_label = QtWidgets.QLabel(description)
    description_label.setObjectName("workspaceDescription")
    description_label.setWordWrap(True)
    layout.addWidget(title_label)
    layout.addWidget(description_label)
    return widget

def _scrollable(content):
    area = QtWidgets.QScrollArea()
    area.setObjectName("workspaceScroll")
    area.setWidgetResizable(True)
    area.setFrameShape(QtWidgets.QFrame.NoFrame)
    area.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarAsNeeded)
    area.setVerticalScrollBarPolicy(QtCore.Qt.ScrollBarAsNeeded)
    area.setWidget(content)
    return area


class AnalyticsWorkspace(QtCore.QObject):
    """Own the three analytical pages and their synchronized selection."""

    residueSelected = QtCore.pyqtSignal(str)
    representativeSelected = QtCore.pyqtSignal(str)
    artifactChanged = QtCore.pyqtSignal(object)

    def __init__(self, parent=None):
        super().__init__(parent)
        self._result: RunResult = make_result()
        self._source_result: RunResult = self._result
        self._tasks = AnalysisTaskRunner(self)
        self._comparison: RunResult | None = None
        self._source_comparison: RunResult | None = None
        self._ligand_group: str | None = None
        self._comparison_ligand_group: str | None = None
        self._series: ObservationSeries | None = None
        self._comparison_series: ObservationSeries | None = None
        self._mode = "docking"
        self.bootstrap_seed = 2026
        self.confidence_level = 0.95
        self._selected_residue = ""
        self._fingerprint_matrix = fingerprint_matrix(self._result)
        self._fingerprint_similarity = fingerprint_similarity(
            self._fingerprint_matrix
        )
        self.state_analysis = interaction_state_analysis(
            self._fingerprint_matrix,
            AnalysisContext(mode=self._mode),
        )
        self.residue_page = self._build_residue_page()
        self.fingerprint_page = self._build_fingerprint_page()
        self.compare_page = self._build_compare_page()
        self.scroll_areas = (
            self.residue_page,
            self.fingerprint_page,
            self.compare_page,
        )
        self.refresh()

    def _build_residue_page(self):
        content = QtWidgets.QWidget()
        content.setMinimumWidth(780)
        layout = QtWidgets.QVBoxLayout(content)
        layout.setContentsMargins(24, 20, 24, 24)
        layout.setSpacing(16)
        layout.addWidget(
            _heading(
                "Interaction profile by residue",
                "Each channel counts once per pose or saved frame. Select a "
                "residue to keep charts, evidence and tables in context.",
            )
        )
        controls = QtWidgets.QHBoxLayout()
        controls.addWidget(QtWidgets.QLabel("DockLens Lens:"))
        self.residue_selector = QtWidgets.QComboBox()
        self.residue_selector.setMinimumWidth(180)
        self.residue_selector.currentTextChanged.connect(self.select_residue)
        controls.addWidget(self.residue_selector)
        self.residue_metric = QtWidgets.QLabel("No interaction evidence loaded.")
        self.residue_metric.setObjectName("metricStrip")
        controls.addWidget(self.residue_metric)
        controls.addStretch(1)
        layout.addLayout(controls)
        self.residue_panel = ChartPanel(minimum_height=390)
        self.residue_panel.setObjectName("analysisField")
        layout.addWidget(self.residue_panel)
        self.residue_barcode_panel = ChartPanel(minimum_height=360)
        self.residue_barcode_panel.setObjectName("analysisField")
        layout.addWidget(self.residue_barcode_panel)
        return _scrollable(content)

    def _build_fingerprint_page(self):
        content = QtWidgets.QWidget()
        content.setMinimumWidth(800)
        layout = QtWidgets.QVBoxLayout(content)
        layout.setContentsMargins(24, 20, 24, 24)
        layout.setSpacing(16)
        layout.addWidget(
            _heading(
                "Interaction fingerprints",
                "Compare pose or frame states with binary residue × "
                "interaction features and Jaccard/Tanimoto similarity.",
            )
        )
        self.fingerprint_status = QtWidgets.QLabel()
        self.fingerprint_status.setObjectName("metricStrip")
        self.fingerprint_status.setWordWrap(True)
        export_row = QtWidgets.QHBoxLayout()
        export_row.addWidget(self.fingerprint_status, 1)
        export_row.addWidget(QtWidgets.QLabel("Export view:"))
        self.fingerprint_export_combo = QtWidgets.QComboBox()
        for label, value in (
            ("Interaction fingerprint", "fingerprint"),
            ("Similarity matrix", "similarity"),
            ("State populations", "population"),
            ("State timeline / pose order", "timeline"),
            ("State transitions", "transitions"),
            ("Confidence intervals", "confidence"),
        ):
            self.fingerprint_export_combo.addItem(label, value)
        export_row.addWidget(self.fingerprint_export_combo)
        layout.addLayout(export_row)
        self.fingerprint_panel = ChartPanel(minimum_height=390)
        self.fingerprint_panel.setObjectName("analysisField")
        layout.addWidget(self.fingerprint_panel)
        lower = QtWidgets.QSplitter(QtCore.Qt.Horizontal)
        self.similarity_panel = ChartPanel(minimum_height=360)
        self.similarity_panel.setObjectName("analysisField")
        lower.addWidget(self.similarity_panel)
        cluster_box = QtWidgets.QWidget()
        cluster_layout = QtWidgets.QVBoxLayout(cluster_box)
        cluster_title = QtWidgets.QLabel("Threshold clusters and medoids")
        cluster_title.setObjectName("sectionTitle")
        cluster_layout.addWidget(cluster_title)
        self.cluster_table = QtWidgets.QTableWidget(0, 4)
        self.cluster_table.setHorizontalHeaderLabels(
            ("Cluster", "Size", "Medoid", "Mean similarity")
        )
        self.cluster_table.setEditTriggers(
            QtWidgets.QAbstractItemView.NoEditTriggers
        )
        self.cluster_table.setSelectionBehavior(
            QtWidgets.QAbstractItemView.SelectRows
        )
        self.cluster_table.horizontalHeader().setStretchLastSection(True)
        cluster_layout.addWidget(self.cluster_table)
        lower.addWidget(cluster_box)
        lower.setSizes((520, 360))
        layout.addWidget(lower)
        layout.addWidget(self._build_dynamic_analysis())
        return _scrollable(content)

    def _build_dynamic_analysis(self):
        return build_dynamic_analysis(self, ChartPanel)

    def _build_compare_page(self):
        content = QtWidgets.QWidget()
        content.setMinimumWidth(800)
        layout = QtWidgets.QVBoxLayout(content)
        layout.setContentsMargins(24, 20, 24, 24)
        layout.setSpacing(16)
        layout.addWidget(
            _heading(
                "Compare interaction evidence",
                "System B minus System A uses independent denominators. "
                "Docking→MD retention remains a separate, explicit analysis.",
            )
        )
        self.compare_status = QtWidgets.QLabel(
            "Load a comparison dataset to activate differential analyses."
        )
        self.compare_status.setObjectName("metricStrip")
        self.compare_status.setWordWrap(True)
        mode_row = QtWidgets.QHBoxLayout()
        mode_row.addWidget(self.compare_status, 1)
        self.compare_analysis_combo = QtWidgets.QComboBox()
        self.compare_analysis_combo.addItem(
            "System A vs B (same evidence mode)", "difference"
        )
        self.compare_analysis_combo.addItem(
            "Docking A → MD B retention", "retention"
        )
        self.compare_analysis_combo.currentIndexChanged.connect(
            self.refresh_compare
        )
        mode_row.addWidget(self.compare_analysis_combo)
        mode_row.addWidget(QtWidgets.QLabel("A role:"))
        self.system_a_role = QtWidgets.QComboBox()
        self.system_a_role.addItem("Docking", "docking")
        self.system_a_role.addItem("MD frames", "md")
        self.system_a_role.currentIndexChanged.connect(
            self._comparison_roles_changed
        )
        mode_row.addWidget(self.system_a_role)
        mode_row.addWidget(QtWidgets.QLabel("B role:"))
        self.system_b_role = QtWidgets.QComboBox()
        self.system_b_role.addItem("Docking", "docking")
        self.system_b_role.addItem("MD frames", "md")
        self.system_b_role.currentIndexChanged.connect(
            self._comparison_roles_changed
        )
        mode_row.addWidget(self.system_b_role)
        self.load_comparison_trajectory_map_button = QtWidgets.QPushButton(
            "Load B trajectory map"
        )
        self.load_comparison_trajectory_map_button.clicked.connect(
            lambda: self._load_trajectory_map(comparison=True)
        )
        mode_row.addWidget(self.load_comparison_trajectory_map_button)
        self.compare_confidence_button = QtWidgets.QPushButton(
            "Compute Δ confidence intervals"
        )
        self.compare_confidence_button.clicked.connect(
            self._compute_comparison_uncertainty
        )
        mode_row.addWidget(self.compare_confidence_button)
        layout.addLayout(mode_row)
        self.compare_panel = ChartPanel(minimum_height=420)
        self.compare_panel.setObjectName("analysisField")
        layout.addWidget(self.compare_panel)
        self.retention_panel = ChartPanel(minimum_height=330)
        self.retention_panel.setObjectName("analysisField")
        layout.addWidget(self.retention_panel)
        self.compare_uncertainty_panel = ChartPanel(minimum_height=380)
        self.compare_uncertainty_panel.setObjectName("analysisField")
        self.compare_uncertainty_panel.setVisible(False)
        layout.addWidget(self.compare_uncertainty_panel)
        return _scrollable(content)

    @property
    def selected_residue(self):
        return self._selected_residue

    def set_result(self, result: RunResult | None):
        self._source_result = result or make_result()
        valid_groups = {
            group.key for group in ligand_groups(self._source_result)
        }
        if self._ligand_group not in valid_groups:
            self._ligand_group = None
        self._result = subset_run_result(
            self._source_result, self._ligand_group
        )
        self.refresh()

    def set_comparison(self, result: RunResult | None):
        self._source_comparison = result
        if result is None:
            self._comparison_series = None
            self._comparison_ligand_group = None
            self._comparison = None
        else:
            valid_groups = {group.key for group in ligand_groups(result)}
            if self._comparison_ligand_group not in valid_groups:
                self._comparison_ligand_group = None
            self._comparison = subset_run_result(
                result, self._comparison_ligand_group
            )
        self.refresh_compare()

    def set_ligand_group(
        self,
        group_key: str | None,
        *,
        comparison: bool = False,
    ):
        if comparison:
            if self._source_comparison is None:
                self._comparison_ligand_group = None
                self._comparison = None
            else:
                self._comparison = subset_run_result(
                    self._source_comparison, group_key
                )
                self._comparison_ligand_group = group_key
            self.refresh_compare()
            return
        self._result = subset_run_result(self._source_result, group_key)
        self._ligand_group = group_key
        self.refresh()

    def set_observation_series(
        self,
        series: ObservationSeries | None,
        *,
        comparison: bool = False,
        refresh: bool = True,
        result: RunResult | None = None,
    ):
        if series is not None and series.mode != "md":
            raise ValueError("trajectory maps must use MD observation series")
        if series is not None:
            target_result = (
                result
                if result is not None
                else (
                    self._source_comparison
                    if comparison
                    else self._source_result
                )
            )
            matrix = fingerprint_matrix(target_result)
            if set(series.observation_ids) != set(matrix.index):
                raise ValueError(
                    "trajectory map must match the analyzed observation IDs"
                )
        if comparison:
            self._comparison_series = series
        else:
            self._series = series
        if refresh:
            if comparison:
                self.refresh_compare()
            else:
                self._refresh_states()

    def clear_observation_series(self):
        self._series = None
        self._comparison_series = None

    def _load_trajectory_map(self, *, comparison: bool):
        dataset_mode = (
            self.system_b_role.currentData() if comparison else self._mode
        )
        if dataset_mode != "md":
            QtWidgets.QMessageBox.information(
                None,
                "Trajectory map is for MD",
                "Change this dataset role to MD frames before loading a map.",
            )
            return
        path, _ = QtWidgets.QFileDialog.getOpenFileName(
            None,
            "Load trajectory map",
            "",
            "CSV trajectory map (*.csv)",
        )
        if not path:
            return
        try:
            source = QtCore.QFileInfo(path)
            if source.size() > 32 * 1024 * 1024:
                raise ValueError("trajectory map exceeds 32 MB")
            matrix = fingerprint_matrix(
                self._source_comparison
                if comparison
                else self._source_result
            )
            frame = pd.read_csv(path)
            if len(frame) > 2_000_000:
                raise ValueError("trajectory map has too many rows")
            series = observation_series_from_dataframe(
                frame,
                expected_ids=matrix.index,
                default_time_step_ns=self.time_step_spin.value(),
            )
            self.set_observation_series(
                series,
                comparison=comparison,
            )
        except (OSError, UnicodeError, ValueError, pd.errors.ParserError):
            QtWidgets.QMessageBox.warning(
                None,
                "Trajectory map could not be loaded",
                "Use a CSV containing exactly one row per analyzed "
                "observation and the columns observation_id, replica_id, "
                "frame_index and optional time_ns.",
            )

    def set_mode(self, mode: str):
        AnalysisContext(mode=mode)
        self._mode = mode
        for selector in (self.system_a_role, self.system_b_role):
            selector.blockSignals(True)
            selector.setCurrentIndex(selector.findData(mode))
            selector.blockSignals(False)
        self.time_step_spin.setEnabled(mode == "md")
        self.load_trajectory_map_button.setEnabled(mode == "md")
        self.refresh()

    def _active_md_series(self, *, comparison=False):
        result = self._comparison if comparison else self._result
        if result is None:
            return None
        explicit = (
            self._comparison_series if comparison else self._series
        )
        observation_ids = fingerprint_matrix(result).index
        if explicit is not None:
            return subset_observation_series(explicit, observation_ids)
        return default_md_series_for_result(
            result,
            time_step_ns=self.time_step_spin.value(),
        )

    def select_residue(self, residue: str):
        residue = str(residue or "")
        if residue == self._selected_residue:
            return
        self._selected_residue = residue
        self.residueSelected.emit(residue)

    def refresh(self):
        self._tasks.invalidate()
        prevalence = residue_type_prevalence(self._result)
        self._fingerprint_matrix = fingerprint_matrix(self._result)
        similarity_matrix = self._fingerprint_matrix
        if len(similarity_matrix.index) > _MAX_SIMILARITY_OBSERVATIONS:
            last = len(similarity_matrix.index) - 1
            indices = tuple(
                round(
                    position * last
                    / (_MAX_SIMILARITY_OBSERVATIONS - 1)
                )
                for position in range(_MAX_SIMILARITY_OBSERVATIONS)
            )
            similarity_matrix = similarity_matrix.iloc[list(indices)]
            self.fingerprint_status.setText(
                f"Barcode: {len(self._fingerprint_matrix.index)} observations · "
                f"similarity/clustering: {_MAX_SIMILARITY_OBSERVATIONS} "
                "evenly sampled observations"
            )
        else:
            self.fingerprint_status.setText(
                f"{len(self._fingerprint_matrix.index)} observations · "
                "similarity and clustering use the complete matrix"
            )
        self._similarity_source_matrix = similarity_matrix
        self._fingerprint_similarity = fingerprint_similarity(
            self._similarity_source_matrix
        )
        residue_artifact = build_residue_chart(
            self._result, mode=self._mode, prevalence=prevalence
        )
        barcode_artifact = build_fingerprint_chart(
            self._result,
            mode=self._mode,
            matrix=self._fingerprint_matrix,
        )
        self.residue_panel.set_artifact(residue_artifact)
        self.residue_barcode_panel.set_artifact(barcode_artifact)
        self.fingerprint_panel.set_artifact(
            build_fingerprint_chart(
                self._result,
                mode=self._mode,
                matrix=self._fingerprint_matrix,
            )
        )
        self.similarity_panel.set_artifact(
            build_similarity_chart(
                self._result,
                matrix=self._similarity_source_matrix,
                similarity=self._fingerprint_similarity,
            )
        )
        self._refresh_residue_selector()
        self._refresh_clusters()
        self._refresh_states()
        self.refresh_compare()
        self.artifactChanged.emit(residue_artifact)

    def _refresh_residue_selector(self):
        residues = sorted(
            set(self.residue_panel.artifact.data["receptor_residue"])
            if not self.residue_panel.artifact.data.empty
            else ()
        )
        current = self._selected_residue
        self.residue_selector.blockSignals(True)
        self.residue_selector.clear()
        self.residue_selector.addItems(residues)
        if current in residues:
            self.residue_selector.setCurrentText(current)
        self.residue_selector.blockSignals(False)
        selected = current if current in residues else (residues[0] if residues else "")
        if selected != self._selected_residue:
            self._selected_residue = selected
            self.residueSelected.emit(selected)
        data = self.residue_panel.artifact.data
        total = self.residue_panel.artifact.metadata["total_observations"]
        measure = "poses" if self._mode == "docking" else "saved frames"
        self.residue_metric.setText(
            f"{len(residues)} residues · {total} {measure} · "
            f"{len(data)} residue/type channels"
        )

    def _refresh_clusters(self):
        matrix_artifact = self.fingerprint_panel.artifact
        if matrix_artifact.data.empty:
            clusters = ()
        else:
            clusters = fingerprint_clusters(
                self._similarity_source_matrix,
                # Clustering follows the same disclosed sample as similarity.
                threshold=0.65,
                similarity=self._fingerprint_similarity,
            )
        self.cluster_table.setRowCount(len(clusters))
        for row, cluster in enumerate(clusters):
            values = (
                cluster.cluster_id,
                len(cluster.members),
                cluster.medoid,
                f"{cluster.mean_similarity:.3f}",
            )
            for column, value in enumerate(values):
                self.cluster_table.setItem(
                    row, column, QtWidgets.QTableWidgetItem(str(value))
                )

    def _refresh_states(self):
        if not hasattr(self, "state_threshold_spin"):
            return
        context = AnalysisContext(
            mode=self._mode,
            time_step_ns=self.time_step_spin.value(),
        )
        self.state_analysis = interaction_state_analysis(
            self._fingerprint_matrix,
            context,
            threshold=self.state_threshold_spin.value(),
            max_training_observations=_MAX_SIMILARITY_OBSERVATIONS,
            series=(
                self._active_md_series()
                if self._mode == "md"
                else None
            ),
        )
        self.state_population_panel.set_artifact(
            build_state_population_chart(self.state_analysis)
        )
        self.state_timeline_panel.set_artifact(
            build_state_timeline_chart(self.state_analysis)
        )
        temporal = self._mode == "md"
        self.dynamic_group_box.setTitle(
            "Dynamic interaction states" if temporal else "Pose families"
        )
        self.dynamic_tabs.setTabText(1, "Timeline" if temporal else "Pose order")
        header = self.state_table.horizontalHeaderItem(5)
        if header is not None:
            header.setText("Mean dwell" if temporal else "Temporal metric")
        self.dynamic_tabs.setTabEnabled(2, temporal)
        self.dynamic_tabs.setTabEnabled(3, temporal)
        self.compute_uncertainty_button.setEnabled(
            temporal and len(self._fingerprint_matrix.index) > 0
        )
        for data in ("transitions", "confidence"):
            item = self.fingerprint_export_combo.model().item(
                self.fingerprint_export_combo.findData(data)
            )
            if item is not None:
                item.setEnabled(temporal)
        if (
            not temporal
            and self.fingerprint_export_combo.currentData()
            in {"transitions", "confidence"}
        ):
            self.fingerprint_export_combo.setCurrentIndex(
                self.fingerprint_export_combo.findData("fingerprint")
            )
        if temporal:
            self.transition_panel.set_artifact(
                build_transition_chart(self.state_analysis)
            )
        self.uncertainty_panel.set_artifact(
            build_uncertainty_chart(pd.DataFrame())
        )
        self._refresh_state_table()
        term = (
            "pose families"
            if self._mode == "docking"
            else "interaction states"
        )
        sampling = (
            f"trained on {self.state_analysis.training_observations} evenly "
            f"sampled observations; {len(self.state_analysis.outlier_observations)} "
            "outliers"
            if self.state_analysis.sampled
            else "all observations used"
        )
        temporal_note = (
            "Transitions are descriptive and preserve frame gaps and replica "
            "boundaries."
            if temporal
            else "Docking groups are pose families; no temporal interpretation."
        )
        if temporal:
            active_series = self.state_analysis.series
            replicas = len(
                {point.replica_id for point in active_series.points}
            )
            gaps = len(active_series.points) - replicas - len(
                active_series.transition_pairs()
            )
            map_note = (
                f"Explicit trajectory map: {replicas} "
                f"{'replica' if replicas == 1 else 'replicas'}, "
                f"{max(gaps, 0)} frame gap(s)."
                if self._series is not None
                else f"Automatic file boundaries: {replicas} "
                f"{'series' if replicas == 1 else 'series'}; load a "
                "trajectory map to declare replicas or gaps explicitly."
            )
        else:
            map_note = ""
        self.state_status.setText(
            f"{len(self.state_analysis.states)} {term} · complete-link "
            f"threshold {self.state_analysis.threshold:.2f} · {sampling}. "
            + temporal_note
            + (" " + map_note if map_note else "")
        )

    def _refresh_state_table(self):
        frame = state_summary_frame(self.state_analysis)
        self.state_table.setRowCount(len(frame))
        for row_index, row in enumerate(frame.itertuples(index=False)):
            features = "; ".join(
                f"{residue} · {kind}"
                for residue, kind in row.characteristic_features
            )
            mean_dwell = (
                "not applicable to docking"
                if row.mean_dwell_observations is None
                else f"{row.mean_dwell_observations:.2f} frames"
            )
            values = (
                row.state_id,
                row.population_count,
                f"{row.population_pct:.1f}",
                row.representative,
                f"{row.mean_similarity:.3f}",
                mean_dwell,
                features or "No consensus contact",
            )
            for column, value in enumerate(values):
                self.state_table.setItem(
                    row_index,
                    column,
                    QtWidgets.QTableWidgetItem(str(value)),
                )
        self.open_representative_button.setEnabled(False)

    def _state_selection_changed(self):
        self.open_representative_button.setEnabled(
            bool(self.state_table.selectedItems())
        )

    def _emit_representative(self):
        rows = self.state_table.selectionModel().selectedRows()
        if not rows:
            return
        item = self.state_table.item(rows[0].row(), 3)
        if item is not None and item.text():
            self.representativeSelected.emit(item.text())

    def _compute_uncertainty(self):
        if self._mode != "md" or self._fingerprint_matrix.empty:
            return
        self.compute_uncertainty_button.setEnabled(False)
        self.compute_uncertainty_button.setText("Computing…")
        self._tasks.start(
            block_bootstrap_occupancy,
            self._fingerprint_matrix.copy(deep=True),
            AnalysisContext(
                mode="md",
                time_step_ns=self.time_step_spin.value(),
            ),
            iterations=int(self.bootstrap_iterations_combo.currentData()),
            block_size=(
                self.bootstrap_block_size_spin.value() or None
            ),
            seed=self.bootstrap_seed,
            confidence_level=self.confidence_level,
            series=self._active_md_series(),
            on_success=self._uncertainty_ready,
            on_error=self._uncertainty_failed,
            on_settled=self._uncertainty_settled,
        )

    def _uncertainty_ready(self, intervals):
        self.uncertainty_panel.set_artifact(
            build_uncertainty_chart(intervals)
        )
        self.dynamic_tabs.setCurrentIndex(3)

    def _uncertainty_failed(self, _error_name):
        self.state_status.setText(
            "Confidence intervals could not be computed for this selection."
        )

    def _uncertainty_settled(self):
        self.compute_uncertainty_button.setText(
            "Compute confidence intervals"
        )
        self.compute_uncertainty_button.setEnabled(self._mode == "md")

    def refresh_compare(self):
        self._tasks.invalidate()
        comparison = self._comparison or make_result()
        compare_artifact = build_comparison_chart(
            self._result, comparison, mode=self._mode
        )
        self.compare_panel.set_artifact(compare_artifact)
        retention_artifact = build_retention_chart(
            self._result, comparison
        )
        self.retention_panel.set_artifact(retention_artifact)
        roles_allow_retention = (
            self._comparison is not None
            and self.system_a_role.currentData() == "docking"
            and self.system_b_role.currentData() == "md"
        )
        self.compare_confidence_button.setEnabled(
            self._comparison is not None
            and self._mode == "md"
            and self.compare_analysis_combo.currentData() == "difference"
        )
        self.load_comparison_trajectory_map_button.setEnabled(
            self._comparison is not None
            and self.system_b_role.currentData() == "md"
        )
        self.compare_uncertainty_panel.setVisible(False)
        retention_item = self.compare_analysis_combo.model().item(1)
        if retention_item is not None:
            retention_item.setEnabled(roles_allow_retention)
        if (
            self.compare_analysis_combo.currentData() == "retention"
            and not roles_allow_retention
        ):
            self.compare_analysis_combo.blockSignals(True)
            self.compare_analysis_combo.setCurrentIndex(0)
            self.compare_analysis_combo.blockSignals(False)
        retention_mode = (
            self.compare_analysis_combo.currentData() == "retention"
        )
        self.compare_panel.setVisible(not retention_mode)
        self.retention_panel.setVisible(retention_mode)
        if self._comparison is None:
            self.compare_status.setText(
                "Load a comparison dataset to activate differential analyses."
            )
        elif retention_mode:
            self.compare_status.setText(
                "System A is interpreted as docking poses; System B as saved "
                "MD frames. Retained threshold: 50% occupancy."
            )
        else:
            a_count = compare_artifact.metadata["system_a_observations"]
            b_count = compare_artifact.metadata["system_b_observations"]
            measure = (
                "docking prevalence"
                if self._mode == "docking"
                else "saved-frame occupancy"
            )
            self.compare_status.setText(
                f"System A: {a_count} observations · System B: "
                f"{b_count} observations · {measure} · Δ = B − A"
            )

    def _comparison_roles_changed(self):
        self.refresh_compare()

    def _compute_comparison_uncertainty(self):
        if self._comparison is None or self._mode != "md":
            return
        self.compare_confidence_button.setEnabled(False)
        self.compare_confidence_button.setText("Computing Δ…")
        self._tasks.start(
            block_bootstrap_difference,
            fingerprint_matrix(self._result).copy(deep=True),
            fingerprint_matrix(self._comparison).copy(deep=True),
            AnalysisContext(
                mode="md",
                time_step_ns=self.time_step_spin.value(),
            ),
            iterations=int(self.bootstrap_iterations_combo.currentData()),
            block_size_a=(
                self.bootstrap_block_size_spin.value() or None
            ),
            block_size_b=(
                self.bootstrap_block_size_spin.value() or None
            ),
            seed=self.bootstrap_seed,
            confidence_level=self.confidence_level,
            series_a=self._active_md_series(),
            series_b=self._active_md_series(comparison=True),
            on_success=self._comparison_uncertainty_ready,
            on_error=self._comparison_uncertainty_failed,
            on_settled=self._comparison_uncertainty_settled,
        )

    def _comparison_uncertainty_ready(self, intervals):
        if (
            self._comparison is None
            or self._mode != "md"
            or self.compare_analysis_combo.currentData() != "difference"
        ):
            return
        self.compare_uncertainty_panel.set_artifact(
            build_difference_uncertainty_chart(intervals)
        )
        self.compare_uncertainty_panel.setVisible(True)

    def _comparison_uncertainty_failed(self, _error_name):
        self.compare_status.setText(
            "Comparison confidence intervals could not be computed."
        )

    def _comparison_uncertainty_settled(self):
        self.compare_confidence_button.setText(
            "Compute Δ confidence intervals"
        )
        self.compare_confidence_button.setEnabled(
            self._comparison is not None
            and self._mode == "md"
            and self.compare_analysis_combo.currentData() == "difference"
        )

    def current_artifact(self, workspace_index):
        if workspace_index == 0:
            return self.residue_panel.artifact
        if workspace_index == 1:
            panels = {
                "fingerprint": self.fingerprint_panel,
                "similarity": self.similarity_panel,
                "population": self.state_population_panel,
                "timeline": self.state_timeline_panel,
                "transitions": self.transition_panel,
                "confidence": self.uncertainty_panel,
            }
            panel = panels.get(
                self.fingerprint_export_combo.currentData(),
                self.fingerprint_panel,
            )
            return panel.artifact or self.fingerprint_panel.artifact
        if workspace_index == 2:
            if (
                self.compare_uncertainty_panel.isVisible()
                and self._mode == "md"
                and self.compare_analysis_combo.currentData() == "difference"
            ):
                return self.compare_uncertainty_panel.artifact
            if self.compare_analysis_combo.currentData() == "retention":
                return self.retention_panel.artifact
            return self.compare_panel.artifact
        return None

    def dispose(self):
        """Release every native Matplotlib canvas before Qt application exit."""
        self._tasks.wait_for_done()
        for panel in (
            self.residue_panel,
            self.residue_barcode_panel,
            self.fingerprint_panel,
            self.similarity_panel,
            self.compare_panel,
            self.retention_panel,
            self.state_population_panel,
            self.state_timeline_panel,
            self.transition_panel,
            self.uncertainty_panel,
            self.compare_uncertainty_panel,
        ):
            panel._dispose_canvas()

    def evidence_text(self, residue):
        if not residue:
            return "No residue selected."
        data = self.residue_panel.artifact.data
        rows = data[data["receptor_residue"] == residue]
        if rows.empty:
            return f"{residue}: no evidence in the current profile."
        metric = "frequency" if self._mode == "docking" else "occupancy"
        lines = [f"{residue} · {metric} by interaction channel"]
        for row in rows.sort_values(
            "prevalence_pct", ascending=False
        ).itertuples():
            lines.append(
                f"{row.interaction_type}: {row.prevalence_pct:.1f}% "
                f"({row.observation_count}/{row.total_observations})"
            )
        if self._mode == "md":
            metrics = [
                item
                for item in episode_statistics(
                    self._result,
                    AnalysisContext(mode="md"),
                    series=self._active_md_series(),
                )
                if item.receptor_residue == residue
            ]
            if metrics:
                longest = max(
                    metrics, key=lambda item: item.longest_episode_observations
                )
                lines.append(
                    f"Longest episode: {longest.longest_episode_observations} "
                    "saved frames"
                )
        return "\n".join(lines)


__all__ = ["AnalyticsWorkspace", "ChartPanel"]
