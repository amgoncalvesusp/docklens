"""Qt widgets for the DockLens 1.0 analytical workspaces."""

from __future__ import annotations

from PyQt5 import QtCore, QtWidgets
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg

from .analytics import (
    AnalysisContext,
    episode_statistics,
    fingerprint_clusters,
    fingerprint_matrix,
    fingerprint_similarity,
    residue_type_prevalence,
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
    artifactChanged = QtCore.pyqtSignal(object)

    def __init__(self, parent=None):
        super().__init__(parent)
        self._result: RunResult = make_result()
        self._comparison: RunResult | None = None
        self._mode = "docking"
        self._selected_residue = ""
        self._fingerprint_matrix = fingerprint_matrix(self._result)
        self._fingerprint_similarity = fingerprint_similarity(
            self._fingerprint_matrix
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
        layout.addWidget(self.fingerprint_status)
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
        return _scrollable(content)

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
        layout.addLayout(mode_row)
        self.compare_panel = ChartPanel(minimum_height=420)
        self.compare_panel.setObjectName("analysisField")
        layout.addWidget(self.compare_panel)
        self.retention_panel = ChartPanel(minimum_height=330)
        self.retention_panel.setObjectName("analysisField")
        layout.addWidget(self.retention_panel)
        return _scrollable(content)

    @property
    def selected_residue(self):
        return self._selected_residue

    def set_result(self, result: RunResult | None):
        self._result = result or make_result()
        self.refresh()

    def set_comparison(self, result: RunResult | None):
        self._comparison = result
        self.refresh_compare()

    def set_mode(self, mode: str):
        AnalysisContext(mode=mode)
        self._mode = mode
        for selector in (self.system_a_role, self.system_b_role):
            selector.blockSignals(True)
            selector.setCurrentIndex(selector.findData(mode))
            selector.blockSignals(False)
        self.refresh()

    def select_residue(self, residue: str):
        residue = str(residue or "")
        if residue == self._selected_residue:
            return
        self._selected_residue = residue
        self.residueSelected.emit(residue)

    def refresh(self):
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

    def refresh_compare(self):
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

    def current_artifact(self, workspace_index):
        if workspace_index == 0:
            return self.residue_panel.artifact
        if workspace_index == 1:
            return self.fingerprint_panel.artifact
        if workspace_index == 2:
            if self.compare_analysis_combo.currentData() == "retention":
                return self.retention_panel.artifact
            return self.compare_panel.artifact
        return None

    def dispose(self):
        """Release every native Matplotlib canvas before Qt application exit."""
        for panel in (
            self.residue_panel,
            self.residue_barcode_panel,
            self.fingerprint_panel,
            self.similarity_panel,
            self.compare_panel,
            self.retention_panel,
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
                    self._result, AnalysisContext(mode="md")
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
