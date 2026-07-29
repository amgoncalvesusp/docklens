"""Declarative controls for pose-family and dynamic-state analysis."""

from __future__ import annotations

from PyQt5 import QtCore, QtWidgets


def build_dynamic_analysis(owner, panel_type):
    """Build the advanced fingerprint section owned by AnalyticsWorkspace."""
    box = QtWidgets.QGroupBox("Pose families")
    owner.dynamic_group_box = box
    layout = QtWidgets.QVBoxLayout(box)
    controls = QtWidgets.QHBoxLayout()
    controls.addWidget(QtWidgets.QLabel("Similarity threshold:"))
    owner.state_threshold_spin = QtWidgets.QDoubleSpinBox()
    owner.state_threshold_spin.setRange(0.0, 1.0)
    owner.state_threshold_spin.setSingleStep(0.05)
    owner.state_threshold_spin.setDecimals(2)
    owner.state_threshold_spin.setValue(0.65)
    owner.state_threshold_spin.setToolTip(
        "Complete-link Jaccard/Tanimoto threshold; every pair in a "
        "family must meet this value."
    )
    owner.state_threshold_spin.valueChanged.connect(owner._refresh_states)
    controls.addWidget(owner.state_threshold_spin)
    controls.addWidget(QtWidgets.QLabel("Saved-frame step (ns):"))
    owner.time_step_spin = QtWidgets.QDoubleSpinBox()
    owner.time_step_spin.setRange(0.000001, 1_000_000.0)
    owner.time_step_spin.setDecimals(6)
    owner.time_step_spin.setValue(1.0)
    owner.time_step_spin.valueChanged.connect(owner._refresh_states)
    controls.addWidget(owner.time_step_spin)
    owner.load_trajectory_map_button = QtWidgets.QPushButton(
        "Load trajectory map"
    )
    owner.load_trajectory_map_button.setToolTip(
        "CSV columns: observation_id, replica_id, frame_index and optional "
        "time_ns. Row order defines the trajectory order."
    )
    owner.load_trajectory_map_button.clicked.connect(
        lambda: owner._load_trajectory_map(comparison=False)
    )
    controls.addWidget(owner.load_trajectory_map_button)
    controls.addWidget(QtWidgets.QLabel("Bootstrap:"))
    owner.bootstrap_iterations_combo = QtWidgets.QComboBox()
    for value in (250, 500, 1000, 2000):
        owner.bootstrap_iterations_combo.addItem(f"{value} resamples", value)
    owner.bootstrap_iterations_combo.setCurrentIndex(1)
    controls.addWidget(owner.bootstrap_iterations_combo)
    controls.addWidget(QtWidgets.QLabel("Block:"))
    owner.bootstrap_block_size_spin = QtWidgets.QSpinBox()
    owner.bootstrap_block_size_spin.setRange(0, 1_000_000)
    owner.bootstrap_block_size_spin.setSpecialValueText("Auto")
    owner.bootstrap_block_size_spin.setToolTip(
        "Number of consecutive saved frames per circular bootstrap block. "
        "Auto uses the square-root heuristic and records the resolved value."
    )
    controls.addWidget(owner.bootstrap_block_size_spin)
    owner.compute_uncertainty_button = QtWidgets.QPushButton(
        "Compute confidence intervals"
    )
    owner.compute_uncertainty_button.clicked.connect(
        owner._compute_uncertainty
    )
    controls.addWidget(owner.compute_uncertainty_button)
    controls.addStretch(1)
    layout.addLayout(controls)
    owner.state_status = QtWidgets.QLabel()
    owner.state_status.setObjectName("metricStrip")
    owner.state_status.setWordWrap(True)
    layout.addWidget(owner.state_status)
    owner.dynamic_tabs = QtWidgets.QTabWidget()
    owner.state_population_panel = panel_type(minimum_height=340)
    owner.state_timeline_panel = panel_type(minimum_height=270)
    owner.transition_panel = panel_type(minimum_height=360)
    owner.uncertainty_panel = panel_type(minimum_height=380)
    for title, panel in (
        ("Population", owner.state_population_panel),
        ("Pose order", owner.state_timeline_panel),
        ("Transitions", owner.transition_panel),
        ("Confidence", owner.uncertainty_panel),
    ):
        panel.setObjectName("analysisField")
        owner.dynamic_tabs.addTab(panel, title)
    layout.addWidget(owner.dynamic_tabs)
    state_title = QtWidgets.QLabel("Representatives and defining contacts")
    state_title.setObjectName("sectionTitle")
    layout.addWidget(state_title)
    owner.state_table = QtWidgets.QTableWidget(0, 7)
    owner.state_table.setHorizontalHeaderLabels(
        (
            "Group",
            "Population",
            "Population %",
            "Representative",
            "Mean similarity",
            "Temporal metric",
            "Characteristic interactions",
        )
    )
    owner.state_table.setEditTriggers(
        QtWidgets.QAbstractItemView.NoEditTriggers
    )
    owner.state_table.setSelectionBehavior(
        QtWidgets.QAbstractItemView.SelectRows
    )
    owner.state_table.setSelectionMode(
        QtWidgets.QAbstractItemView.SingleSelection
    )
    owner.state_table.horizontalHeader().setStretchLastSection(True)
    owner.state_table.itemSelectionChanged.connect(
        owner._state_selection_changed
    )
    owner.state_table.itemDoubleClicked.connect(
        lambda _item: owner._emit_representative()
    )
    layout.addWidget(owner.state_table)
    owner.open_representative_button = QtWidgets.QPushButton(
        "Open representative structure"
    )
    owner.open_representative_button.setEnabled(False)
    owner.open_representative_button.clicked.connect(
        owner._emit_representative
    )
    layout.addWidget(
        owner.open_representative_button,
        alignment=QtCore.Qt.AlignLeft,
    )
    return box


__all__ = ["build_dynamic_analysis"]
