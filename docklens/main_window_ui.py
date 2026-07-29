"""Declarative construction of the DockLens 1.0 desktop shell."""

from __future__ import annotations

import pandas as pd
from PyQt5 import QtCore, QtGui, QtWidgets

from .analytics_widgets import AnalyticsWorkspace
from .interaction_core import VALID_TYPES
from .table_models import DataFrameModel, MultiFilterProxy


def _header(window, resource_path):
    bar = QtWidgets.QHBoxLayout()
    bar.setContentsMargins(16, 8, 16, 8)
    logo = QtWidgets.QLabel()
    pixmap = QtGui.QPixmap(resource_path("docklens_icon.png"))
    if not pixmap.isNull():
        logo.setPixmap(
            pixmap.scaledToHeight(36, QtCore.Qt.SmoothTransformation)
        )
    bar.addWidget(logo)
    titles = QtWidgets.QVBoxLayout()
    title = QtWidgets.QLabel("DockLens 1.0")
    title.setObjectName("title")
    title_palette = title.palette()
    title_palette.setColor(QtGui.QPalette.WindowText, QtGui.QColor("#FFFFFF"))
    title.setPalette(title_palette)
    window.brand_title = title
    subtitle = QtWidgets.QLabel("Reciprocal-space interaction atlas")
    subtitle.setObjectName("subtitle")
    subtitle_palette = subtitle.palette()
    subtitle_palette.setColor(
        QtGui.QPalette.WindowText, QtGui.QColor("#DCE7EE")
    )
    subtitle.setPalette(subtitle_palette)
    window.brand_subtitle = subtitle
    titles.addWidget(title)
    titles.addWidget(subtitle)
    bar.addLayout(titles)
    bar.addStretch(1)
    window.dataset_context = QtWidgets.QLabel("No dataset loaded")
    window.dataset_context.setObjectName("datasetContext")
    context_palette = window.dataset_context.palette()
    context_palette.setColor(
        QtGui.QPalette.WindowText, QtGui.QColor("#DCE7EE")
    )
    window.dataset_context.setPalette(context_palette)
    bar.addWidget(window.dataset_context)
    return bar


def _toolbar(window):
    layout = QtWidgets.QVBoxLayout()
    action_row = QtWidgets.QHBoxLayout()
    open_files = QtWidgets.QPushButton("Open files")
    open_folder = QtWidgets.QPushButton("Open folder")
    run = QtWidgets.QPushButton("Run detection")
    run.setObjectName("primary")
    window.run_detection_button = run
    reset = QtWidgets.QPushButton("Reset")
    open_files.clicked.connect(window._open_files)
    open_folder.clicked.connect(window._open_folder)
    run.clicked.connect(window._run)
    reset.clicked.connect(window._reset)
    for button in (open_files, open_folder, run, reset):
        action_row.addWidget(button)
    action_row.addSpacing(16)
    action_row.addWidget(QtWidgets.QLabel("Criteria:"))
    window.preset_combo = QtWidgets.QComboBox()
    window.preset_combo.addItem("PLIP (default)", "plip")
    window.preset_combo.addItem(
        "DS-calibrated beta (explicit-H geometry)", "dsv"
    )
    action_row.addWidget(window.preset_combo)
    action_row.addSpacing(12)
    action_row.addWidget(QtWidgets.QLabel("Evidence:"))
    window.mode_combo = QtWidgets.QComboBox()
    window.mode_combo.addItem("Docking poses", "docking")
    window.mode_combo.addItem("Molecular dynamics frames", "md")
    window.mode_combo.currentIndexChanged.connect(window._mode_changed)
    action_row.addWidget(window.mode_combo)
    action_row.addSpacing(12)
    action_row.addWidget(QtWidgets.QLabel("Profile:"))
    window.analysis_combo = QtWidgets.QComboBox()
    window.analysis_combo.addItem("Complete", "complete")
    window.analysis_combo.addItem(
        "Discovery Studio-like", "ds_like"
    )
    window.analysis_combo.currentIndexChanged.connect(window._refresh_tables)
    action_row.addWidget(window.analysis_combo)
    action_row.addStretch(1)
    comparison = QtWidgets.QPushButton("Load system B")
    comparison.clicked.connect(window._load_comparison)
    action_row.addWidget(comparison)
    layout.addLayout(action_row)

    scope_row = QtWidgets.QHBoxLayout()
    scope_row.addWidget(QtWidgets.QLabel("Chart scope — ligand/file A:"))
    window.primary_ligand_combo = QtWidgets.QComboBox()
    window.primary_ligand_combo.setMinimumContentsLength(18)
    window.primary_ligand_combo.setSizeAdjustPolicy(
        QtWidgets.QComboBox.AdjustToMinimumContentsLengthWithIcon
    )
    window.primary_ligand_combo.addItem(
        "All ligands / uploaded files", None
    )
    window.primary_ligand_combo.setEnabled(False)
    window.primary_ligand_combo.currentIndexChanged.connect(
        window._chart_scope_changed
    )
    scope_row.addWidget(window.primary_ligand_combo)
    scope_row.addWidget(QtWidgets.QLabel("System B:"))
    window.comparison_ligand_combo = QtWidgets.QComboBox()
    window.comparison_ligand_combo.setMinimumContentsLength(16)
    window.comparison_ligand_combo.setSizeAdjustPolicy(
        QtWidgets.QComboBox.AdjustToMinimumContentsLengthWithIcon
    )
    window.comparison_ligand_combo.addItem(
        "All ligands / uploaded files", None
    )
    window.comparison_ligand_combo.setEnabled(False)
    window.comparison_ligand_combo.currentIndexChanged.connect(
        window._chart_scope_changed
    )
    scope_row.addWidget(window.comparison_ligand_combo)
    window.chart_scope_status = QtWidgets.QLabel(
        "Charts use all loaded observations."
    )
    window.chart_scope_status.setObjectName("workspaceDescription")
    scope_row.addWidget(window.chart_scope_status, 1)
    window.export_figure_button = QtWidgets.QPushButton("Export figure")
    window.export_figure_button.clicked.connect(window._export_figure)
    scope_row.addWidget(window.export_figure_button)
    export_csv = QtWidgets.QPushButton("CSV")
    export_xlsx = QtWidgets.QPushButton("XLSX")
    export_csv.clicked.connect(window._export_csv)
    export_xlsx.clicked.connect(window._export_xlsx)
    scope_row.addWidget(export_csv)
    scope_row.addWidget(export_xlsx)
    layout.addLayout(scope_row)
    return layout


def _navigation(window):
    rail = QtWidgets.QWidget()
    rail.setObjectName("navigationRail")
    rail.setMinimumWidth(132)
    rail.setMaximumWidth(172)
    layout = QtWidgets.QVBoxLayout(rail)
    layout.setContentsMargins(8, 16, 8, 12)
    layout.setSpacing(6)
    group = QtWidgets.QButtonGroup(window)
    group.setExclusive(True)
    window.workspace_buttons = []
    for index, label in enumerate(
        ("Residues", "Fingerprint", "Compare", "Tables")
    ):
        button = QtWidgets.QPushButton(label)
        button.setObjectName("navButton")
        button.setCheckable(True)
        button.setChecked(index == 0)
        button.setToolTip(f"Open {label}")
        button.clicked.connect(
            lambda _checked=False, target=index: (
                window.workspace_stack.setCurrentIndex(target)
            )
        )
        group.addButton(button)
        window.workspace_buttons.append(button)
        layout.addWidget(button)
    layout.addStretch(1)
    window.open_project_button = QtWidgets.QPushButton("Open project")
    window.open_project_button.setObjectName("navButton")
    window.open_project_button.clicked.connect(window._open_project)
    layout.addWidget(window.open_project_button)
    window.save_project_button = QtWidgets.QPushButton("Save project")
    window.save_project_button.setObjectName("navButton")
    window.save_project_button.setEnabled(False)
    window.save_project_button.clicked.connect(window._save_project)
    layout.addWidget(window.save_project_button)
    about = QtWidgets.QPushButton("About")
    about.setObjectName("navButton")
    about.clicked.connect(window._about)
    layout.addWidget(about)
    layout.activate()
    rail.setMinimumHeight(layout.sizeHint().height())
    window.nav_rail = rail
    area = QtWidgets.QScrollArea()
    area.setObjectName("navigationScroll")
    area.setWidgetResizable(True)
    area.setFrameShape(QtWidgets.QFrame.NoFrame)
    area.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff)
    area.setVerticalScrollBarPolicy(QtCore.Qt.ScrollBarAsNeeded)
    area.setWidget(rail)
    window.nav_scroll_area = area
    return area


def _filters_card(window):
    card = QtWidgets.QGroupBox("Filters & key residues")
    layout = QtWidgets.QGridLayout(card)
    layout.addWidget(QtWidgets.QLabel("Key residues:"), 0, 0)
    window.key_edit = QtWidgets.QLineEdit()
    window.key_edit.setPlaceholderText(
        "e.g. SER70; LYS73; GLU166 (space, comma, semicolon or line break)"
    )
    window.key_edit.editingFinished.connect(window._key_text_changed)
    layout.addWidget(window.key_edit, 0, 1, 1, 2)
    layout.addWidget(QtWidgets.QLabel("Search:"), 0, 3)
    window.search_edit = QtWidgets.QLineEdit()
    window.search_edit.setPlaceholderText("residue / file / ligand")
    window.search_edit.textChanged.connect(window._apply_filters)
    layout.addWidget(window.search_edit, 0, 4)
    window.key_only_cb = QtWidgets.QCheckBox("Key residues only")
    window.key_only_cb.stateChanged.connect(window._apply_filters)
    layout.addWidget(window.key_only_cb, 0, 5)
    layout.addWidget(
        QtWidgets.QLabel("Pick key residues from detected protein:"),
        1,
        0,
        1,
        2,
    )
    window.res_filter = QtWidgets.QLineEdit()
    window.res_filter.setPlaceholderText("filter list...")
    window.res_filter.textChanged.connect(window._filter_residue_list)
    layout.addWidget(window.res_filter, 1, 3, 1, 3)
    window.res_list = QtWidgets.QListWidget()
    window.res_list.setMaximumHeight(140)
    window.res_list.itemChanged.connect(window._residue_checks_changed)
    layout.addWidget(window.res_list, 2, 0, 1, 6)
    window.key_status = QtWidgets.QLabel("No key residues configured.")
    window.key_status.setWordWrap(True)
    window.key_status.setObjectName("workspaceDescription")
    layout.addWidget(window.key_status, 3, 0, 1, 6)
    return card


def _type_card(window):
    card = QtWidgets.QGroupBox("Interaction types (filter Detail table)")
    grid = QtWidgets.QGridLayout(card)
    window.type_boxes = {}
    for index, interaction_type in enumerate(VALID_TYPES):
        checkbox = QtWidgets.QCheckBox(interaction_type)
        checkbox.setChecked(True)
        checkbox.stateChanged.connect(window._apply_filters)
        window.type_boxes[interaction_type] = checkbox
        grid.addWidget(checkbox, index // 6, index % 6)
    return card


def _tables(window):
    window.tabs = QtWidgets.QTabWidget()
    window.summary_view = QtWidgets.QTableView()
    window.coverage_view = QtWidgets.QTableView()
    window.detail_view = QtWidgets.QTableView()
    for view in (
        window.summary_view,
        window.coverage_view,
        window.detail_view,
    ):
        view.setSortingEnabled(True)
        view.setAlternatingRowColors(True)
        view.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectRows)
        view.horizontalHeader().setStretchLastSection(True)
    window.tabs.addTab(window.summary_view, "Summary")
    window.tabs.addTab(window.coverage_view, "Key Residue Coverage")
    window.tabs.addTab(window.detail_view, "Detail")
    return window.tabs


def _tables_workspace(window):
    content = QtWidgets.QWidget()
    content.setMinimumWidth(760)
    layout = QtWidgets.QVBoxLayout(content)
    layout.setContentsMargins(24, 20, 24, 24)
    layout.setSpacing(14)
    heading = QtWidgets.QLabel("Auditable interaction tables")
    heading.setObjectName("workspaceTitle")
    description = QtWidgets.QLabel(
        "Raw atom-pair evidence, key-residue coverage and pose summaries "
        "remain available for both Complete and Discovery Studio-like views."
    )
    description.setObjectName("workspaceDescription")
    description.setWordWrap(True)
    layout.addWidget(heading)
    layout.addWidget(description)
    layout.addWidget(_filters_card(window))
    layout.addWidget(_type_card(window))
    layout.addWidget(_tables(window), 1)
    area = QtWidgets.QScrollArea()
    area.setObjectName("workspaceScroll")
    area.setWidgetResizable(True)
    area.setFrameShape(QtWidgets.QFrame.NoFrame)
    area.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarAsNeeded)
    area.setVerticalScrollBarPolicy(QtCore.Qt.ScrollBarAsNeeded)
    area.setWidget(content)
    return area


def _lens_panel(window):
    panel = QtWidgets.QWidget()
    panel.setObjectName("lensPanel")
    panel.setMinimumWidth(220)
    panel.setMaximumWidth(330)
    layout = QtWidgets.QVBoxLayout(panel)
    layout.setContentsMargins(16, 18, 16, 16)
    title = QtWidgets.QLabel("DockLens Lens")
    title.setObjectName("sectionTitle")
    layout.addWidget(title)
    subtitle = QtWidgets.QLabel(
        "Synchronized residue evidence for the active analytical profile."
    )
    subtitle.setObjectName("workspaceDescription")
    subtitle.setWordWrap(True)
    layout.addWidget(subtitle)
    layout.addSpacing(12)
    window.lens_residue = QtWidgets.QLabel("No selection")
    window.lens_residue.setObjectName("workspaceTitle")
    window.lens_residue.setWordWrap(True)
    layout.addWidget(window.lens_residue)
    window.lens_evidence = QtWidgets.QLabel(
        "Run detection to inspect evidence."
    )
    window.lens_evidence.setTextInteractionFlags(
        QtCore.Qt.TextSelectableByMouse
    )
    window.lens_evidence.setAlignment(QtCore.Qt.AlignTop)
    window.lens_evidence.setWordWrap(True)
    layout.addWidget(window.lens_evidence)
    layout.addStretch(1)
    window.lens_profile = QtWidgets.QLabel("Profile: Complete")
    window.lens_profile.setObjectName("metricStrip")
    window.lens_profile.setWordWrap(True)
    layout.addWidget(window.lens_profile)
    window.lens_panel = panel
    return panel


def _init_models(window):
    window.summary_model = DataFrameModel(pd.DataFrame())
    window.coverage_model = DataFrameModel(pd.DataFrame())
    window.detail_model = DataFrameModel(
        pd.DataFrame(), colour_type_col="interaction_type"
    )
    window.summary_proxy = MultiFilterProxy()
    window.summary_proxy.setSourceModel(window.summary_model)
    window.coverage_proxy = MultiFilterProxy()
    window.coverage_proxy.setSourceModel(window.coverage_model)
    window.detail_proxy = MultiFilterProxy()
    window.detail_proxy.setSourceModel(window.detail_model)
    window.summary_view.setModel(window.summary_proxy)
    window.coverage_view.setModel(window.coverage_proxy)
    window.detail_view.setModel(window.detail_proxy)


def build_main_window_ui(window, resource_path):
    """Populate a ``MainWindow`` while keeping behavioral methods separate."""
    central = QtWidgets.QWidget()
    window.setCentralWidget(central)
    root = QtWidgets.QVBoxLayout(central)
    root.setContentsMargins(0, 0, 0, 0)
    root.setSpacing(0)
    top_bar = QtWidgets.QWidget()
    top_bar.setObjectName("topBar")
    top_bar.setLayout(_header(window, resource_path))
    root.addWidget(top_bar)
    context_bar = QtWidgets.QWidget()
    context_layout = _toolbar(window)
    context_layout.setContentsMargins(14, 8, 14, 8)
    context_bar.setLayout(context_layout)
    context_bar.setMinimumWidth(context_bar.sizeHint().width())
    context_scroll = QtWidgets.QScrollArea()
    context_scroll.setObjectName("contextScroll")
    context_scroll.setWidgetResizable(True)
    context_scroll.setFrameShape(QtWidgets.QFrame.NoFrame)
    context_scroll.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarAsNeeded)
    context_scroll.setVerticalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff)
    context_scroll.setWidget(context_bar)
    context_scroll.setFixedHeight(
        context_bar.sizeHint().height()
        + context_scroll.style().pixelMetric(
            QtWidgets.QStyle.PM_ScrollBarExtent
        )
        + 2
    )
    window.context_scroll_area = context_scroll
    root.addWidget(context_scroll)
    window.analytics_workspace = AnalyticsWorkspace(window)
    window.analytics_workspace.residueSelected.connect(window._update_lens)
    window.analytics_workspace.representativeSelected.connect(
        window._open_representative
    )
    window.workspace_stack = QtWidgets.QStackedWidget()
    window.workspace_stack.addWidget(window.analytics_workspace.residue_page)
    window.workspace_stack.addWidget(
        window.analytics_workspace.fingerprint_page
    )
    window.workspace_stack.addWidget(window.analytics_workspace.compare_page)
    tables_page = _tables_workspace(window)
    window.workspace_stack.addWidget(tables_page)
    window.workspace_stack.currentChanged.connect(window._workspace_changed)
    window.workspace_scroll_areas = list(
        window.analytics_workspace.scroll_areas
    ) + [tables_page]
    body = QtWidgets.QSplitter(QtCore.Qt.Horizontal)
    body.setChildrenCollapsible(False)
    body.addWidget(_navigation(window))
    body.addWidget(window.workspace_stack)
    body.addWidget(_lens_panel(window))
    body.setStretchFactor(0, 0)
    body.setStretchFactor(1, 1)
    body.setStretchFactor(2, 0)
    body.setSizes((150, 850, 260))
    window.body_splitter = body
    root.addWidget(body, 1)
    window.status = window.statusBar()
    window.status.showMessage(
        "Open a file, list or folder, then Run detection."
    )
    _init_models(window)


__all__ = ["build_main_window_ui"]
