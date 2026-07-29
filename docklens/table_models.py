"""Qt table models shared by the DockLens auditable-table workspace."""

from __future__ import annotations

from PyQt5 import QtCore, QtGui

from .interaction_core import INTERACTION_COLORS, color_hex


class DataFrameModel(QtCore.QAbstractTableModel):
    """Table model over a pandas DataFrame with type-aware sort and colour."""

    def __init__(self, dataframe, colour_type_col=None):
        super().__init__()
        self._df = dataframe.reset_index(drop=True)
        self._colour_type_col = colour_type_col

    def set_dataframe(self, dataframe):
        self.beginResetModel()
        self._df = dataframe.reset_index(drop=True)
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
        column = self._df.columns[index.column()]
        value = self._df.iat[index.row(), index.column()]
        if role in (QtCore.Qt.DisplayRole, QtCore.Qt.EditRole):
            if value is None or (
                isinstance(value, float) and value != value
            ):
                return ""
            return str(value)
        if role == QtCore.Qt.UserRole:
            return value
        if role == QtCore.Qt.BackgroundRole and self._colour_type_col == column:
            if value in INTERACTION_COLORS:
                return QtGui.QBrush(QtGui.QColor(color_hex(value)))
        if role == QtCore.Qt.ForegroundRole and self._colour_type_col == column:
            if value in INTERACTION_COLORS:
                return QtGui.QBrush(QtGui.QColor("white"))
        return None


class MultiFilterProxy(QtCore.QSortFilterProxyModel):
    """Sort numerically where possible and filter table evidence."""

    def __init__(self):
        super().__init__()
        self.type_filter = None
        self.text_filter = ""
        self.key_only = False
        self.setSortRole(QtCore.Qt.UserRole)

    def _column(self, name):
        model = self.sourceModel()
        for column in range(model.columnCount()):
            if model.headerData(column, QtCore.Qt.Horizontal) == name:
                return column
        return -1

    def lessThan(self, left, right):
        left_value = self.sourceModel().data(left, QtCore.Qt.UserRole)
        right_value = self.sourceModel().data(right, QtCore.Qt.UserRole)
        try:
            return float(left_value) < float(right_value)
        except (TypeError, ValueError):
            return str(left_value) < str(right_value)

    def filterAcceptsRow(self, row, parent):
        model = self.sourceModel()

        def cell(name):
            column = self._column(name)
            if column < 0:
                return ""
            return model.data(model.index(row, column)) or ""

        type_column = self._column("interaction_type")
        if self.type_filter is not None and type_column >= 0:
            if cell("interaction_type") not in self.type_filter:
                return False
        if self.key_only and self._column("is_key_residue") >= 0:
            if cell("is_key_residue") not in ("True", "true", "1"):
                return False
        if self.text_filter:
            needle = self.text_filter.lower()
            haystack = " ".join(
                cell(name)
                for name in ("receptor_residue", "source_file", "ligand_id")
            ).lower()
            if needle not in haystack:
                return False
        return True


__all__ = ["DataFrameModel", "MultiFilterProxy"]
