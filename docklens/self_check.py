"""Non-visual packaged-application smoke check used by release builds."""

from __future__ import annotations

import os
import sys
import tempfile
from pathlib import Path

from openpyxl import load_workbook

from . import __version__
from .export import export_xlsx
from .results import AnalysisParameters, make_result


def run_self_check():
    parameters = AnalysisParameters(app_version=__version__)
    result = make_result(parameters=parameters)
    with tempfile.TemporaryDirectory(prefix="docklens-self-check-") as directory:
        output = Path(export_xlsx(result, Path(directory) / "self-check.xlsx"))
        workbook = load_workbook(output, read_only=True)
        expected = ["Summary", "Residue Matrix", "Detail", "Parameters", "Input QC"]
        actual = workbook.sheetnames
        workbook.close()
        if actual != expected:
            raise RuntimeError("Unexpected workbook schema: %r" % actual)

    # Exercise the packaged Qt runtime and real main window last. A frozen,
    # windowed bootloader is terminated immediately after successful creation;
    # allowing Qt to tear itself down can leave the smoke-test child alive.
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    from PyQt5 import QtWidgets

    from .main_window import MainWindow

    existing_app = QtWidgets.QApplication.instance()
    application = existing_app or QtWidgets.QApplication(["DockLens", "--self-check"])
    window = MainWindow()
    if getattr(sys, "frozen", False):
        os._exit(0)
    window.close()
    window.deleteLater()
    application.processEvents()
    if existing_app is None:
        application.quit()
    return 0
