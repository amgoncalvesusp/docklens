"""Frozen-app entry point for PyInstaller (imports the docklens package)."""

import sys

from docklens.main_window import launch

if __name__ == "__main__":
    sys.exit(launch())
