"""Frozen-app entry point for PyInstaller (imports the docklens package)."""

import sys

from docklens.app import main

if __name__ == "__main__":
    sys.exit(main())
