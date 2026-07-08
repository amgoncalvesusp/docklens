"""
app.py — entry point for the DockLens desktop application.

Run with:  python -m docklens.app
"""

from __future__ import annotations

import sys

from .main_window import launch


def main():
    return launch()


if __name__ == "__main__":
    sys.exit(main())
