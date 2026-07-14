"""
app.py — entry point for the DockLens desktop application.

Run with:  python -m docklens.app
"""

from __future__ import annotations

import sys


def main():
    if "--self-check" in sys.argv:
        from .self_check import run_self_check

        return run_self_check()
    from .main_window import launch

    return launch()


if __name__ == "__main__":
    sys.exit(main())
