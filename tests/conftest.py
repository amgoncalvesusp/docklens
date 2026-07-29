from __future__ import annotations

import os
from pathlib import Path

import pytest


# GUI tests must not initialize the native Windows platform integration inside
# an automated, non-interactive session. Set this before pytest-qt creates its
# QApplication so teardown follows the same headless path as CI/self-check.
os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

FIXTURE_DIR = Path(__file__).parent / "fixtures"


@pytest.fixture
def fixture_path():
    """Return absolute paths rooted in the repository-owned fixture folder."""

    def resolve(name: str) -> str:
        path = (FIXTURE_DIR / name).resolve()
        assert path.is_file(), f"Missing test fixture: {path}"
        return str(path)

    return resolve
