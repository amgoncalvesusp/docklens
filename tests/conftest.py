from __future__ import annotations

from pathlib import Path

import pytest


FIXTURE_DIR = Path(__file__).parent / "fixtures"


@pytest.fixture
def fixture_path():
    """Return absolute paths rooted in the repository-owned fixture folder."""

    def resolve(name: str) -> str:
        path = (FIXTURE_DIR / name).resolve()
        assert path.is_file(), f"Missing test fixture: {path}"
        return str(path)

    return resolve
