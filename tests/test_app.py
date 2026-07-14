"""Entry-point and packaged smoke-check regression tests."""

from __future__ import annotations

import sys

import pytest

from docklens import app, main_window, self_check


def test_app_routes_self_check(monkeypatch):
    monkeypatch.setattr(sys, "argv", ["docklens", "--self-check"])
    monkeypatch.setattr(self_check, "run_self_check", lambda: 17)

    assert app.main() == 17


def test_app_routes_desktop_launch(monkeypatch):
    monkeypatch.setattr(sys, "argv", ["docklens"])
    monkeypatch.setattr(main_window, "launch", lambda: 23)

    assert app.main() == 23


def test_packaged_self_check_creates_and_reopens_workbook():
    assert self_check.run_self_check() == 0


def test_frozen_self_check_forces_clean_bootloader_exit(monkeypatch):
    monkeypatch.setattr(self_check.sys, "frozen", True, raising=False)
    monkeypatch.setattr(
        self_check.os,
        "_exit",
        lambda code: (_ for _ in ()).throw(SystemExit(code)),
    )

    with pytest.raises(SystemExit) as stopped:
        self_check.run_self_check()

    assert stopped.value.code == 0
