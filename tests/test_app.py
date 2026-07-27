"""Entry-point and packaged smoke-check regression tests."""

from __future__ import annotations

import sys
from types import SimpleNamespace

from docklens import app, main_window, self_check


def test_app_routes_self_check(monkeypatch):
    monkeypatch.setattr(sys, "argv", ["docklens", "--self-check"])
    monkeypatch.setattr(self_check, "run_self_check", lambda: 17)

    assert app.main() == 17


def test_app_routes_desktop_launch(monkeypatch):
    monkeypatch.setattr(sys, "argv", ["docklens"])
    monkeypatch.setattr(main_window, "launch", lambda: 23)

    assert app.main() == 23


def test_app_routes_vinalab_manifest_launch(monkeypatch, tmp_path):
    manifest = tmp_path / "launch.json"
    manifest.write_text("{}", encoding="utf-8")
    loaded = object()
    monkeypatch.setattr(sys, "argv", ["docklens", "--manifest", str(manifest)])
    monkeypatch.setattr("docklens.integration_manifest.load_manifest", lambda path: loaded)
    monkeypatch.setattr(main_window, "launch", lambda launch_manifest=None: 31 if launch_manifest is loaded else 0)

    assert app.main() == 31


def test_app_checks_vinalab_manifest_without_opening_gui(monkeypatch, tmp_path, capsys):
    manifest = tmp_path / "launch.json"
    manifest.write_text("{}", encoding="utf-8")
    loaded = SimpleNamespace(
        receptor_path=tmp_path / "receptor.pdb",
        poses_path=tmp_path / "poses.pdbqt",
        key_residues=("SER1A",),
        hbond_preset="plip",
    )
    analyzed = SimpleNamespace(summaries=(1, 2), details=(1, 2, 3), input_qc=())
    monkeypatch.setattr(sys, "argv", ["docklens", "--check-manifest", str(manifest)])
    monkeypatch.setattr("docklens.integration_manifest.load_manifest", lambda path: loaded)
    monkeypatch.setattr("docklens.batch_runner.run_paired", lambda *args, **kwargs: analyzed)

    assert app.main() == 0
    assert '"poses": 2' in capsys.readouterr().out


def test_app_contains_unexpected_manifest_analysis_failure(monkeypatch, tmp_path, capsys):
    manifest = tmp_path / "launch.json"
    manifest.write_text("{}", encoding="utf-8")
    loaded = SimpleNamespace(
        receptor_path=tmp_path / "receptor.pdb",
        poses_path=tmp_path / "poses.pdbqt",
        key_residues=(),
        hbond_preset="plip",
    )
    monkeypatch.setattr(sys, "argv", ["docklens", "--check-manifest", str(manifest)])
    monkeypatch.setattr("docklens.integration_manifest.load_manifest", lambda path: loaded)
    monkeypatch.setattr(
        "docklens.batch_runner.run_paired",
        lambda *args, **kwargs: (_ for _ in ()).throw(RuntimeError("internal detail")),
    )

    assert app.main() == 2
    captured = capsys.readouterr()
    assert "analysis_failed" in captured.err
    assert "internal detail" not in captured.err


def test_manifest_error_is_safe_in_windowed_build_without_stderr(monkeypatch, tmp_path):
    manifest = tmp_path / "launch.json"
    monkeypatch.setattr(sys, "argv", ["docklens", "--manifest", str(manifest)])
    monkeypatch.setattr(sys, "stderr", None)
    monkeypatch.setattr(sys, "stdout", None)

    assert app.main() == 2


def test_packaged_self_check_creates_and_reopens_workbook():
    assert self_check.run_self_check() == 0


def test_self_check_returns_through_bootloader(monkeypatch):
    monkeypatch.setattr(
        self_check.os,
        "_exit",
        lambda code: (_ for _ in ()).throw(AssertionError("os._exit must not bypass the bootloader")),
    )

    assert self_check.run_self_check() == 0
