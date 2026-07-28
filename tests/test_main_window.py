"""Offscreen UI integration tests for the v0.3 workflow."""

from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

from PyQt5 import QtWidgets

from docklens import batch_runner as br
from docklens.main_window import MainWindow
from docklens.results import Detail, Endpoint, Summary, make_result


def _endpoint(side, atom_name="X"):
    return Endpoint(
        side=side,
        kind="atom",
        atom_name=atom_name,
        atom_serials=(1,),
        resname="LIG" if side == "ligand" else "LYS",
        resseq="1" if side == "ligand" else "76",
        chain="" if side == "ligand" else "A",
    )


def _detail(interaction_id, interaction_type, distance):
    return Detail(
        ligand_id="PEP",
        source_file="pose.mol2",
        interaction_type=interaction_type,
        subtype="",
        ligand=_endpoint("ligand", "C1"),
        receptor=_endpoint("receptor", "NZ"),
        distance_A=distance,
        source_id="source-1",
        pose_id="pose-1",
        interaction_id=interaction_id,
        pose=1,
    )


def test_window_runs_filters_and_recomputes_keys(qtbot, fixture_path):
    window = MainWindow()
    qtbot.addWidget(window)
    window._files = [fixture_path("minimal_complex.pdb")]

    window._run()

    assert window._result.summaries
    residue = window._result.details[0].receptor_residue
    window.key_edit.setText(residue)
    window._key_text_changed()
    assert window._result.key_residues == frozenset({residue})
    assert window.detail_proxy.rowCount() > 0
    assert window.coverage_model.rowCount() == len(window._result.summaries)
    assert "matched" in window.key_status.text().lower()


def test_window_normalizes_and_reports_invalid_key_residue_input(qtbot, fixture_path):
    window = MainWindow()
    qtbot.addWidget(window)
    window._files = [fixture_path("minimal_complex.pdb")]
    window._run()
    residue = sorted(window._result.receptor_residues)[0].rstrip("_")

    window.key_edit.setText(f"{residue}; GLU999; not-a-residue")
    window._key_text_changed()

    assert window.key_edit.text() == " ".join(sorted((residue, "GLU999")))
    status = window.key_status.text().lower()
    assert "1 of 2" in status
    assert "unmatched: glu999" in status
    assert "invalid: not-a-residue" in status


def test_window_exports_selected_scope(qtbot, fixture_path, tmp_path, monkeypatch):
    window = MainWindow()
    qtbot.addWidget(window)
    window._result = br.run([fixture_path("minimal_complex.pdb")])
    output = tmp_path / "ui-export.xlsx"
    monkeypatch.setattr(
        QtWidgets.QFileDialog,
        "getSaveFileName",
        lambda *args, **kwargs: (str(output), "Excel (*.xlsx)"),
    )
    monkeypatch.setattr(
        window,
        "_choose_export_filter",
        lambda include_matrix_mode: br.ExportFilter(
            scope="filtered", interaction_types=frozenset({"hbond"})
        ),
    )
    monkeypatch.setattr(QtWidgets.QMessageBox, "information", lambda *args: None)

    window._export_xlsx()

    assert Path(output).is_file()


def test_window_analysis_profile_filters_visible_tables(qtbot):
    window = MainWindow()
    qtbot.addWidget(window)
    window._result = make_result(
        details=(
            _detail("i1", "hbond", 3.0),
            _detail("i2", "alkyl", 3.7),
            _detail("i3", "saltbridge", 4.3),
        ),
        summaries=(
            Summary(
                ligand_id="PEP",
                source_file="pose.mol2",
                sol=None,
                pose=1,
                docking_score=-8.0,
                n_total_interactions=3,
                n_key_residue_interactions=0,
                counts={"hbond": 1, "alkyl": 1, "saltbridge": 1},
                source_id="source-1",
                pose_id="pose-1",
            ),
        ),
    )

    window._refresh_tables()
    assert window.detail_proxy.rowCount() == 3

    window.analysis_combo.setCurrentIndex(window.analysis_combo.findData("ds_like"))

    assert window.detail_proxy.rowCount() == 2
    assert window.summary_model.rowCount() == 1


def test_window_loads_explicit_vinalab_pair(qtbot, fixture_path):
    window = MainWindow()
    qtbot.addWidget(window)
    manifest = SimpleNamespace(
        receptor_path=Path(fixture_path("minimal_complex.pdb")),
        poses_path=Path(fixture_path("two_poses_sol3.pdbqt")),
        hbond_preset="plip",
        key_residues=("SER1A",),
        result_path=None,
    )

    window.load_manifest(manifest)

    assert len(window._result.summaries) == 2
    assert all(
        item.resolution_method == "paired-manifest" for item in window._result.summaries
    )
    assert "DockingHub" in window.status.currentMessage()


def test_manifest_rerun_keeps_type_checkboxes_as_visual_filters(
    qtbot, monkeypatch, tmp_path
):
    window = MainWindow()
    qtbot.addWidget(window)
    manifest = SimpleNamespace(
        receptor_path=tmp_path / "receptor.pdb",
        poses_path=tmp_path / "poses.pdbqt",
        hbond_preset="plip",
        key_residues=(),
        result_path=None,
    )
    captured = {}

    def fake_run_paired(*args, **kwargs):
        captured.update(kwargs)
        return make_result()

    monkeypatch.setattr(br, "run_paired", fake_run_paired)
    window._launch_manifest = manifest
    window._files = [str(manifest.receptor_path), str(manifest.poses_path)]
    window.type_boxes["alkyl"].setChecked(False)

    window._run()

    assert "types" not in captured


def test_window_keeps_running_when_manifest_pair_cannot_be_parsed(
    qtbot, monkeypatch, tmp_path
):
    window = MainWindow()
    qtbot.addWidget(window)
    manifest = SimpleNamespace(
        receptor_path=tmp_path / "receptor.pdb",
        poses_path=tmp_path / "poses.pdbqt",
        hbond_preset="plip",
        key_residues=(),
        result_path=None,
    )
    monkeypatch.setattr(
        br,
        "run_paired",
        lambda *args, **kwargs: (_ for _ in ()).throw(RuntimeError("bad pair")),
    )
    shown = []
    monkeypatch.setattr(
        QtWidgets.QMessageBox, "critical", lambda *args: shown.append(args)
    )

    window.load_manifest(manifest)

    assert window._result is None
    assert shown
    assert "could not" in window.status.currentMessage().lower()


def test_window_writes_vinalab_roundtrip_after_manifest_analysis(
    qtbot, fixture_path, monkeypatch, tmp_path
):
    window = MainWindow()
    qtbot.addWidget(window)
    manifest = SimpleNamespace(
        receptor_path=Path(fixture_path("minimal_complex.pdb")),
        poses_path=Path(fixture_path("two_poses_sol3.pdbqt")),
        hbond_preset="plip",
        key_residues=(),
        result_path=tmp_path / "docklens-result.json",
    )
    captured = []
    monkeypatch.setattr(
        "docklens.main_window.write_integration_result",
        lambda contract, result: captured.append((contract, result)),
    )

    window.load_manifest(manifest)

    assert captured == [(manifest, window._result)]
    window.key_edit.setText("SER1A")
    window._key_text_changed()
    assert len(captured) == 2
    assert captured[-1] == (manifest, window._result)
