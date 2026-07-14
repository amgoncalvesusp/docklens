"""Pure transformations from immutable analysis results to export views."""

from __future__ import annotations

from dataclasses import replace

import pandas as pd

from .interaction_core import VALID_TYPES
from .results import ExportFilter, RunResult


SUMMARY_BASE_COLS = [
    "ligand_id",
    "source_file",
    "sol",
    "pose",
    "docking_score",
    "n_total_interactions",
    "n_key_residue_interactions",
]
SUMMARY_META_COLS = ["source_path", "source_id", "pose_id", "resolution_method"]
DETAIL_BASE_COLS = [
    "ligand_id",
    "source_file",
    "interaction_type",
    "subtype",
    "ligand_atom",
    "receptor_residue",
    "receptor_atom",
    "distance_A",
    "is_key_residue",
]
DETAIL_META_COLS = [
    "source_id",
    "pose_id",
    "interaction_id",
    "sol",
    "pose",
    "docking_score",
    "resolution_method",
    "source_path",
    "ligand_residue",
    "ligand_resname",
    "ligand_resseq",
    "ligand_chain",
    "ligand_atom_serials",
    "ligand_element",
    "ligand_role",
    "receptor_resname",
    "receptor_resseq",
    "receptor_chain",
    "receptor_atom_serials",
    "receptor_element",
    "receptor_role",
    "water_residue",
    "water_atom",
    "water_atom_serials",
    "water_element",
    "receptor_water_distance_A",
    "ligand_water_distance_A",
    "water_angle_deg",
]


def _matches_text(detail, needle):
    if not needle:
        return True
    haystack = " ".join(
        (detail.receptor_residue, detail.source_file, detail.ligand_id)
    ).lower()
    return needle in haystack


def build_export_view(result: RunResult, export_filter: ExportFilter) -> RunResult:
    if export_filter.scope == "all":
        return result
    needle = export_filter.text.strip().lower()
    details = tuple(
        detail
        for detail in result.details
        if (
            export_filter.interaction_types is None
            or detail.interaction_type in export_filter.interaction_types
        )
        and (not export_filter.key_only or detail.is_key_residue)
        and _matches_text(detail, needle)
    )
    by_pose = {}
    for detail in details:
        entry = by_pose.setdefault(
            detail.pose_id,
            {"total": 0, "key": 0, "counts": {kind: 0 for kind in VALID_TYPES}},
        )
        entry["total"] += 1
        entry["key"] += int(detail.is_key_residue)
        entry["counts"][detail.interaction_type] += 1
    empty = {"total": 0, "key": 0, "counts": {kind: 0 for kind in VALID_TYPES}}
    summaries = tuple(
        replace(
            summary,
            n_total_interactions=by_pose.get(summary.pose_id, empty)["total"],
            n_key_residue_interactions=by_pose.get(summary.pose_id, empty)["key"],
            counts=by_pose.get(summary.pose_id, empty)["counts"],
        )
        for summary in result.summaries
    )
    return replace(result, details=details, summaries=summaries)


def summary_dataframe(result: RunResult) -> pd.DataFrame:
    rows = []
    for summary in result.summaries:
        row = {
            "ligand_id": summary.ligand_id,
            "source_file": summary.source_file,
            "sol": summary.sol,
            "pose": summary.pose,
            "docking_score": summary.docking_score,
            "n_total_interactions": summary.n_total_interactions,
            "n_key_residue_interactions": summary.n_key_residue_interactions,
            "source_path": summary.source_path,
            "source_id": summary.source_id,
            "pose_id": summary.pose_id,
            "resolution_method": summary.resolution_method,
        }
        for kind in VALID_TYPES:
            row[kind] = summary.counts.get(kind, 0)
        rows.append(row)
    return pd.DataFrame(
        rows, columns=SUMMARY_BASE_COLS + SUMMARY_META_COLS + list(VALID_TYPES)
    )


def _serials(endpoint):
    return ";".join(str(value) for value in endpoint.atom_serials)


def detail_dataframe(result: RunResult) -> pd.DataFrame:
    rows = []
    for detail in result.details:
        water = detail.water
        rows.append(
            {
                "ligand_id": detail.ligand_id,
                "source_file": detail.source_file,
                "interaction_type": detail.interaction_type,
                "subtype": detail.subtype,
                "ligand_atom": detail.ligand_atom,
                "receptor_residue": detail.receptor_residue,
                "receptor_atom": detail.receptor_atom,
                "distance_A": detail.distance_A,
                "is_key_residue": detail.is_key_residue,
                "source_id": detail.source_id,
                "pose_id": detail.pose_id,
                "interaction_id": detail.interaction_id,
                "sol": detail.sol,
                "pose": detail.pose,
                "docking_score": detail.docking_score,
                "resolution_method": detail.resolution_method,
                "source_path": detail.source_path,
                "ligand_residue": detail.ligand.residue_id,
                "ligand_resname": detail.ligand.resname,
                "ligand_resseq": detail.ligand.resseq,
                "ligand_chain": detail.ligand.chain,
                "ligand_atom_serials": _serials(detail.ligand),
                "ligand_element": detail.ligand.element,
                "ligand_role": detail.ligand.role,
                "receptor_resname": detail.receptor.resname,
                "receptor_resseq": detail.receptor.resseq,
                "receptor_chain": detail.receptor.chain,
                "receptor_atom_serials": _serials(detail.receptor),
                "receptor_element": detail.receptor.element,
                "receptor_role": detail.receptor.role,
                "water_residue": water.residue_id if water else "",
                "water_atom": water.atom_name if water else "",
                "water_atom_serials": _serials(water) if water else "",
                "water_element": water.element if water else "",
                "receptor_water_distance_A": detail.receptor_water_distance_A,
                "ligand_water_distance_A": detail.ligand_water_distance_A,
                "water_angle_deg": detail.water_angle_deg,
            }
        )
    return pd.DataFrame(rows, columns=DETAIL_BASE_COLS + DETAIL_META_COLS)


def residue_matrix_dataframe(result: RunResult, mode="count") -> pd.DataFrame:
    if mode not in {"count", "presence"}:
        raise ValueError("mode must be 'count' or 'presence'")
    identity = [
        ("Identity", "ligand_id"),
        ("Identity", "source_file"),
        ("Identity", "source_path"),
        ("Identity", "source_id"),
        ("Identity", "pose_id"),
        ("Identity", "sol"),
        ("Identity", "pose"),
        ("Identity", "docking_score"),
    ]
    observed = sorted(
        {(detail.receptor_residue, detail.interaction_type) for detail in result.details},
        key=lambda value: (value[0], VALID_TYPES.index(value[1])),
    )
    columns = pd.MultiIndex.from_tuples(identity + observed)
    counts = {}
    for detail in result.details:
        key = (detail.pose_id, detail.receptor_residue, detail.interaction_type)
        counts[key] = counts.get(key, 0) + 1
    rows = []
    for summary in result.summaries:
        row = [
            summary.ligand_id,
            summary.source_file,
            summary.source_path,
            summary.source_id,
            summary.pose_id,
            summary.sol,
            summary.pose,
            summary.docking_score,
        ]
        for residue, kind in observed:
            value = counts.get((summary.pose_id, residue, kind), 0)
            row.append(int(value > 0) if mode == "presence" else value)
        rows.append(row)
    return pd.DataFrame(rows, columns=columns)


def parameters_dataframe(result: RunResult, export_filter=None) -> pd.DataFrame:
    export_filter = export_filter or ExportFilter()
    params = result.parameters
    rows = [
        ("schema_version", params.schema_version),
        ("app_version", params.app_version),
        ("started_at", params.started_at),
        ("hbond_preset", params.hbond_preset),
        ("interaction_types", ", ".join(params.interaction_types)),
        ("key_residues", ", ".join(params.key_residues)),
        ("counting_unit", params.counting_unit),
        ("export_scope", export_filter.scope),
        ("matrix_mode", export_filter.matrix_mode),
        ("filter_text", export_filter.text),
        ("filter_key_only", export_filter.key_only),
        (
            "filter_interaction_types",
            "" if export_filter.interaction_types is None else ", ".join(sorted(export_filter.interaction_types)),
        ),
    ]
    rows.extend(("cutoff.%s" % key, value) for key, value in params.cutoffs)
    return pd.DataFrame(rows, columns=["parameter", "value"])


def input_qc_dataframe(result: RunResult) -> pd.DataFrame:
    rows = []
    for record in result.input_qc:
        rows.append(
            {
                "source_id": record.source_id,
                "source_file": record.source_file,
                "source_path": record.source_path,
                "pose_id": record.pose_id,
                "status": record.status,
                "code": record.code,
                "message": record.message,
                "format": record.format,
                "poses_found": record.poses_found,
                "poses_processed": record.poses_processed,
                "resolution_method": record.resolution_method,
                "receptor_atoms": record.receptor_atoms,
                "ligand_atoms": record.ligand_atoms,
                "water_atoms": record.water_atoms,
                "warnings": "; ".join(record.warnings),
            }
        )
    return pd.DataFrame(
        rows,
        columns=[
            "source_id",
            "source_file",
            "source_path",
            "pose_id",
            "status",
            "code",
            "message",
            "format",
            "poses_found",
            "poses_processed",
            "resolution_method",
            "receptor_atoms",
            "ligand_atoms",
            "water_atoms",
            "warnings",
        ],
    )
