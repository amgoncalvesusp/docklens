"""
export.py — write the Summary and Detail tables to CSV / XLSX.

CSV : two files ('<prefix>_summary.csv' and '<prefix>_detail.csv').
XLSX: one workbook with two sheets ('Summary', 'Detail'); interaction-type
      cells/columns are shaded with the Okabe-Ito palette for visual parity with
      the PyMOL plugin.
"""

from __future__ import annotations

import pandas as pd

from .interaction_core import INTERACTION_COLORS, VALID_TYPES, color_hex

SUMMARY_BASE_COLS = [
    "ligand_id",
    "source_file",
    "sol",
    "pose",
    "docking_score",
    "n_total_interactions",
    "n_key_residue_interactions",
]
DETAIL_COLS = [
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


def summary_dataframe(result) -> pd.DataFrame:
    rows = []
    for s in result.summaries:
        row = {
            "ligand_id": s.ligand_id,
            "source_file": s.source_file,
            "sol": s.sol,
            "pose": s.pose,
            "docking_score": s.docking_score,
            "n_total_interactions": s.n_total_interactions,
            "n_key_residue_interactions": s.n_key_residue_interactions,
        }
        for t in VALID_TYPES:
            row[t] = s.counts.get(t, 0)
        rows.append(row)
    cols = SUMMARY_BASE_COLS + list(VALID_TYPES)
    return pd.DataFrame(rows, columns=cols)


def detail_dataframe(result) -> pd.DataFrame:
    rows = []
    for d in result.details:
        rows.append(
            {
                "ligand_id": d.ligand_id,
                "source_file": d.source_file,
                "interaction_type": d.interaction_type,
                "subtype": d.subtype,
                "ligand_atom": d.ligand_atom,
                "receptor_residue": d.receptor_residue,
                "receptor_atom": d.receptor_atom,
                "distance_A": d.distance,
                "is_key_residue": d.is_key_residue,
            }
        )
    return pd.DataFrame(rows, columns=DETAIL_COLS)


def export_csv(result, path_prefix) -> list:
    """Write two CSV files; return the written paths."""
    if path_prefix.lower().endswith(".csv"):
        path_prefix = path_prefix[:-4]
    s_path = "%s_summary.csv" % path_prefix
    d_path = "%s_detail.csv" % path_prefix
    summary_dataframe(result).to_csv(s_path, index=False)
    detail_dataframe(result).to_csv(d_path, index=False)
    return [s_path, d_path]


def _fill(hex_color):
    from openpyxl.styles import PatternFill

    argb = "FF" + hex_color.lstrip("#").upper()
    return PatternFill(start_color=argb, end_color=argb, fill_type="solid")


def export_xlsx(result, path) -> str:
    """Write a two-sheet .xlsx with Okabe-Ito shading on interaction types."""
    if not path.lower().endswith(".xlsx"):
        path = path + ".xlsx"
    df_summary = summary_dataframe(result)
    df_detail = detail_dataframe(result)

    with pd.ExcelWriter(path, engine="openpyxl") as writer:
        df_summary.to_excel(writer, sheet_name="Summary", index=False)
        df_detail.to_excel(writer, sheet_name="Detail", index=False)
        wb = writer.book

        # Summary: shade the per-type count column headers.
        ws_s = wb["Summary"]
        header = {c.value: c.column for c in ws_s[1]}
        for t in VALID_TYPES:
            col = header.get(t)
            if col:
                ws_s.cell(row=1, column=col).fill = _fill(color_hex(t))

        # Detail: shade each interaction_type cell.
        ws_d = wb["Detail"]
        d_header = {c.value: c.column for c in ws_d[1]}
        type_col = d_header.get("interaction_type")
        if type_col:
            for r in range(2, ws_d.max_row + 1):
                cell = ws_d.cell(row=r, column=type_col)
                t = cell.value
                if t in INTERACTION_COLORS:
                    cell.fill = _fill(color_hex(t))
    return path
