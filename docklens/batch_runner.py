"""
batch_runner.py — drive parsing, entity resolution and detection over inputs.

Input modes: a single file, a list of files, or a folder (recursive scan of
.mol2/.pdb/.pdbqt). Produces two tables:
  * Summary  — one row per ligand/pose (counts, score, key-residue count).
  * Detail   — one row per interaction.

Key-residue membership is recomputed cheaply (``recompute_key``) without redoing
the geometric detection.
"""

from __future__ import annotations

import os
from dataclasses import dataclass, field

from .entity_resolver import resolve, resolve_manual, set_sides
from .interaction_core import (
    VALID_TYPES,
    compute_interactions,
    endpoint_name,
    endpoint_resid,
    endpoint_side,
)
from .parser_mol2 import parse_mol2
from .parser_pdb import parse_pdb
from .parser_pdbqt import parse_pdbqt

_PARSERS = {".mol2": parse_mol2, ".pdb": parse_pdb, ".pdbqt": parse_pdbqt}
SUPPORTED_EXT = tuple(_PARSERS)


@dataclass
class Detail:
    ligand_id: str
    source_file: str
    interaction_type: str
    subtype: str
    ligand_atom: str
    receptor_residue: str
    receptor_atom: str
    distance: float
    is_key_residue: bool = False
    _res_nochain: str = ""  # internal: resn+resi for chain-agnostic key matching


@dataclass
class Summary:
    ligand_id: str
    source_file: str
    sol: object
    pose: int
    docking_score: object
    n_total_interactions: int
    n_key_residue_interactions: int
    counts: dict = field(default_factory=dict)


@dataclass
class Pending:
    """A pose whose ligand/receptor split needs user confirmation (fallback)."""

    pose: object
    resolution: object
    preview: str
    source_file: str


@dataclass
class RunResult:
    details: list = field(default_factory=list)
    summaries: list = field(default_factory=list)
    pending: list = field(default_factory=list)  # needs confirmation
    key_residues: set = field(default_factory=set)


# ---------------------------------------------------------------------------
# File gathering + parsing
# ---------------------------------------------------------------------------


def gather_files(paths):
    """Expand a path or list of paths/folders into a flat list of input files."""
    if isinstance(paths, str):
        paths = [paths]
    out = []
    for p in paths:
        if os.path.isdir(p):
            for root, _dirs, files in os.walk(p):
                for f in files:
                    if f.lower().endswith(SUPPORTED_EXT):
                        out.append(os.path.join(root, f))
        elif os.path.isfile(p) and p.lower().endswith(SUPPORTED_EXT):
            out.append(p)
    return out


def parse_file(path):
    """Parse any supported file into a list of ParsedPose."""
    ext = os.path.splitext(path)[1].lower()
    parser = _PARSERS.get(ext)
    if parser is None:
        raise ValueError("Unsupported file type: %s" % path)
    return parser(path)


# ---------------------------------------------------------------------------
# Key-residue matching
# ---------------------------------------------------------------------------


def normalize_key_residues(items):
    """Normalise a user key-residue list to an upper-case set."""
    if isinstance(items, str):
        items = items.replace(",", " ").split()
    return {str(x).strip().upper() for x in items if str(x).strip()}


def _is_key(res_tag, res_nochain, key_set):
    if not key_set:
        return False
    return res_tag.upper() in key_set or res_nochain.upper() in key_set


# ---------------------------------------------------------------------------
# Detection driver
# ---------------------------------------------------------------------------


def _detail_from_interaction(it, ligand_id, source_file, key_set):
    a, b = it["a_obj"], it["b_obj"]
    if endpoint_side(a) == "ligand":
        lig, rec = a, b
    elif endpoint_side(b) == "ligand":
        lig, rec = b, a
    else:
        lig, rec = a, b  # water-bridge receptor leg (no ligand endpoint)
    rec_res = endpoint_resid(rec)
    rec_obj = rec.atoms[0] if hasattr(rec, "atoms") else rec
    res_nochain = "%s%s" % (rec_obj.resn, rec_obj.resi)
    return Detail(
        ligand_id=ligand_id,
        source_file=os.path.basename(source_file),
        interaction_type=it["type"],
        subtype=it.get("subtype", ""),
        ligand_atom=endpoint_name(lig),
        receptor_residue=rec_res,
        receptor_atom=endpoint_name(rec),
        distance=round(float(it["dist"]), 2),
        is_key_residue=_is_key(rec_res, res_nochain, key_set),
        _res_nochain=res_nochain,
    )


def _summarize(ligand_id, source_file, sol, pose_no, score, details, key_set):
    counts = {t: 0 for t in VALID_TYPES}
    n_key = 0
    for d in details:
        counts[d.interaction_type] = counts.get(d.interaction_type, 0) + 1
        if d.is_key_residue:
            n_key += 1
    return Summary(
        ligand_id=ligand_id,
        source_file=os.path.basename(source_file),
        sol=sol,
        pose=pose_no,
        docking_score=score,
        n_total_interactions=len(details),
        n_key_residue_interactions=n_key,
        counts=counts,
    )


def run(
    paths, types=None, key_residues=None, confirm_fallback=False, manual_overrides=None
):
    """Run detection over the given inputs.

    types            list of interaction types (default: all).
    key_residues     iterable / string of key residue tags.
    confirm_fallback if True, fallback resolutions are run anyway (headless);
                     if False, they are collected in RunResult.pending instead.
    manual_overrides {source_file: set_of_ligand_serials} to force a split.
    """
    key_set = normalize_key_residues(key_residues or [])
    manual_overrides = manual_overrides or {}
    result = RunResult(key_residues=key_set)

    for path in gather_files(paths):
        stem = os.path.splitext(os.path.basename(path))[0]
        try:
            poses = parse_file(path)
        except Exception as exc:  # noqa: BLE001 - report, keep going
            print("[batch] failed to parse %s: %s" % (path, exc))
            continue

        for pose in poses:
            if path in manual_overrides:
                resolutions = [resolve_manual(pose, manual_overrides[path])]
            else:
                resolutions = resolve(pose, import_stem=stem)

            for res in resolutions:
                if res.needs_confirmation and not confirm_fallback:
                    result.pending.append(Pending(pose, res, res.preview, path))
                    continue
                if not res.ligand_atoms:
                    continue
                set_sides(res)
                inters = compute_interactions(
                    res.receptor_atoms, res.ligand_atoms, res.waters, types
                )
                details = [
                    _detail_from_interaction(it, res.ligand_id, path, key_set)
                    for it in inters
                ]
                pose_no = pose.pose_index + 1
                summ = _summarize(
                    res.ligand_id, path, pose.sol, pose_no, pose.score, details, key_set
                )
                result.details.extend(details)
                result.summaries.append(summ)

    return result


def recompute_key(result: RunResult, key_residues):
    """Re-evaluate key-residue flags/counts WITHOUT redoing detection."""
    key_set = normalize_key_residues(key_residues or [])
    result.key_residues = key_set
    for d in result.details:
        d.is_key_residue = _is_key(d.receptor_residue, d._res_nochain, key_set)
    per_ligand = {}
    for d in result.details:
        if d.is_key_residue:
            per_ligand[(d.ligand_id, d.source_file)] = (
                per_ligand.get((d.ligand_id, d.source_file), 0) + 1
            )
    for s in result.summaries:
        s.n_key_residue_interactions = per_ligand.get((s.ligand_id, s.source_file), 0)
    return result
