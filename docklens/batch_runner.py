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
from dataclasses import dataclass
from datetime import datetime, timezone

from .entity_resolver import resolve, resolve_manual, set_sides
from .interaction_core import (
    HBOND_PRESETS,
    VALID_TYPES,
    compute_interactions,
    cutoffs_for_preset,
    endpoint_name,
    endpoint_resid,
    endpoint_side,
)
from . import __version__
from .parser_mol2 import parse_mol2
from .parser_pdb import parse_pdb
from .parser_pdbqt import parse_pdbqt
from .results import (
    AnalysisParameters,
    Detail,
    Endpoint,
    ExportFilter as ExportFilter,
    InputQC,
    RunResult,
    Summary,
    make_result,
    with_key_residues,
)

_PARSERS = {".mol2": parse_mol2, ".pdb": parse_pdb, ".pdbqt": parse_pdbqt}
SUPPORTED_EXT = tuple(_PARSERS)
DEFAULT_MAX_FILE_SIZE_BYTES = 256 * 1024 * 1024


@dataclass(frozen=True)
class Pending:
    """A pose whose ligand/receptor split needs user confirmation (fallback)."""

    pose: object
    resolution: object
    preview: str
    source_file: str
    source_id: str = ""
    pose_id: str = ""
    resolution_index: int = 1


@dataclass(frozen=True)
class _InputCandidate:
    """One supported source or one manifest-only input problem."""

    path: str
    code: str = ""
    status: str = "success"
    message: str = ""


# ---------------------------------------------------------------------------
# File gathering + parsing
# ---------------------------------------------------------------------------


def _coerce_paths(paths):
    if isinstance(paths, (str, os.PathLike)):
        paths = [paths]
    if paths is None:
        raise ValueError("paths must contain at least one input")
    try:
        normalized = [os.path.abspath(os.fspath(path)) for path in paths]
    except TypeError as exc:
        raise ValueError("paths must contain filesystem paths") from exc
    if not normalized:
        raise ValueError("paths must contain at least one input")
    return normalized


def _gather_inputs(paths):
    """Expand inputs while retaining direct-input and empty-folder failures."""
    candidates = []
    for path in _coerce_paths(paths):
        if os.path.isdir(path):
            discovered = []
            for root, _dirs, files in os.walk(path):
                for name in files:
                    if name.lower().endswith(SUPPORTED_EXT):
                        discovered.append(os.path.join(root, name))
            if discovered:
                candidates.extend(_InputCandidate(item) for item in discovered)
            else:
                candidates.append(
                    _InputCandidate(
                        path,
                        code="no_supported_files",
                        status="warning",
                        message="Directory contains no supported structure files.",
                    )
                )
        elif not os.path.exists(path):
            candidates.append(
                _InputCandidate(
                    path,
                    code="missing_input",
                    status="error",
                    message="Input path does not exist.",
                )
            )
        elif not os.path.isfile(path) or not path.lower().endswith(SUPPORTED_EXT):
            candidates.append(
                _InputCandidate(
                    path,
                    code="unsupported_input",
                    status="error",
                    message="Input is not a supported MOL2, PDB, or PDBQT file.",
                )
            )
        else:
            candidates.append(_InputCandidate(path))

    unique = {}
    for candidate in candidates:
        unique.setdefault(os.path.normcase(candidate.path), candidate)
    return sorted(unique.values(), key=lambda item: item.path.lower())


def gather_files(paths):
    """Return the supported files discovered from the supplied paths."""
    return [candidate.path for candidate in _gather_inputs(paths) if not candidate.code]


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


def _normalize_types(types):
    if types is None:
        return tuple(VALID_TYPES)
    if isinstance(types, str):
        types = types.replace(",", " ").split()
    normalized = []
    for item in types:
        kind = str(item).strip().lower()
        if kind and kind not in normalized:
            normalized.append(kind)
    if not normalized:
        raise ValueError("At least one interaction type is required")
    unknown = [kind for kind in normalized if kind not in VALID_TYPES]
    if unknown:
        raise ValueError("Unknown interaction type: %s" % ", ".join(unknown))
    return tuple(normalized)


def _normalize_preset(value):
    preset = str(value).strip().lower()
    if preset not in HBOND_PRESETS:
        raise ValueError("Unknown H-bond preset: %s" % preset)
    return preset


def _safe_exception_message(action, exc):
    """Describe a failure without copying attacker-controlled exception text."""
    return "%s (%s)." % (action, type(exc).__name__)


def _safe_qc_text(value, limit=500):
    """Make resolver diagnostics single-line, bounded and spreadsheet-safe."""
    text = " ".join(str(value or "").split())
    text = "".join(char for char in text if char.isprintable())
    if len(text) > limit:
        text = text[: limit - 1] + "…"
    if text.lstrip().startswith(("=", "+", "-", "@")):
        text = "'" + text
    return text


def _is_key(res_tag, res_nochain, key_set):
    if not key_set:
        return False
    return res_tag.upper() in key_set or res_nochain.upper() in key_set


# ---------------------------------------------------------------------------
# Detection driver
# ---------------------------------------------------------------------------


def _endpoint(obj, role=""):
    atoms = tuple(obj.atoms) if hasattr(obj, "atoms") else (obj,)
    first = atoms[0]
    serials = tuple(sorted(a.serial for a in atoms if a.serial is not None))
    return Endpoint(
        side=endpoint_side(obj) or "",
        kind="ring" if hasattr(obj, "atoms") else "atom",
        atom_name=endpoint_name(obj),
        atom_serials=serials,
        resname=first.resn,
        resseq=first.resi,
        chain=first.chain,
        element=first.elem,
        role=role or "",
    )


def _detail_from_interaction(
    it,
    ligand_id,
    source_file,
    key_set,
    *,
    source_id="S000000",
    pose_id="S000000:P0001:R001",
    interaction_index=1,
    pose_no=1,
    sol=None,
    score=None,
    resolution_method="",
):
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
    if lig is a:
        lig_role, rec_role = it.get("a_role", ""), it.get("b_role", "")
    else:
        lig_role, rec_role = it.get("b_role", ""), it.get("a_role", "")
    water_obj = it.get("water_obj")
    return Detail(
        ligand_id=ligand_id,
        source_file=os.path.basename(source_file),
        interaction_type=it["type"],
        subtype=it.get("subtype", ""),
        ligand=_endpoint(lig, lig_role),
        receptor=_endpoint(rec, rec_role),
        distance_A=float(it["dist"]) if it.get("dist") is not None else None,
        source_id=source_id,
        pose_id=pose_id,
        interaction_id="%s:I%06d" % (pose_id, interaction_index),
        pose=pose_no,
        sol=sol,
        docking_score=score,
        source_path=os.path.abspath(source_file),
        resolution_method=resolution_method,
        is_key_residue=_is_key(rec_res, res_nochain, key_set),
        water=_endpoint(water_obj, "bridge") if water_obj is not None else None,
        receptor_water_distance_A=it.get("receptor_water_distance"),
        ligand_water_distance_A=it.get("ligand_water_distance"),
        water_angle_deg=it.get("water_angle"),
    )


def _summarize(
    ligand_id,
    source_file,
    sol,
    pose_no,
    score,
    details,
    key_set,
    *,
    source_id="",
    pose_id="",
    resolution_method="",
):
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
        source_id=source_id,
        pose_id=pose_id,
        source_path=os.path.abspath(source_file),
        resolution_method=resolution_method,
    )


def run(
    paths,
    types=None,
    key_residues=None,
    confirm_fallback=False,
    manual_overrides=None,
    hbond_preset="plip",
    max_file_size_bytes=DEFAULT_MAX_FILE_SIZE_BYTES,
):
    """Run detection over the given inputs.

    types            list of interaction types (default: all).
    key_residues     iterable / string of key residue tags.
    confirm_fallback if True, fallback resolutions are run anyway (headless);
                     if False, they are collected in RunResult.pending instead.
    manual_overrides {source_file: set_of_ligand_serials} to force a split.
    hbond_preset     'plip' (default, permissive) or 'dsv' (Discovery-Studio-like,
                     stricter H-bond distance/angle).
    """
    if isinstance(max_file_size_bytes, bool) or not isinstance(
        max_file_size_bytes, int
    ) or max_file_size_bytes <= 0:
        raise ValueError("max_file_size_bytes must be a positive integer")
    requested_types = _normalize_types(types)
    hbond_preset = _normalize_preset(hbond_preset)
    key_set = normalize_key_residues(key_residues or [])
    manual_overrides = manual_overrides or {}
    details_out = []
    summaries_out = []
    pending_out = []
    qc_out = []
    receptor_residues = set()
    effective_cutoffs = cutoffs_for_preset(hbond_preset)
    parameters = AnalysisParameters(
        app_version=__version__,
        started_at=datetime.now(timezone.utc).isoformat(),
        hbond_preset=hbond_preset,
        cutoffs=tuple(
            sorted((key, float(value)) for key, value in effective_cutoffs.items())
        ),
        interaction_types=requested_types,
        key_residues=tuple(sorted(key_set)),
    )

    for source_number, candidate in enumerate(_gather_inputs(paths), 1):
        path = candidate.path
        source_id = "S%06d" % source_number
        stem = os.path.splitext(os.path.basename(path))[0]
        fmt = os.path.splitext(path)[1].lstrip(".").lower()
        if candidate.code:
            qc_out.append(
                InputQC(
                    source_id=source_id,
                    source_file=os.path.basename(path),
                    source_path=path,
                    status=candidate.status,
                    code=candidate.code,
                    message=candidate.message,
                    format=fmt,
                )
            )
            continue
        try:
            file_size = os.path.getsize(path)
        except OSError as exc:
            qc_out.append(
                InputQC(
                    source_id=source_id,
                    source_file=os.path.basename(path),
                    source_path=path,
                    status="error",
                    code="file_stat_error",
                    message=_safe_exception_message("Could not inspect input file", exc),
                    format=fmt,
                )
            )
            continue
        if file_size > max_file_size_bytes:
            qc_out.append(
                InputQC(
                    source_id=source_id,
                    source_file=os.path.basename(path),
                    source_path=path,
                    status="error",
                    code="file_too_large",
                    message="Input exceeds the configured size limit (%d bytes)."
                    % max_file_size_bytes,
                    format=fmt,
                )
            )
            continue
        try:
            poses = parse_file(path)
        except Exception as exc:  # noqa: BLE001 - report, keep going
            qc_out.append(
                InputQC(
                    source_id=source_id,
                    source_file=os.path.basename(path),
                    source_path=os.path.abspath(path),
                    status="error",
                    code="parse_error",
                    message=_safe_exception_message("Could not parse input file", exc),
                    format=fmt,
                )
            )
            continue
        if not poses:
            qc_out.append(
                InputQC(
                    source_id=source_id,
                    source_file=os.path.basename(path),
                    source_path=os.path.abspath(path),
                    status="warning",
                    code="no_poses",
                    message="No structural poses were found.",
                    format=fmt,
                )
            )
            continue

        for pose in poses:
            pose_no = pose.pose_index + 1
            try:
                if path in manual_overrides:
                    resolutions = [resolve_manual(pose, manual_overrides[path])]
                else:
                    resolutions = resolve(pose, import_stem=stem)
            except Exception as exc:  # noqa: BLE001 - isolate malformed poses
                qc_out.append(
                    InputQC(
                        source_id=source_id,
                        source_file=os.path.basename(path),
                        source_path=path,
                        pose_id="%s:P%04d:R000" % (source_id, pose_no),
                        status="error",
                        code="pose_error",
                        message=_safe_exception_message("Could not process pose", exc),
                        format=fmt,
                        poses_found=len(poses),
                    )
                )
                continue

            for resolution_index, res in enumerate(resolutions, 1):
                pose_id = "%s:P%04d:R%03d" % (
                    source_id,
                    pose_no,
                    resolution_index,
                )
                if res.needs_confirmation and not confirm_fallback:
                    pending_out.append(
                        Pending(
                            pose,
                            res,
                            res.preview,
                            path,
                            source_id,
                            pose_id,
                            resolution_index,
                        )
                    )
                    qc_out.append(
                        InputQC(
                            source_id=source_id,
                            source_file=os.path.basename(path),
                            source_path=os.path.abspath(path),
                            pose_id=pose_id,
                            status="pending",
                            code="confirmation_required",
                            message=_safe_qc_text(res.preview),
                            format=fmt,
                            poses_found=len(poses),
                            resolution_method=res.method,
                            receptor_atoms=len(res.receptor_atoms),
                            ligand_atoms=len(res.ligand_atoms),
                            water_atoms=len(res.waters),
                            warnings=tuple(
                                _safe_qc_text(item) for item in res.warnings
                            ),
                        )
                    )
                    continue
                if not res.ligand_atoms:
                    qc_out.append(
                        InputQC(
                            source_id=source_id,
                            source_file=os.path.basename(path),
                            source_path=os.path.abspath(path),
                            pose_id=pose_id,
                            status="warning",
                            code="no_ligand",
                            message="No ligand atoms were identified.",
                            format=fmt,
                            poses_found=len(poses),
                            resolution_method=res.method,
                            receptor_atoms=len(res.receptor_atoms),
                            water_atoms=len(res.waters),
                        )
                    )
                    continue
                try:
                    set_sides(res)
                    residue_tags = {a.res_tag() for a in res.receptor_atoms}
                    inters = compute_interactions(
                        res.receptor_atoms,
                        res.ligand_atoms,
                        res.waters,
                        requested_types,
                        cutoffs=effective_cutoffs,
                    )
                    details = [
                        _detail_from_interaction(
                            it,
                            res.ligand_id,
                            path,
                            key_set,
                            source_id=source_id,
                            pose_id=pose_id,
                            interaction_index=index,
                            pose_no=pose_no,
                            sol=pose.sol,
                            score=pose.score,
                            resolution_method=res.method,
                        )
                        for index, it in enumerate(inters, 1)
                    ]
                    summ = _summarize(
                        res.ligand_id,
                        path,
                        pose.sol,
                        pose_no,
                        pose.score,
                        details,
                        key_set,
                        source_id=source_id,
                        pose_id=pose_id,
                        resolution_method=res.method,
                    )
                except Exception as exc:  # noqa: BLE001 - isolate one resolution
                    qc_out.append(
                        InputQC(
                            source_id=source_id,
                            source_file=os.path.basename(path),
                            source_path=path,
                            pose_id=pose_id,
                            status="error",
                            code="pose_error",
                            message=_safe_exception_message(
                                "Could not process pose", exc
                            ),
                            format=fmt,
                            poses_found=len(poses),
                            resolution_method=_safe_qc_text(
                                getattr(res, "method", "")
                            ),
                        )
                    )
                    continue
                receptor_residues.update(residue_tags)
                details_out.extend(details)
                summaries_out.append(summ)
                qc_out.append(
                    InputQC(
                        source_id=source_id,
                        source_file=os.path.basename(path),
                        source_path=os.path.abspath(path),
                        pose_id=pose_id,
                        status="warning" if res.warnings else "success",
                        code="resolution_warning" if res.warnings else "",
                        message=_safe_qc_text("; ".join(res.warnings)),
                        format=fmt,
                        poses_found=len(poses),
                        poses_processed=1,
                        resolution_method=res.method,
                        receptor_atoms=len(res.receptor_atoms),
                        ligand_atoms=len(res.ligand_atoms),
                        water_atoms=len(res.waters),
                        warnings=tuple(_safe_qc_text(item) for item in res.warnings),
                    )
                )

    return make_result(
        details=details_out,
        summaries=summaries_out,
        pending=pending_out,
        key_residues=key_set,
        receptor_residues=receptor_residues,
        input_qc=qc_out,
        parameters=parameters,
    )


def recompute_key(result: RunResult, key_residues):
    """Re-evaluate key-residue flags/counts WITHOUT redoing detection."""
    return with_key_residues(result, key_residues or [])
