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

from .entity_resolver import (
    Resolution,
    _split_waters,
    resolve,
    resolve_manual,
    set_sides,
)
from .interaction_core import (
    HBOND_PRESETS,
    VALID_TYPES,
    compute_interactions,
    cutoffs_for_preset,
    default_types_for_profile,
    endpoint_name,
    endpoint_resid,
    endpoint_side,
)
from . import __version__
from .input_plan import InputJob, InputPlan, run_shared_receptor_plan
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
from .residue_keys import normalize_key_residues as _normalize_key_residues

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
    return set(_normalize_key_residues(items))


def _normalize_types(types, profile="plip"):
    if types is None:
        return tuple(default_types_for_profile(profile))
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
        raise ValueError("Unknown scientific profile: %s" % preset)
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
    source_pose_label="",
    source_molecule_name="",
    score_type="",
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
    hydrogen_obj = it.get("hydrogen_obj")
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
        source_pose_label=source_pose_label,
        source_molecule_name=source_molecule_name,
        score_type=score_type,
        is_key_residue=_is_key(rec_res, res_nochain, key_set),
        water=_endpoint(water_obj, "bridge") if water_obj is not None else None,
        receptor_water_distance_A=it.get("receptor_water_distance"),
        ligand_water_distance_A=it.get("ligand_water_distance"),
        water_angle_deg=it.get("water_angle"),
        chemistry_basis=it.get("chemistry_basis", ""),
        chemistry_confidence=it.get("confidence", ""),
        hydrogen_atom=hydrogen_obj.name if hydrogen_obj is not None else "",
        hydrogen_atom_serial=(
            hydrogen_obj.serial if hydrogen_obj is not None else None
        ),
        hydrogen_acceptor_distance_A=it.get("hydrogen_acceptor_distance"),
        donor_hydrogen_acceptor_angle_deg=it.get("donor_hydrogen_acceptor_angle"),
        hydrogen_acceptor_base_angle_deg=it.get("hydrogen_acceptor_base_angle"),
        theta_deg=it.get("theta"),
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
    source_pose_label="",
    source_molecule_name="",
    score_type="",
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
        source_pose_label=source_pose_label,
        source_molecule_name=source_molecule_name,
        score_type=score_type,
    )


def _plan_from_gathered_inputs(paths):
    candidates = _gather_inputs(paths)
    return InputPlan(
        InputJob(
            kind="combined",
            source_path=candidate.path,
            group_id="combined_%d" % index,
            input_status=candidate.status,
            input_code=candidate.code,
            input_message=candidate.message,
        )
        for index, candidate in enumerate(candidates, 1)
    )


def _validate_plan_size(max_file_size_bytes):
    if (
        isinstance(max_file_size_bytes, bool)
        or not isinstance(max_file_size_bytes, int)
        or max_file_size_bytes <= 0
    ):
        raise ValueError("max_file_size_bytes must be a positive integer")


def _parameters(types, key_set, preset):
    effective_cutoffs = cutoffs_for_preset(preset)
    return AnalysisParameters(
        app_version=__version__,
        started_at=datetime.now(timezone.utc).isoformat(),
        hbond_preset=preset,
        cutoffs=tuple(
            sorted((key, float(value)) for key, value in effective_cutoffs.items())
        ),
        interaction_types=types,
        key_residues=tuple(sorted(key_set)),
    )


def _qc_common(job, source_id, path, fmt, **extra):
    return dict(
        source_id=source_id,
        source_file=os.path.basename(path),
        source_path=os.path.abspath(path),
        format=fmt,
        input_mode=job.kind,
        receptor_source_file=job.receptor_source_file,
        receptor_source_path=job.receptor_path or "",
        group_id=job.group_id,
        **extra,
    )


def _source_error(job, source_id, max_file_size_bytes):
    path = job.source_path
    fmt = os.path.splitext(path)[1].lstrip(".").lower()
    common = _qc_common(job, source_id, path, fmt)
    if job.input_code:
        return InputQC(
            **common,
            status=job.input_status,
            code=job.input_code,
            message=_safe_qc_text(job.input_message),
        )
    if not os.path.exists(path):
        return InputQC(
            **common,
            status="error",
            code="missing_input",
            message="Input path does not exist.",
        )
    if not os.path.isfile(path) or not path.lower().endswith(SUPPORTED_EXT):
        return InputQC(
            **common,
            status="error",
            code="unsupported_input",
            message="Input is not a supported MOL2, PDB, or PDBQT file.",
        )
    try:
        size = os.path.getsize(path)
    except OSError as exc:
        return InputQC(
            **common,
            status="error",
            code="file_stat_error",
            message=_safe_exception_message("Could not inspect input file", exc),
        )
    if size > max_file_size_bytes:
        return InputQC(
            **common,
            status="error",
            code="file_too_large",
            message="Input exceeds the configured size limit (%d bytes)."
            % max_file_size_bytes,
        )
    return None


def _load_receptor(job, max_file_size_bytes, cache):
    path = job.receptor_path
    key = os.path.normcase(path or "")
    if key in cache:
        return cache[key]
    if not path or not os.path.isfile(path) or not path.lower().endswith(SUPPORTED_EXT):
        value = (
            None,
            None,
            "receptor_missing",
            "Receptor is not a supported structure file.",
        )
        cache[key] = value
        return value
    try:
        if os.path.getsize(path) > max_file_size_bytes:
            raise OverflowError
        poses = parse_file(path)
        if not poses:
            value = (
                None,
                None,
                "receptor_no_models",
                "Receptor contains no structural model.",
            )
        elif len(poses) != 1:
            value = (
                None,
                None,
                "receptor_multiple_models",
                "Paired receptor must contain exactly one structural model.",
            )
        else:
            atoms, waters = _split_waters(poses[0].atoms)
            value = (atoms, waters, "", "")
    except OverflowError:
        value = (
            None,
            None,
            "receptor_file_too_large",
            "Receptor exceeds the configured size limit (%d bytes)."
            % max_file_size_bytes,
        )
    except Exception as exc:  # noqa: BLE001 - isolate dependent jobs
        value = (
            None,
            None,
            "receptor_parse_error",
            _safe_exception_message("Could not parse receptor", exc),
        )
    cache[key] = value
    return value


def _run_resolution(
    *,
    pose,
    resolution,
    path,
    source_id,
    pose_id,
    pose_no,
    ligand_id,
    source_pose_label,
    source_molecule_name,
    score_type,
    key_set,
    requested_types,
    effective_cutoffs,
    hbond_preset,
):
    set_sides(resolution)
    interactions = compute_interactions(
        resolution.receptor_atoms,
        resolution.ligand_atoms,
        resolution.waters,
        requested_types,
        cutoffs=effective_cutoffs,
        chemistry_profile=hbond_preset,
    )
    details = [
        _detail_from_interaction(
            item,
            ligand_id,
            path,
            key_set,
            source_id=source_id,
            pose_id=pose_id,
            interaction_index=index,
            pose_no=pose_no,
            sol=pose.sol,
            score=pose.score,
            resolution_method=resolution.method,
            source_pose_label=source_pose_label,
            source_molecule_name=source_molecule_name,
            score_type=score_type,
        )
        for index, item in enumerate(interactions, 1)
    ]
    summary = _summarize(
        ligand_id,
        path,
        pose.sol,
        pose_no,
        pose.score,
        details,
        key_set,
        source_id=source_id,
        pose_id=pose_id,
        resolution_method=resolution.method,
        source_pose_label=source_pose_label,
        source_molecule_name=source_molecule_name,
        score_type=score_type,
    )
    return summary, details


def run_plan(
    input_plan,
    types=None,
    key_residues=None,
    confirm_fallback=False,
    manual_overrides=None,
    hbond_preset="plip",
    max_file_size_bytes=DEFAULT_MAX_FILE_SIZE_BYTES,
):
    """Run all jobs in one immutable plan through the generic core."""

    if not isinstance(input_plan, InputPlan):
        raise TypeError("input_plan must be an InputPlan")
    _validate_plan_size(max_file_size_bytes)
    preset = _normalize_preset(hbond_preset)
    requested_types = _normalize_types(types, preset)
    key_set = normalize_key_residues(key_residues or [])
    manual_overrides = manual_overrides or {}
    effective_cutoffs = cutoffs_for_preset(preset)
    details_out, summaries_out, pending_out, qc_out = [], [], [], []
    receptor_residues = set()
    receptor_cache = {}
    parameters = _parameters(requested_types, key_set, preset)

    for source_number, job in enumerate(input_plan.jobs, 1):
        source_id = "S%06d" % source_number
        path = job.source_path
        fmt = os.path.splitext(path)[1].lstrip(".").lower()
        source_error = _source_error(job, source_id, max_file_size_bytes)
        if source_error is not None:
            qc_out.append(source_error)
            continue

        receptor_atoms = receptor_waters = None
        if job.kind == "paired":
            receptor_atoms, receptor_waters, code, message = _load_receptor(
                job, max_file_size_bytes, receptor_cache
            )
            if code:
                qc_out.append(
                    InputQC(
                        **_qc_common(
                            job,
                            source_id,
                            path,
                            fmt,
                            status="error",
                            code=code,
                            message=message,
                        )
                    )
                )
                continue
            receptor_residues.update(atom.res_tag() for atom in receptor_atoms)

        try:
            poses = parse_file(path)
        except Exception as exc:  # noqa: BLE001 - isolate one source file
            qc_out.append(
                InputQC(
                    **_qc_common(
                        job,
                        source_id,
                        path,
                        fmt,
                        status="error",
                        code="parse_error",
                        message=_safe_exception_message("Could not parse input file", exc),
                    )
                )
            )
            continue
        if not poses:
            qc_out.append(
                InputQC(
                    **_qc_common(
                        job,
                        source_id,
                        path,
                        fmt,
                        status="warning",
                        code="no_poses",
                        message="No structural poses were found.",
                    )
                )
            )
            continue

        stem = os.path.splitext(os.path.basename(path))[0]
        for pose in poses:
            pose_no = pose.pose_index + 1
            pose_id_prefix = "%s:P%04d" % (source_id, pose_no)
            source_molecule_name = getattr(pose, "molecule_name", "") or ""
            source_pose_label = getattr(pose, "source_pose_label", "") or ""
            score_type = getattr(pose, "score_type", "") or ""
            try:
                if job.kind == "paired":
                    ligand_id = (
                        getattr(pose, "ligand_id_hint", "")
                        or source_molecule_name
                        or "%s_pose_%04d" % (stem, pose_no)
                    )
                    resolutions = [
                        Resolution(
                            receptor_atoms,
                            list(pose.atoms),
                            receptor_waters,
                            ligand_id,
                            job.resolution_method or "paired-plan",
                        )
                    ]
                else:
                    override = manual_overrides.get(path)
                    resolutions = (
                        [resolve_manual(pose, override)]
                        if override is not None
                        else resolve(pose, import_stem=stem)
                    )
            except Exception as exc:  # noqa: BLE001 - isolate malformed poses
                qc_out.append(
                    InputQC(
                        **_qc_common(
                            job,
                            source_id,
                            path,
                            fmt,
                            pose_id=pose_id_prefix + ":R000",
                            status="error",
                            code="pose_error",
                            message=_safe_exception_message("Could not process pose", exc),
                            poses_found=len(poses),
                            source_pose_label=source_pose_label,
                            source_molecule_name=source_molecule_name,
                            docking_score=pose.score,
                            score_type=score_type,
                        )
                    )
                )
                continue

            for resolution_index, resolution in enumerate(resolutions, 1):
                pose_id = "%s:R%03d" % (pose_id_prefix, resolution_index)
                common = _qc_common(
                    job,
                    source_id,
                    path,
                    fmt,
                    pose_id=pose_id,
                    poses_found=len(poses),
                    resolution_method=resolution.method,
                    receptor_atoms=len(resolution.receptor_atoms),
                    ligand_atoms=len(resolution.ligand_atoms),
                    water_atoms=len(resolution.waters),
                    warnings=tuple(_safe_qc_text(item) for item in resolution.warnings),
                    ligand_id=(
                        ligand_id
                        if job.kind == "paired"
                        else resolution.ligand_id
                    ),
                    source_pose_label=source_pose_label,
                    source_molecule_name=source_molecule_name,
                    docking_score=pose.score,
                    score_type=score_type,
                )
                if resolution.needs_confirmation and not confirm_fallback:
                    pending_out.append(
                        Pending(
                            pose,
                            resolution,
                            resolution.preview,
                            path,
                            source_id,
                            pose_id,
                            resolution_index,
                        )
                    )
                    qc_out.append(
                        InputQC(
                            **common,
                            status="pending",
                            code="confirmation_required",
                            message=_safe_qc_text(resolution.preview),
                        )
                    )
                    continue
                if not resolution.ligand_atoms:
                    qc_out.append(
                        InputQC(
                            **common,
                            status="warning",
                            code="no_ligand",
                            message="No ligand atoms were identified.",
                        )
                    )
                    continue
                try:
                    summary, details = _run_resolution(
                        pose=pose,
                        resolution=resolution,
                        path=path,
                        source_id=source_id,
                        pose_id=pose_id,
                        pose_no=pose_no,
                        ligand_id=(
                            ligand_id
                            if job.kind == "paired"
                            else resolution.ligand_id
                        ),
                        source_pose_label=source_pose_label,
                        source_molecule_name=source_molecule_name,
                        score_type=score_type,
                        key_set=key_set,
                        requested_types=requested_types,
                        effective_cutoffs=effective_cutoffs,
                        hbond_preset=preset,
                    )
                except Exception as exc:  # noqa: BLE001 - isolate one resolution
                    qc_out.append(
                        InputQC(
                            **common,
                            status="error",
                            code="pose_error",
                            message=_safe_exception_message("Could not process pose", exc),
                        )
                    )
                    continue
                details_out.extend(details)
                summaries_out.append(summary)
                receptor_residues.update(
                    atom.res_tag() for atom in resolution.receptor_atoms
                )
                qc_out.append(
                    InputQC(
                        **common,
                        status="warning" if resolution.warnings else "success",
                        code="resolution_warning" if resolution.warnings else "",
                        message=_safe_qc_text("; ".join(resolution.warnings)),
                        poses_processed=1,
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


def run(
    paths,
    types=None,
    key_residues=None,
    confirm_fallback=False,
    manual_overrides=None,
    hbond_preset="plip",
    max_file_size_bytes=DEFAULT_MAX_FILE_SIZE_BYTES,
):
    """Backward-compatible wrapper that treats each input as combined."""

    return run_plan(
        _plan_from_gathered_inputs(paths),
        types=types,
        key_residues=key_residues,
        confirm_fallback=confirm_fallback,
        manual_overrides=manual_overrides,
        hbond_preset=hbond_preset,
        max_file_size_bytes=max_file_size_bytes,
    )


def run_shared_receptor(
    receptor_path,
    ligand_paths,
    types=None,
    key_residues=None,
    confirm_fallback=False,
    manual_overrides=None,
    hbond_preset="plip",
    max_file_size_bytes=DEFAULT_MAX_FILE_SIZE_BYTES,
):
    """Analyze multiple ligand/pose files against one cached receptor."""

    return run_plan(
        run_shared_receptor_plan(receptor_path, ligand_paths),
        types=types,
        key_residues=key_residues,
        confirm_fallback=confirm_fallback,
        manual_overrides=manual_overrides,
        hbond_preset=hbond_preset,
        max_file_size_bytes=max_file_size_bytes,
    )


def run_paired(
    receptor_path,
    poses_path,
    types=None,
    key_residues=None,
    hbond_preset="plip",
    max_file_size_bytes=DEFAULT_MAX_FILE_SIZE_BYTES,
):
    """Compatibility wrapper for the DockingHub receptor/poses contract."""

    receptor_path, poses_path = _coerce_paths([receptor_path, poses_path])
    plan = InputPlan(
        (
            InputJob(
                kind="paired",
                source_path=poses_path,
                receptor_path=receptor_path,
                group_id="dockhub",
                resolution_method="paired-manifest",
            ),
        )
    )
    return run_plan(
        plan,
        types=types,
        key_residues=key_residues,
        hbond_preset=hbond_preset,
        max_file_size_bytes=max_file_size_bytes,
    )


def recompute_key(result: RunResult, key_residues):
    """Re-evaluate key-residue flags/counts WITHOUT redoing detection."""
    return with_key_residues(result, key_residues or [])
