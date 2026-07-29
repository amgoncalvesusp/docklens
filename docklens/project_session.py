"""Safe, reproducible DockLens project containers.

``.docklens`` files are ZIP containers with a manifest, a methods record and
separate JSON members for completed ``RunResult`` datasets. Inputs are
referenced by absolute path and SHA-256 rather than copied into the container.
"""

from __future__ import annotations

from dataclasses import dataclass, fields
import hashlib
import json
import math
import os
from pathlib import Path, PurePosixPath
import re
import tempfile
from typing import Any, Mapping
import zipfile

from .observation_series import ObservationPoint, ObservationSeries
from .ligand_selection import ligand_groups
from .results import (
    AnalysisParameters,
    Detail,
    Endpoint,
    InputQC,
    RunResult,
    Summary,
)

PROJECT_SCHEMA_VERSION = "3"
SUPPORTED_PROJECT_SCHEMA_VERSIONS = frozenset({"1", "2", "3"})
MAX_PROJECT_BYTES = 5 * 1024 * 1024
MAX_UNCOMPRESSED_BYTES = 32 * 1024 * 1024
MAX_ARCHIVE_ENTRIES = 8
MAX_JSON_DEPTH = 40
_FIXED_ENTRIES = frozenset({"manifest.json", "methods.txt"})


def is_local_filesystem_path(value: str) -> bool:
    """Reject network/URI forms without touching the referenced location."""
    path = str(value)
    return bool(path) and not (
        path.startswith(("\\\\", "//"))
        or re.match(r"^[A-Za-z][A-Za-z0-9+.-]*://", path)
    )


class ProjectIntegrityError(ValueError):
    """Raised when a project member differs from its signed manifest metadata."""


@dataclass(frozen=True)
class ProjectInput:
    """Provenance for one external input without embedding its contents."""

    path: str
    sha256: str
    size_bytes: int
    modified_ns: int

    def __post_init__(self) -> None:
        path = str(self.path)
        if not is_local_filesystem_path(path):
            raise ValueError("project inputs must use local filesystem paths")
        digest = str(self.sha256).lower()
        if not _valid_sha(digest):
            raise ValueError("project input SHA-256 is invalid")
        if int(self.size_bytes) < 0 or int(self.modified_ns) < 0:
            raise ValueError("project input metadata must be non-negative")
        object.__setattr__(self, "path", path)
        object.__setattr__(self, "sha256", digest)
        object.__setattr__(self, "size_bytes", int(self.size_bytes))
        object.__setattr__(self, "modified_ns", int(self.modified_ns))

    @property
    def source_path(self) -> str:
        """Compatibility alias used by UI code."""

        return self.path


@dataclass(frozen=True)
class ProjectDataset:
    """One analyzed dataset and, when available, its cached immutable result."""

    label: str
    mode: str
    time_step_ns: float | None
    inputs: tuple[ProjectInput, ...] = ()
    result: RunResult | None = None
    observation_series: ObservationSeries | None = None

    def __post_init__(self) -> None:
        if self.mode not in {"docking", "md"}:
            raise ValueError("dataset mode must be 'docking' or 'md'")
        if self.time_step_ns is not None:
            value = float(self.time_step_ns)
            if not math.isfinite(value) or value <= 0:
                raise ValueError("time_step_ns must be a positive finite number")
            object.__setattr__(self, "time_step_ns", value)
        object.__setattr__(self, "inputs", tuple(self.inputs))
        if (
            self.observation_series is not None
            and self.observation_series.mode != self.mode
        ):
            raise ValueError(
                "observation series mode must match the dataset mode"
            )
        if self.result is not None and self.observation_series is not None:
            result_ids = {summary.pose_id for summary in self.result.summaries}
            if set(self.observation_series.observation_ids) != result_ids:
                raise ValueError(
                    "observation series must match the cached result IDs"
                )


@dataclass(frozen=True)
class ProjectState:
    """Settings needed to reproduce or resume a DockLens analysis."""

    app_version: str
    analysis_profile: str
    hbond_preset: str
    key_residues: tuple[str, ...]
    selected_types: tuple[str, ...]
    active_workspace: str
    selected_residue: str
    state_threshold: float
    bootstrap_iterations: int
    primary: ProjectDataset
    comparison: ProjectDataset | None = None
    bootstrap_block_size: int | None = None
    bootstrap_seed: int = 2026
    confidence_level: float = 0.95
    primary_ligand_group: str | None = None
    comparison_ligand_group: str | None = None
    observation_label_mode: str = "ligand"
    heatmap_group_by: str = "source"
    heatmap_feature_level: str = "residue_type"
    heatmap_top_n: int | None = 40

    def __post_init__(self) -> None:
        threshold = float(self.state_threshold)
        if not math.isfinite(threshold) or not 0 <= threshold <= 1:
            raise ValueError("state_threshold must be between 0 and 1")
        iterations = int(self.bootstrap_iterations)
        if iterations < 0 or iterations > 1_000_000:
            raise ValueError("bootstrap_iterations is outside the supported range")
        object.__setattr__(self, "state_threshold", threshold)
        object.__setattr__(self, "bootstrap_iterations", iterations)
        block_size = self.bootstrap_block_size
        if block_size is not None:
            block_size = int(block_size)
            if block_size < 1:
                raise ValueError(
                    "bootstrap_block_size must be at least one"
                )
        seed = int(self.bootstrap_seed)
        if seed < 0:
            raise ValueError("bootstrap_seed must be non-negative")
        confidence = float(self.confidence_level)
        if not math.isfinite(confidence) or not 0 < confidence < 1:
            raise ValueError("confidence_level must be between zero and one")
        object.__setattr__(self, "bootstrap_block_size", block_size)
        object.__setattr__(self, "bootstrap_seed", seed)
        object.__setattr__(self, "confidence_level", confidence)
        for name in ("primary_ligand_group", "comparison_ligand_group"):
            value = getattr(self, name)
            if value is not None and not str(value):
                raise ValueError(f"{name} must be non-empty or None")
            object.__setattr__(
                self, name, None if value is None else str(value)
            )
        if self.observation_label_mode not in {"ligand", "file", "index"}:
            raise ValueError("observation_label_mode is not supported")
        if self.heatmap_group_by not in {"source", "observation"}:
            raise ValueError("heatmap_group_by is not supported")
        if self.heatmap_feature_level not in {"residue_type", "residue"}:
            raise ValueError("heatmap_feature_level is not supported")
        if self.heatmap_top_n is not None:
            top_n = int(self.heatmap_top_n)
            if top_n < 1 or top_n > 10_000:
                raise ValueError("heatmap_top_n is outside the supported range")
            object.__setattr__(self, "heatmap_top_n", top_n)
        object.__setattr__(self, "key_residues", tuple(self.key_residues))
        object.__setattr__(self, "selected_types", tuple(self.selected_types))


def build_project_input(source: str | os.PathLike[str]) -> ProjectInput:
    """Hash a regular source file and return immutable provenance."""

    path = Path(source).expanduser().resolve(strict=True)
    if not path.is_file():
        raise ValueError("project input must be a regular file")
    stat = path.stat()
    return ProjectInput(
        path=str(path),
        sha256=_sha256(path),
        size_bytes=stat.st_size,
        modified_ns=stat.st_mtime_ns,
    )


def validate_project_inputs(project: ProjectState) -> tuple[str, ...]:
    """Report missing or changed sources without modifying project state."""

    messages: list[str] = []
    datasets = (project.primary,) + (
        (project.comparison,) if project.comparison is not None else ()
    )
    for dataset in datasets:
        for source in dataset.inputs:
            path = Path(source.path)
            if not path.is_file():
                messages.append("Source is missing: %s" % source.path)
                continue
            try:
                stat = path.stat()
                changed = (
                    stat.st_size != source.size_bytes
                    or _sha256(path) != source.sha256
                )
            except OSError:
                messages.append("Source could not be read: %s" % source.path)
                continue
            if changed:
                messages.append("Source has changed: %s" % source.path)
    return tuple(messages)


def methods_summary(project: ProjectState) -> str:
    """Return the methods disclosure stored beside and inside a project."""

    profile = (
        "Discovery Studio-like"
        if project.analysis_profile == "ds_like"
        else project.analysis_profile.replace("_", " ").title()
    )
    comparison = (
        "\nComparison dataset: %s (%s)."
        % (project.comparison.label, project.comparison.mode)
        if project.comparison is not None
        else ""
    )
    block_size = (
        f"{project.bootstrap_block_size} saved frames"
        if project.bootstrap_block_size is not None
        else "automatic square-root heuristic"
    )
    time_axis = (
        f"{project.primary.time_step_ns:g} ns between saved frames"
        if project.primary.mode == "md"
        and project.primary.time_step_ns is not None
        else "not temporal"
    )
    axis_source = _dataset_axis_summary(project.primary)
    comparison_axis = (
        "\nSystem B trajectory mapping: "
        + _dataset_axis_summary(project.comparison)
        + "."
        if project.comparison is not None
        else ""
    )
    primary_scope = _ligand_scope_summary(
        project.primary, project.primary_ligand_group
    )
    comparison_scope = (
        "\nSystem B chart scope: "
        + _ligand_scope_summary(
            project.comparison, project.comparison_ligand_group
        )
        + "."
        if project.comparison is not None
        else ""
    )
    label_name = {
        "ligand": "Ligand name",
        "file": "Uploaded file name",
        "index": "Pose / frame index",
    }[project.observation_label_mode]
    heatmap_rows = (
        "ligand/uploaded-file groups"
        if project.heatmap_group_by == "source"
        else "individual observations"
    )
    heatmap_features = (
        "residue × interaction type"
        if project.heatmap_feature_level == "residue_type"
        else "residue with any interaction"
    )
    heatmap_limit = (
        "all features"
        if project.heatmap_top_n is None
        else f"top {project.heatmap_top_n} features"
    )
    return (
        "DockLens reproducible analysis\n"
        "Application version: {version}\n"
        "Analysis profile: {profile}\n"
        "Hydrogen-bond preset: {preset}\n"
        "Primary dataset: {label} ({mode}).{comparison}\n"
        "Primary observation axis: {time_axis}.\n"
        "Trajectory mapping: {axis_source}.{comparison_axis}\n"
        "Primary chart scope: {primary_scope}.{comparison_scope}\n"
        "Observation labels: {label_name}.\n"
        "Heatmap: {heatmap_rows}; {heatmap_features}; {heatmap_limit}.\n"
        "Interaction counting: one presence per observation, receptor residue, "
        "and interaction type.\n"
        "Fingerprint similarity: Jaccard/Tanimoto coefficient.\n"
        "Interaction-state clustering: complete-link Jaccard/Tanimoto, "
        "threshold {threshold:.3f}.\n"
        "Uncertainty: circular moving-block bootstrap with {iterations} "
        "iterations, block length {block_size}, seed {seed}, and "
        "{confidence:.1f}% confidence; temporally adjacent frames are not "
        "assumed independent.\n"
        "Key residues: {residues}\n"
        "Selected interaction types: {types}\n"
    ).format(
        version=project.app_version,
        profile=profile,
        preset=project.hbond_preset,
        label=project.primary.label,
        mode=project.primary.mode,
        comparison=comparison,
        time_axis=time_axis,
        axis_source=axis_source,
        comparison_axis=comparison_axis,
        primary_scope=primary_scope,
        comparison_scope=comparison_scope,
        label_name=label_name,
        heatmap_rows=heatmap_rows,
        heatmap_features=heatmap_features,
        heatmap_limit=heatmap_limit,
        threshold=project.state_threshold,
        iterations=project.bootstrap_iterations,
        block_size=block_size,
        seed=project.bootstrap_seed,
        confidence=100.0 * project.confidence_level,
        residues=", ".join(project.key_residues) or "none",
        types=", ".join(project.selected_types) or "all",
    )


def _dataset_axis_summary(dataset: ProjectDataset) -> str:
    if dataset.mode != "md":
        return "not temporal"
    if dataset.observation_series is None:
        return "implicit single contiguous series"
    replicas = len(
        {point.replica_id for point in dataset.observation_series.points}
    )
    return (
        f"explicit trajectory map with {replicas} "
        f"{'replica' if replicas == 1 else 'replicas'}"
    )


def _ligand_scope_summary(
    dataset: ProjectDataset,
    group_key: str | None,
) -> str:
    if group_key is None:
        return "all ligand/uploaded-file groups, weighted by observation count"
    if dataset.result is None:
        return f"saved group {group_key}"
    group = next(
        (
            item
            for item in ligand_groups(dataset.result)
            if item.key == group_key
        ),
        None,
    )
    return group.label if group is not None else f"unavailable group {group_key}"


def save_project(
    project: ProjectState, destination: str | os.PathLike[str]
) -> tuple[Path, Path]:
    """Atomically write a safe project container and methods sidecar."""

    if not isinstance(project, ProjectState):
        raise TypeError("project must be a ProjectState")
    destination_path = Path(destination).expanduser().resolve()
    destination_path.parent.mkdir(parents=True, exist_ok=True)
    methods_path = destination_path.with_name(
        "%s_methods.txt" % destination_path.stem
    )
    methods = methods_summary(project)
    result_entries: dict[str, bytes] = {}
    project_payload = _project_to_payload(project, result_entries)
    manifest = {
        "schema_version": PROJECT_SCHEMA_VERSION,
        "methods_sha256": _sha256_bytes(methods.encode("utf-8")),
        "project": project_payload,
    }
    encoded = json.dumps(
        manifest,
        ensure_ascii=False,
        indent=2,
        sort_keys=True,
        allow_nan=False,
    ).encode("utf-8")
    if len(encoded) > MAX_UNCOMPRESSED_BYTES:
        raise ValueError("project is too large")

    entries = {
        "manifest.json": encoded,
        "methods.txt": methods.encode("utf-8"),
        **result_entries,
    }
    _atomic_zip_write(destination_path, entries)
    _atomic_write(methods_path, methods.encode("utf-8"))
    return destination_path, methods_path


def load_project(source: str | os.PathLike[str]) -> ProjectState:
    """Load and validate a project without extracting archive entries."""

    path = Path(source).expanduser().resolve(strict=True)
    size = path.stat().st_size
    if size > MAX_PROJECT_BYTES:
        raise ValueError("project file is too large")
    if not zipfile.is_zipfile(path):
        # Read-only legacy JSON diagnostics keep older projects intelligible.
        # New projects are always written as ZIP containers.
        return _load_legacy_json(path)

    try:
        with zipfile.ZipFile(path, "r") as archive:
            infos = archive.infolist()
            _validate_archive_infos(infos)
            manifest_bytes = _read_bounded(
                archive, "manifest.json", MAX_UNCOMPRESSED_BYTES
            )
            document = _parse_manifest(manifest_bytes)
            schema_version = str(document.get("schema_version", ""))
            if schema_version not in SUPPORTED_PROJECT_SCHEMA_VERSIONS:
                raise ValueError("unsupported project schema")
            names = {info.filename for info in infos}
            if "methods.txt" not in names:
                raise ValueError("project is missing required entries")
            methods_bytes = _read_bounded(archive, "methods.txt", 1024 * 1024)
            methods_hash = document.get("methods_sha256")
            if not isinstance(methods_hash, str) or not _valid_sha(methods_hash):
                raise ValueError("project methods hash is invalid")
            if _sha256_bytes(methods_bytes) != methods_hash.lower():
                raise ValueError("project methods integrity hash does not match")
            referenced: set[str] = set()

            def load_result(reference: Mapping[str, Any]) -> RunResult:
                entry = _text(reference, "entry")
                expected_hash = _sha_text(reference, "sha256")
                expected_size = _integer(reference, "size_bytes")
                if not _is_result_entry(entry) or entry not in names:
                    raise ProjectIntegrityError(
                        "cached result entry is missing or unsafe"
                    )
                encoded = _read_bounded(archive, entry, MAX_UNCOMPRESSED_BYTES)
                if len(encoded) != expected_size:
                    raise ProjectIntegrityError(
                        "cached result integrity size does not match"
                    )
                if _sha256_bytes(encoded) != expected_hash:
                    raise ProjectIntegrityError(
                        "cached result integrity hash does not match"
                    )
                referenced.add(entry)
                result_document = _parse_json_object(
                    encoded, "cached analysis result"
                )
                return _result_from_payload(result_document)

            project = _decode_document(document, load_result)
            unreferenced = {
                name for name in names if _is_result_entry(name)
            } - referenced
            if unreferenced:
                raise ValueError("project contains an unreferenced result entry")
    except (OSError, zipfile.BadZipFile, RuntimeError) as exc:
        raise ValueError("invalid DockLens project container") from exc

    expected = methods_summary(project).encode("utf-8")
    if schema_version == PROJECT_SCHEMA_VERSION and methods_bytes != expected:
        raise ValueError("project methods record does not match its manifest")
    return project


def _load_legacy_json(path: Path) -> ProjectState:
    try:
        document = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise ValueError("invalid project document") from exc
    if not isinstance(document, dict):
        raise ValueError("project document must be an object")
    schema_version = str(document.get("schema_version", ""))
    if schema_version not in SUPPORTED_PROJECT_SCHEMA_VERSIONS:
        raise ValueError("unsupported project schema")
    return _decode_document(document)


def _parse_manifest(encoded: bytes) -> Mapping[str, Any]:
    return _parse_json_object(encoded, "project manifest")


def _parse_json_object(encoded: bytes, description: str) -> Mapping[str, Any]:
    try:
        document = json.loads(encoded.decode("utf-8"))
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise ValueError("invalid %s" % description) from exc
    if not isinstance(document, dict):
        raise ValueError("%s must be an object" % description)
    _check_json_depth(document)
    return document


def _decode_document(
    document: Mapping[str, Any],
    result_loader=None,
) -> ProjectState:
    if (
        str(document.get("schema_version", ""))
        not in SUPPORTED_PROJECT_SCHEMA_VERSIONS
    ):
        raise ValueError("unsupported project schema")
    payload = document.get("project")
    if not isinstance(payload, dict):
        raise ValueError("project manifest is missing project state")
    try:
        primary = _dataset_from_payload(
            _mapping(payload, "primary"), result_loader
        )
        comparison_payload = payload.get("comparison")
        comparison = (
            _dataset_from_payload(comparison_payload, result_loader)
            if isinstance(comparison_payload, dict)
            else None
        )
        return ProjectState(
            app_version=_text(payload, "app_version"),
            analysis_profile=_text(payload, "analysis_profile"),
            hbond_preset=_text(payload, "hbond_preset"),
            key_residues=_text_tuple(payload, "key_residues"),
            selected_types=_text_tuple(payload, "selected_types"),
            active_workspace=_text(payload, "active_workspace"),
            selected_residue=_text(payload, "selected_residue"),
            state_threshold=_number(payload, "state_threshold"),
            bootstrap_iterations=_integer(payload, "bootstrap_iterations"),
            primary=primary,
            comparison=comparison,
            bootstrap_block_size=(
                None
                if payload.get("bootstrap_block_size") is None
                else _strict_int(payload["bootstrap_block_size"])
            ),
            bootstrap_seed=_strict_int(payload.get("bootstrap_seed", 2026)),
            confidence_level=_number_value(
                payload.get("confidence_level", 0.95),
                "confidence_level",
            ),
            primary_ligand_group=(
                None
                if payload.get("primary_ligand_group") is None
                else _text(payload, "primary_ligand_group")
            ),
            comparison_ligand_group=(
                None
                if payload.get("comparison_ligand_group") is None
                else _text(payload, "comparison_ligand_group")
            ),
            observation_label_mode=str(
                payload.get("observation_label_mode", "ligand")
            ),
            heatmap_group_by=str(
                payload.get("heatmap_group_by", "source")
            ),
            heatmap_feature_level=str(
                payload.get("heatmap_feature_level", "residue_type")
            ),
            heatmap_top_n=(
                None
                if payload.get("heatmap_top_n", 40) is None
                else _strict_int(payload.get("heatmap_top_n", 40))
            ),
        )
    except ProjectIntegrityError:
        raise
    except (KeyError, TypeError, ValueError) as exc:
        raise ValueError("invalid project state") from exc


def _project_to_payload(
    project: ProjectState, result_entries: dict[str, bytes]
) -> dict[str, Any]:
    used_names: set[str] = set()
    primary_entry = _result_entry_name(project.primary.label, used_names)
    comparison_entry = (
        _result_entry_name(project.comparison.label, used_names)
        if project.comparison is not None
        else None
    )
    return {
        "app_version": project.app_version,
        "analysis_profile": project.analysis_profile,
        "hbond_preset": project.hbond_preset,
        "key_residues": list(project.key_residues),
        "selected_types": list(project.selected_types),
        "active_workspace": project.active_workspace,
        "selected_residue": project.selected_residue,
        "state_threshold": project.state_threshold,
        "bootstrap_iterations": project.bootstrap_iterations,
        "primary": _dataset_to_payload(
            project.primary, primary_entry, result_entries
        ),
        "comparison": (
            _dataset_to_payload(
                project.comparison, comparison_entry, result_entries
            )
            if project.comparison is not None
            else None
        ),
        "bootstrap_block_size": project.bootstrap_block_size,
        "bootstrap_seed": project.bootstrap_seed,
        "confidence_level": project.confidence_level,
        "primary_ligand_group": project.primary_ligand_group,
        "comparison_ligand_group": project.comparison_ligand_group,
        "observation_label_mode": project.observation_label_mode,
        "heatmap_group_by": project.heatmap_group_by,
        "heatmap_feature_level": project.heatmap_feature_level,
        "heatmap_top_n": project.heatmap_top_n,
    }


def _dataset_to_payload(
    dataset: ProjectDataset,
    result_entry: str | None,
    result_entries: dict[str, bytes],
) -> dict[str, Any]:
    result_reference = None
    if dataset.result is not None:
        if result_entry is None:
            raise ValueError("cached result entry name is missing")
        encoded = json.dumps(
            _result_to_payload(dataset.result),
            ensure_ascii=False,
            indent=2,
            sort_keys=True,
            allow_nan=False,
        ).encode("utf-8")
        if len(encoded) > MAX_UNCOMPRESSED_BYTES:
            raise ValueError("cached analysis result is too large")
        result_entries[result_entry] = encoded
        result_reference = {
            "entry": result_entry,
            "sha256": _sha256_bytes(encoded),
            "size_bytes": len(encoded),
        }
    return {
        "label": dataset.label,
        "mode": dataset.mode,
        "time_step_ns": dataset.time_step_ns,
        "inputs": [
            {
                "path": item.path,
                "sha256": item.sha256,
                "size_bytes": item.size_bytes,
                "modified_ns": item.modified_ns,
            }
            for item in dataset.inputs
        ],
        "result": result_reference,
        "observation_series": (
            {
                "mode": dataset.observation_series.mode,
                "time_step_ns": dataset.observation_series.time_step_ns,
                "points": [
                    {
                        "observation_id": point.observation_id,
                        "ordinal": point.ordinal,
                        "frame_index": point.frame_index,
                        "time_ns": point.time_ns,
                        "replica_id": point.replica_id,
                    }
                    for point in dataset.observation_series.points
                ],
            }
            if dataset.observation_series is not None
            else None
        ),
    }


def _dataset_from_payload(
    payload: Mapping[str, Any], result_loader=None
) -> ProjectDataset:
    inputs_payload = payload.get("inputs", [])
    if not isinstance(inputs_payload, list):
        raise ValueError("dataset inputs must be a list")
    inputs = tuple(
        ProjectInput(
            path=_text(item, "path"),
            sha256=_sha_text(item, "sha256"),
            size_bytes=_integer(item, "size_bytes"),
            modified_ns=_integer(item, "modified_ns"),
        )
        for item in inputs_payload
        if isinstance(item, dict)
    )
    if len(inputs) != len(inputs_payload):
        raise ValueError("invalid dataset input")
    result_payload = payload.get("result")
    if result_payload is not None and not isinstance(result_payload, dict):
        raise ValueError("dataset result reference must be an object")
    if result_payload is not None and result_loader is None:
        raise ValueError("cached result cannot be loaded from this document")
    series_payload = payload.get("observation_series")
    if series_payload is not None and not isinstance(series_payload, dict):
        raise ValueError("observation_series must be an object")
    observation_series = None
    if isinstance(series_payload, dict):
        point_payloads = _list(series_payload, "points")
        observation_series = ObservationSeries(
            mode=_text(series_payload, "mode"),
            time_step_ns=(
                None
                if series_payload.get("time_step_ns") is None
                else _number(series_payload, "time_step_ns")
            ),
            points=tuple(
                ObservationPoint(
                    observation_id=_text(item, "observation_id"),
                    ordinal=_integer(item, "ordinal"),
                    frame_index=(
                        None
                        if item.get("frame_index") is None
                        else _integer(item, "frame_index")
                    ),
                    time_ns=(
                        None
                        if item.get("time_ns") is None
                        else _number(item, "time_ns")
                    ),
                    replica_id=_text(item, "replica_id"),
                )
                for item in point_payloads
                if isinstance(item, dict)
            ),
        )
        if len(observation_series.points) != len(point_payloads):
            raise ValueError("invalid observation-series point")
    return ProjectDataset(
        label=_text(payload, "label"),
        mode=_text(payload, "mode"),
        time_step_ns=(
            None
            if payload.get("time_step_ns") is None
            else _number(payload, "time_step_ns")
        ),
        inputs=inputs,
        result=(
            result_loader(result_payload)
            if isinstance(result_payload, dict)
            else None
        ),
        observation_series=observation_series,
    )


def _result_to_payload(result: RunResult) -> dict[str, Any]:
    if result.pending:
        raise ValueError(
            "unresolved confirmations must be completed before saving a cached result"
        )
    return {
        "details": [_dataclass_payload(item) for item in result.details],
        "summaries": [_dataclass_payload(item) for item in result.summaries],
        "pending": [],
        "key_residues": sorted(result.key_residues),
        "receptor_residues": sorted(result.receptor_residues),
        "input_qc": [_dataclass_payload(item) for item in result.input_qc],
        "parameters": _dataclass_payload(result.parameters),
    }


def _result_from_payload(payload: Mapping[str, Any]) -> RunResult:
    try:
        details = tuple(_detail_from_payload(item) for item in _list(payload, "details"))
        summaries = tuple(
            _summary_from_payload(item) for item in _list(payload, "summaries")
        )
        if _list(payload, "pending"):
            raise ValueError("cached result contains unresolved confirmations")
        qc = tuple(_construct(InputQC, item) for item in _list(payload, "input_qc"))
        parameters = _construct(
            AnalysisParameters, _mapping(payload, "parameters")
        )
        return RunResult(
            details=details,
            summaries=summaries,
            pending=(),
            key_residues=frozenset(_text_tuple(payload, "key_residues")),
            receptor_residues=frozenset(
                _text_tuple(payload, "receptor_residues")
            ),
            input_qc=qc,
            parameters=parameters,
        )
    except (KeyError, TypeError, ValueError) as exc:
        raise ValueError("invalid cached analysis result") from exc


def _detail_from_payload(payload: Mapping[str, Any]) -> Detail:
    values = dict(payload)
    values["ligand"] = _construct(Endpoint, _mapping(payload, "ligand"))
    values["receptor"] = _construct(Endpoint, _mapping(payload, "receptor"))
    water = payload.get("water")
    values["water"] = _construct(Endpoint, water) if isinstance(water, dict) else None
    return _construct(Detail, values)


def _summary_from_payload(payload: Mapping[str, Any]) -> Summary:
    values = dict(payload)
    counts = values.get("counts")
    if not isinstance(counts, dict):
        raise ValueError("summary counts must be an object")
    values["counts"] = {
        str(key): _strict_int(value) for key, value in counts.items()
    }
    return _construct(Summary, values)


def _dataclass_payload(instance: Any) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for item in fields(instance):
        value = getattr(instance, item.name)
        if hasattr(value, "__dataclass_fields__"):
            result[item.name] = _dataclass_payload(value)
        elif isinstance(value, Mapping):
            result[item.name] = dict(value)
        elif isinstance(value, (tuple, list, frozenset, set)):
            result[item.name] = list(value)
        else:
            result[item.name] = value
    return result


def _construct(cls, payload: Mapping[str, Any]):
    if not isinstance(payload, Mapping):
        raise ValueError("record must be an object")
    allowed = {item.name for item in fields(cls)}
    if not set(payload).issubset(allowed):
        raise ValueError("record contains unknown fields")
    values = dict(payload)
    for name in ("atom_serials", "warnings", "cutoffs", "interaction_types", "key_residues"):
        if name in values:
            values[name] = tuple(values[name])
    if "cutoffs" in values:
        values["cutoffs"] = tuple(
            (str(name), float(value)) for name, value in values["cutoffs"]
        )
    return cls(**values)


def _result_entry_name(label: str, used: set[str]) -> str:
    slug = re.sub(r"[^a-z0-9]+", "-", label.strip().lower()).strip("-")
    if not slug:
        slug = "dataset"
    slug = slug[:64].rstrip("-")
    candidate = slug
    suffix = 2
    while candidate in used:
        candidate = "%s-%d" % (slug[:58].rstrip("-"), suffix)
        suffix += 1
    used.add(candidate)
    return "datasets/%s/run_result.json" % candidate


def _is_result_entry(name: str) -> bool:
    parts = PurePosixPath(name).parts
    return (
        len(parts) == 3
        and parts[0] == "datasets"
        and bool(parts[1])
        and parts[2] == "run_result.json"
    )


def _validate_archive_infos(infos: list[zipfile.ZipInfo]) -> None:
    if len(infos) > MAX_ARCHIVE_ENTRIES:
        raise ValueError("project contains too many entries")
    names: set[str] = set()
    total = 0
    for info in infos:
        name = info.filename
        pure = PurePosixPath(name)
        if (
            not name
            or "\\" in name
            or pure.is_absolute()
            or ".." in pure.parts
            or ":" in pure.parts[0]
        ):
            raise ValueError("project contains an unsafe entry name")
        if name in names:
            raise ValueError("project contains duplicate entries")
        if name not in _FIXED_ENTRIES and not _is_result_entry(name):
            raise ValueError("project contains an unknown entry")
        if info.is_dir() or ((info.external_attr >> 16) & 0o170000) == 0o120000:
            raise ValueError("project contains an unsupported entry")
        if info.file_size < 0 or info.file_size > MAX_UNCOMPRESSED_BYTES:
            raise ValueError("project entry is too large")
        total += info.file_size
        if total > MAX_UNCOMPRESSED_BYTES:
            raise ValueError("project contents are too large")
        names.add(name)
    if "manifest.json" not in names:
        raise ValueError("project is missing required entries")


def _read_bounded(
    archive: zipfile.ZipFile, name: str, maximum: int
) -> bytes:
    with archive.open(name, "r") as handle:
        data = handle.read(maximum + 1)
    if len(data) > maximum:
        raise ValueError("project entry is too large")
    return data


def _atomic_zip_write(destination: Path, entries: Mapping[str, bytes]) -> None:
    temporary: Path | None = None
    try:
        with tempfile.NamedTemporaryFile(
            dir=destination.parent, delete=False
        ) as handle:
            temporary = Path(handle.name)
        with zipfile.ZipFile(
            temporary,
            "w",
            compression=zipfile.ZIP_DEFLATED,
            compresslevel=9,
        ) as archive:
            for name, payload in entries.items():
                archive.writestr(name, payload)
        if temporary.stat().st_size > MAX_PROJECT_BYTES:
            raise ValueError("project file is too large")
        with temporary.open("r+b") as handle:
            os.fsync(handle.fileno())
        os.replace(temporary, destination)
        temporary = None
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)


def _atomic_write(destination: Path, encoded: bytes) -> None:
    temporary: Path | None = None
    try:
        with tempfile.NamedTemporaryFile(
            dir=destination.parent, delete=False
        ) as handle:
            temporary = Path(handle.name)
            handle.write(encoded)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, destination)
        temporary = None
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def _valid_sha(value: str) -> bool:
    normalized = value.lower()
    return len(normalized) == 64 and all(
        character in "0123456789abcdef" for character in normalized
    )


def _check_json_depth(value: Any, depth: int = 0) -> None:
    if depth > MAX_JSON_DEPTH:
        raise ValueError("project manifest is too deeply nested")
    if isinstance(value, dict):
        for key, item in value.items():
            if not isinstance(key, str):
                raise ValueError("project manifest keys must be text")
            _check_json_depth(item, depth + 1)
    elif isinstance(value, list):
        for item in value:
            _check_json_depth(item, depth + 1)


def _mapping(payload: Mapping[str, Any], name: str) -> Mapping[str, Any]:
    value = payload[name]
    if not isinstance(value, dict):
        raise ValueError("%s must be an object" % name)
    return value


def _list(payload: Mapping[str, Any], name: str) -> list[Any]:
    value = payload[name]
    if not isinstance(value, list):
        raise ValueError("%s must be a list" % name)
    return value


def _text(payload: Mapping[str, Any], name: str) -> str:
    value = payload[name]
    if not isinstance(value, str):
        raise ValueError("%s must be text" % name)
    return value


def _sha_text(payload: Mapping[str, Any], name: str) -> str:
    value = _text(payload, name).lower()
    if not _valid_sha(value):
        raise ValueError("%s must be a SHA-256 digest" % name)
    return value


def _text_tuple(payload: Mapping[str, Any], name: str) -> tuple[str, ...]:
    value = payload[name]
    if not isinstance(value, list) or any(not isinstance(item, str) for item in value):
        raise ValueError("%s must be a text list" % name)
    return tuple(value)


def _number(payload: Mapping[str, Any], name: str) -> float:
    value = payload[name]
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError("%s must be numeric" % name)
    number = float(value)
    if not math.isfinite(number):
        raise ValueError("%s must be finite" % name)
    return number


def _number_value(value: Any, name: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError("%s must be numeric" % name)
    number = float(value)
    if not math.isfinite(number):
        raise ValueError("%s must be finite" % name)
    return number


def _strict_int(value: Any) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        raise ValueError("value must be an integer")
    return value


def _integer(payload: Mapping[str, Any], name: str) -> int:
    return _strict_int(payload[name])
