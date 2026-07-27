"""Validated launch contract for pairing one receptor with docking poses."""

from __future__ import annotations

from dataclasses import dataclass
import hashlib
import hmac
import json
from pathlib import Path
import re

from .interaction_core import HBOND_PRESETS


SCHEMA = "docklens-launch-v1"
MAX_MANIFEST_BYTES = 1024 * 1024
MAX_STRUCTURE_BYTES = 256 * 1024 * 1024
SUPPORTED_SUFFIXES = {".pdb", ".pdbqt", ".mol2"}
_HASH = re.compile(r"^[0-9a-f]{64}$")


class ManifestError(ValueError):
    def __init__(self, code: str, message: str) -> None:
        super().__init__(message)
        self.code = code


@dataclass(frozen=True)
class LaunchManifest:
    source_path: Path
    project_root: Path
    run_id: str
    receptor_path: Path
    receptor_sha256: str
    poses_path: Path
    poses_sha256: str
    hbond_preset: str
    key_residues: tuple[str, ...]
    result_path: Path | None = None


def load_manifest(path: Path | str) -> LaunchManifest:
    source = Path(path).expanduser().resolve()
    if not source.is_file():
        raise ManifestError("manifest_not_found", "Launch manifest was not found.")
    if source.stat().st_size > MAX_MANIFEST_BYTES:
        raise ManifestError("manifest_too_large", "Launch manifest exceeds the size limit.")
    try:
        payload = json.loads(source.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        raise ManifestError("manifest_invalid", "Launch manifest is not valid UTF-8 JSON.") from exc
    if not isinstance(payload, dict) or payload.get("schema") != SCHEMA:
        raise ManifestError("schema_unsupported", "Launch manifest schema is unsupported.")
    project_value = payload.get("project_root")
    if not isinstance(project_value, str) or not project_value:
        raise ManifestError("project_root_required", "A project root is required.")
    project_root = Path(project_value).expanduser().resolve()
    if not project_root.is_dir():
        raise ManifestError("project_root_invalid", "The project root is not available.")
    _require_within(source, project_root, "manifest_outside_project")
    receptor_path, receptor_hash = _validated_input(payload.get("receptor"), project_root)
    poses_path, poses_hash = _validated_input(payload.get("poses"), project_root)
    options = payload.get("options", {})
    if not isinstance(options, dict):
        raise ManifestError("options_invalid", "Manifest options must be an object.")
    preset = options.get("hbond_preset", "plip")
    if preset not in HBOND_PRESETS:
        raise ManifestError("preset_invalid", "The H-bond preset is unsupported.")
    key_values = options.get("key_residues", [])
    if not isinstance(key_values, list) or not all(isinstance(item, str) for item in key_values):
        raise ManifestError("key_residues_invalid", "Key residues must be a list of strings.")
    if len(key_values) > 10_000:
        raise ManifestError("key_residues_too_many", "Too many key residues were supplied.")
    run_id = payload.get("run_id", "")
    if not isinstance(run_id, str) or not run_id.strip() or len(run_id) > 200:
        raise ManifestError("run_id_invalid", "Run ID is invalid.")
    result_path = _validated_output(payload.get("result"), project_root, source)
    return LaunchManifest(
        source,
        project_root,
        run_id,
        receptor_path,
        receptor_hash,
        poses_path,
        poses_hash,
        preset,
        tuple(item.strip().upper() for item in key_values if item.strip()),
        result_path,
    )


def _validated_input(value, project_root: Path) -> tuple[Path, str]:
    if not isinstance(value, dict):
        raise ManifestError("input_required", "Receptor and poses inputs are required.")
    raw_path = value.get("path")
    expected_hash = value.get("sha256")
    if not isinstance(raw_path, str) or not raw_path:
        raise ManifestError("input_path_required", "Input path is required.")
    if not isinstance(expected_hash, str) or not _HASH.fullmatch(expected_hash.lower()):
        raise ManifestError("hash_invalid", "Input SHA-256 is invalid.")
    candidate = Path(raw_path)
    resolved = (project_root / candidate).resolve() if not candidate.is_absolute() else candidate.resolve()
    _require_within(resolved, project_root, "path_outside_project")
    if not resolved.is_file():
        raise ManifestError("input_not_found", "A manifest input was not found.")
    if resolved.suffix.lower() not in SUPPORTED_SUFFIXES:
        raise ManifestError("input_format_unsupported", "A manifest input format is unsupported.")
    if resolved.stat().st_size > MAX_STRUCTURE_BYTES:
        raise ManifestError("input_too_large", "A manifest input exceeds the size limit.")
    actual_hash = _sha256(resolved)
    if not hmac.compare_digest(actual_hash, expected_hash.lower()):
        raise ManifestError("hash_mismatch", "A manifest input no longer matches its SHA-256.")
    return resolved, actual_hash


def _validated_output(value, project_root: Path, source_path: Path) -> Path | None:
    if value is None:
        return None
    if not isinstance(value, dict) or value.get("schema") != "docklens-result-v1":
        raise ManifestError("result_contract_invalid", "The result contract is invalid.")
    raw_path = value.get("path")
    if not isinstance(raw_path, str) or not raw_path:
        raise ManifestError("result_path_required", "A result output path is required.")
    candidate = Path(raw_path)
    resolved = (project_root / candidate).resolve() if not candidate.is_absolute() else candidate.resolve()
    reports_root = (project_root / "reports").resolve()
    if resolved == source_path:
        raise ManifestError("result_path_collision", "The result output collides with the launch manifest.")
    if not reports_root.is_dir():
        raise ManifestError("result_parent_invalid", "The reports directory is unavailable.")
    _require_within(resolved, reports_root, "result_outside_reports")
    if (
        resolved.suffix.lower() != ".json"
        or not resolved.name.startswith("docklens_result_")
        or resolved.is_dir()
    ):
        raise ManifestError("result_path_invalid", "The result output must be a JSON file.")
    if not resolved.parent.is_dir():
        raise ManifestError("result_parent_invalid", "The result output directory is unavailable.")
    return resolved


def _require_within(path: Path, root: Path, code: str) -> None:
    if path != root and root not in path.parents:
        raise ManifestError(code, "A manifest path is outside the project root.")


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()
