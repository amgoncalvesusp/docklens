"""Atomic DockLens result contract consumed by DockingHub."""

from __future__ import annotations

from collections import Counter
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import tempfile

from . import __version__
from .integration_manifest import LaunchManifest, ManifestError

RESULT_SCHEMA = "docklens-result-v1"
MAX_RESULT_BYTES = 128 * 1024 * 1024


def write_integration_result(manifest: LaunchManifest, result) -> Path | None:
    """Serialize a paired analysis without mutating the immutable result."""
    destination = manifest.result_path
    if destination is None:
        return None
    type_counts = Counter(detail.interaction_type for detail in result.details)
    payload = {
        "schema": RESULT_SCHEMA,
        "run_id": manifest.run_id,
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "generated_by": {"application": "DockLens", "version": __version__},
        "launch_manifest_sha256": _sha256(manifest.source_path),
        "inputs": {
            "receptor_sha256": manifest.receptor_sha256,
            "poses_sha256": manifest.poses_sha256,
        },
        "parameters": {
            "hbond_preset": result.parameters.hbond_preset,
            "key_residues": list(result.parameters.key_residues),
            "interaction_types": list(result.parameters.interaction_types),
            "counting_unit": result.parameters.counting_unit,
        },
        "totals": {
            "poses": len(result.summaries),
            "interactions": len(result.details),
            "errors": sum(record.status == "error" for record in result.input_qc),
            "by_type": dict(sorted(type_counts.items())),
        },
        "poses": [_summary_payload(item) for item in result.summaries],
        "interactions": [_detail_payload(item) for item in result.details],
        "input_qc": [_qc_payload(item) for item in result.input_qc],
    }
    encoded = json.dumps(payload, ensure_ascii=False, indent=2, sort_keys=True).encode("utf-8")
    if len(encoded) > MAX_RESULT_BYTES:
        raise ManifestError("result_too_large", "The DockLens result exceeds the size limit.")
    _atomic_write(destination, encoded)
    return destination


def _summary_payload(item) -> dict[str, object]:
    return {
        "pose_id": item.pose_id,
        "pose": item.pose,
        "ligand_id": item.ligand_id,
        "docking_score": item.docking_score,
        "total_interactions": item.n_total_interactions,
        "key_residue_interactions": item.n_key_residue_interactions,
        "counts": dict(sorted(item.counts.items())),
        "resolution_method": item.resolution_method,
    }


def _detail_payload(item) -> dict[str, object]:
    return {
        "interaction_id": item.interaction_id,
        "pose_id": item.pose_id,
        "pose": item.pose,
        "type": item.interaction_type,
        "subtype": item.subtype,
        "receptor_residue": item.receptor_residue,
        "receptor_atom": item.receptor_atom,
        "ligand_atom": item.ligand_atom,
        "distance_A": item.distance_A,
        "is_key_residue": item.is_key_residue,
        "chemistry_basis": item.chemistry_basis,
        "chemistry_confidence": item.chemistry_confidence,
        "hydrogen_atom": item.hydrogen_atom,
        "hydrogen_atom_serial": item.hydrogen_atom_serial,
        "hydrogen_acceptor_distance_A": item.hydrogen_acceptor_distance_A,
        "donor_hydrogen_acceptor_angle_deg": (
            item.donor_hydrogen_acceptor_angle_deg
        ),
        "hydrogen_acceptor_base_angle_deg": (
            item.hydrogen_acceptor_base_angle_deg
        ),
        "theta_deg": item.theta_deg,
    }


def _qc_payload(item) -> dict[str, object]:
    return {
        "pose_id": item.pose_id,
        "status": item.status,
        "code": item.code,
        "message": item.message,
        "format": item.format,
        "resolution_method": item.resolution_method,
        "receptor_atoms": item.receptor_atoms,
        "ligand_atoms": item.ligand_atoms,
        "water_atoms": item.water_atoms,
    }


def _atomic_write(destination: Path, encoded: bytes) -> None:
    destination = destination.resolve()
    with tempfile.NamedTemporaryFile(dir=destination.parent, delete=False) as handle:
        temporary = Path(handle.name)
        handle.write(encoded)
        handle.flush()
        os.fsync(handle.fileno())
    try:
        os.replace(temporary, destination)
    finally:
        temporary.unlink(missing_ok=True)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()
