"""Source-native molecule and pose identity parsing.

The parser deliberately recognizes only a small set of unambiguous pose-label
patterns.  A pipe-delimited name is not treated as GOLD metadata merely because
it contains a pipe.
"""

from __future__ import annotations

from dataclasses import dataclass
import re


_SOURCE_POSE_RE = re.compile(
    r"(?:dock|pose|sol|solution)[_-]?\d+\Z", re.IGNORECASE
)


@dataclass(frozen=True)
class SourcePoseIdentity:
    molecule_name: str = ""
    ligand_id: str = ""
    source_pose_label: str = ""


def parse_source_pose_identity(molecule_name: str | None) -> SourcePoseIdentity:
    """Return identity fields without rewriting the source molecule name."""

    original = "" if molecule_name is None else str(molecule_name)
    fields = tuple(part.strip() for part in original.split("|") if part.strip())
    if not fields:
        return SourcePoseIdentity(molecule_name=original)

    pose_label = next(
        (field for field in reversed(fields) if _SOURCE_POSE_RE.fullmatch(field)),
        "",
    )
    return SourcePoseIdentity(
        molecule_name=original,
        ligand_id=fields[0],
        source_pose_label=pose_label,
    )


__all__ = ["SourcePoseIdentity", "parse_source_pose_identity"]
