"""Immutable chart scopes for ligands grouped by uploaded source file."""

from __future__ import annotations

from dataclasses import dataclass, replace
from typing import Iterable

from .observation_series import ObservationPoint, ObservationSeries
from .results import RunResult


@dataclass(frozen=True)
class LigandGroup:
    """One uploaded source and all analytical observations resolved from it."""

    key: str
    source_id: str
    source_file: str
    source_path: str
    ligand_ids: tuple[str, ...]
    observation_ids: tuple[str, ...]

    @property
    def observation_count(self) -> int:
        return len(self.observation_ids)

    @property
    def label(self) -> str:
        if len(self.ligand_ids) == 1:
            ligands = self.ligand_ids[0]
        elif self.ligand_ids:
            ligands = f"{len(self.ligand_ids)} detected ligand labels"
        else:
            ligands = "unlabelled ligand"
        source = self.source_file or self.source_id or "uploaded source"
        identity = f" [{self.source_id}]" if self.source_id else ""
        return (
            f"{source} · {ligands} "
            f"({self.observation_count} observations){identity}"
        )


def _source_key(item) -> str:
    return str(
        getattr(item, "source_id", "")
        or getattr(item, "source_path", "")
        or getattr(item, "source_file", "")
    )


def ligand_groups(result: RunResult) -> tuple[LigandGroup, ...]:
    """Catalog uploaded ligand sources without splitting multipose files."""

    grouped: dict[str, dict[str, object]] = {}
    for summary in result.summaries:
        key = _source_key(summary)
        if not key:
            key = "unidentified-source"
        entry = grouped.setdefault(
            key,
            {
                "source_id": summary.source_id,
                "source_file": summary.source_file,
                "source_path": summary.source_path,
                "ligand_ids": [],
                "observation_ids": [],
            },
        )
        ligand_ids = entry["ligand_ids"]
        observation_ids = entry["observation_ids"]
        if summary.ligand_id not in ligand_ids:
            ligand_ids.append(summary.ligand_id)
        if summary.pose_id not in observation_ids:
            observation_ids.append(summary.pose_id)
    return tuple(
        LigandGroup(
            key=key,
            source_id=str(entry["source_id"]),
            source_file=str(entry["source_file"]),
            source_path=str(entry["source_path"]),
            ligand_ids=tuple(entry["ligand_ids"]),
            observation_ids=tuple(entry["observation_ids"]),
        )
        for key, entry in grouped.items()
    )


def subset_run_result(
    result: RunResult,
    group_key: str | None,
) -> RunResult:
    """Return a source-scoped result while retaining empty observations."""

    if group_key is None:
        return result
    groups = {group.key: group for group in ligand_groups(result)}
    if group_key not in groups:
        raise ValueError("unknown ligand/uploaded-file group")
    allowed_ids = frozenset(groups[group_key].observation_ids)
    summaries = tuple(
        item for item in result.summaries if item.pose_id in allowed_ids
    )
    details = tuple(
        item for item in result.details if item.pose_id in allowed_ids
    )
    pending = tuple(
        item
        for item in result.pending
        if getattr(item, "pose_id", "") in allowed_ids
    )
    qc = tuple(
        item for item in result.input_qc if _source_key(item) == group_key
    )
    return replace(
        result,
        summaries=summaries,
        details=details,
        pending=pending,
        input_qc=qc,
    )


def subset_observation_series(
    series: ObservationSeries,
    observation_ids: Iterable[str],
) -> ObservationSeries:
    """Restrict an explicit series without renumbering frames or time."""

    allowed = frozenset(str(value) for value in observation_ids)
    known = frozenset(series.observation_ids)
    if not allowed.issubset(known):
        raise ValueError("observation IDs are not present in the series")
    return ObservationSeries(
        mode=series.mode,
        points=tuple(
            point for point in series.points if point.observation_id in allowed
        ),
        time_step_ns=series.time_step_ns,
    )


def default_md_series_for_result(
    result: RunResult,
    *,
    time_step_ns: float,
) -> ObservationSeries:
    """Build a disclosed MD fallback with one replica per uploaded source."""

    source_by_observation = {
        item.pose_id: (_source_key(item) or "dataset-1")
        for item in result.summaries
    }
    ordered_ids = list(source_by_observation)
    for detail in result.details:
        if detail.pose_id not in source_by_observation:
            source_by_observation[detail.pose_id] = (
                _source_key(detail) or "dataset-1"
            )
            ordered_ids.append(detail.pose_id)
    frame_by_replica: dict[str, int] = {}
    points = []
    for ordinal, observation_id in enumerate(ordered_ids):
        replica_id = source_by_observation[observation_id]
        frame_index = frame_by_replica.get(replica_id, 0)
        frame_by_replica[replica_id] = frame_index + 1
        points.append(
            ObservationPoint(
                observation_id=observation_id,
                ordinal=ordinal,
                frame_index=frame_index,
                time_ns=frame_index * float(time_step_ns),
                replica_id=replica_id,
            )
        )
    return ObservationSeries(
        mode="md",
        points=tuple(points),
        time_step_ns=float(time_step_ns),
    )


__all__ = [
    "LigandGroup",
    "default_md_series_for_result",
    "ligand_groups",
    "subset_observation_series",
    "subset_run_result",
]
