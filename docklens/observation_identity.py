"""Human-readable observation labels separated from scientific identifiers."""

from __future__ import annotations

from collections import Counter
from pathlib import Path
import re
from types import MappingProxyType
from typing import TYPE_CHECKING, Mapping

from .ligand_selection import ligand_groups
from .results import RunResult

if TYPE_CHECKING:
    from .observation_series import ObservationSeries


VALID_OBSERVATION_LABEL_MODES = frozenset({"ligand", "file", "index"})
_GENERIC_LIGAND_LABELS = frozenset({"", "LIG", "RES", "RES1", "UNK", "UNL"})
_POSE_SUFFIX = re.compile(r"(?i)(?:[_\-\s]*pose[_\-\s]*\d+)$")


def _ordered_observation_ids(result: RunResult) -> tuple[str, ...]:
    ordered: list[str] = []
    seen: set[str] = set()
    for item in tuple(result.summaries) + tuple(result.details):
        observation_id = str(item.pose_id)
        if observation_id and observation_id not in seen:
            seen.add(observation_id)
            ordered.append(observation_id)
    return tuple(ordered)


def _summary_by_observation(result: RunResult) -> dict[str, object]:
    return {
        str(summary.pose_id): summary
        for summary in result.summaries
        if summary.pose_id
    }


def _series_points(series: ObservationSeries | None) -> dict[str, object]:
    if series is None:
        return {}
    return {point.observation_id: point for point in series.points}


def _source_filename(summary) -> str:
    return str(summary.source_file or summary.source_id or "Uploaded file")


def _ligand_name(summary) -> str:
    value = str(summary.ligand_id or "").strip()
    if value.upper() in _GENERIC_LIGAND_LABELS:
        return _source_filename(summary)
    return value


def _observation_suffix(
    *,
    mode: str,
    ordinal: int,
    summary,
    point,
    include_replica: bool,
) -> str:
    if mode == "docking":
        pose = getattr(summary, "pose", None)
        return f"Pose {pose if pose is not None else ordinal + 1}"
    frame = getattr(point, "frame_index", None)
    text = f"Frame {frame if frame is not None else ordinal}"
    replica = str(getattr(point, "replica_id", "") or "")
    if include_replica and replica:
        text += f" · {replica}"
    return text


def _make_unique(
    labels: list[str],
    observation_ids: tuple[str, ...],
) -> tuple[str, ...]:
    counts = Counter(labels)
    seen: Counter[str] = Counter()
    unique: list[str] = []
    for label, observation_id in zip(labels, observation_ids):
        if counts[label] == 1:
            unique.append(label)
            continue
        seen[label] += 1
        suffix = observation_id.split(":", 1)[0] or str(seen[label])
        candidate = f"{label} · {suffix}"
        if candidate in unique:
            candidate = f"{candidate} [{seen[label]}]"
        unique.append(candidate)
    return tuple(unique)


def observation_labels(
    result: RunResult,
    *,
    mode: str = "docking",
    label_mode: str = "ligand",
    series: ObservationSeries | None = None,
) -> Mapping[str, str]:
    """Map stable observation IDs to unique, human-first display labels."""

    if mode not in {"docking", "md"}:
        raise ValueError("mode must be 'docking' or 'md'")
    if label_mode not in VALID_OBSERVATION_LABEL_MODES:
        raise ValueError("label_mode is not supported")
    if series is not None and series.mode != "md":
        raise ValueError("series must use MD mode")

    observation_ids = _ordered_observation_ids(result)
    summaries = _summary_by_observation(result)
    points = _series_points(series)
    replicas = {
        str(point.replica_id)
        for point in points.values()
        if getattr(point, "replica_id", "")
    }
    include_replica = mode == "md" and len(replicas) > 1

    bases: list[str] = []
    suffixes: list[str] = []
    for ordinal, observation_id in enumerate(observation_ids):
        summary = summaries.get(observation_id)
        point = points.get(observation_id)
        suffix = _observation_suffix(
            mode=mode,
            ordinal=ordinal,
            summary=summary,
            point=point,
            include_replica=include_replica,
        )
        suffixes.append(suffix)
        if label_mode == "index":
            if mode == "docking":
                bases.append(f"Pose {ordinal + 1}")
            else:
                bases.append(suffix)
        elif summary is None:
            bases.append(observation_id)
        elif label_mode == "file":
            bases.append(_source_filename(summary))
        else:
            bases.append(_ligand_name(summary))

    if label_mode != "index":
        counts = Counter(bases)
        bases = [
            f"{base} · {suffix}" if counts[base] > 1 else base
            for base, suffix in zip(bases, suffixes)
        ]
    unique = _make_unique(bases, observation_ids)
    return MappingProxyType(dict(zip(observation_ids, unique)))


def _group_ligand_label(group) -> str:
    labels = tuple(
        value.strip() for value in group.ligand_ids if str(value).strip()
    )
    normalized = tuple(_POSE_SUFFIX.sub("", value) for value in labels)
    distinct = tuple(dict.fromkeys(normalized))
    if len(distinct) == 1 and distinct[0].upper() not in _GENERIC_LIGAND_LABELS:
        return distinct[0]
    if len(labels) == 1 and labels[0].upper() not in _GENERIC_LIGAND_LABELS:
        return labels[0]
    return Path(group.source_file or group.source_id).stem or "Ligand"


def source_group_labels(
    result: RunResult,
    *,
    label_mode: str = "ligand",
) -> Mapping[str, str]:
    """Return deterministic labels for uploaded ligand/file groups."""

    if label_mode not in VALID_OBSERVATION_LABEL_MODES:
        raise ValueError("label_mode is not supported")
    groups = ligand_groups(result)
    bases: list[str] = []
    for ordinal, group in enumerate(groups, 1):
        if label_mode == "file":
            bases.append(group.source_file or group.source_id or f"File {ordinal}")
        elif label_mode == "index":
            bases.append(f"Ligand/file {ordinal}")
        else:
            bases.append(_group_ligand_label(group))
    counts = Counter(bases)
    labels = [
        (
            f"{base} · {group.source_id or ordinal}"
            if counts[base] > 1
            else base
        )
        for ordinal, (base, group) in enumerate(zip(bases, groups), 1)
    ]
    return MappingProxyType(
        {group.key: label for group, label in zip(groups, labels)}
    )


__all__ = [
    "VALID_OBSERVATION_LABEL_MODES",
    "observation_labels",
    "source_group_labels",
]
