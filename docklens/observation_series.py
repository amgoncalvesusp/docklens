"""Explicit observation axes for docking poses and molecular-dynamics frames."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable

import pandas as pd


_MODES = frozenset({"docking", "md"})


@dataclass(frozen=True)
class ObservationPoint:
    """One ordered observation, optionally located in a trajectory replica."""

    observation_id: str
    ordinal: int
    frame_index: int | None = None
    time_ns: float | None = None
    replica_id: str = "replica-1"

    def __post_init__(self):
        if not str(self.observation_id):
            raise ValueError("observation_id must not be empty")
        if self.ordinal < 0:
            raise ValueError("ordinal must be non-negative")
        if self.frame_index is not None and self.frame_index < 0:
            raise ValueError("frame_index must be non-negative")
        if self.time_ns is not None and self.time_ns < 0:
            raise ValueError("time_ns must be non-negative")
        if not str(self.replica_id):
            raise ValueError("replica_id must not be empty")


@dataclass(frozen=True)
class ObservationSeries:
    """Immutable ordering and replica boundaries for an analytical dataset."""

    mode: str
    points: tuple[ObservationPoint, ...]
    time_step_ns: float | None = None

    def __post_init__(self):
        if self.mode not in _MODES:
            raise ValueError("mode must be 'docking' or 'md'")
        if self.time_step_ns is not None and self.time_step_ns <= 0:
            raise ValueError("time_step_ns must be greater than zero")
        points = tuple(self.points)
        object.__setattr__(self, "points", points)
        identifiers = tuple(point.observation_id for point in points)
        if len(set(identifiers)) != len(identifiers):
            raise ValueError("observation IDs must be unique")
        ordinals = tuple(point.ordinal for point in points)
        if len(set(ordinals)) != len(ordinals):
            raise ValueError("ordinals must be unique")
        if ordinals != tuple(sorted(ordinals)):
            raise ValueError("points must be ordered by ordinal")
        previous_by_replica: dict[str, ObservationPoint] = {}
        for point in points:
            previous = previous_by_replica.get(point.replica_id)
            if previous is not None:
                if (
                    previous.frame_index is not None
                    and point.frame_index is not None
                    and point.frame_index <= previous.frame_index
                ):
                    raise ValueError(
                        "frame_index must increase within each replica"
                    )
                if (
                    previous.time_ns is not None
                    and point.time_ns is not None
                    and point.time_ns <= previous.time_ns
                ):
                    raise ValueError("time_ns must increase within each replica")
            previous_by_replica[point.replica_id] = point

    @property
    def observation_ids(self) -> tuple[str, ...]:
        return tuple(point.observation_id for point in self.points)

    def point_map(self) -> dict[str, ObservationPoint]:
        return {point.observation_id: point for point in self.points}

    def transition_pairs(
        self,
        *,
        lag: int = 1,
        require_consecutive_frames: bool = True,
    ) -> tuple[tuple[str, str], ...]:
        """Return observed transitions without crossing replica boundaries."""
        if self.mode != "md":
            raise ValueError("temporal transitions are available only for MD")
        if lag < 1:
            raise ValueError("lag must be at least one")
        grouped: dict[str, list[ObservationPoint]] = {}
        replica_order = []
        for point in self.points:
            if point.replica_id not in grouped:
                grouped[point.replica_id] = []
                replica_order.append(point.replica_id)
            grouped[point.replica_id].append(point)
        pairs = []
        for replica_id in replica_order:
            replica = grouped[replica_id]
            for index in range(len(replica) - lag):
                left = replica[index]
                right = replica[index + lag]
                if (
                    require_consecutive_frames
                    and left.frame_index is not None
                    and right.frame_index is not None
                    and right.frame_index - left.frame_index != lag
                ):
                    continue
                pairs.append((left.observation_id, right.observation_id))
        return tuple(pairs)


def default_observation_series(
    observation_ids: Iterable[str],
    *,
    mode: str,
    time_step_ns: float = 1.0,
) -> ObservationSeries:
    """Create the disclosed single-series fallback used by the desktop UI."""
    identifiers = tuple(str(value) for value in observation_ids)
    points = tuple(
        ObservationPoint(
            observation_id=identifier,
            ordinal=index,
            frame_index=(index if mode == "md" else None),
            time_ns=(index * time_step_ns if mode == "md" else None),
        )
        for index, identifier in enumerate(identifiers)
    )
    return ObservationSeries(
        mode=mode,
        points=points,
        time_step_ns=(time_step_ns if mode == "md" else None),
    )


def observation_series_from_dataframe(
    frame: pd.DataFrame,
    *,
    expected_ids: Iterable[str],
    default_time_step_ns: float,
) -> ObservationSeries:
    """Build an MD axis from a validated user-supplied trajectory map.

    Row order is the analytical order. ``time_ns`` is optional; when absent,
    it is derived from ``frame_index`` and the disclosed saved-frame step.
    """

    required = {"observation_id", "replica_id", "frame_index"}
    missing_columns = required.difference(frame.columns)
    if missing_columns:
        raise ValueError(
            "trajectory map is missing columns: "
            + ", ".join(sorted(missing_columns))
        )
    if not pd.api.types.is_numeric_dtype(type(default_time_step_ns)):
        raise ValueError("default_time_step_ns must be numeric")
    step = float(default_time_step_ns)
    if not pd.notna(step) or step <= 0:
        raise ValueError("default_time_step_ns must be greater than zero")

    identifiers = tuple(str(value).strip() for value in frame["observation_id"])
    expected = tuple(str(value) for value in expected_ids)
    if (
        len(identifiers) != len(expected)
        or len(set(identifiers)) != len(identifiers)
        or set(identifiers) != set(expected)
    ):
        raise ValueError(
            "trajectory map must contain exactly the analyzed observation IDs"
        )
    replica_ids = tuple(str(value).strip() for value in frame["replica_id"])
    if any(not value for value in identifiers + replica_ids):
        raise ValueError("observation_id and replica_id must not be empty")

    numeric_frames = pd.to_numeric(frame["frame_index"], errors="coerce")
    if (
        numeric_frames.isna().any()
        or (numeric_frames < 0).any()
        or (numeric_frames % 1 != 0).any()
    ):
        raise ValueError("frame_index values must be non-negative integers")
    frame_indices = tuple(int(value) for value in numeric_frames)

    if "time_ns" in frame.columns:
        numeric_times = pd.to_numeric(frame["time_ns"], errors="coerce")
        if numeric_times.isna().any() or (numeric_times < 0).any():
            raise ValueError("time_ns values must be non-negative numbers")
        times = tuple(float(value) for value in numeric_times)
    else:
        times = tuple(index * step for index in frame_indices)

    points = tuple(
        ObservationPoint(
            observation_id=identifier,
            ordinal=ordinal,
            frame_index=frame_indices[ordinal],
            time_ns=times[ordinal],
            replica_id=replica_ids[ordinal],
        )
        for ordinal, identifier in enumerate(identifiers)
    )
    return ObservationSeries(mode="md", points=points, time_step_ns=step)


__all__ = [
    "ObservationPoint",
    "ObservationSeries",
    "default_observation_series",
    "observation_series_from_dataframe",
]
