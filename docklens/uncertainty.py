"""Uncertainty estimates for interaction occupancies in ordered MD frames.

The implementation uses a circular moving-block bootstrap so that consecutive
saved frames are resampled together.  This preserves short-range temporal
dependence better than treating every frame as an independent observation.
"""

from __future__ import annotations

from math import ceil, sqrt
from numbers import Integral, Real
from typing import Final

import numpy as np
import pandas as pd

from .analytics import AnalysisContext
from .observation_series import ObservationSeries


_METHOD: Final = "circular moving-block bootstrap"
_MIN_ITERATIONS: Final = 100
_MIN_RELIABLE_FRAMES: Final = 8

_OCCUPANCY_COLUMNS: Final = [
    "receptor_residue",
    "interaction_type",
    "occupancy_pct",
    "ci_low_pct",
    "ci_high_pct",
    "observation_count",
    "total_observations",
    "confidence_level_pct",
    "block_size",
    "iterations",
    "seed",
    "method",
    "replica_count",
    "insufficient_data",
    "warning",
]

_DIFFERENCE_COLUMNS: Final = [
    "receptor_residue",
    "interaction_type",
    "occupancy_a_pct",
    "occupancy_b_pct",
    "delta_pct_points",
    "ci_low_pct_points",
    "ci_high_pct_points",
    "observation_count_a",
    "observation_count_b",
    "total_observations_a",
    "total_observations_b",
    "confidence_level_pct",
    "block_size_a",
    "block_size_b",
    "iterations",
    "seed",
    "method",
    "delta_definition",
    "replica_count_a",
    "replica_count_b",
    "insufficient_data",
    "warning",
]


def _validate_context(context: AnalysisContext) -> None:
    if not isinstance(context, AnalysisContext) or context.mode != "md":
        raise ValueError("block bootstrap uncertainty is available only for MD")


def _validate_iterations(iterations: int) -> int:
    if (
        isinstance(iterations, bool)
        or not isinstance(iterations, Integral)
        or iterations < _MIN_ITERATIONS
    ):
        raise ValueError(
            f"iterations must be an integer of at least {_MIN_ITERATIONS}"
        )
    return int(iterations)


def _validate_seed(seed: int) -> int:
    if isinstance(seed, bool) or not isinstance(seed, Integral) or seed < 0:
        raise ValueError("seed must be a non-negative integer")
    return int(seed)


def _validate_confidence_level(confidence_level: float) -> float:
    if (
        isinstance(confidence_level, bool)
        or not isinstance(confidence_level, Real)
        or not 0.0 < float(confidence_level) < 1.0
    ):
        raise ValueError("confidence_level must be between zero and one")
    return float(confidence_level)


def _validated_matrix(
    matrix: pd.DataFrame,
) -> tuple[np.ndarray, tuple[tuple[str, str], ...]]:
    if not isinstance(matrix, pd.DataFrame):
        raise TypeError("matrix must be a pandas DataFrame")
    if matrix.empty and len(matrix.index) == 0:
        raise ValueError("matrix must contain at least one saved MD frame")
    if not matrix.index.is_unique:
        raise ValueError("matrix frame identifiers must be unique")
    if not matrix.columns.is_unique:
        raise ValueError("matrix interaction features must be unique")
    if not isinstance(matrix.columns, pd.MultiIndex) or matrix.columns.nlevels != 2:
        raise ValueError(
            "matrix columns must identify receptor residue and interaction type"
        )
    if matrix.isna().to_numpy().any():
        raise ValueError("matrix must not contain missing interaction states")

    raw = matrix.to_numpy(copy=True)
    if raw.size and not np.isin(raw, (False, True, 0, 1)).all():
        raise ValueError("matrix interaction states must be binary")
    values = raw.astype(bool, copy=False)
    features = tuple((str(residue), str(kind)) for residue, kind in matrix.columns)
    if len(set(features)) != len(features):
        raise ValueError(
            "matrix interaction features must remain unique as tidy labels"
        )
    return values, features


def _resolved_block_size(block_size: int | None, max_block_size: int) -> int:
    value = (
        max(1, int(round(sqrt(max_block_size))))
        if block_size is None
        else block_size
    )
    if (
        isinstance(value, bool)
        or not isinstance(value, Integral)
        or not 1 <= value <= max_block_size
    ):
        raise ValueError(
            "block_size must be an integer between one and the number of "
            "frames in the shortest replica"
        )
    return int(value)


def _data_warning(
    frame_count: int,
    block_size: int,
    replica_indices: tuple[np.ndarray, ...],
) -> tuple[bool, str]:
    warnings = []
    if frame_count < _MIN_RELIABLE_FRAMES:
        warnings.append(
            f"Trajectory has fewer than {_MIN_RELIABLE_FRAMES} saved frames; "
            "the interval is exploratory."
        )
    effective_blocks = sum(
        ceil(len(replica) / block_size) for replica in replica_indices
    )
    if effective_blocks < 2:
        warnings.append(
            "Block size leaves fewer than two resampled blocks; "
            "the interval may understate uncertainty."
        )
    return bool(warnings), " ".join(warnings)


def _replica_indices(
    matrix: pd.DataFrame,
    series: ObservationSeries | None,
) -> tuple[np.ndarray, ...]:
    if series is None:
        return (np.arange(len(matrix.index), dtype=np.intp),)
    if not isinstance(series, ObservationSeries) or series.mode != "md":
        raise ValueError("series must be an MD ObservationSeries")
    matrix_ids = tuple(matrix.index)
    if set(series.observation_ids) != set(matrix_ids):
        raise ValueError("observation series IDs must match the matrix index")
    positions = {observation_id: index for index, observation_id in enumerate(matrix_ids)}
    grouped: dict[str, list[int]] = {}
    for point in series.points:
        grouped.setdefault(point.replica_id, []).append(
            positions[point.observation_id]
        )
    if not grouped:
        raise ValueError("observation series must contain at least one replica")
    return tuple(np.asarray(indices, dtype=np.intp) for indices in grouped.values())


def _bootstrap_means(
    values: np.ndarray,
    *,
    iterations: int,
    block_size: int,
    replica_indices: tuple[np.ndarray, ...],
    rng: np.random.Generator,
) -> np.ndarray:
    """Return bootstrap means without materializing a 3-D sample array."""
    frame_count, feature_count = values.shape
    offsets = np.arange(block_size, dtype=np.intp)
    replica_starts = tuple(
        rng.integers(
            0,
            len(replica),
            size=(iterations, ceil(len(replica) / block_size)),
            dtype=np.intp,
        )
        for replica in replica_indices
    )
    means = np.empty((iterations, feature_count), dtype=float)
    numeric_values = values.astype(float, copy=False)
    for iteration in range(iterations):
        weights = np.zeros(frame_count, dtype=np.intp)
        for replica, starts in zip(replica_indices, replica_starts):
            local_indices = (
                starts[iteration, :, np.newaxis] + offsets[np.newaxis, :]
            ) % len(replica)
            sampled_indices = replica[local_indices.reshape(-1)[: len(replica)]]
            weights += np.bincount(sampled_indices, minlength=frame_count)
        means[iteration] = weights @ numeric_values / frame_count
    return means


def _interval(
    samples: np.ndarray, confidence_level: float
) -> tuple[np.ndarray, np.ndarray]:
    tail = (1.0 - confidence_level) / 2.0
    bounds = np.quantile(samples, (tail, 1.0 - tail), axis=0)
    return bounds[0], bounds[1]


def block_bootstrap_occupancy(
    matrix: pd.DataFrame,
    context: AnalysisContext,
    *,
    iterations: int = 2000,
    block_size: int | None = None,
    seed: int = 2026,
    confidence_level: float = 0.95,
    series: ObservationSeries | None = None,
) -> pd.DataFrame:
    """Estimate occupancy confidence intervals for ordered saved MD frames.

    Each output row is one receptor-residue/interaction-type feature.  The
    point estimate is calculated from the original frames, while percentile
    confidence limits are calculated from circular moving-block resamples.
    """
    _validate_context(context)
    validated_iterations = _validate_iterations(iterations)
    validated_seed = _validate_seed(seed)
    validated_level = _validate_confidence_level(confidence_level)
    values, features = _validated_matrix(matrix)
    frame_count = values.shape[0]
    replicas = _replica_indices(matrix, series)
    resolved_block = _resolved_block_size(
        block_size, min(len(replica) for replica in replicas)
    )
    insufficient, warning = _data_warning(
        frame_count, resolved_block, replicas
    )

    if not features:
        return pd.DataFrame(columns=_OCCUPANCY_COLUMNS)

    samples = _bootstrap_means(
        values,
        iterations=validated_iterations,
        block_size=resolved_block,
        replica_indices=replicas,
        rng=np.random.default_rng(validated_seed),
    )
    low, high = _interval(samples, validated_level)
    occupancy = values.mean(axis=0)
    counts = values.sum(axis=0)
    rows = [
        {
            "receptor_residue": residue,
            "interaction_type": kind,
            "occupancy_pct": float(100.0 * occupancy[index]),
            "ci_low_pct": float(100.0 * low[index]),
            "ci_high_pct": float(100.0 * high[index]),
            "observation_count": int(counts[index]),
            "total_observations": frame_count,
            "confidence_level_pct": 100.0 * validated_level,
            "block_size": resolved_block,
            "iterations": validated_iterations,
            "seed": validated_seed,
            "method": _METHOD,
            "replica_count": len(replicas),
            "insufficient_data": insufficient,
            "warning": warning,
        }
        for index, (residue, kind) in enumerate(features)
    ]
    return pd.DataFrame(rows, columns=_OCCUPANCY_COLUMNS)


def block_bootstrap_difference(
    matrix_a: pd.DataFrame,
    matrix_b: pd.DataFrame,
    context: AnalysisContext,
    *,
    iterations: int = 2000,
    block_size_a: int | None = None,
    block_size_b: int | None = None,
    seed: int = 2026,
    confidence_level: float = 0.95,
    series_a: ObservationSeries | None = None,
    series_b: ObservationSeries | None = None,
) -> pd.DataFrame:
    """Estimate the occupancy difference B - A for independent MD systems."""
    _validate_context(context)
    validated_iterations = _validate_iterations(iterations)
    validated_seed = _validate_seed(seed)
    validated_level = _validate_confidence_level(confidence_level)
    values_a, features_a = _validated_matrix(matrix_a)
    values_b, features_b = _validated_matrix(matrix_b)
    replicas_a = _replica_indices(matrix_a, series_a)
    replicas_b = _replica_indices(matrix_b, series_b)
    resolved_a = _resolved_block_size(
        block_size_a, min(len(replica) for replica in replicas_a)
    )
    resolved_b = _resolved_block_size(
        block_size_b, min(len(replica) for replica in replicas_b)
    )

    features = tuple(sorted(set(features_a) | set(features_b)))
    if not features:
        return pd.DataFrame(columns=_DIFFERENCE_COLUMNS)

    index_a = {feature: index for index, feature in enumerate(features_a)}
    index_b = {feature: index for index, feature in enumerate(features_b)}
    aligned_a = np.zeros((values_a.shape[0], len(features)), dtype=bool)
    aligned_b = np.zeros((values_b.shape[0], len(features)), dtype=bool)
    for target, feature in enumerate(features):
        if feature in index_a:
            aligned_a[:, target] = values_a[:, index_a[feature]]
        if feature in index_b:
            aligned_b[:, target] = values_b[:, index_b[feature]]

    child_seeds = np.random.SeedSequence(validated_seed).spawn(2)
    samples_a = _bootstrap_means(
        aligned_a,
        iterations=validated_iterations,
        block_size=resolved_a,
        replica_indices=replicas_a,
        rng=np.random.default_rng(child_seeds[0]),
    )
    samples_b = _bootstrap_means(
        aligned_b,
        iterations=validated_iterations,
        block_size=resolved_b,
        replica_indices=replicas_b,
        rng=np.random.default_rng(child_seeds[1]),
    )
    difference_samples = samples_b - samples_a
    low, high = _interval(difference_samples, validated_level)
    occupancy_a = aligned_a.mean(axis=0)
    occupancy_b = aligned_b.mean(axis=0)
    counts_a = aligned_a.sum(axis=0)
    counts_b = aligned_b.sum(axis=0)
    insufficient_a, warning_a = _data_warning(
        values_a.shape[0], resolved_a, replicas_a
    )
    insufficient_b, warning_b = _data_warning(
        values_b.shape[0], resolved_b, replicas_b
    )
    warnings = " ".join(
        value
        for value in (
            f"System A: {warning_a}" if warning_a else "",
            f"System B: {warning_b}" if warning_b else "",
        )
        if value
    )

    rows = [
        {
            "receptor_residue": residue,
            "interaction_type": kind,
            "occupancy_a_pct": float(100.0 * occupancy_a[index]),
            "occupancy_b_pct": float(100.0 * occupancy_b[index]),
            "delta_pct_points": float(
                100.0 * (occupancy_b[index] - occupancy_a[index])
            ),
            "ci_low_pct_points": float(100.0 * low[index]),
            "ci_high_pct_points": float(100.0 * high[index]),
            "observation_count_a": int(counts_a[index]),
            "observation_count_b": int(counts_b[index]),
            "total_observations_a": values_a.shape[0],
            "total_observations_b": values_b.shape[0],
            "confidence_level_pct": 100.0 * validated_level,
            "block_size_a": resolved_a,
            "block_size_b": resolved_b,
            "iterations": validated_iterations,
            "seed": validated_seed,
            "method": _METHOD,
            "delta_definition": "B - A",
            "replica_count_a": len(replicas_a),
            "replica_count_b": len(replicas_b),
            "insufficient_data": insufficient_a or insufficient_b,
            "warning": warnings,
        }
        for index, (residue, kind) in enumerate(features)
    ]
    return pd.DataFrame(rows, columns=_DIFFERENCE_COLUMNS)


__all__ = [
    "block_bootstrap_difference",
    "block_bootstrap_occupancy",
]
