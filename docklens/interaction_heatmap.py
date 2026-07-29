"""Sparse-first interaction heatmaps for ligands, poses and saved frames."""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from pathlib import Path
import re
from types import MappingProxyType
from typing import TYPE_CHECKING, Mapping

import pandas as pd

from .observation_identity import (
    VALID_OBSERVATION_LABEL_MODES,
    observation_labels,
)
from .results import RunResult

if TYPE_CHECKING:
    from .observation_series import ObservationSeries


MAX_HEATMAP_ROWS = 400
MAX_HEATMAP_CELLS = 120_000
_GENERIC_LIGAND_LABELS = frozenset({"", "LIG", "RES", "RES1", "UNK", "UNL"})
_POSE_SUFFIX = re.compile(r"(?i)(?:[_\-\s]*pose[_\-\s]*\d+)$")
_CELL_COLUMNS = (
    "row_id",
    "row_label",
    "source_id",
    "ligand_id",
    "observation_id",
    "receptor_residue",
    "interaction_type",
    "present_observations",
    "total_observations",
    "value_pct",
)


@dataclass(frozen=True)
class InteractionHeatmap:
    """Matrix plus tidy, auditable rows used by a heatmap figure."""

    matrix: pd.DataFrame
    cell_data: pd.DataFrame
    metadata: Mapping[str, object]

    def __post_init__(self) -> None:
        object.__setattr__(self, "matrix", self.matrix.copy(deep=True))
        object.__setattr__(self, "cell_data", self.cell_data.copy(deep=True))
        object.__setattr__(
            self,
            "metadata",
            MappingProxyType(dict(self.metadata)),
        )


@dataclass(frozen=True)
class _HeatmapRow:
    row_id: str
    row_label: str
    source_id: str
    ligand_id: str
    observation_ids: tuple[str, ...]


def _source_key(item) -> str:
    return str(
        getattr(item, "source_id", "")
        or getattr(item, "source_path", "")
        or getattr(item, "source_file", "")
        or "unidentified-source"
    )


def _ordered_observation_ids(result: RunResult) -> tuple[str, ...]:
    ordered: list[str] = []
    seen: set[str] = set()
    for item in tuple(result.summaries) + tuple(result.details):
        observation_id = str(item.pose_id)
        if observation_id and observation_id not in seen:
            seen.add(observation_id)
            ordered.append(observation_id)
    return tuple(ordered)


def _feature(detail, feature_level: str):
    residue = str(detail.receptor_residue)
    if feature_level == "residue":
        return residue
    return residue, str(detail.interaction_type)


def _presence_catalog(
    result: RunResult,
    feature_level: str,
) -> tuple[tuple[str, ...], dict[str, set[object]], tuple[object, ...]]:
    observation_ids = _ordered_observation_ids(result)
    presence: dict[str, set[object]] = {
        observation_id: set() for observation_id in observation_ids
    }
    features: set[object] = set()
    for detail in result.details:
        observation_id = str(detail.pose_id)
        if observation_id not in presence:
            presence[observation_id] = set()
        value = _feature(detail, feature_level)
        presence[observation_id].add(value)
        features.add(value)
    return (
        tuple(presence),
        presence,
        tuple(sorted(features, key=_feature_text)),
    )


def _canonical_ligand(summary) -> str:
    raw = str(summary.ligand_id or "").strip()
    normalized = _POSE_SUFFIX.sub("", raw).strip()
    if normalized.upper() not in _GENERIC_LIGAND_LABELS:
        return normalized
    return Path(str(summary.source_file or summary.source_id)).stem or "Ligand"


def _unique_row_labels(rows: list[dict[str, object]]) -> tuple[_HeatmapRow, ...]:
    counts = Counter(str(row["base_label"]) for row in rows)
    used: set[str] = set()
    output: list[_HeatmapRow] = []
    for row in rows:
        base = str(row["base_label"])
        label = base
        if counts[base] > 1:
            label = f"{base} · {row['source_id'] or row['row_id']}"
        root = label
        suffix = 2
        while label in used:
            label = f"{root} [{suffix}]"
            suffix += 1
        used.add(label)
        output.append(
            _HeatmapRow(
                row_id=str(row["row_id"]),
                row_label=label,
                source_id=str(row["source_id"]),
                ligand_id=str(row["ligand_id"]),
                observation_ids=tuple(row["observation_ids"]),
            )
        )
    return tuple(output)


def _ligand_rows(
    result: RunResult,
    label_mode: str,
) -> tuple[_HeatmapRow, ...]:
    grouped: dict[tuple[str, str], dict[str, object]] = {}
    source_ligands: dict[str, set[str]] = {}
    for summary in result.summaries:
        source_key = _source_key(summary)
        ligand = _canonical_ligand(summary)
        key = source_key, ligand
        source_ligands.setdefault(source_key, set()).add(ligand)
        entry = grouped.setdefault(
            key,
            {
                "source_id": str(summary.source_id),
                "source_file": str(summary.source_file),
                "ligand_id": ligand,
                "observation_ids": [],
            },
        )
        observations = entry["observation_ids"]
        if summary.pose_id not in observations:
            observations.append(str(summary.pose_id))

    rows: list[dict[str, object]] = []
    for ordinal, ((source_key, ligand), entry) in enumerate(
        grouped.items(), 1
    ):
        if label_mode == "index":
            base_label = f"Ligand/file {ordinal}"
        elif label_mode == "file":
            base_label = str(entry["source_file"] or source_key)
            if len(source_ligands[source_key]) > 1:
                base_label = f"{base_label} · {ligand}"
        else:
            base_label = ligand
        rows.append(
            {
                "row_id": f"{source_key}::{ligand}",
                "base_label": base_label,
                "source_id": entry["source_id"],
                "ligand_id": ligand,
                "observation_ids": tuple(entry["observation_ids"]),
            }
        )
    return _unique_row_labels(rows)


def _validate_series(
    observation_ids: tuple[str, ...],
    mode: str,
    series: ObservationSeries | None,
) -> tuple[str, ...]:
    if series is None:
        return observation_ids
    if mode != "md" or series.mode != "md":
        raise ValueError("series is only valid for MD heatmaps")
    if set(series.observation_ids) != set(observation_ids):
        raise ValueError("series must match all heatmap observation IDs")
    return tuple(series.observation_ids)


def _observation_rows(
    result: RunResult,
    observation_ids: tuple[str, ...],
    *,
    label_mode: str,
    mode: str,
    series: ObservationSeries | None,
) -> tuple[_HeatmapRow, ...]:
    ordered_ids = _validate_series(observation_ids, mode, series)
    labels = observation_labels(
        result,
        mode=mode,
        label_mode=label_mode,
        series=series,
    )
    summaries = {item.pose_id: item for item in result.summaries}
    return tuple(
        _HeatmapRow(
            row_id=observation_id,
            row_label=labels[observation_id],
            source_id=str(
                getattr(summaries.get(observation_id), "source_id", "") or ""
            ),
            ligand_id=str(
                getattr(summaries.get(observation_id), "ligand_id", "") or ""
            ),
            observation_ids=(observation_id,),
        )
        for observation_id in ordered_ids
    )


def _feature_text(feature) -> str:
    if isinstance(feature, tuple):
        return " · ".join(str(value) for value in feature)
    return str(feature)


def _rank_features(
    features: tuple[object, ...],
    rows: tuple[_HeatmapRow, ...],
    presence: Mapping[str, set[object]],
    *,
    group_by: str,
) -> tuple[object, ...]:
    scores: dict[object, tuple[float, float, int]] = {}
    if group_by == "source":
        maximum: Counter[object] = Counter()
        normalized_sum: Counter[object] = Counter()
        group_presence: Counter[object] = Counter()
        for row in rows:
            denominator = len(row.observation_ids)
            if not denominator:
                continue
            counts: Counter[object] = Counter()
            for observation_id in row.observation_ids:
                counts.update(presence[observation_id])
            for feature, count in counts.items():
                normalized = count / denominator
                maximum[feature] = max(maximum[feature], normalized)
                normalized_sum[feature] += normalized
                group_presence[feature] += 1
        scores = {
            feature: (
                float(maximum[feature]),
                (
                    float(normalized_sum[feature]) / len(rows)
                    if rows
                    else 0.0
                ),
                int(group_presence[feature]),
            )
            for feature in features
        }
    else:
        counts: Counter[object] = Counter()
        for observed in presence.values():
            counts.update(observed)
        scores = {
            feature: (float(counts[feature]), 0.0, 0)
            for feature in features
        }

    return tuple(
        sorted(
            features,
            key=lambda feature: (
                *(-value for value in scores[feature]),
                _feature_text(feature),
            ),
        )
    )


def _even_positions(total: int, limit: int) -> tuple[int, ...]:
    if total <= limit:
        return tuple(range(total))
    if limit <= 1:
        return (0,)
    return tuple(
        round(position * (total - 1) / (limit - 1))
        for position in range(limit)
    )


def _apply_limits(
    rows: tuple[_HeatmapRow, ...],
    ranked_features: tuple[object, ...],
    *,
    group_by: str,
    top_n: int | None,
) -> tuple[
    tuple[_HeatmapRow, ...],
    tuple[object, ...],
    bool,
    bool,
    int,
]:
    user_feature_count = min(
        len(ranked_features),
        top_n if top_n is not None else len(ranked_features),
    )
    row_limit = MAX_HEATMAP_ROWS
    if group_by == "observation" and user_feature_count:
        row_limit = min(
            row_limit,
            max(1, MAX_HEATMAP_CELLS // user_feature_count),
        )
    positions = _even_positions(len(rows), row_limit)
    selected_rows = tuple(rows[position] for position in positions)
    rows_sampled = len(selected_rows) < len(rows)
    hard_feature_count = max(
        1,
        MAX_HEATMAP_CELLS // max(1, len(selected_rows)),
    )
    selected_feature_count = min(user_feature_count, hard_feature_count)
    selected_features = ranked_features[:selected_feature_count]
    hard_feature_limit = selected_feature_count < user_feature_count
    return (
        selected_rows,
        selected_features,
        rows_sampled,
        hard_feature_limit,
        user_feature_count,
    )


def _materialize_matrix(
    rows: tuple[_HeatmapRow, ...],
    features: tuple[object, ...],
    presence: Mapping[str, set[object]],
    *,
    group_by: str,
) -> pd.DataFrame:
    values: list[list[float]] = []
    for row in rows:
        denominator = len(row.observation_ids)
        values.append(
            [
                (
                    100.0
                    * sum(
                        feature in presence[observation_id]
                        for observation_id in row.observation_ids
                    )
                    / denominator
                    if denominator
                    else 0.0
                )
                for feature in features
            ]
        )
    return pd.DataFrame(
        values,
        index=[row.row_label for row in rows],
        columns=features,
        dtype=float,
    )


def _tidy_cells(
    matrix: pd.DataFrame,
    rows: tuple[_HeatmapRow, ...],
    *,
    feature_level: str,
) -> pd.DataFrame:
    records: list[dict[str, object]] = []
    for row in rows:
        denominator = len(row.observation_ids)
        observation_id = (
            row.observation_ids[0]
            if len(row.observation_ids) == 1
            else ""
        )
        for feature in matrix.columns:
            value = float(matrix.loc[row.row_label, feature])
            if feature_level == "residue_type":
                residue, interaction_type = feature
            else:
                residue, interaction_type = feature, ""
            records.append(
                {
                    "row_id": row.row_id,
                    "row_label": row.row_label,
                    "source_id": row.source_id,
                    "ligand_id": row.ligand_id,
                    "observation_id": observation_id,
                    "receptor_residue": str(residue),
                    "interaction_type": str(interaction_type),
                    "present_observations": round(
                        denominator * value / 100.0
                    ),
                    "total_observations": denominator,
                    "value_pct": value,
                }
            )
    return pd.DataFrame(records, columns=_CELL_COLUMNS)


def build_interaction_heatmap_data(
    result: RunResult,
    *,
    group_by: str = "source",
    feature_level: str = "residue_type",
    label_mode: str = "ligand",
    mode: str = "docking",
    series: ObservationSeries | None = None,
    top_n: int | None = 40,
) -> InteractionHeatmap:
    """Aggregate semantic contacts without materializing an unbounded matrix."""

    if group_by not in {"source", "observation"}:
        raise ValueError("group_by must be 'source' or 'observation'")
    if feature_level not in {"residue_type", "residue"}:
        raise ValueError(
            "feature_level must be 'residue_type' or 'residue'"
        )
    if mode not in {"docking", "md"}:
        raise ValueError("mode must be 'docking' or 'md'")
    if label_mode not in VALID_OBSERVATION_LABEL_MODES:
        raise ValueError("label_mode is not supported")
    if top_n is not None and (
        isinstance(top_n, bool) or not isinstance(top_n, int) or top_n < 1
    ):
        raise ValueError("top_n must be a positive integer or None")

    observation_ids, presence, features = _presence_catalog(
        result, feature_level
    )
    rows = (
        _ligand_rows(result, label_mode)
        if group_by == "source"
        else _observation_rows(
            result,
            observation_ids,
            label_mode=label_mode,
            mode=mode,
            series=series,
        )
    )
    ranked_features = _rank_features(
        features,
        rows,
        presence,
        group_by=group_by,
    )
    (
        selected_rows,
        selected_features,
        rows_sampled,
        hard_feature_limit,
        features_after_user_limit,
    ) = _apply_limits(
        rows,
        ranked_features,
        group_by=group_by,
        top_n=top_n,
    )
    matrix = _materialize_matrix(
        selected_rows,
        selected_features,
        presence,
        group_by=group_by,
    )
    tidy = _tidy_cells(
        matrix,
        selected_rows,
        feature_level=feature_level,
    )
    feature_counting = (
        "binary presence per observation, receptor residue and interaction type"
        if feature_level == "residue_type"
        else "binary any-interaction presence per observation and residue"
    )
    return InteractionHeatmap(
        matrix=matrix,
        cell_data=tidy,
        metadata={
            "mode": mode,
            "group_by": group_by,
            "feature_level": feature_level,
            "label_mode": label_mode,
            "top_n": top_n,
            "features_before_limit": len(features),
            "features_after_user_limit": features_after_user_limit,
            "features_displayed": len(matrix.columns),
            "rows_before_limit": len(rows),
            "rows_displayed": len(matrix.index),
            "rows_sampled": rows_sampled,
            "hard_feature_limit_applied": hard_feature_limit,
            "hard_cell_limit": MAX_HEATMAP_CELLS,
            "total_observations": len(observation_ids),
            "row_denominators": {
                row.row_label: len(row.observation_ids)
                for row in selected_rows
            },
            "value_unit": "percent",
            "feature_counting": feature_counting,
            "feature_ranking": (
                "maximum normalized row prevalence, then mean row prevalence"
                if group_by == "source"
                else "observation prevalence"
            ),
            "zero_contact_observations_included": True,
            "row_order": (
                "uploaded ligand/source order"
                if group_by == "source"
                else (
                    "explicit trajectory-map order"
                    if series is not None
                    else "observation order"
                )
            ),
        },
    )


__all__ = [
    "InteractionHeatmap",
    "MAX_HEATMAP_CELLS",
    "MAX_HEATMAP_ROWS",
    "build_interaction_heatmap_data",
]
