"""Matplotlib figure builders backed by the auditable analytics layer."""

from __future__ import annotations

from dataclasses import dataclass
from types import MappingProxyType
from typing import Mapping

import numpy as np
import pandas as pd
from matplotlib.colors import to_rgba
from matplotlib.figure import Figure

from .analytics import (
    differential_prevalence,
    docking_md_retention,
    fingerprint_matrix,
    fingerprint_similarity,
    residue_type_prevalence,
)
from .interaction_core import color_hex
from .interaction_heatmap import build_interaction_heatmap_data
from .results import RunResult


_INK = "#12202F"
_MUTED = "#607080"
_GRID = "#DCE3E8"
_PAPER = "#FFFFFF"


@dataclass(frozen=True)
class ChartArtifact:
    """A figure and the exact tidy rows used to draw it."""

    kind: str
    figure: Figure
    data: pd.DataFrame
    metadata: Mapping[str, object]

    def __post_init__(self):
        object.__setattr__(self, "data", self.data.copy(deep=True))
        object.__setattr__(
            self, "metadata", MappingProxyType(dict(self.metadata))
        )


def _figure(width=10.0, height=5.4) -> Figure:
    figure = Figure(figsize=(width, height), facecolor=_PAPER, layout="constrained")
    return figure


def _style_axis(axis):
    axis.set_facecolor(_PAPER)
    axis.tick_params(colors=_MUTED, labelsize=9)
    axis.xaxis.label.set_color(_INK)
    axis.yaxis.label.set_color(_INK)
    axis.title.set_color(_INK)
    axis.grid(axis="x", color=_GRID, linewidth=0.7, alpha=0.8)
    axis.set_axisbelow(True)
    axis.spines["top"].set_visible(False)
    axis.spines["right"].set_visible(False)
    axis.spines["left"].set_color(_GRID)
    axis.spines["bottom"].set_color(_GRID)


def _empty_axis(axis, message):
    axis.text(
        0.5,
        0.5,
        message,
        ha="center",
        va="center",
        color=_MUTED,
        transform=axis.transAxes,
    )
    axis.set_xticks([])
    axis.set_yticks([])


def build_residue_chart(
    result: RunResult,
    *,
    mode: str = "docking",
    top_n: int = 20,
    prevalence: pd.DataFrame | None = None,
) -> ChartArtifact:
    """Build stacked residue/type bars using consolidated observations."""
    if mode not in {"docking", "md"}:
        raise ValueError("mode must be 'docking' or 'md'")
    if top_n < 1:
        raise ValueError("top_n must be at least one")
    data = (
        prevalence.copy(deep=True)
        if prevalence is not None
        else residue_type_prevalence(result)
    )
    figure = _figure()
    axis = figure.add_subplot(111)
    _style_axis(axis)
    label = (
        "Docking frequency (% of poses)"
        if mode == "docking"
        else "MD occupancy (% of saved frames)"
    )
    if data.empty:
        _empty_axis(axis, "No interaction evidence for this selection")
    else:
        order = (
            data.groupby("receptor_residue")["prevalence_pct"]
            .sum()
            .sort_values(ascending=False)
            .head(top_n)
            .index.tolist()
        )
        selected = data[data["receptor_residue"].isin(order)]
        pivot = selected.pivot_table(
            index="receptor_residue",
            columns="interaction_type",
            values="prevalence_pct",
            fill_value=0.0,
            aggfunc="sum",
        ).reindex(order)
        y_positions = np.arange(len(pivot))
        left = np.zeros(len(pivot), dtype=float)
        for kind in pivot.columns:
            values = pivot[kind].to_numpy(dtype=float)
            axis.barh(
                y_positions,
                values,
                left=left,
                height=0.66,
                label=kind,
                color=color_hex(kind),
                edgecolor=_PAPER,
                linewidth=0.5,
            )
            left += values
        axis.set_yticks(y_positions, labels=pivot.index)
        axis.invert_yaxis()
        axis.legend(
            loc="upper center",
            bbox_to_anchor=(0.5, -0.14),
            ncols=min(4, max(1, len(pivot.columns))),
            frameon=False,
            fontsize=8,
        )
    axis.set_xlabel(label)
    axis.set_title("Interaction profile by receptor residue", loc="left")
    total = len({summary.pose_id for summary in result.summaries})
    return ChartArtifact(
        kind="residue-prevalence",
        figure=figure,
        data=data,
        metadata={
            "mode": mode,
            "total_observations": total,
            "denominator": "poses" if mode == "docking" else "saved frames",
            "counting_unit": (
                "observation × receptor residue × interaction type"
            ),
            "stacking_note": (
                "Interaction channels are not mutually exclusive; stacked "
                "values may exceed 100%."
            ),
        },
    )


def _fingerprint_long(
    matrix: pd.DataFrame,
    labels: Mapping[str, str] | None = None,
) -> pd.DataFrame:
    rows = []
    observation_label = "observation_id"
    for observation_id, values in matrix.iterrows():
        for (residue, kind), present in values.items():
            rows.append(
                {
                    observation_label: observation_id,
                    "observation_label": (
                        labels.get(str(observation_id), str(observation_id))
                        if labels is not None
                        else str(observation_id)
                    ),
                    "receptor_residue": residue,
                    "interaction_type": kind,
                    "present": bool(present),
                }
            )
    return pd.DataFrame(
        rows,
        columns=[
            observation_label,
            "observation_label",
            "receptor_residue",
            "interaction_type",
            "present",
        ],
    )


def build_fingerprint_chart(
    result: RunResult,
    *,
    mode: str = "docking",
    matrix: pd.DataFrame | None = None,
    observation_labels: Mapping[str, str] | None = None,
) -> ChartArtifact:
    """Build an observation × interaction-feature barcode."""
    if mode not in {"docking", "md"}:
        raise ValueError("mode must be 'docking' or 'md'")
    matrix = (
        matrix.copy(deep=True)
        if matrix is not None
        else fingerprint_matrix(result)
    )
    data = _fingerprint_long(matrix, observation_labels)
    figure = _figure(width=10.8, height=5.7)
    axis = figure.add_subplot(111)
    _style_axis(axis)
    if matrix.empty or matrix.shape[1] == 0:
        _empty_axis(axis, "No fingerprint features for this selection")
    else:
        values = matrix.to_numpy(dtype=bool)
        pixels = np.empty((values.shape[0], values.shape[1], 4), dtype=float)
        pixels[:] = to_rgba("#EEF2F4")
        for column_index, (_residue, kind) in enumerate(matrix.columns):
            pixels[values[:, column_index], column_index] = to_rgba(
                color_hex(kind)
            )
        axis.imshow(
            pixels,
            aspect="auto",
            interpolation="nearest",
        )
        feature_labels = [
            f"{residue}\n{kind}" for residue, kind in matrix.columns
        ]
        axis.set_xticks(np.arange(len(feature_labels)), labels=feature_labels)
        axis.tick_params(axis="x", labelrotation=90, labelsize=7)
        display_labels = [
            (
                observation_labels.get(str(value), str(value))
                if observation_labels is not None
                else str(value)
            )
            for value in matrix.index
        ]
        axis.set_yticks(
            np.arange(len(matrix.index)),
            labels=display_labels,
        )
        axis.grid(False)
    axis.set_xlabel("Receptor residue × interaction type")
    axis.set_ylabel("Pose" if mode == "docking" else "Saved frame")
    axis.set_title("Interaction fingerprint barcode", loc="left")
    return ChartArtifact(
        kind="interaction-fingerprint",
        figure=figure,
        data=data,
        metadata={
            "mode": mode,
            "total_observations": len(matrix.index),
            "total_features": len(matrix.columns),
            "counting_unit": "binary presence per observation × feature",
        },
    )


def build_similarity_chart(
    result: RunResult,
    *,
    matrix: pd.DataFrame | None = None,
    similarity: pd.DataFrame | None = None,
    observation_labels: Mapping[str, str] | None = None,
) -> ChartArtifact:
    """Build a Tanimoto/Jaccard similarity heatmap for observations."""
    matrix = (
        matrix.copy(deep=True)
        if matrix is not None
        else fingerprint_matrix(result)
    )
    similarity = (
        similarity.copy(deep=True)
        if similarity is not None
        else fingerprint_similarity(matrix)
    )
    rows = [
        {
            "observation_a": left,
            "observation_b": right,
            "observation_a_label": (
                observation_labels.get(str(left), str(left))
                if observation_labels is not None
                else str(left)
            ),
            "observation_b_label": (
                observation_labels.get(str(right), str(right))
                if observation_labels is not None
                else str(right)
            ),
            "tanimoto_similarity": float(similarity.loc[left, right]),
        }
        for left in similarity.index
        for right in similarity.columns
    ]
    data = pd.DataFrame(
        rows,
        columns=[
            "observation_a",
            "observation_b",
            "observation_a_label",
            "observation_b_label",
            "tanimoto_similarity",
        ],
    )
    figure = _figure(width=7.2, height=6.2)
    axis = figure.add_subplot(111)
    _style_axis(axis)
    if similarity.empty:
        _empty_axis(axis, "No observations available for similarity analysis")
    else:
        image = axis.imshow(
            similarity.to_numpy(dtype=float),
            cmap="Blues",
            vmin=0,
            vmax=1,
            aspect="equal",
            interpolation="nearest",
        )
        labels = [
            (
                observation_labels.get(str(value), str(value))
                if observation_labels is not None
                else str(value)
            )
            for value in similarity.index
        ]
        axis.set_xticks(np.arange(len(labels)), labels=labels)
        axis.set_yticks(np.arange(len(labels)), labels=labels)
        axis.tick_params(axis="x", labelrotation=90, labelsize=7)
        axis.tick_params(axis="y", labelsize=7)
        axis.grid(False)
        figure.colorbar(image, ax=axis, label="Tanimoto similarity")
    axis.set_xlabel("Observation")
    axis.set_ylabel("Observation")
    axis.set_title("Fingerprint similarity", loc="left")
    return ChartArtifact(
        kind="fingerprint-similarity",
        figure=figure,
        data=data,
        metadata={
            "metric": "Jaccard/Tanimoto",
            "empty_fingerprint_convention": "empty/empty = 1",
            "total_observations": len(similarity.index),
        },
    )


def build_interaction_heatmap_chart(
    result: RunResult,
    *,
    group_by: str = "source",
    feature_level: str = "residue_type",
    label_mode: str = "ligand",
    mode: str = "docking",
    series=None,
    top_n: int | None = 40,
) -> ChartArtifact:
    """Build a normalized ligand/file or observation interaction heatmap."""

    heatmap = build_interaction_heatmap_data(
        result,
        group_by=group_by,
        feature_level=feature_level,
        label_mode=label_mode,
        mode=mode,
        series=series,
        top_n=top_n,
    )
    height = min(12.0, max(4.8, 0.26 * len(heatmap.matrix.index) + 2.8))
    figure = _figure(width=11.2, height=height)
    axis = figure.add_subplot(111)
    _style_axis(axis)
    if heatmap.matrix.empty or heatmap.matrix.shape[1] == 0:
        _empty_axis(axis, "No interaction features for this heatmap")
    else:
        image = axis.imshow(
            heatmap.matrix.to_numpy(dtype=float),
            cmap="YlGnBu",
            vmin=0,
            vmax=100,
            aspect="auto",
            interpolation="nearest",
        )
        feature_labels = [_feature_label(value) for value in heatmap.matrix.columns]
        axis.set_xticks(
            np.arange(len(feature_labels)),
            labels=feature_labels,
        )
        axis.tick_params(axis="x", labelrotation=90, labelsize=7)
        row_labels = heatmap.matrix.index.tolist()
        tick_step = max(1, int(np.ceil(len(row_labels) / 60)))
        tick_positions = np.arange(0, len(row_labels), tick_step)
        axis.set_yticks(
            tick_positions,
            labels=[row_labels[index] for index in tick_positions],
        )
        axis.tick_params(axis="y", labelsize=8)
        axis.grid(False)
        figure.colorbar(
            image,
            ax=axis,
            label="Interaction frequency / occupancy (%)",
        )
    axis.set_xlabel(
        "Receptor residue × interaction type"
        if feature_level == "residue_type"
        else "Receptor residue (any interaction)"
    )
    axis.set_ylabel(
        "Ligand / uploaded file"
        if group_by == "source"
        else ("Pose" if mode == "docking" else "Saved frame")
    )
    axis.set_title("Interaction comparison heatmap", loc="left")
    return ChartArtifact(
        kind="interaction-comparison-heatmap",
        figure=figure,
        data=heatmap.cell_data,
        metadata=heatmap.metadata,
    )


def _feature_label(feature) -> str:
    if isinstance(feature, tuple):
        return "\n".join(str(value) for value in feature)
    return str(feature)


def build_comparison_chart(
    system_a: RunResult, system_b: RunResult, *, mode: str = "docking"
) -> ChartArtifact:
    """Build signed B-minus-A prevalence/occupancy bars."""
    if mode not in {"docking", "md"}:
        raise ValueError("mode must be 'docking' or 'md'")
    data = differential_prevalence(system_a, system_b)
    figure = _figure(width=10.0, height=5.8)
    axis = figure.add_subplot(111)
    _style_axis(axis)
    if data.empty:
        _empty_axis(axis, "Load two systems to compare interaction evidence")
    else:
        ordered = data.assign(
            magnitude=data["delta_pct_points"].abs()
        ).sort_values(["magnitude", "receptor_residue"], ascending=[False, True])
        ordered = ordered.head(24).iloc[::-1]
        labels = [
            f"{row.receptor_residue} · {row.interaction_type}"
            for row in ordered.itertuples()
        ]
        values = ordered["delta_pct_points"].to_numpy(dtype=float)
        colors = np.where(values >= 0, "#0072B2", "#D55E00")
        axis.barh(
            np.arange(len(values)),
            values,
            color=colors,
            height=0.68,
            edgecolor=_PAPER,
            linewidth=0.5,
        )
        axis.axvline(0, color=_INK, linewidth=0.9)
        axis.set_yticks(np.arange(len(values)), labels=labels)
    measure = "prevalence" if mode == "docking" else "occupancy"
    axis.set_xlabel(f"Δ {measure} (B − A, percentage points)")
    axis.set_title("Differential interaction evidence", loc="left")
    return ChartArtifact(
        kind="differential-prevalence",
        figure=figure,
        data=data,
        metadata={
            "mode": mode,
            "system_a_observations": len(observation_ids(system_a)),
            "system_b_observations": len(observation_ids(system_b)),
            "delta_definition": "B − A",
            "counting_unit": (
                "observation × receptor residue × interaction type"
            ),
        },
    )


def build_retention_chart(
    docking: RunResult,
    md: RunResult,
    *,
    retained_threshold_pct: float = 50.0,
) -> ChartArtifact:
    """Build docking-to-MD retention category counts."""
    data = docking_md_retention(
        docking, md, retained_threshold_pct=retained_threshold_pct
    )
    categories = ("retained", "intermittent", "lost", "gained")
    colours = {
        "retained": "#009E73",
        "intermittent": "#CC79A7",
        "lost": "#D55E00",
        "gained": "#0072B2",
    }
    counts = {
        category: int((data["category"] == category).sum())
        for category in categories
    }
    figure = _figure(width=8.4, height=4.8)
    axis = figure.add_subplot(111)
    _style_axis(axis)
    if data.empty:
        _empty_axis(axis, "Load docking and MD evidence to assess retention")
    else:
        axis.barh(
            np.arange(len(categories)),
            [counts[category] for category in categories],
            color=[colours[category] for category in categories],
            height=0.62,
        )
        axis.set_yticks(
            np.arange(len(categories)),
            labels=[category.capitalize() for category in categories],
        )
        axis.invert_yaxis()
        for index, category in enumerate(categories):
            axis.text(
                counts[category] + 0.05,
                index,
                str(counts[category]),
                va="center",
                color=_INK,
                fontsize=9,
            )
    axis.set_xlabel("Residue × interaction-type features")
    axis.set_title("Docking → MD retention", loc="left")
    return ChartArtifact(
        kind="docking-md-retention",
        figure=figure,
        data=data,
        metadata={
            "retained_threshold_pct": retained_threshold_pct,
            "docking_observations": len(observation_ids(docking)),
            "md_observations": len(observation_ids(md)),
            "category_definition": (
                "retained ≥ threshold; intermittent >0 and < threshold; "
                "lost = 0; gained absent in docking"
            ),
        },
    )


# Local import avoids obscuring the chart-building dependencies above.
from .analytics import observation_ids  # noqa: E402


__all__ = [
    "ChartArtifact",
    "build_comparison_chart",
    "build_fingerprint_chart",
    "build_interaction_heatmap_chart",
    "build_retention_chart",
    "build_residue_chart",
    "build_similarity_chart",
]
