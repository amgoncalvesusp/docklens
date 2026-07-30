"""Figures for pose families, MD interaction states and uncertainty."""

from __future__ import annotations

import numpy as np
import pandas as pd
from matplotlib.colors import ListedColormap

from .dynamic_states import (
    InteractionStateAnalysis,
    state_assignment_frame,
    state_summary_frame,
    state_transition_frame,
)
from .plotting import ChartArtifact, _empty_axis, _figure, _style_axis


_STATE_COLOURS = (
    "#1B4965",
    "#5FA8D3",
    "#7A5195",
    "#BC5090",
    "#4C956C",
    "#B07D62",
    "#577590",
    "#F28E2B",
    "#8A817C",
    "#648FFF",
)
_OUTLIER_COLOUR = "#AEBBC4"
_INK = "#12202F"


def _state_colour_map(state_ids):
    identifiers = tuple(dict.fromkeys(state_ids))
    return {
        state_id: (
            _OUTLIER_COLOUR
            if state_id == "OUTLIER"
            else _STATE_COLOURS[index % len(_STATE_COLOURS)]
        )
        for index, state_id in enumerate(identifiers)
    }


def build_state_timeline_chart(
    analysis: InteractionStateAnalysis,
) -> ChartArtifact:
    """Build an ordered MD-state timeline or docking-family ribbon."""
    data = state_assignment_frame(analysis)
    figure = _figure(width=10.4, height=3.4)
    axis = figure.add_subplot(111)
    _style_axis(axis)
    if data.empty:
        _empty_axis(axis, "No observations available for state analysis")
    else:
        states = list(dict.fromkeys(data["state_id"].tolist()))
        state_index = {state: index for index, state in enumerate(states)}
        colours = _state_colour_map(states)
        values = np.array(
            [[state_index[state] for state in data["state_id"]]], dtype=float
        )
        cmap = ListedColormap([colours[state] for state in states])
        if analysis.mode == "md":
            x_values = data["time_ns"].to_numpy(dtype=float)
            step = analysis.series.time_step_ns or analysis.context.time_step_ns
            extent = (
                float(x_values[0]),
                float(x_values[-1] + step),
                0,
                1,
            )
            axis.imshow(
                values,
                aspect="auto",
                interpolation="nearest",
                cmap=cmap,
                vmin=-0.5,
                vmax=max(0.5, len(states) - 0.5),
                extent=extent,
            )
            axis.set_xlabel("Simulation time (ns)")
            title = "Interaction-state timeline"
        else:
            axis.imshow(
                values,
                aspect="auto",
                interpolation="nearest",
                cmap=cmap,
                vmin=-0.5,
                vmax=max(0.5, len(states) - 0.5),
                extent=(0, len(data), 0, 1),
            )
            axis.set_xlabel("Docking pose order (not time)")
            title = "Pose-family membership"
        axis.set_yticks([])
        axis.grid(False)
        handles = [
            axis.scatter([], [], marker="s", color=colours[state], label=state)
            for state in states
        ]
        axis.legend(
            handles=handles,
            loc="upper center",
            bbox_to_anchor=(0.5, -0.24),
            ncols=min(6, len(states)),
            frameon=False,
        )
        axis.set_title(title, loc="left")
    return ChartArtifact(
        kind="interaction-state-timeline",
        figure=figure,
        data=data,
        metadata={
            "mode": analysis.mode,
            "time_step_ns": (
                analysis.series.time_step_ns or analysis.context.time_step_ns
            ),
            "clustering_method": analysis.method,
            "similarity_threshold": analysis.threshold,
            "sampled": analysis.sampled,
            "training_observations": analysis.training_observations,
            "total_observations": analysis.total_observations,
            "outlier_observations": len(analysis.outlier_observations),
        },
    )


def build_state_population_chart(
    analysis: InteractionStateAnalysis,
) -> ChartArtifact:
    """Build population bars with modality-safe pose-family/state wording."""
    data = state_summary_frame(analysis)
    figure = _figure(width=8.4, height=4.8)
    axis = figure.add_subplot(111)
    _style_axis(axis)
    if data.empty:
        _empty_axis(axis, "No pose families or interaction states available")
    else:
        colours = _state_colour_map(data["state_id"].tolist())
        positions = np.arange(len(data))
        axis.barh(
            positions,
            data["population_pct"].to_numpy(dtype=float),
            color=[colours[state] for state in data["state_id"]],
            height=0.64,
        )
        labels = [
            f"{row.state_id} · representative {row.representative}"
            for row in data.itertuples()
        ]
        axis.set_yticks(positions, labels=labels)
        axis.invert_yaxis()
    if analysis.mode == "docking":
        title = "Pose-family population"
        label = "Docking frequency (% of poses)"
        kind = "pose-family-population"
    else:
        title = "Interaction-state population"
        label = "MD occupancy (% of saved frames)"
        kind = "interaction-state-population"
    axis.set_title(title, loc="left")
    axis.set_xlabel(label)
    return ChartArtifact(
        kind=kind,
        figure=figure,
        data=data,
        metadata={
            "mode": analysis.mode,
            "clustering_method": analysis.method,
            "similarity_threshold": analysis.threshold,
            "population_denominator": analysis.total_observations,
        },
    )


def build_transition_chart(
    analysis: InteractionStateAnalysis,
    *,
    lag: int = 1,
) -> ChartArtifact:
    """Build a descriptive observed-transition probability heatmap."""
    data = state_transition_frame(analysis, lag=lag)
    figure = _figure(width=7.2, height=5.8)
    axis = figure.add_subplot(111)
    _style_axis(axis)
    states = [
        state.state_id for state in analysis.states
    ]
    if data.empty:
        _empty_axis(axis, "No observed transitions at the selected lag")
    else:
        states = list(
            dict.fromkeys(
                states
                + data["from_state"].tolist()
                + data["to_state"].tolist()
            )
        )
        matrix = data.pivot_table(
            index="from_state",
            columns="to_state",
            values="transition_probability_pct",
            aggfunc="sum",
        ).reindex(index=states, columns=states)
        image = axis.imshow(
            matrix.to_numpy(dtype=float),
            cmap="Blues",
            vmin=0,
            vmax=100,
            interpolation="nearest",
        )
        axis.set_xticks(np.arange(len(states)), labels=states)
        axis.set_yticks(np.arange(len(states)), labels=states)
        axis.grid(False)
        figure.colorbar(image, ax=axis, label="Observed probability (%)")
        for row_index, source in enumerate(states):
            for column_index, target in enumerate(states):
                matches = data[
                    (data["from_state"] == source)
                    & (data["to_state"] == target)
                ]
                if not matches.empty:
                    row = matches.iloc[0]
                    axis.text(
                        column_index,
                        row_index,
                        f"{int(row.transition_count)}",
                        ha="center",
                        va="center",
                        color=_INK,
                        fontsize=8,
                    )
    axis.set_xlabel("To state")
    axis.set_ylabel("From state")
    axis.set_title("Observed state transitions", loc="left")
    return ChartArtifact(
        kind="observed-state-transitions",
        figure=figure,
        data=data,
        metadata={
            "mode": analysis.mode,
            "lag_observations": lag,
            "interpretation": (
                "descriptive observed transitions; not a validated Markov model"
            ),
            "replica_boundaries_preserved": True,
            "frame_gaps_excluded": True,
        },
    )


def build_uncertainty_chart(
    intervals: pd.DataFrame,
    *,
    top_n: int = 24,
) -> ChartArtifact:
    """Build saved-frame occupancy estimates with block-bootstrap intervals."""
    if top_n < 1:
        raise ValueError("top_n must be at least one")
    data = intervals.copy(deep=True)
    figure = _figure(width=9.2, height=5.8)
    axis = figure.add_subplot(111)
    _style_axis(axis)
    if data.empty:
        _empty_axis(axis, "Compute MD confidence intervals to view uncertainty")
    else:
        ordered = data.sort_values(
            ["occupancy_pct", "receptor_residue", "interaction_type"],
            ascending=[False, True, True],
        ).head(top_n).iloc[::-1]
        estimates = ordered["occupancy_pct"].to_numpy(dtype=float)
        lower = estimates - ordered["ci_low_pct"].to_numpy(dtype=float)
        upper = ordered["ci_high_pct"].to_numpy(dtype=float) - estimates
        positions = np.arange(len(ordered))
        labels = [
            f"{row.receptor_residue} · {row.interaction_type}"
            for row in ordered.itertuples()
        ]
        axis.errorbar(
            estimates,
            positions,
            xerr=np.vstack((lower, upper)),
            fmt="o",
            color="#0072B2",
            ecolor="#5B6976",
            capsize=3,
            markersize=5,
        )
        axis.set_yticks(positions, labels=labels)
        axis.set_xlim(0, 100)
    axis.set_xlabel("MD occupancy with confidence interval (%)")
    axis.set_title("Interaction occupancy uncertainty", loc="left")
    first = data.iloc[0] if not data.empty else None
    return ChartArtifact(
        kind="md-occupancy-confidence-intervals",
        figure=figure,
        data=data,
        metadata={
            "bootstrap_method": (
                first["method"]
                if first is not None
                else "circular moving-block bootstrap"
            ),
            "iterations": (
                int(first["iterations"]) if first is not None else None
            ),
            "block_size": (
                int(first["block_size"]) if first is not None else None
            ),
            "confidence_level": (
                float(first["confidence_level_pct"]) / 100.0
                if first is not None
                else None
            ),
            "multiplicity_adjustment": "none",
        },
    )


def build_difference_uncertainty_chart(
    intervals: pd.DataFrame,
    *,
    top_n: int = 24,
) -> ChartArtifact:
    """Build independent-system B-minus-A occupancy intervals."""
    if top_n < 1:
        raise ValueError("top_n must be at least one")
    data = intervals.copy(deep=True)
    figure = _figure(width=9.4, height=5.8)
    axis = figure.add_subplot(111)
    _style_axis(axis)
    if data.empty:
        _empty_axis(
            axis,
            "Compute two-system MD intervals to compare uncertainty",
        )
    else:
        ordered = data.assign(
            magnitude=data["delta_pct_points"].abs()
        ).sort_values(
            ["magnitude", "receptor_residue", "interaction_type"],
            ascending=[False, True, True],
        ).head(top_n).iloc[::-1]
        estimates = ordered["delta_pct_points"].to_numpy(dtype=float)
        lower = estimates - ordered["ci_low_pct_points"].to_numpy(dtype=float)
        upper = (
            ordered["ci_high_pct_points"].to_numpy(dtype=float) - estimates
        )
        positions = np.arange(len(ordered))
        labels = [
            f"{row.receptor_residue} · {row.interaction_type}"
            for row in ordered.itertuples()
        ]
        axis.errorbar(
            estimates,
            positions,
            xerr=np.vstack((lower, upper)),
            fmt="o",
            color="#0072B2",
            ecolor="#5B6976",
            capsize=3,
            markersize=5,
        )
        axis.axvline(0, color=_INK, linewidth=0.9)
        axis.set_yticks(positions, labels=labels)
    axis.set_xlabel("MD occupancy difference B - A (percentage points)")
    axis.set_title("Differential occupancy uncertainty", loc="left")
    first = data.iloc[0] if not data.empty else None
    return ChartArtifact(
        kind="md-differential-occupancy-confidence-intervals",
        figure=figure,
        data=data,
        metadata={
            "bootstrap_method": (
                first["method"]
                if first is not None
                else "circular moving-block bootstrap"
            ),
            "iterations": (
                int(first["iterations"]) if first is not None else None
            ),
            "block_size_a": (
                int(first["block_size_a"]) if first is not None else None
            ),
            "block_size_b": (
                int(first["block_size_b"]) if first is not None else None
            ),
            "confidence_level": (
                float(first["confidence_level_pct"]) / 100.0
                if first is not None
                else None
            ),
            "design": "independent B - A",
            "multiplicity_adjustment": "none",
        },
    )


__all__ = [
    "build_difference_uncertainty_chart",
    "build_state_population_chart",
    "build_state_timeline_chart",
    "build_transition_chart",
    "build_uncertainty_chart",
]
