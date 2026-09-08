"""Deterministic ligand ranking for charts only; scientific tables stay complete."""

from __future__ import annotations

from dataclasses import dataclass, replace
from pathlib import Path
from types import MappingProxyType
from typing import Mapping

from .results import RunResult


MAX_CHART_LIGANDS = 100
MAX_CHART_OBSERVATIONS = 100
RANKING_CRITERIA = {
    "interaction_count": "Interaction count",
    "fingerprint_features": "Fingerprint features",
}


def ligand_identity(item) -> tuple[str, str]:
    """Use canonical compound ID, never the raw vendor name of each pose."""
    source = str(item.source_id or item.source_path or item.source_file or "unidentified-source")
    ligand = str(item.ligand_id or "").strip()
    if not ligand:
        ligand = str(item.source_molecule_name or Path(item.source_file).stem or "Ligand")
    return source, ligand


def subset_chart_observations(result: RunResult, ids) -> RunResult:
    allowed = frozenset(ids)
    return replace(
        result,
        summaries=tuple(item for item in result.summaries if item.pose_id in allowed),
        details=tuple(item for item in result.details if item.pose_id in allowed),
        pending=tuple(item for item in result.pending if getattr(item, "pose_id", "") in allowed),
    )


@dataclass(frozen=True)
class ChartSelection:
    result: RunResult
    metadata: Mapping[str, object]
    display_observation_ids: tuple[str, ...]

    @property
    def note(self) -> str:
        criterion = ("mean interactions/pose" if self.metadata["ranking_criterion"] == "interaction_count"
                     else "mean distinct fingerprint features/pose")
        return (f"Top {self.metadata['ligands_displayed']} of {self.metadata['ligands_before_limit']} "
                f"ligands by {criterion} (limit 100). "
                "Interaction evidence ranking, not binding affinity. Tables retain all ligands.")

    def observation_result(self) -> RunResult:
        return subset_chart_observations(self.result, self.display_observation_ids)


def select_chart_ligands(result: RunResult, criterion: str = "interaction_count") -> ChartSelection:
    """Rank mean per observation, preserving zero-contact poses and source identity."""
    if criterion not in RANKING_CRITERIA:
        raise ValueError("unknown chart ranking criterion")
    observations = {item.pose_id: item for item in result.summaries}
    features: dict[str, set[tuple[str, str]]] = {}
    detail_counts: dict[str, int] = {}
    for item in result.details:
        observations.setdefault(item.pose_id, item)
        detail_counts[item.pose_id] = detail_counts.get(item.pose_id, 0) + 1
        if criterion == "fingerprint_features":
            features.setdefault(item.pose_id, set()).add((item.receptor_residue, item.interaction_type))
    groups: dict[tuple[str, str], list[str]] = {}
    scores: dict[str, float] = {}
    for pose_id, item in observations.items():
        groups.setdefault(ligand_identity(item), []).append(pose_id)
        scores[pose_id] = (float(len(features.get(pose_id, ()))) if criterion == "fingerprint_features"
                          else float(getattr(item, "n_total_interactions", detail_counts.get(pose_id, 0))))
    ranked = sorted(groups, key=lambda key: (-sum(scores[p] for p in groups[key]) / len(groups[key]), key))
    selected_keys = ranked[:MAX_CHART_LIGANDS]
    allowed = {pose for key in selected_keys for pose in groups[key]}
    ordered = tuple(pose for pose in observations if pose in allowed)
    # Cover every ranked ligand before assigning remaining display slots to poses.
    representatives = {min(groups[key], key=lambda pose: (-scores[pose], pose)) for key in selected_keys}
    display = set(representatives)
    for pose in ordered:
        if len(display) >= MAX_CHART_OBSERVATIONS:
            break
        display.add(pose)
    shown = tuple(pose for pose in ordered if pose in display)
    ranking = tuple({"source_id": key[0], "ligand_id": key[1],
                     "mean_per_observation": sum(scores[p] for p in groups[key]) / len(groups[key]),
                     "observations": len(groups[key])} for key in selected_keys)
    return ChartSelection(
        subset_chart_observations(result, allowed),
        MappingProxyType({
            "ranking_criterion": criterion,
            "ranking_aggregation": "arithmetic mean per pose or saved frame, including zero-contact observations",
            "ranking_identity": "source identity + canonical ligand_id; raw molecule name only as fallback",
            "ranking_tie_break": "source identity then canonical ligand ID, ascending",
            "ligand_limit": MAX_CHART_LIGANDS,
            "ligands_before_limit": len(groups), "ligands_displayed": len(selected_keys),
            "ligand_ranking": ranking,
            "selected_observations": len(ordered),
            "display_observation_limit": MAX_CHART_OBSERVATIONS,
            "display_observation_ids": shown,
            "observation_display_policy": "highest-ranked pose per selected ligand, then remaining slots in original observation order; ties by observation ID",
            "affinity_inference": False,
        }),
        shown,
    )
