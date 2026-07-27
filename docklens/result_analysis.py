"""Pure, immutable analyses derived from DockLens interaction results.

The metrics in this module deliberately keep atomic interaction-pair counts
separate from residue coverage.  This makes rankings auditable: repeated
contacts with one residue cannot be mistaken for coverage of several key
residues.
"""

from __future__ import annotations

from dataclasses import dataclass
from types import MappingProxyType
from typing import Iterable, Mapping, Sequence, Tuple

from .residue_keys import parse_key_residues
from .results import Detail, RunResult


@dataclass(frozen=True)
class PoseKeyResidueMetrics:
    """Key-residue evidence for one pose."""

    pose_id: str
    raw_key_pair_count: int
    distinct_key_residues: Tuple[str, ...]
    configured_key_count: int
    coverage: float
    conventional_hbond_residues: Tuple[str, ...]
    interaction_types: Tuple[str, ...]

    @property
    def distinct_key_residue_count(self) -> int:
        return len(self.distinct_key_residues)

    @property
    def conventional_hbond_residue_count(self) -> int:
        return len(self.conventional_hbond_residues)

    @property
    def interaction_type_diversity(self) -> int:
        return len(self.interaction_types)


@dataclass(frozen=True)
class PoseResidueCounts:
    """Atomic interaction-pair counts for one pose and each configured key."""

    pose_id: str
    counts: Mapping[str, int]

    def __post_init__(self) -> None:
        object.__setattr__(self, "counts", MappingProxyType(dict(self.counts)))


@dataclass(frozen=True)
class KeyResidueMatrix:
    """A deterministic pose-by-key-residue count matrix."""

    residues: Tuple[str, ...]
    rows: Tuple[PoseResidueCounts, ...]


def _normalise_key(value: object) -> str:
    return "".join(str(value).strip().upper().split())


def _normalise_keys(values: Iterable[object]) -> Tuple[str, ...]:
    return parse_key_residues(values).keys


def configured_key_residues(result: RunResult) -> Tuple[str, ...]:
    """Return the authoritative, normalised key-residue denominator.

    ``RunResult.key_residues`` is authoritative after a result has been
    recomputed.  Parameter keys are a compatibility fallback for results
    loaded before that recomputation.
    """

    values: Iterable[object]
    if result.key_residues:
        values = result.key_residues
    else:
        values = result.parameters.key_residues
    return _normalise_keys(values)


def detail_key_residues(
    detail: Detail,
    key_residues: Sequence[str],
) -> Tuple[str, ...]:
    """Return configured keys matched by one interaction detail.

    A chain-qualified key (``SER70A``) matches only that residue.  A
    chainless key (``SER70``) matches the same residue on any chain.
    Returning configured identifiers, rather than inferred identifiers,
    keeps matrix columns aligned with the coverage denominator.
    """

    exact = _normalise_key(detail.receptor_residue)
    chainless = _normalise_key(detail._res_nochain)
    keys = _normalise_keys(key_residues)
    exact_matches = tuple(
        key
        for key in keys
        if key == exact
    )
    if exact_matches:
        return exact_matches
    return tuple(key for key in keys if key == chainless)


def _pose_ids(result: RunResult) -> Tuple[str, ...]:
    ordered = []
    seen = set()
    for pose_id in (
        *(summary.pose_id for summary in result.summaries),
        *(detail.pose_id for detail in result.details),
    ):
        if pose_id not in seen:
            seen.add(pose_id)
            ordered.append(pose_id)
    return tuple(ordered)


def analyze_key_residues(
    result: RunResult,
) -> Tuple[PoseKeyResidueMetrics, ...]:
    """Calculate auditable key-residue metrics for every known pose.

    Coverage is the fraction of configured keys matched by at least one
    interaction.  It is defined as zero when no keys are configured.
    ``carbon_hbond`` is intentionally excluded from the conventional H-bond
    residue metric.
    """

    keys = configured_key_residues(result)
    details_by_pose: dict[str, list[Detail]] = {}
    for detail in result.details:
        details_by_pose.setdefault(detail.pose_id, []).append(detail)

    metrics = []
    for pose_id in _pose_ids(result):
        raw_pairs = 0
        matched_residues = set()
        hbond_residues = set()
        interaction_types = set()

        for detail in details_by_pose.get(pose_id, ()):
            matched = detail_key_residues(detail, keys)
            if not matched:
                continue
            raw_pairs += 1
            matched_residues.update(matched)
            interaction_types.add(detail.interaction_type)
            if detail.interaction_type == "hbond":
                hbond_residues.update(matched)

        coverage = len(matched_residues) / len(keys) if keys else 0.0
        metrics.append(
            PoseKeyResidueMetrics(
                pose_id=pose_id,
                raw_key_pair_count=raw_pairs,
                distinct_key_residues=tuple(sorted(matched_residues)),
                configured_key_count=len(keys),
                coverage=coverage,
                conventional_hbond_residues=tuple(sorted(hbond_residues)),
                interaction_types=tuple(sorted(interaction_types)),
            )
        )

    return tuple(metrics)


def build_key_residue_matrix(result: RunResult) -> KeyResidueMatrix:
    """Build an immutable matrix of raw pairs by pose and configured key."""

    keys = configured_key_residues(result)
    counts = {
        pose_id: {key: 0 for key in keys}
        for pose_id in _pose_ids(result)
    }
    for detail in result.details:
        pose_counts = counts.setdefault(
            detail.pose_id,
            {key: 0 for key in keys},
        )
        for key in detail_key_residues(detail, keys):
            pose_counts[key] += 1

    rows = tuple(
        PoseResidueCounts(pose_id=pose_id, counts=counts[pose_id])
        for pose_id in _pose_ids(result)
    )
    return KeyResidueMatrix(residues=keys, rows=rows)
