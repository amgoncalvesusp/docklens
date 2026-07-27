"""Immutable result contracts shared by the runner, UI and exporters."""

from __future__ import annotations

from dataclasses import dataclass, field, replace
from types import MappingProxyType
from typing import Mapping, Optional, Tuple

from .analysis_profiles import VALID_ANALYSIS_PROFILES, normalize_analysis_profile
from .residue_keys import normalize_key_residues


@dataclass(frozen=True)
class Endpoint:
    side: str
    kind: str
    atom_name: str
    atom_serials: Tuple[int, ...]
    resname: str
    resseq: str
    chain: str = ""
    element: str = ""
    role: str = ""

    @property
    def residue_id(self) -> str:
        return "%s%s%s" % (self.resname, self.resseq, self.chain)


@dataclass(frozen=True)
class Detail:
    ligand_id: str
    source_file: str
    interaction_type: str
    subtype: str
    ligand: Endpoint
    receptor: Endpoint
    distance_A: Optional[float]
    source_id: str
    pose_id: str
    interaction_id: str
    pose: int
    sol: Optional[int] = None
    docking_score: Optional[float] = None
    source_path: str = ""
    resolution_method: str = ""
    is_key_residue: bool = False
    water: Optional[Endpoint] = None
    receptor_water_distance_A: Optional[float] = None
    ligand_water_distance_A: Optional[float] = None
    water_angle_deg: Optional[float] = None
    chemistry_basis: str = ""
    chemistry_confidence: str = ""
    hydrogen_atom: str = ""
    hydrogen_atom_serial: Optional[int] = None
    hydrogen_acceptor_distance_A: Optional[float] = None
    donor_hydrogen_acceptor_angle_deg: Optional[float] = None
    hydrogen_acceptor_base_angle_deg: Optional[float] = None
    theta_deg: Optional[float] = None

    @property
    def ligand_atom(self) -> str:
        return self.ligand.atom_name

    @property
    def receptor_residue(self) -> str:
        return self.receptor.residue_id

    @property
    def receptor_atom(self) -> str:
        return self.receptor.atom_name

    @property
    def distance(self) -> Optional[float]:
        return self.distance_A

    @property
    def _res_nochain(self) -> str:
        return "%s%s" % (self.receptor.resname, self.receptor.resseq)


@dataclass(frozen=True)
class Summary:
    ligand_id: str
    source_file: str
    sol: Optional[int]
    pose: int
    docking_score: Optional[float]
    n_total_interactions: int
    n_key_residue_interactions: int
    counts: Mapping[str, int]
    source_id: str = ""
    pose_id: str = ""
    source_path: str = ""
    resolution_method: str = ""

    def __post_init__(self):
        object.__setattr__(self, "counts", MappingProxyType(dict(self.counts)))


@dataclass(frozen=True)
class InputQC:
    source_id: str
    source_file: str
    source_path: str = ""
    pose_id: str = ""
    status: str = "success"
    code: str = ""
    message: str = ""
    format: str = ""
    poses_found: int = 0
    poses_processed: int = 0
    resolution_method: str = ""
    receptor_atoms: int = 0
    ligand_atoms: int = 0
    water_atoms: int = 0
    warnings: Tuple[str, ...] = ()


@dataclass(frozen=True)
class AnalysisParameters:
    schema_version: str = "3"
    app_version: str = ""
    started_at: str = ""
    hbond_preset: str = "plip"
    cutoffs: Tuple[Tuple[str, float], ...] = ()
    interaction_types: Tuple[str, ...] = ()
    key_residues: Tuple[str, ...] = ()
    counting_unit: str = "semantic_interaction"


@dataclass(frozen=True)
class ExportFilter:
    scope: str = "all"
    interaction_types: Optional[frozenset] = None
    text: str = ""
    key_only: bool = False
    matrix_mode: str = "count"
    analysis_profile: str = "complete"

    def __post_init__(self):
        if self.scope not in {"all", "filtered"}:
            raise ValueError("scope must be 'all' or 'filtered'")
        if self.matrix_mode not in {"count", "presence"}:
            raise ValueError("matrix_mode must be 'count' or 'presence'")
        if self.analysis_profile not in VALID_ANALYSIS_PROFILES:
            raise ValueError(
                "analysis_profile must be one of %s"
                % ", ".join(VALID_ANALYSIS_PROFILES)
            )
        object.__setattr__(
            self, "analysis_profile", normalize_analysis_profile(self.analysis_profile)
        )


@dataclass(frozen=True)
class RunResult:
    details: Tuple[Detail, ...] = ()
    summaries: Tuple[Summary, ...] = ()
    pending: tuple = ()
    key_residues: frozenset = frozenset()
    receptor_residues: frozenset = frozenset()
    input_qc: Tuple[InputQC, ...] = ()
    parameters: AnalysisParameters = field(default_factory=AnalysisParameters)


def make_result(
    *,
    details=(),
    summaries=(),
    pending=(),
    key_residues=(),
    receptor_residues=(),
    input_qc=(),
    parameters=None,
) -> RunResult:
    return RunResult(
        details=tuple(details),
        summaries=tuple(summaries),
        pending=tuple(pending),
        key_residues=frozenset(key_residues),
        receptor_residues=frozenset(receptor_residues),
        input_qc=tuple(input_qc),
        parameters=parameters or AnalysisParameters(),
    )


def _normalise_keys(items) -> frozenset:
    return normalize_key_residues(items)


def with_key_residues(result: RunResult, key_residues) -> RunResult:
    keys = _normalise_keys(key_residues)
    details = tuple(
        replace(
            detail,
            is_key_residue=(
                detail.receptor_residue.upper() in keys
                or detail._res_nochain.upper() in keys
            ),
        )
        for detail in result.details
    )
    per_pose = {}
    for detail in details:
        if detail.is_key_residue:
            per_pose[detail.pose_id] = per_pose.get(detail.pose_id, 0) + 1
    summaries = tuple(
        replace(
            summary,
            n_key_residue_interactions=per_pose.get(summary.pose_id, 0),
        )
        for summary in result.summaries
    )
    parameters = replace(result.parameters, key_residues=tuple(sorted(keys)))
    return replace(
        result,
        details=details,
        summaries=summaries,
        key_residues=keys,
        parameters=parameters,
    )
