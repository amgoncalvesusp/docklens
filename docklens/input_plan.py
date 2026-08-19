"""Immutable input plans shared by the UI, runner and project persistence.

DockLens deliberately reduces all user-facing loading workflows to two input
kinds: ``combined`` files containing receptor and ligand material, and
``paired`` jobs where ``receptor_path`` is shared by one ligand/pose source.
The runner can therefore remain generic while the UI is free to organize
files into receptor groups.
"""

from __future__ import annotations

from dataclasses import dataclass
import os
from pathlib import Path


SUPPORTED_SUFFIXES = (".mol2", ".pdb", ".pdbqt")
VALID_KINDS = frozenset({"paired", "combined"})


def _absolute(value: str | os.PathLike[str]) -> str:
    return os.path.abspath(os.fspath(value))


@dataclass(frozen=True)
class InputJob:
    """One source file and its optional external receptor provenance."""

    kind: str
    source_path: str
    receptor_path: str | None = None
    group_id: str = ""
    resolution_method: str = ""
    input_status: str = "success"
    input_code: str = ""
    input_message: str = ""

    def __post_init__(self) -> None:
        kind = str(self.kind).strip().lower()
        if kind not in VALID_KINDS:
            raise ValueError("kind must be 'paired' or 'combined'")
        source_path = _absolute(self.source_path)
        if kind == "paired":
            if not self.receptor_path:
                raise ValueError("paired jobs require receptor_path")
            receptor_path = _absolute(self.receptor_path)
        else:
            if self.receptor_path is not None:
                raise ValueError("combined jobs cannot define receptor_path")
            receptor_path = None
        group_id = str(self.group_id or "").strip()
        object.__setattr__(self, "kind", kind)
        object.__setattr__(self, "source_path", source_path)
        object.__setattr__(self, "receptor_path", receptor_path)
        object.__setattr__(self, "group_id", group_id)
        object.__setattr__(self, "resolution_method", str(self.resolution_method or ""))
        object.__setattr__(self, "input_status", str(self.input_status or "success"))
        object.__setattr__(self, "input_code", str(self.input_code or ""))
        object.__setattr__(self, "input_message", str(self.input_message or ""))

    @property
    def receptor_source_file(self) -> str:
        return os.path.basename(self.receptor_path) if self.receptor_path else ""


@dataclass(frozen=True)
class InputPlan:
    """Ordered, deduplicated collection of jobs to analyze."""

    jobs: tuple[InputJob, ...] = ()

    def __post_init__(self) -> None:
        normalized = []
        seen = set()
        for job in self.jobs:
            if not isinstance(job, InputJob):
                raise TypeError("InputPlan jobs must be InputJob instances")
            identity = (
                job.kind,
                os.path.normcase(job.source_path),
                os.path.normcase(job.receptor_path or ""),
                job.group_id,
            )
            if identity in seen:
                continue
            seen.add(identity)
            normalized.append(job)
        object.__setattr__(self, "jobs", tuple(normalized))

    def paths(self) -> tuple[str, ...]:
        """Return all declared source and receptor paths without duplicates."""

        values = []
        seen = set()
        for job in self.jobs:
            for value in (job.receptor_path, job.source_path):
                if value is None:
                    continue
                key = os.path.normcase(value)
                if key not in seen:
                    seen.add(key)
                    values.append(value)
        return tuple(values)


def plan_from_paths(paths) -> InputPlan:
    """Create a combined plan from files or recursively scanned folders."""

    if isinstance(paths, (str, os.PathLike)):
        paths = [paths]
    if paths is None:
        raise ValueError("paths must contain at least one input")
    expanded = []
    for value in paths:
        path = Path(value).expanduser()
        if path.is_dir():
            expanded.extend(
                item
                for item in sorted(path.rglob("*"), key=lambda item: str(item).casefold())
                if item.is_file() and item.suffix.lower() in SUPPORTED_SUFFIXES
            )
        else:
            expanded.append(path)
    if not expanded:
        raise ValueError("paths must contain at least one input")
    jobs = tuple(
        InputJob(
            kind="combined",
            source_path=str(path),
            group_id="combined_%d" % index,
        )
        for index, path in enumerate(expanded, 1)
    )
    return InputPlan(jobs)


def run_shared_receptor_plan(
    receptor_path: str | os.PathLike[str], ligand_paths
) -> InputPlan:
    """Build paired jobs that share one receptor without duplicating it."""

    if isinstance(ligand_paths, (str, os.PathLike)):
        ligand_paths = [ligand_paths]
    if ligand_paths is None:
        raise ValueError("ligand_paths must contain at least one input")
    receptor = _absolute(receptor_path)
    jobs = tuple(
        InputJob(
            kind="paired",
            source_path=os.fspath(path),
            receptor_path=receptor,
            group_id="receptor_1",
        )
        for path in ligand_paths
    )
    if not jobs:
        raise ValueError("ligand_paths must contain at least one input")
    return InputPlan(jobs)


__all__ = [
    "InputJob",
    "InputPlan",
    "SUPPORTED_SUFFIXES",
    "VALID_KINDS",
    "plan_from_paths",
    "run_shared_receptor_plan",
]
