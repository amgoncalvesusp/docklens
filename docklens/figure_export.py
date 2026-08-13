"""Atomic export of publication figures with data and provenance sidecars."""

from __future__ import annotations

import json
import os
import re
import shutil
import tempfile
from collections.abc import Callable, Mapping
from pathlib import Path

import pandas as pd

from . import __version__
from .plotting import ChartArtifact


PUBLICATION_TIFF_DPI = 600
_TIFF_COMPRESSION = "tiff_lzw"
_FORMATS = frozenset({"png", "svg", "pdf", "tiff"})
_FORMULA_PREFIXES = ("=", "+", "-", "@")
_CHART_NAME_PATTERN = re.compile(r"[A-Za-z0-9][A-Za-z0-9._-]*\Z")


def _safe_frame(frame: pd.DataFrame) -> pd.DataFrame:
    result = frame.copy(deep=True)

    def safe(value):
        if not isinstance(value, str):
            return value
        significant = value.lstrip(" \t\r\n")
        return "'" + value if significant.startswith(_FORMULA_PREFIXES) else value

    for column in result.columns:
        result[column] = result[column].map(safe)
    return result


def _atomic_figure(artifact, output, *, file_format, dpi):
    output.parent.mkdir(parents=True, exist_ok=True)
    handle = tempfile.NamedTemporaryFile(
        suffix=f".{file_format}",
        prefix=".docklens-figure-",
        dir=output.parent,
        delete=False,
    )
    temporary = Path(handle.name)
    handle.close()
    try:
        options = {
            "format": file_format,
            "dpi": dpi,
            "bbox_inches": "tight",
            "facecolor": artifact.figure.get_facecolor(),
        }
        if file_format == "tiff":
            options["pil_kwargs"] = {"compression": _TIFF_COMPRESSION}
        artifact.figure.savefig(temporary, **options)
        os.replace(temporary, output)
    except Exception:
        temporary.unlink(missing_ok=True)
        raise


def _atomic_text(output: Path, text: str):
    output.parent.mkdir(parents=True, exist_ok=True)
    handle = tempfile.NamedTemporaryFile(
        mode="w",
        suffix=output.suffix,
        prefix=".docklens-sidecar-",
        dir=output.parent,
        encoding="utf-8",
        newline="",
        delete=False,
    )
    temporary = Path(handle.name)
    try:
        handle.write(text)
        handle.close()
        os.replace(temporary, output)
    except Exception:
        handle.close()
        temporary.unlink(missing_ok=True)
        raise


def _bundle_paths(prefix: Path, formats: tuple[str, ...]) -> tuple[Path, ...]:
    return tuple(
        [*(prefix.with_suffix(f".{file_format}") for file_format in formats),
         prefix.with_name(prefix.name + "_data.csv"),
         prefix.with_name(prefix.name + "_manifest.json")]
    )


def _write_figure_bundle(
    artifact: ChartArtifact,
    prefix: Path,
    *,
    formats: tuple[str, ...],
    dpi=300,
    extra_metadata=None,
) -> tuple[Path, ...]:
    """Write one bundle into a staging directory."""
    for file_format in formats:
        output = prefix.with_suffix(f".{file_format}")
        _atomic_figure(
            artifact, output, file_format=file_format, dpi=int(dpi)
        )
    data_path = prefix.with_name(prefix.name + "_data.csv")
    csv_text = _safe_frame(artifact.data).to_csv(index=False)
    _atomic_text(data_path, csv_text)
    metadata = {
        "schema": "docklens-figure-bundle-v1",
        "docklens_version": __version__,
        "kind": artifact.kind,
        "dpi": int(dpi),
        "formats": list(formats),
        **(
            {"tiff_compression": _TIFF_COMPRESSION}
            if "tiff" in formats
            else {}
        ),
        **dict(artifact.metadata),
        **dict(extra_metadata or {}),
    }
    manifest_path = prefix.with_name(prefix.name + "_manifest.json")
    _atomic_text(
        manifest_path,
        json.dumps(metadata, ensure_ascii=False, indent=2, sort_keys=True),
    )
    return _bundle_paths(prefix, formats)


def _reserve_backup(destination: Path) -> Path:
    handle = tempfile.NamedTemporaryFile(
        prefix=".docklens-backup-",
        dir=destination.parent,
        delete=False,
    )
    backup = Path(handle.name)
    handle.close()
    backup.unlink(missing_ok=True)
    return backup


def _commit_files(staged: tuple[Path, ...], destinations: tuple[Path, ...]):
    if len(staged) != len(destinations):
        raise ValueError("staged and destination file counts must match")
    backups = []
    installed = []
    try:
        for source, destination in zip(staged, destinations):
            destination.parent.mkdir(parents=True, exist_ok=True)
            backup = None
            if destination.exists():
                backup = _reserve_backup(destination)
                os.replace(destination, backup)
            backups.append((backup, destination))
            os.replace(source, destination)
            installed.append(destination)
    except Exception:
        for destination in installed:
            destination.unlink(missing_ok=True)
        for backup, destination in reversed(backups):
            if backup is not None and backup.exists():
                os.replace(backup, destination)
        raise
    else:
        for backup, _destination in backups:
            if backup is not None:
                backup.unlink(missing_ok=True)


def _normalize_formats(formats) -> tuple[str, ...]:
    requested_formats = (str(value).lower().lstrip(".") for value in formats)
    normalized_formats = tuple(
        dict.fromkeys(
            "tiff" if value == "tif" else value for value in requested_formats
        )
    )
    if not normalized_formats or any(
        value not in _FORMATS for value in normalized_formats
    ):
        raise ValueError("formats must contain only png, svg, pdf or tiff")
    return normalized_formats


def export_figure_bundle(
    artifact: ChartArtifact,
    path_prefix,
    *,
    formats=("png",),
    dpi=300,
    extra_metadata=None,
) -> tuple[str, ...]:
    """Export figure(s), exact chart rows and a reproducibility manifest."""
    normalized_formats = _normalize_formats(formats)
    if dpi < 1:
        raise ValueError("dpi must be greater than zero")
    prefix = Path(os.fspath(path_prefix))
    if prefix.suffix.lower().lstrip(".") in _FORMATS | {"tif", "csv", "json"}:
        prefix = prefix.with_suffix("")
    prefix.parent.mkdir(parents=True, exist_ok=True)
    staging = Path(
        tempfile.mkdtemp(prefix=f".docklens-bundle-{prefix.name}-", dir=prefix.parent)
    )
    try:
        staged_prefix = staging / prefix.name
        staged_files = _write_figure_bundle(
            artifact,
            staged_prefix,
            formats=normalized_formats,
            dpi=dpi,
            extra_metadata=extra_metadata,
        )
        destinations = _bundle_paths(prefix, normalized_formats)
        _commit_files(staged_files, destinations)
        return tuple(str(path) for path in destinations)
    finally:
        shutil.rmtree(staging, ignore_errors=True)


def export_tiff_collection(
    artifacts: Mapping[str, ChartArtifact],
    output_directory,
    *,
    dpi=PUBLICATION_TIFF_DPI,
    metadata_factory: Callable[[ChartArtifact], Mapping] | None = None,
) -> tuple[str, ...]:
    """Export every named chart as an auditable LZW TIFF bundle."""
    items = tuple(artifacts.items())
    if not items:
        raise ValueError("artifacts must contain at least one chart")
    if dpi < 1:
        raise ValueError("dpi must be greater than zero")
    if metadata_factory is not None and not callable(metadata_factory):
        raise TypeError("metadata_factory must be callable")
    for name, artifact in items:
        if not _CHART_NAME_PATTERN.fullmatch(str(name)):
            raise ValueError(f"invalid chart name: {name!r}")
        if not isinstance(artifact, ChartArtifact):
            raise TypeError(f"chart {name!r} must be a ChartArtifact")

    directory = Path(os.fspath(output_directory))
    directory.mkdir(parents=True, exist_ok=True)
    staging = Path(
        tempfile.mkdtemp(prefix=".docklens-tiff-collection-", dir=directory)
    )
    staged_files = []
    destinations = []
    figures = []
    try:
        for name, artifact in items:
            metadata = metadata_factory(artifact) if metadata_factory else None
            staged_prefix = staging / str(name)
            bundle = _write_figure_bundle(
                artifact,
                staged_prefix,
                formats=("tiff",),
                dpi=int(dpi),
                extra_metadata=metadata,
            )
            staged_files.extend(bundle)
            destinations.extend(
                directory / path.relative_to(staging) for path in bundle
            )
            figures.append(
                {
                    "name": str(name),
                    "kind": artifact.kind,
                    "files": [path.name for path in bundle],
                }
            )

        staged_manifest = staging / "docklens_tiff_collection_manifest.json"
        _atomic_text(
            staged_manifest,
            json.dumps(
                {
                    "schema": "docklens-tiff-collection-v1",
                    "docklens_version": __version__,
                    "dpi": int(dpi),
                    "format": "tiff",
                    "compression": _TIFF_COMPRESSION,
                    "figures": figures,
                },
                ensure_ascii=False,
                indent=2,
                sort_keys=True,
            ),
        )
        staged_files.append(staged_manifest)
        destinations.append(directory / staged_manifest.name)
        _commit_files(tuple(staged_files), tuple(destinations))
        return tuple(str(path) for path in destinations)
    finally:
        shutil.rmtree(staging, ignore_errors=True)


__all__ = [
    "PUBLICATION_TIFF_DPI",
    "export_figure_bundle",
    "export_tiff_collection",
]
