"""Atomic export of publication figures with data and provenance sidecars."""

from __future__ import annotations

import json
import os
import tempfile
from pathlib import Path

import pandas as pd

from . import __version__
from .plotting import ChartArtifact


_FORMATS = frozenset({"png", "svg", "pdf"})
_FORMULA_PREFIXES = ("=", "+", "-", "@")


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
        artifact.figure.savefig(
            temporary,
            format=file_format,
            dpi=dpi,
            bbox_inches="tight",
            facecolor=artifact.figure.get_facecolor(),
        )
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


def export_figure_bundle(
    artifact: ChartArtifact,
    path_prefix,
    *,
    formats=("png",),
    dpi=300,
    extra_metadata=None,
) -> tuple[str, ...]:
    """Export figure(s), exact chart rows and a reproducibility manifest."""
    normalized_formats = tuple(
        dict.fromkeys(str(value).lower().lstrip(".") for value in formats)
    )
    if not normalized_formats or any(
        value not in _FORMATS for value in normalized_formats
    ):
        raise ValueError("formats must contain only png, svg or pdf")
    if dpi < 1:
        raise ValueError("dpi must be greater than zero")
    prefix = Path(os.fspath(path_prefix))
    if prefix.suffix.lower().lstrip(".") in _FORMATS | {"csv", "json"}:
        prefix = prefix.with_suffix("")
    written = []
    for file_format in normalized_formats:
        output = prefix.with_suffix(f".{file_format}")
        _atomic_figure(
            artifact, output, file_format=file_format, dpi=int(dpi)
        )
        written.append(str(output))
    data_path = prefix.with_name(prefix.name + "_data.csv")
    csv_text = _safe_frame(artifact.data).to_csv(index=False)
    _atomic_text(data_path, csv_text)
    written.append(str(data_path))
    metadata = {
        "schema": "docklens-figure-bundle-v1",
        "docklens_version": __version__,
        "kind": artifact.kind,
        "dpi": int(dpi),
        "formats": list(normalized_formats),
        **dict(artifact.metadata),
        **dict(extra_metadata or {}),
    }
    manifest_path = prefix.with_name(prefix.name + "_manifest.json")
    _atomic_text(
        manifest_path,
        json.dumps(metadata, ensure_ascii=False, indent=2, sort_keys=True),
    )
    written.append(str(manifest_path))
    return tuple(written)


__all__ = ["export_figure_bundle"]
