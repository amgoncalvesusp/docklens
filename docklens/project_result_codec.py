"""Incremental result storage for DockLens project schema 5.

Each NDJSON chunk is bounded and the reader processes raw chunk bytes one at a
time. The cached ``RunResult`` remains immutable at the public boundary, but
the container no longer requires a monolithic JSON representation of it.
"""

from __future__ import annotations

from dataclasses import fields
import hashlib
import json
import re
from typing import Any, Mapping

from .results import AnalysisParameters, Detail, InputQC, RunResult, Summary
from .results import Endpoint


RESULT_FORMAT = "ndjson-v1"
MAX_RESULT_CHUNK_BYTES = 16 * 1024 * 1024
MAX_PROJECT_BYTES = 128 * 1024 * 1024
MAX_TOTAL_UNCOMPRESSED_BYTES = 512 * 1024 * 1024
MAX_ARCHIVE_ENTRIES = 128
MAX_JSON_DEPTH = 40

_CHUNK_RE = re.compile(
    r"^results/system-(?:a|b)/(?:summaries|details|input-qc)/part-\d{4}\.ndjson$"
)


class ProjectLimitError(ValueError):
    """The project cannot be represented within DockLens safety limits."""


class ProjectSerializationError(ValueError):
    """A result could not be encoded or decoded as a project payload."""


class ResultIntegrityError(ValueError):
    """A result chunk differs from its manifest declaration."""


def is_chunk_entry(name: str) -> bool:
    return bool(_CHUNK_RE.fullmatch(str(name)))


def encode_result(
    result: RunResult,
    dataset_key: str,
    entries: dict[str, bytes],
) -> dict[str, Any]:
    """Encode result records as bounded NDJSON chunks and return a reference."""

    if dataset_key not in {"system-a", "system-b"}:
        raise ProjectSerializationError("invalid project dataset key")
    if result.pending:
        raise ProjectSerializationError(
            "unresolved confirmations must be completed before saving a cached result"
        )

    categories = {
        "summaries": result.summaries,
        "details": result.details,
        "input-qc": result.input_qc,
    }
    chunk_refs = {
        category: _encode_category(records, dataset_key, category, entries)
        for category, records in categories.items()
    }
    return {
        "format": RESULT_FORMAT,
        "dataset_key": dataset_key,
        "chunks": chunk_refs,
        "pending": [],
        "key_residues": sorted(result.key_residues),
        "receptor_residues": sorted(result.receptor_residues),
        "parameters": _dataclass_payload(result.parameters),
    }


def decode_result(
    archive,
    reference: Mapping[str, Any],
    names: set[str],
) -> tuple[RunResult, set[str]]:
    """Read and verify result chunks one line at a time."""

    if not isinstance(reference, Mapping):
        raise ProjectSerializationError("cached result reference must be an object")
    allowed = {
        "format",
        "dataset_key",
        "chunks",
        "pending",
        "key_residues",
        "receptor_residues",
        "parameters",
    }
    if not set(reference).issubset(allowed):
        raise ProjectSerializationError("cached result reference contains unknown fields")
    if reference.get("format") != RESULT_FORMAT:
        raise ProjectSerializationError("unsupported cached result format")
    if reference.get("dataset_key") not in {"system-a", "system-b"}:
        raise ProjectSerializationError("cached result dataset key is invalid")
    pending = reference.get("pending", [])
    if not isinstance(pending, list) or pending:
        raise ProjectSerializationError("cached result contains unresolved confirmations")
    key_residues = _text_list(reference, "key_residues")
    receptor_residues = _text_list(reference, "receptor_residues")
    parameters = _construct(AnalysisParameters, _mapping(reference, "parameters"))
    chunks = reference.get("chunks")
    if not isinstance(chunks, dict):
        raise ProjectSerializationError("cached result chunks must be an object")
    if set(chunks) != {"summaries", "details", "input-qc"}:
        raise ProjectSerializationError("cached result chunk categories are invalid")

    records: dict[str, list[Any]] = {
        "summaries": [],
        "details": [],
        "input-qc": [],
    }
    referenced: set[str] = set()
    for category in ("summaries", "details", "input-qc"):
        category_refs = chunks[category]
        if not isinstance(category_refs, list):
            raise ProjectSerializationError("cached result chunk list is invalid")
        for chunk_ref in category_refs:
            if not isinstance(chunk_ref, Mapping):
                raise ProjectSerializationError("cached result chunk metadata is invalid")
            if not set(chunk_ref).issubset(
                {"entry", "sha256", "size_bytes", "records"}
            ):
                raise ProjectSerializationError(
                    "cached result chunk metadata contains unknown fields"
                )
            entry = _text(chunk_ref, "entry")
            if (
                not is_chunk_entry(entry)
                or entry not in names
                or entry in referenced
                or not entry.startswith("results/system-")
                or not entry.startswith(
                    "results/%s/" % reference["dataset_key"]
                )
            ):
                raise ResultIntegrityError("cached result entry is missing or unsafe")
            expected_hash = _sha_text(chunk_ref, "sha256")
            expected_size = _integer(chunk_ref, "size_bytes")
            expected_records = _integer(chunk_ref, "records")
            if expected_size < 0 or expected_records < 0:
                raise ProjectSerializationError("cached result chunk metadata is negative")
            info = archive.getinfo(entry)
            if info.file_size > MAX_RESULT_CHUNK_BYTES:
                raise ProjectLimitError("cached result chunk is too large")
            decoded, actual_size, actual_hash, count = _read_chunk(
                archive, entry, category
            )
            if actual_size != expected_size:
                raise ResultIntegrityError("cached result integrity size does not match")
            if actual_hash != expected_hash:
                raise ResultIntegrityError("cached result integrity hash does not match")
            if count != expected_records:
                raise ResultIntegrityError("cached result integrity record count does not match")
            records[category].extend(decoded)
            referenced.add(entry)

    try:
        result = RunResult(
            details=tuple(records["details"]),
            summaries=tuple(records["summaries"]),
            pending=(),
            key_residues=frozenset(key_residues),
            receptor_residues=frozenset(receptor_residues),
            input_qc=tuple(records["input-qc"]),
            parameters=parameters,
        )
    except (KeyError, TypeError, ValueError) as exc:
        raise ProjectSerializationError("invalid cached analysis result") from exc
    return result, referenced


def _encode_category(records, dataset_key, category, entries):
    chunks = []
    buffer = bytearray()
    record_count = 0
    digest = hashlib.sha256()

    def flush() -> None:
        nonlocal buffer, record_count, digest
        if not buffer:
            return
        entry = (
            f"results/{dataset_key}/{category}/part-{len(chunks) + 1:04d}.ndjson"
        )
        payload = bytes(buffer)
        entries[entry] = payload
        chunks.append(
            {
                "entry": entry,
                "sha256": digest.hexdigest(),
                "size_bytes": len(payload),
                "records": record_count,
            }
        )
        buffer = bytearray()
        record_count = 0
        digest = hashlib.sha256()

    for record in records:
        try:
            line = (
                json.dumps(
                    _dataclass_payload(record),
                    ensure_ascii=False,
                    sort_keys=True,
                    separators=(",", ":"),
                    allow_nan=False,
                ).encode("utf-8")
                + b"\n"
            )
        except (TypeError, ValueError) as exc:
            raise ProjectSerializationError("could not serialize analysis result") from exc
        if len(line) > MAX_RESULT_CHUNK_BYTES:
            raise ProjectLimitError("one cached result record exceeds the safe chunk limit")
        if buffer and len(buffer) + len(line) > MAX_RESULT_CHUNK_BYTES:
            flush()
        buffer.extend(line)
        digest.update(line)
        record_count += 1
    flush()
    return chunks


def _read_chunk(archive, entry, category):
    digest = hashlib.sha256()
    decoded = []
    total = 0
    count = 0
    try:
        with archive.open(entry, "r") as handle:
            for raw_line in handle:
                total += len(raw_line)
                if total > MAX_RESULT_CHUNK_BYTES:
                    raise ProjectLimitError("cached result chunk is too large")
                digest.update(raw_line)
                if not raw_line.strip():
                    raise ResultIntegrityError(
                        "cached result integrity content is invalid"
                    )
                document = json.loads(raw_line.decode("utf-8"))
                if not isinstance(document, dict):
                    raise ProjectSerializationError("cached result record must be an object")
                _check_json_depth(document)
                decoded.append(_decode_record(category, document))
                count += 1
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise ResultIntegrityError(
            "cached result integrity content is invalid"
        ) from exc
    except ProjectSerializationError as exc:
        raise ResultIntegrityError(
            "cached result integrity content is invalid"
        ) from exc
    except ValueError:
        raise
    except OSError as exc:
        raise ProjectSerializationError("could not read cached result chunk") from exc
    return decoded, total, digest.hexdigest(), count


def _dataclass_payload(instance: Any) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for item in fields(instance):
        value = getattr(instance, item.name)
        if hasattr(value, "__dataclass_fields__"):
            result[item.name] = _dataclass_payload(value)
        elif isinstance(value, Mapping):
            result[item.name] = dict(value)
        elif isinstance(value, (tuple, list, frozenset, set)):
            result[item.name] = list(value)
        else:
            result[item.name] = value
    return result


def _decode_record(category: str, payload: Mapping[str, Any]):
    try:
        if category == "details":
            return _detail_from_payload(payload)
        if category == "summaries":
            return _summary_from_payload(payload)
        if category == "input-qc":
            return _construct(InputQC, payload)
    except (KeyError, TypeError, ValueError) as exc:
        raise ProjectSerializationError("invalid cached result record") from exc
    raise ProjectSerializationError("cached result category is invalid")


def _detail_from_payload(payload):
    values = dict(payload)
    values["ligand"] = _construct(Endpoint, _mapping(payload, "ligand"))
    values["receptor"] = _construct(Endpoint, _mapping(payload, "receptor"))
    water = payload.get("water")
    values["water"] = _construct(Endpoint, water) if isinstance(water, dict) else None
    return _construct(Detail, values)


def _summary_from_payload(payload):
    values = dict(payload)
    counts = values.get("counts")
    if not isinstance(counts, dict):
        raise ValueError("summary counts must be an object")
    values["counts"] = {
        str(key): _strict_int(value) for key, value in counts.items()
    }
    return _construct(Summary, values)


def _construct(cls, payload):
    if not isinstance(payload, Mapping):
        raise ValueError("record must be an object")
    allowed = {item.name for item in fields(cls)}
    if not set(payload).issubset(allowed):
        raise ValueError("record contains unknown fields")
    values = dict(payload)
    for name in (
        "atom_serials",
        "warnings",
        "cutoffs",
        "interaction_types",
        "key_residues",
    ):
        if name in values:
            values[name] = tuple(values[name])
    if "cutoffs" in values:
        values["cutoffs"] = tuple(
            (str(name), float(value)) for name, value in values["cutoffs"]
        )
    return cls(**values)


def _mapping(payload, name):
    value = payload[name]
    if not isinstance(value, dict):
        raise ValueError("%s must be an object" % name)
    return value


def _text(payload, name):
    value = payload[name]
    if not isinstance(value, str):
        raise ValueError("%s must be text" % name)
    return value


def _text_list(payload, name):
    value = payload[name]
    if not isinstance(value, list) or any(not isinstance(item, str) for item in value):
        raise ValueError("%s must be a text list" % name)
    return tuple(value)


def _sha_text(payload, name):
    value = _text(payload, name).lower()
    if not _valid_sha(value):
        raise ValueError("%s must be a SHA-256 digest" % name)
    return value


def _integer(payload, name):
    value = payload[name]
    if isinstance(value, bool) or not isinstance(value, int):
        raise ValueError("%s must be an integer" % name)
    return value


def _strict_int(value):
    if isinstance(value, bool) or not isinstance(value, int):
        raise ValueError("value must be an integer")
    return value


def _valid_sha(value):
    return len(value) == 64 and all(char in "0123456789abcdef" for char in value)


def _sha256_bytes(payload):
    return hashlib.sha256(payload).hexdigest()


def _check_json_depth(value: Any, depth: int = 0) -> None:
    if depth > MAX_JSON_DEPTH:
        raise ProjectSerializationError("cached result record is too deeply nested")
    if isinstance(value, dict):
        for key, item in value.items():
            if not isinstance(key, str):
                raise ProjectSerializationError("cached result record keys must be text")
            _check_json_depth(item, depth + 1)
    elif isinstance(value, list):
        for item in value:
            _check_json_depth(item, depth + 1)


__all__ = [
    "MAX_ARCHIVE_ENTRIES",
    "MAX_PROJECT_BYTES",
    "MAX_RESULT_CHUNK_BYTES",
    "MAX_TOTAL_UNCOMPRESSED_BYTES",
    "ProjectLimitError",
    "ProjectSerializationError",
    "RESULT_FORMAT",
    "ResultIntegrityError",
    "decode_result",
    "encode_result",
    "is_chunk_entry",
]
