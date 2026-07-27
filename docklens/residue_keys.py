"""Canonical parsing and matching for user-selected receptor residues."""

from __future__ import annotations

import re
from dataclasses import dataclass
from typing import Iterable, Tuple


_RESIDUE_KEY = re.compile(r"^[A-Z]{1,4}-?\d+[A-Z]?$")
_RESIDUE_NAME = re.compile(r"^[A-Z]{1,4}$")
_RESIDUE_NUMBER = re.compile(r"^-?\d+[A-Z]?$")
_CHAINED_RESIDUE = re.compile(r"^([A-Z]{1,4}-?\d+)([A-Z])$")


@dataclass(frozen=True)
class KeyResidueParse:
    keys: Tuple[str, ...] = ()
    invalid: Tuple[str, ...] = ()


@dataclass(frozen=True)
class KeyResidueMatch:
    matched_keys: Tuple[str, ...] = ()
    matched_residues: Tuple[str, ...] = ()
    ambiguous_keys: Tuple[str, ...] = ()
    unmatched_keys: Tuple[str, ...] = ()


def _tokenize(items) -> Tuple[str, ...]:
    if isinstance(items, str):
        text = items
    else:
        text = " ".join(str(item) for item in (items or ()))
    text = re.sub(r"[,;\r\n\t]+", " ", text.upper())
    return tuple(part for part in text.split() if part)


def parse_key_residues(items) -> KeyResidueParse:
    """Parse common residue-list notation without silently keeping punctuation."""
    tokens = _tokenize(items)
    keys = set()
    invalid = []
    index = 0
    while index < len(tokens):
        token = tokens[index]
        if (
            _RESIDUE_NAME.fullmatch(token)
            and index + 1 < len(tokens)
            and _RESIDUE_NUMBER.fullmatch(tokens[index + 1])
        ):
            token = token + tokens[index + 1]
            index += 2
        else:
            index += 1
        if _RESIDUE_KEY.fullmatch(token):
            keys.add(token)
        else:
            invalid.append(token)
    return KeyResidueParse(
        keys=tuple(sorted(keys)),
        invalid=tuple(sorted(set(invalid))),
    )


def normalize_key_residues(items) -> frozenset:
    return frozenset(parse_key_residues(items).keys)


def _without_chain(residue: str) -> str:
    match = _CHAINED_RESIDUE.fullmatch(residue)
    return match.group(1) if match else residue


def match_key_residues(
    keys: Iterable[str],
    receptor_residues: Iterable[str],
) -> KeyResidueMatch:
    """Resolve configured keys against concrete receptor residues and chains."""
    normalized_keys = normalize_key_residues(keys)
    residues = tuple(sorted({str(value).strip().upper() for value in receptor_residues}))
    matched_keys = []
    matched_residues = set()
    ambiguous = []
    unmatched = []
    for key in sorted(normalized_keys):
        concrete = tuple(
            residue
            for residue in residues
            if residue == key or _without_chain(residue) == key
        )
        if not concrete:
            unmatched.append(key)
            continue
        matched_keys.append(key)
        matched_residues.update(concrete)
        if len(concrete) > 1 and key not in concrete:
            ambiguous.append(key)
    return KeyResidueMatch(
        matched_keys=tuple(matched_keys),
        matched_residues=tuple(sorted(matched_residues)),
        ambiguous_keys=tuple(ambiguous),
        unmatched_keys=tuple(unmatched),
    )
