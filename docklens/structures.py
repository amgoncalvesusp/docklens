"""
structures.py — shared parser data container + distance-based bond inference.

Glue between the four format parsers and the entity resolver. Holds the
``ParsedPose`` schema every parser returns, plus ``infer_bonds`` (covalent-radii
+ tolerance) used when a format has no explicit connectivity (PDB without
CONECT, PDBQT).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional

import numpy as np

from .interaction_core import Atom  # noqa: F401  (re-exported for parsers)

# Covalent radii (Angstrom), Cordero et al. 2008 (subset relevant to bio/docking).
COVALENT_RADII = {
    "H": 0.31,
    "C": 0.76,
    "N": 0.71,
    "O": 0.66,
    "F": 0.57,
    "P": 1.07,
    "S": 1.05,
    "Cl": 1.02,
    "Br": 1.20,
    "I": 1.39,
    "B": 0.84,
    "Si": 1.11,
    "Se": 1.20,
    "Na": 1.66,
    "K": 2.03,
    "Mg": 1.41,
    "Ca": 1.76,
    "Mn": 1.39,
    "Fe": 1.32,
    "Co": 1.26,
    "Ni": 1.24,
    "Cu": 1.32,
    "Zn": 1.22,
    "Cd": 1.44,
    "Hg": 1.32,
}
_DEFAULT_RADIUS = 0.77
# Bond if dist <= r1 + r2 + tolerance. 0.45 A is a common, forgiving value.
BOND_TOLERANCE = 0.45


@dataclass
class ParsedPose:
    """One structural pose parsed from a file.

    atoms          list[Atom] in file order (neighbors filled if bonds known)
    atom_by_serial serial -> Atom (serial = file atom id, 1-based)
    ccdc_ligand    set of serials flagged CCDC_LIGAND (mol2), or None
    ccdc_receptor  set of serials flagged CCDC_AMINOACID (mol2), or None
    group_subst    subst_id -> substructure type ('RESIDUE'/'GROUP'), or {}
    is_hetatm      serial -> bool (PDB/PDBQT HETATM flag), or {}
    fmt            'mol2' | 'pdb' | 'pdbqt'
    source_file    path
    pose_index     0-based index within a multi-MODEL file
    score          docking score if present in the file, else None
    sol            'sol<N>' number parsed from the filename, else None
    """

    atoms: list
    atom_by_serial: dict
    fmt: str
    source_file: str
    ccdc_ligand: Optional[set] = None
    ccdc_receptor: Optional[set] = None
    group_subst: dict = field(default_factory=dict)
    is_hetatm: dict = field(default_factory=dict)
    pose_index: int = 0
    score: Optional[float] = None
    sol: Optional[int] = None


def normalize_element(sym: str) -> str:
    """Normalise an element symbol to standard capitalisation ('CL' -> 'Cl')."""
    sym = (sym or "").strip()
    if not sym:
        return ""
    if len(sym) == 1:
        return sym.upper()
    return sym[0].upper() + sym[1:].lower()


def infer_bonds(atoms, tol: float = BOND_TOLERANCE) -> None:
    """Fill `atom.neighbors` by distance using covalent radii + tolerance.

    Uses a uniform spatial grid (cell size = max bond length) so it stays close
    to O(N) instead of O(N^2) for large receptors.
    """
    if not atoms:
        return
    coords = np.array([a.coord for a in atoms], dtype=float)
    radii = np.array([COVALENT_RADII.get(a.elem, _DEFAULT_RADIUS) for a in atoms])
    max_r = float(radii.max()) if len(radii) else _DEFAULT_RADIUS
    cell = max(2.0 * max_r + tol, 1.0)  # any bond fits within one cell step

    mins = coords.min(axis=0)
    cells = np.floor((coords - mins) / cell).astype(int)
    grid = {}
    for i, c in enumerate(map(tuple, cells)):
        grid.setdefault(c, []).append(i)

    neighbor_offsets = [
        (dx, dy, dz) for dx in (-1, 0, 1) for dy in (-1, 0, 1) for dz in (-1, 0, 1)
    ]
    for i in range(len(atoms)):
        ci = tuple(cells[i])
        ai = atoms[i]
        for off in neighbor_offsets:
            key = (ci[0] + off[0], ci[1] + off[1], ci[2] + off[2])
            for j in grid.get(key, ()):
                if j <= i:
                    continue
                aj = atoms[j]
                d = float(np.linalg.norm(coords[i] - coords[j]))
                if d < 0.4:  # overlapping / duplicate atoms: not a real bond
                    continue
                if d <= radii[i] + radii[j] + tol:
                    ai.neighbors.append(aj)
                    aj.neighbors.append(ai)


def parse_sol(filename: str):
    """Return the integer N from a 'sol<N>' token in the filename, or None."""
    import re

    m = re.search(r"sol(\d+)", filename, flags=re.IGNORECASE)
    return int(m.group(1)) if m else None
