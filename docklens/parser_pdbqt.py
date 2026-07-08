"""
parser_pdbqt.py — AutoDock PDBQT reader.

Same column layout as PDB for coordinates/names, plus the AutoDock atom type in
columns 78-79 (used for the element). Each MODEL/ENDMDL block is a separate pose.
Docking score is read best-effort from the file only:
  * AutoDock Vina : 'REMARK VINA RESULT:  <score> ...'
  * AutoDock4     : '... Estimated Free Energy of Binding = <score> ...'
Score is left as None when absent (never invented). Bonds are inferred by
distance (PDBQT has no explicit connectivity).
"""

from __future__ import annotations

import re

from .interaction_core import Atom
from .structures import ParsedPose, infer_bonds, normalize_element, parse_sol

# AutoDock atom type -> element.
_AD_TYPE_ELEM = {
    "A": "C",  # aromatic carbon
    "C": "C",
    "N": "N",
    "NA": "N",
    "NS": "N",
    "O": "O",
    "OA": "O",
    "OS": "O",
    "S": "S",
    "SA": "S",
    "H": "H",
    "HD": "H",
    "HS": "H",
    "F": "F",
    "Cl": "Cl",
    "CL": "Cl",
    "Br": "Br",
    "BR": "Br",
    "I": "I",
    "P": "P",
    "Mg": "Mg",
    "MG": "Mg",
    "Zn": "Zn",
    "ZN": "Zn",
    "Mn": "Mn",
    "MN": "Mn",
    "Ca": "Ca",
    "CA": "Ca",
    "Fe": "Fe",
    "FE": "Fe",
    "K": "K",
    "Na": "Na",
}

_VINA = re.compile(r"REMARK\s+VINA\s+RESULT:\s*(-?\d+\.?\d*)", re.IGNORECASE)
_AD4 = re.compile(
    r"Estimated Free Energy of Binding\s*=?\s*(-?\d+\.?\d*)", re.IGNORECASE
)


def _element_from_adtype(adtype, name):
    adtype = (adtype or "").strip()
    if adtype in _AD_TYPE_ELEM:
        return _AD_TYPE_ELEM[adtype]
    if adtype:
        return normalize_element(adtype[:2] if len(adtype) > 1 else adtype)
    letters = "".join(ch for ch in name if ch.isalpha())
    return letters[:1].upper() if letters else ""


def _parse_atom_line(line, idx):
    serial = int(line[6:11]) if line[6:11].strip() else idx + 1
    name = line[12:16].strip()
    resn = line[17:20].strip()
    chain = line[21:22].strip()
    resi = line[22:26].strip()
    x = float(line[30:38])
    y = float(line[38:46])
    z = float(line[46:54])
    adtype = line[77:79].strip() if len(line) >= 79 else line[76:].strip()
    elem = _element_from_adtype(adtype, name)
    return Atom(
        idx=idx,
        elem=elem,
        name=name,
        resn=resn,
        resi=resi,
        chain=chain,
        coord=(x, y, z),
        fcharge=0,  # pdbqt stores partial (Gasteiger) charges only
        serial=serial,
    )


def _finish_pose(atoms, atom_by_serial, source_file, pose_index, score):
    infer_bonds(atoms)
    return ParsedPose(
        atoms=atoms,
        atom_by_serial=atom_by_serial,
        fmt="pdbqt",
        source_file=source_file,
        is_hetatm={a.serial: True for a in atoms},  # a docked ligand is all-HETATM
        pose_index=pose_index,
        score=score,
        sol=parse_sol(source_file),
    )


def parse_pdbqt(path):
    """Parse a .pdbqt file into a list of ParsedPose (one per MODEL, or one)."""
    poses = []
    atoms, atom_by_serial = [], {}
    idx = 0
    score = None

    def flush():
        nonlocal atoms, atom_by_serial, idx, score
        if atoms:
            poses.append(_finish_pose(atoms, atom_by_serial, path, len(poses), score))
        atoms, atom_by_serial, idx, score = [], {}, 0, None

    with open(path, "r", errors="replace") as fh:
        for line in fh:
            rec = line[0:6].strip()
            if rec in ("ATOM", "HETATM"):
                atom = _parse_atom_line(line, idx)
                atoms.append(atom)
                atom_by_serial[atom.serial] = atom
                idx += 1
            elif rec == "ENDMDL":
                flush()
            else:
                m = _VINA.search(line) or _AD4.search(line)
                if m:
                    try:
                        score = float(m.group(1))
                    except ValueError:
                        pass
    flush()
    return poses
