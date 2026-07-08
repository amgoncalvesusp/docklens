"""
parser_pdb.py — PDB reader (ATOM / HETATM, CONECT).

Element comes from columns 77-78 when present, else is guessed from the atom
name. Connectivity uses CONECT records when present; otherwise bonds are
inferred by distance (covalent radii + tolerance) via ``structures.infer_bonds``.
MODEL/ENDMDL blocks are returned as separate poses.
"""

from __future__ import annotations

from .interaction_core import Atom
from .structures import ParsedPose, infer_bonds, normalize_element, parse_sol


def _element_from_name(name):
    """Guess an element from a PDB atom name (fallback when col 77-78 empty)."""
    n = name.strip()
    letters = "".join(ch for ch in n if ch.isalpha())
    if not letters:
        return ""
    two = letters[:2].capitalize()
    if two in {
        "Cl",
        "Br",
        "Fe",
        "Zn",
        "Mg",
        "Mn",
        "Ca",
        "Na",
        "Cu",
        "Ni",
        "Co",
        "Cd",
        "Hg",
        "Se",
    }:
        return two
    return letters[0].upper()


def _parse_charge(field):
    """Parse a PDB charge field like '1+' or '2-' into a signed int."""
    field = (field or "").strip()
    if not field:
        return 0
    sign = -1 if "-" in field else 1
    digits = "".join(ch for ch in field if ch.isdigit())
    return sign * int(digits) if digits else 0


def _parse_atom_line(line, idx):
    record = line[0:6].strip()
    serial = int(line[6:11]) if line[6:11].strip() else idx + 1
    name = line[12:16].strip()
    resn = line[17:20].strip()
    chain = line[21:22].strip()
    resi = line[22:26].strip()
    x = float(line[30:38])
    y = float(line[38:46])
    z = float(line[46:54])
    elem_field = line[76:78].strip() if len(line) >= 78 else ""
    charge_field = line[78:80] if len(line) >= 80 else ""
    elem = normalize_element(elem_field) if elem_field else _element_from_name(name)
    atom = Atom(
        idx=idx,
        elem=elem,
        name=name,
        resn=resn,
        resi=resi,
        chain=chain,
        coord=(x, y, z),
        fcharge=_parse_charge(charge_field),
        serial=serial,
    )
    return atom, (record == "HETATM")


def _finish_pose(atoms, atom_by_serial, is_het, conect, source_file, pose_index):
    if conect:
        for a1, others in conect.items():
            u = atom_by_serial.get(a1)
            if u is None:
                continue
            for a2 in others:
                v = atom_by_serial.get(a2)
                if v is not None and v not in u.neighbors:
                    u.neighbors.append(v)
                    v.neighbors.append(u)
    else:
        infer_bonds(atoms)
    return ParsedPose(
        atoms=atoms,
        atom_by_serial=atom_by_serial,
        fmt="pdb",
        source_file=source_file,
        is_hetatm=is_het,
        pose_index=pose_index,
        sol=parse_sol(source_file),
    )


def parse_pdb(path):
    """Parse a .pdb file into a list of ParsedPose (one per MODEL, or one)."""
    poses = []
    atoms, atom_by_serial, is_het, conect = [], {}, {}, {}
    idx = 0

    def flush():
        nonlocal atoms, atom_by_serial, is_het, conect, idx
        if atoms:
            poses.append(
                _finish_pose(atoms, atom_by_serial, is_het, conect, path, len(poses))
            )
        atoms, atom_by_serial, is_het, conect, idx = [], {}, {}, {}, 0

    with open(path, "r", errors="replace") as fh:
        for line in fh:
            rec = line[0:6].strip()
            if rec in ("ATOM", "HETATM"):
                atom, het = _parse_atom_line(line, idx)
                atoms.append(atom)
                atom_by_serial[atom.serial] = atom
                is_het[atom.serial] = het
                idx += 1
            elif rec == "CONECT":
                nums = line[6:].split()
                if nums:
                    try:
                        base = int(nums[0])
                        conect.setdefault(base, []).extend(int(n) for n in nums[1:])
                    except ValueError:
                        pass
            elif rec == "ENDMDL":
                flush()
    flush()
    return poses
