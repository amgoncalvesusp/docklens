"""
parser_mol2.py — TRIPOS .mol2 reader (Discovery Studio / CCDC GOLD flavour).

Reads @<TRIPOS>ATOM (element from the SYBYL atom type), @<TRIPOS>BOND (explicit
connectivity), @<TRIPOS>SUBSTRUCTURE (per-residue chain + RESIDUE/GROUP type)
and @<TRIPOS>SET (CCDC_LIGAND / CCDC_AMINOACID serial lists). Returns one
``ParsedPose`` per @<TRIPOS>MOLECULE block.
"""

from __future__ import annotations

import re

from .interaction_core import Atom
from .structures import ParsedPose, normalize_element, parse_sol

_RESN_RESI = re.compile(r"^([A-Za-z][A-Za-z0-9]*?)(\d+)$")


def _split_resn_resi(subst_name, subst_id):
    """Split a mol2 subst_name like 'THR30' into ('THR','30').

    Falls back to (subst_name, subst_id) for non-standard names (e.g. the
    peptide-ligand group name '<1>').
    """
    m = _RESN_RESI.match(subst_name or "")
    if m:
        return m.group(1), m.group(2)
    return (subst_name or "LIG"), str(subst_id)


def _element_from_sybyl(sybyl_type):
    """SYBYL type -> element ('C.ar'->'C', 'N.am'->'N', 'Cl'->'Cl')."""
    return normalize_element((sybyl_type or "").split(".")[0])


def _parse_set_block(lines, start_idx):
    """Collect integers from a SET data block starting after its header line.

    The first integer is a count; the rest are atom serials. Continuation lines
    end with a backslash. Returns (serials_set, next_line_index).
    """
    ints = []
    i = start_idx
    while i < len(lines):
        raw = lines[i]
        stripped = raw.rstrip()
        cont = stripped.endswith("\\")
        body = stripped[:-1] if cont else stripped
        for tok in body.split():
            try:
                ints.append(int(tok))
            except ValueError:
                pass
        i += 1
        if not cont:
            break
    serials = set(ints[1:]) if ints else set()  # drop leading count
    return serials, i


def _parse_one_molecule(block_lines, source_file):
    section = None
    atom_rows = []  # (serial, name, x, y, z, sybyl, subst_id, subst_name)
    bonds = []  # (a1, a2)
    chain_by_subst = {}  # subst_id -> chain
    group_subst = {}  # subst_id -> 'RESIDUE'/'GROUP'
    group_info = {}  # subst_id -> raw info string (for GROUP naming)
    ccdc_ligand = None
    ccdc_receptor = None

    i = 0
    while i < len(block_lines):
        line = block_lines[i]
        s = line.strip()
        if s.startswith("@<TRIPOS>"):
            section = s[len("@<TRIPOS>") :].strip().upper()
            i += 1
            continue
        if not s or s.startswith("#"):
            i += 1
            continue

        if section == "ATOM":
            t = s.split()
            if len(t) >= 6:
                serial = int(t[0])
                name = t[1]
                x, y, z = float(t[2]), float(t[3]), float(t[4])
                sybyl = t[5]
                subst_id = int(t[6]) if len(t) >= 7 else 0
                subst_name = t[7] if len(t) >= 8 else ""
                atom_rows.append((serial, name, x, y, z, sybyl, subst_id, subst_name))
            i += 1
            continue

        if section == "BOND":
            t = s.split()
            if len(t) >= 4:
                bonds.append((int(t[1]), int(t[2])))
            i += 1
            continue

        if section == "SUBSTRUCTURE":
            t = s.split()
            # subst_id subst_name root_atom type [dict chain sub_type ...]
            if len(t) >= 4:
                sid = int(t[0])
                stype = t[3].upper()
                group_subst[sid] = stype
                if stype == "RESIDUE" and len(t) >= 6:
                    chain_by_subst[sid] = t[5]
                elif stype == "GROUP":
                    # The GROUP info field can contain spaces (e.g.
                    # 'B1 ZINC218266525_final'); it runs from token 5 up to the
                    # '****' status marker.
                    info_tokens = t[5:]
                    if "****" in info_tokens:
                        info_tokens = info_tokens[: info_tokens.index("****")]
                    group_info[sid] = " ".join(info_tokens)
            i += 1
            continue

        if section == "SET":
            up = s.upper()
            if up.startswith("CCDC_LIGAND"):
                ccdc_ligand, i = _parse_set_block(block_lines, i + 1)
                continue
            if up.startswith("CCDC_AMINOACID"):
                ccdc_receptor, i = _parse_set_block(block_lines, i + 1)
                continue
            i += 1
            continue

        i += 1

    # build Atom objects
    atoms = []
    atom_by_serial = {}
    for idx, (serial, name, x, y, z, sybyl, subst_id, subst_name) in enumerate(
        atom_rows
    ):
        elem = _element_from_sybyl(sybyl)
        resn, resi = _split_resn_resi(subst_name, subst_id)
        chain = chain_by_subst.get(subst_id, "")
        a = Atom(
            idx=idx,
            elem=elem,
            name=name,
            resn=resn,
            resi=resi,
            chain=chain,
            coord=(x, y, z),
            fcharge=0,  # mol2 stores partial charges only; no formal charge
            serial=serial,
            subst_id=subst_id,
        )
        atoms.append(a)
        atom_by_serial[serial] = a

    for a1, a2 in bonds:
        u, v = atom_by_serial.get(a1), atom_by_serial.get(a2)
        if u is not None and v is not None:
            u.neighbors.append(v)
            v.neighbors.append(u)

    pose = ParsedPose(
        atoms=atoms,
        atom_by_serial=atom_by_serial,
        fmt="mol2",
        source_file=source_file,
        ccdc_ligand=ccdc_ligand,
        ccdc_receptor=ccdc_receptor,
        group_subst=group_subst,
        sol=parse_sol(source_file),
    )
    pose.group_info = group_info  # attach raw GROUP info strings (ligand naming)
    return pose


def parse_mol2(path):
    """Parse a .mol2 file into a list of ParsedPose (one per MOLECULE block)."""
    with open(path, "r", errors="replace") as fh:
        text = fh.read()

    # split into MOLECULE blocks, keeping the marker with each block
    marker = "@<TRIPOS>MOLECULE"
    chunks = text.split(marker)
    poses = []
    for chunk in chunks:
        if not chunk.strip():
            continue
        block = (marker + chunk).splitlines()
        pose = _parse_one_molecule(block, path)
        if pose.atoms:
            pose.pose_index = len(poses)
            poses.append(pose)
    return poses
