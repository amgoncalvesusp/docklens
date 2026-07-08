"""
entity_resolver.py — decide which atoms are ligand vs. receptor.

Priority order (NEVER skip a step, NEVER guess silently):
  1. mol2 CCDC tags (CCDC_LIGAND / CCDC_AMINOACID serial lists) — source of truth.
  2. mol2 SUBSTRUCTURE GROUP record  -> that group's atoms are the ligand.
  3. pdb/pdbqt HETATM (excluding water and metal ions) = ligand; ATOM = receptor.
  4. Fallback: smallest connected component / smallest chain by residue count.
     Flagged `needs_confirmation=True` with a preview — the UI MUST confirm.
  5. Manual override (resolve_manual) — the UI selects receptor/ligand explicitly.

Water residues are always split into a separate `waters` list (used by the
water-mediated H-bond detector); metal ions stay on the receptor side.
"""

from __future__ import annotations

from dataclasses import dataclass, field

from .interaction_core import _METALS, _WATER_RESN


@dataclass
class Resolution:
    receptor_atoms: list
    ligand_atoms: list
    waters: list
    ligand_id: str
    method: str  # 'ccdc'|'group'|'hetatm'|'fallback'|'manual'
    needs_confirmation: bool = False
    preview: str = ""
    warnings: list = field(default_factory=list)


def _split_waters(atoms):
    waters, rest = [], []
    for a in atoms:
        if a.resn.upper() in _WATER_RESN:
            waters.append(a)
        else:
            rest.append(a)
    return rest, waters


def set_sides(res: Resolution) -> None:
    """Tag every atom's `.side` for this resolution (call before detection)."""
    for a in res.receptor_atoms:
        a.side = "receptor"
    for a in res.ligand_atoms:
        a.side = "ligand"
    for a in res.waters:
        a.side = "water"


def _ligand_id_from_group_info(info: str, fallback: str) -> str:
    """Extract a readable ligand id from a CCDC GROUP info string.

    e.g. '979|B1 ZINC218266525_final|B1_ZINC218266525|mol2|1|dock47|1'
         -> 'B1_ZINC218266525'
    """
    if not info:
        return fallback
    parts = [p.strip() for p in info.split("|") if p.strip()]
    if len(parts) >= 3:
        return parts[2].replace(" ", "_")
    if len(parts) >= 2:
        return parts[1].replace(" ", "_")
    return fallback


def _ligand_id_from_atoms(ligand_atoms, fallback):
    if not ligand_atoms:
        return fallback
    tag = ligand_atoms[0].res_tag()
    return tag if tag and tag != "_" else fallback


def _connected_components(atoms):
    """Return connected components (by bond graph) as lists of atoms."""
    seen = set()
    comps = []
    for a in atoms:
        if a.idx in seen:
            continue
        stack = [a]
        comp = []
        seen.add(a.idx)
        while stack:
            cur = stack.pop()
            comp.append(cur)
            for nb in cur.neighbors:
                if nb.idx not in seen:
                    seen.add(nb.idx)
                    stack.append(nb)
        comps.append(comp)
    return comps


def _residue_count(atoms):
    return len({a.res_tag() for a in atoms})


def resolve(pose, import_stem=None):
    """Resolve one ParsedPose into a list of Resolution (usually one).

    A mol2 with several GROUP records yields one Resolution per ligand group.
    """
    fallback_id = import_stem or "ligand"
    group_info = getattr(pose, "group_info", {})

    # --- Priority 1: mol2 CCDC tags ---
    if pose.ccdc_ligand:
        lig_serials = pose.ccdc_ligand
        ligand = [a for a in pose.atoms if a.serial in lig_serials]
        if pose.ccdc_receptor:
            receptor = [a for a in pose.atoms if a.serial in pose.ccdc_receptor]
        else:
            receptor = [a for a in pose.atoms if a.serial not in lig_serials]
        receptor, waters = _split_waters(receptor)
        lig_group_ids = [sid for sid, t in pose.group_subst.items() if t == "GROUP"]
        info = group_info.get(lig_group_ids[0], "") if lig_group_ids else ""
        lig_id = _ligand_id_from_group_info(
            info, _ligand_id_from_atoms(ligand, fallback_id)
        )
        return [Resolution(receptor, ligand, waters, lig_id, "ccdc")]

    # --- Priority 2: mol2 SUBSTRUCTURE GROUP ---
    group_ids = [sid for sid, t in pose.group_subst.items() if t == "GROUP"]
    if group_ids:
        out = []
        for gid in group_ids:
            ligand = [a for a in pose.atoms if a.subst_id == gid]
            receptor = [a for a in pose.atoms if a.subst_id != gid]
            receptor, waters = _split_waters(receptor)
            lig_id = _ligand_id_from_group_info(
                group_info.get(gid, ""), _ligand_id_from_atoms(ligand, fallback_id)
            )
            out.append(Resolution(receptor, ligand, waters, lig_id, "group"))
        return out

    # --- Priority 3: pdb/pdbqt HETATM ---
    if pose.is_hetatm:
        ligand, receptor = [], []
        for a in pose.atoms:
            het = pose.is_hetatm.get(a.serial, False)
            is_water = a.resn.upper() in _WATER_RESN
            is_metal = a.elem in _METALS
            if het and not is_water and not is_metal:
                ligand.append(a)
            else:
                receptor.append(a)
        receptor, waters = _split_waters(receptor)
        if ligand:
            lig_id = _ligand_id_from_atoms(ligand, fallback_id)
            return [Resolution(receptor, ligand, waters, lig_id, "hetatm")]

    # --- Priority 4: fallback (smallest component/chain) — needs confirmation ---
    non_water, waters = _split_waters(pose.atoms)
    comps = _connected_components(non_water)
    if len(comps) >= 2:
        comps.sort(key=_residue_count)
        ligand = comps[0]
        receptor = [a for c in comps[1:] for a in c]
        basis = "smallest connected component"
    else:
        by_chain = {}
        for a in non_water:
            by_chain.setdefault(a.chain, []).append(a)
        if len(by_chain) >= 2:
            chains = sorted(by_chain.values(), key=_residue_count)
            ligand = chains[0]
            receptor = [a for c in chains[1:] for a in c]
            basis = "smallest chain"
        else:
            return [
                Resolution(
                    non_water,
                    [],
                    waters,
                    fallback_id,
                    "fallback",
                    needs_confirmation=True,
                    preview="Could not separate ligand from receptor "
                    "automatically. Use manual override.",
                    warnings=["no ligand candidate found"],
                )
            ]

    lig_res = sorted({a.res_tag() for a in ligand})
    preview = (
        "Fallback (%s): ligand candidate = %d atoms, %d residue(s) "
        "[%s ... %s]; chain(s) %s. Confirm before running."
        % (
            basis,
            len(ligand),
            len(lig_res),
            lig_res[0] if lig_res else "?",
            lig_res[-1] if lig_res else "?",
            ",".join(sorted({a.chain or "_" for a in ligand})),
        )
    )
    lig_id = _ligand_id_from_atoms(ligand, fallback_id)
    return [
        Resolution(
            receptor,
            ligand,
            waters,
            lig_id,
            "fallback",
            needs_confirmation=True,
            preview=preview,
        )
    ]


def resolve_manual(pose, ligand_serials):
    """Manual override: the given serials are the ligand, the rest receptor."""
    lig = set(ligand_serials)
    ligand = [a for a in pose.atoms if a.serial in lig]
    receptor = [a for a in pose.atoms if a.serial not in lig]
    receptor, waters = _split_waters(receptor)
    lig_id = _ligand_id_from_atoms(ligand, "ligand")
    return Resolution(receptor, ligand, waters, lig_id, "manual")
