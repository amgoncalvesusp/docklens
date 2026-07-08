"""
interaction_core.py — PyMOL-independent non-covalent interaction detection.

Ported verbatim (geometry + cutoffs unchanged) from the PyMOL plugin
``interactions_plugin.py``. The ONLY differences from the plugin are:

  * ``Atom`` is built from plain fields (not a PyMOL chempy atom) and gains
    parser/UI bookkeeping fields (serial, subst_id, side).
  * ``Atom.res_sele()`` (the only PyMOL-coupled method) is removed; ``res_tag``
    and ``label`` are kept.
  * Interaction dicts carry the endpoint objects (``a_obj``/``b_obj``) instead
    of PyMOL selection strings (``a_sele``/``b_sele``), so the desktop tool can
    resolve ligand vs. receptor side per endpoint. No cutoff or geometric test
    is altered.

Results are therefore checkable against the PyMOL plugin.
"""

from __future__ import annotations

import numpy as np

# ===========================================================================
# Colour palette (Okabe-Ito) — shared with the PyMOL plugin for visual parity
# ===========================================================================

# Okabe & Ito (2008) colour-blind-safe qualitative palette. RGB 0-255.
_OKABE_ITO = {
    "black": (0, 0, 0),
    "orange": (230, 159, 0),
    "skyblue": (86, 180, 233),
    "bluishgreen": (0, 158, 115),
    "yellow": (240, 228, 66),
    "blue": (0, 114, 178),
    "vermillion": (213, 94, 0),
    "reddishpurple": (204, 121, 167),
}

# type -> (okabe-ito colour key, human label)
INTERACTION_COLORS = {
    "hbond": ("skyblue", "Hydrogen bond (conventional)"),
    "carbon_hbond": ("bluishgreen", "Carbon H-bond (weak, C-H...O/N)"),
    "saltbridge": ("vermillion", "Salt bridge / ionic"),
    "pipi": ("reddishpurple", "pi-pi stacking (sandwich/T-shaped)"),
    "pication": ("yellow", "pi-cation"),
    "pialkyl": ("orange", "pi-alkyl"),
    "alkyl": ("blue", "Alkyl-alkyl (hydrophobic)"),
    "halogen": ("black", "Halogen bond"),
    "metal": ("vermillion", "Metal coordination"),
    "water_bridge": ("skyblue", "Water-mediated H-bond"),
    "pi_sulfur": ("reddishpurple", "pi-sulfur"),
    "pi_anion": ("yellow", "pi-anion"),
}


def color_hex(itype: str) -> str:
    """Return the '#RRGGBB' Okabe-Ito colour for an interaction type."""
    okabe_key = INTERACTION_COLORS[itype][0]
    r, g, b = _OKABE_ITO[okabe_key]
    return "#%02X%02X%02X" % (r, g, b)


# ===========================================================================
# Geometric cutoffs — COPIED VERBATIM from interactions_plugin.py (do not edit)
# ===========================================================================
#
# Sources:
#   PLIP    = Salentin et al. 2015, PLIP config.py defaults.
#   Steiner = Steiner, Angew. Chem. Int. Ed. 2002 (weak H-bonds).
#   DS      = Discovery Studio interaction definitions (proprietary; ranges
#             approximate, flagged UNCERTAIN).
#
CUTOFFS = {
    "hbond_dist": 4.1,
    "hbond_angle": 100.0,
    "carbon_hbond_dist": 3.6,  # UNCERTAIN
    "carbon_hbond_angle": 120.0,  # UNCERTAIN
    "saltbridge_dist": 5.5,
    "pipi_dist": 5.5,
    "pipi_offset": 2.0,
    "pipi_angle_dev": 30.0,
    "pication_dist": 6.0,
    "pication_offset": 2.0,
    "pialkyl_dist": 5.0,  # UNCERTAIN
    "alkyl_dist": 4.0,
    "halogen_dist": 4.0,
    "halogen_angle": 135.0,
    "metal_dist": 3.0,
    "water_bridge_min": 2.5,
    "water_bridge_max": 4.1,
    "water_bridge_angle_min": 75.0,
    "water_bridge_angle_max": 140.0,
    "pi_sulfur_dist": 5.3,  # UNCERTAIN
    "pi_anion_dist": 5.0,  # UNCERTAIN
    "pi_anion_offset": 2.0,
}

_CUTOFF_DEFAULTS = dict(CUTOFFS)

VALID_TYPES = list(INTERACTION_COLORS.keys())

# Planarity tolerance for aromatic-ring perception (Angstrom). Verbatim.
_RING_PLANARITY_TOL = 0.15
_RING_ELEMENTS = {"C", "N", "O", "S"}


# ===========================================================================
# Small vector helpers — VERBATIM
# ===========================================================================


def _v(coord):
    return np.asarray(coord, dtype=float)


def _dist(a, b):
    return float(np.linalg.norm(_v(a) - _v(b)))


def _centroid(coords):
    return np.mean(np.asarray(coords, dtype=float), axis=0)


def _plane_normal(coords):
    """Best-fit plane normal via SVD of centred coordinates."""
    pts = np.asarray(coords, dtype=float)
    centred = pts - pts.mean(axis=0)
    _u, _s, vh = np.linalg.svd(centred)
    return vh[2]


def _planar_deviation(coords):
    """Max absolute distance of any atom from the best-fit plane."""
    pts = np.asarray(coords, dtype=float)
    n = _plane_normal(pts)
    return float(np.max(np.abs((pts - pts.mean(axis=0)).dot(n))))


def _angle_at(vertex, p1, p2):
    """Angle p1-vertex-p2 in degrees."""
    a = _v(p1) - _v(vertex)
    b = _v(p2) - _v(vertex)
    na, nb = np.linalg.norm(a), np.linalg.norm(b)
    if na < 1e-6 or nb < 1e-6:
        return 0.0
    cosv = np.clip(a.dot(b) / (na * nb), -1.0, 1.0)
    return float(np.degrees(np.arccos(cosv)))


def _plane_angle(n1, n2):
    """Acute angle (deg) between two plane normals."""
    cosv = abs(np.clip(_v(n1).dot(_v(n2)), -1.0, 1.0))
    return float(np.degrees(np.arccos(cosv)))


def _proj_offset(point, plane_point, normal):
    """Lateral offset of `point` from the axis through `plane_point`."""
    d = _v(point) - _v(plane_point)
    along = d.dot(_v(normal))
    perp = d - along * _v(normal)
    return float(np.linalg.norm(perp))


# ===========================================================================
# Atom / topology — ported (Atom rebuilt from plain fields; res_sele removed)
# ===========================================================================


class Atom(object):
    """Lightweight atom. Geometry-relevant fields match the plugin's Atom;
    `serial`, `subst_id` and `side` are parser/UI bookkeeping."""

    __slots__ = (
        "idx",
        "elem",
        "name",
        "resn",
        "resi",
        "chain",
        "segi",
        "coord",
        "fcharge",
        "neighbors",
        "serial",
        "subst_id",
        "side",
    )

    def __init__(
        self,
        idx,
        elem,
        name,
        resn,
        resi,
        chain="",
        segi="",
        coord=(0.0, 0.0, 0.0),
        fcharge=0,
        serial=None,
        subst_id=None,
    ):
        self.idx = idx
        self.elem = (elem or (name[:1] if name else "")).strip().capitalize()
        self.name = (name or "").strip()
        self.resn = (resn or "").strip()
        self.resi = str(resi).strip()
        self.chain = (chain or "").strip()
        self.segi = (segi or "").strip()
        self.coord = _v(coord)
        try:
            self.fcharge = int(fcharge)
        except Exception:
            self.fcharge = 0
        self.neighbors = []
        self.serial = serial
        self.subst_id = subst_id
        self.side = None  # 'receptor' | 'ligand' | 'water', set before classify

    # --- kept from the plugin ---
    def res_tag(self):
        chain = self.chain if self.chain else "_"
        return "%s%s%s" % (self.resn, self.resi, ("" if chain == "_" else chain))

    def label(self):
        return "%s_%s" % (self.res_tag(), self.name)


def _find_rings(atoms):
    """Topology-based aromatic-ring perception. VERBATIM."""
    adj = {a.idx: [n.idx for n in a.neighbors] for a in atoms}
    by_idx = {a.idx: a for a in atoms}
    rings = set()
    for start in adj:
        stack = [(start, (start,))]
        while stack:
            cur, path = stack.pop()
            for nb in adj[cur]:
                if nb == start and len(path) >= 5:
                    rings.add(frozenset(path))
                elif nb not in path and len(path) < 6:
                    stack.append((nb, path + (nb,)))

    ring_atoms = []
    for ring in rings:
        members = [by_idx[i] for i in ring]
        if not (5 <= len(members) <= 6):
            continue
        if any(a.elem not in _RING_ELEMENTS for a in members):
            continue
        if _planar_deviation([a.coord for a in members]) > _RING_PLANARITY_TOL:
            continue
        ring_atoms.append(members)
    return ring_atoms


class Ring(object):
    """Aromatic ring feature: centroid, plane normal, label. VERBATIM."""

    __slots__ = ("atoms", "centroid", "normal", "tag")

    def __init__(self, members, tag):
        self.atoms = members
        coords = [a.coord for a in members]
        self.centroid = _centroid(coords)
        self.normal = _plane_normal(coords)
        self.tag = tag


def _build_rings(atoms):
    """VERBATIM."""
    rings = _find_rings(atoms)
    per_res = {}
    for members in rings:
        per_res[members[0].res_tag()] = per_res.get(members[0].res_tag(), 0) + 1
    counters = {}
    out = []
    for members in rings:
        res = members[0].res_tag()
        if per_res[res] > 1:
            counters[res] = counters.get(res, 0) + 1
            tag = "%s_ring%d" % (res, counters[res])
        else:
            tag = "%s_ring" % res
        out.append(Ring(members, tag))
    return out


# ===========================================================================
# Chemical feature classification — ported (charged tuples carry the atom obj)
# ===========================================================================

_CATION_RES_ATOMS = {
    "LYS": ["NZ"],
    "ARG": ["NH1", "NH2", "NE"],
    "HIS": ["ND1", "NE2"],
    "HIP": ["ND1", "NE2"],
}
_ANION_RES_ATOMS = {
    "ASP": ["OD1", "OD2"],
    "GLU": ["OE1", "OE2"],
}
_HALOGENS = {"Cl", "Br", "I"}
_HB_ACCEPTOR_ELEMS = {"N", "O", "S", "F"}
_HB_DONOR_ELEMS = {"N", "O"}
_METALS = {"Na", "K", "Mg", "Ca", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Cd", "Hg"}
_WATER_RESN = {"HOH", "WAT", "H2O", "SOL", "TIP", "TIP3", "TIP4", "SPC", "DOD"}


def _h_neighbors(atom):
    return [n for n in atom.neighbors if n.elem == "H"]


def classify(atoms, rings, has_h):
    """Return a dict of feature lists for one molecular side. Ported.

    The only change from the plugin: charged-centre tuples carry the
    representative Atom object as the third element (was a PyMOL selection
    string), so the caller can recover residue/side metadata.
    """
    ring_atom_ids = set(a.idx for r in rings for a in r.atoms)

    donors = []
    carbon_donors = []
    acceptors = []
    cations = []  # (point, label, repr_atom)
    anions = []  # (point, label, repr_atom)
    halogens = []
    alkyl_carbons = []
    metals = []
    sulfurs = []

    # charged centres from formal charge
    for a in atoms:
        if a.fcharge > 0:
            cations.append((a.coord, a.label(), a))
        elif a.fcharge < 0:
            anions.append((a.coord, a.label(), a))

    # protein charged groups (grouped centres)
    grouped_cation = {}
    grouped_anion = {}
    for a in atoms:
        if a.resn in _CATION_RES_ATOMS and a.name in _CATION_RES_ATOMS[a.resn]:
            grouped_cation.setdefault((a.res_tag(), a.resn), []).append(a)
        if a.resn in _ANION_RES_ATOMS and a.name in _ANION_RES_ATOMS[a.resn]:
            grouped_anion.setdefault((a.res_tag(), a.resn), []).append(a)
    for (res, resn), grp in grouped_cation.items():
        pt = _centroid([x.coord for x in grp])
        lbl = "%s_guan" % res if resn == "ARG" else "%s_%s" % (res, grp[0].name)
        cations.append((pt, lbl, grp[0]))
    for (res, _resn), grp in grouped_anion.items():
        pt = _centroid([x.coord for x in grp])
        anions.append((pt, "%s_carboxyl" % res, grp[0]))

    # H-bond donors/acceptors, halogens, alkyl carbons, metals, sulfurs
    for a in atoms:
        if a.elem in _HB_ACCEPTOR_ELEMS:
            if not (a.elem == "N" and a.fcharge > 0):
                acceptors.append(a)
        if a.elem in _HB_DONOR_ELEMS:
            if has_h:
                hs = _h_neighbors(a)
                if hs:
                    donors.append((a, hs))
            else:
                donors.append((a, []))
        if a.elem == "C":
            if has_h:
                hs = _h_neighbors(a)
                if hs:
                    carbon_donors.append((a, hs))
            if a.idx not in ring_atom_ids:
                heavy = [n for n in a.neighbors if n.elem != "H"]
                if heavy and all(n.elem == "C" for n in heavy):
                    alkyl_carbons.append(a)
        if a.elem in _HALOGENS:
            cbonded = [n for n in a.neighbors if n.elem == "C"]
            if cbonded:
                halogens.append((a, cbonded[0]))
        if a.elem in _METALS:
            metals.append(a)
        if a.elem == "S":
            sulfurs.append(a)

    return {
        "donors": donors,
        "carbon_donors": carbon_donors,
        "acceptors": acceptors,
        "cations": cations,
        "anions": anions,
        "halogens": halogens,
        "alkyl": alkyl_carbons,
        "metals": metals,
        "sulfurs": sulfurs,
        "rings": rings,
    }


# ===========================================================================
# Interaction detectors — ported VERBATIM except a_sele/b_sele -> a_obj/b_obj
# ===========================================================================


def _mk(itype, subtype, a_label, b_label, a_point, b_point, a_obj, b_obj):
    return {
        "type": itype,
        "subtype": subtype,
        "a_label": a_label,
        "b_label": b_label,
        "a_point": a_point,
        "b_point": b_point,
        "a_obj": a_obj,
        "b_obj": b_obj,
    }


def _hbond_pairs(feat_a, feat_b, itype, dist_cut, angle_cut, has_h):
    donor_key = "donors" if itype == "hbond" else "carbon_donors"
    out = []
    for donor, hs in feat_a[donor_key]:
        for acc in feat_b["acceptors"]:
            if donor.idx == acc.idx:
                continue
            d = _dist(donor.coord, acc.coord)
            if d > dist_cut:
                continue
            if has_h and hs:
                best = max(_angle_at(h.coord, donor.coord, acc.coord) for h in hs)
                if best < angle_cut:
                    continue
            out.append(
                _mk(
                    itype,
                    "",
                    donor.label(),
                    acc.label(),
                    donor.coord,
                    acc.coord,
                    donor,
                    acc,
                )
            )
    return out


def detect_hbond(fa, fb, has_h):
    c = CUTOFFS
    res = _hbond_pairs(fa, fb, "hbond", c["hbond_dist"], c["hbond_angle"], has_h)
    res += _hbond_pairs(fb, fa, "hbond", c["hbond_dist"], c["hbond_angle"], has_h)
    return res


def detect_carbon_hbond(fa, fb, has_h):
    c = CUTOFFS
    res = _hbond_pairs(
        fa, fb, "carbon_hbond", c["carbon_hbond_dist"], c["carbon_hbond_angle"], has_h
    )
    res += _hbond_pairs(
        fb, fa, "carbon_hbond", c["carbon_hbond_dist"], c["carbon_hbond_angle"], has_h
    )
    return res


def detect_saltbridge(fa, fb):
    cut = CUTOFFS["saltbridge_dist"]
    out = []
    for cats, anis in ((fa["cations"], fb["anions"]), (fb["cations"], fa["anions"])):
        for cpt, clbl, catom in cats:
            for apt, albl, aatom in anis:
                if _dist(cpt, apt) <= cut:
                    out.append(
                        _mk("saltbridge", "", clbl, albl, cpt, apt, catom, aatom)
                    )
    return out


def detect_pipi(fa, fb):
    c = CUTOFFS
    out = []
    for r1 in fa["rings"]:
        for r2 in fb["rings"]:
            if _dist(r1.centroid, r2.centroid) > c["pipi_dist"]:
                continue
            offset = min(
                _proj_offset(r2.centroid, r1.centroid, r1.normal),
                _proj_offset(r1.centroid, r2.centroid, r2.normal),
            )
            if offset > c["pipi_offset"]:
                continue
            ang = _plane_angle(r1.normal, r2.normal)
            dev = c["pipi_angle_dev"]
            if ang <= dev:
                subtype = "sandwich"
            elif ang >= (90.0 - dev):
                subtype = "tshaped"
            else:
                continue
            out.append(
                _mk("pipi", subtype, r1.tag, r2.tag, r1.centroid, r2.centroid, r1, r2)
            )
    return out


def detect_pication(fa, fb):
    c = CUTOFFS
    out = []
    for rings, cats in ((fa["rings"], fb["cations"]), (fb["rings"], fa["cations"])):
        for r in rings:
            for cpt, clbl, catom in cats:
                if _dist(r.centroid, cpt) > c["pication_dist"]:
                    continue
                if _proj_offset(cpt, r.centroid, r.normal) > c["pication_offset"]:
                    continue
                out.append(_mk("pication", "", r.tag, clbl, r.centroid, cpt, r, catom))
    return out


def detect_pialkyl(fa, fb):
    cut = CUTOFFS["pialkyl_dist"]
    out = []
    for rings, alks in ((fa["rings"], fb["alkyl"]), (fb["rings"], fa["alkyl"])):
        for r in rings:
            for a in alks:
                if _dist(r.centroid, a.coord) <= cut:
                    out.append(
                        _mk("pialkyl", "", r.tag, a.label(), r.centroid, a.coord, r, a)
                    )
    return out


def detect_alkyl(fa, fb):
    cut = CUTOFFS["alkyl_dist"]
    out = []
    for a in fa["alkyl"]:
        for b in fb["alkyl"]:
            if _dist(a.coord, b.coord) <= cut:
                out.append(
                    _mk("alkyl", "", a.label(), b.label(), a.coord, b.coord, a, b)
                )
    return out


def detect_halogen(fa, fb):
    c = CUTOFFS
    out = []
    for hals, accs in (
        (fa["halogens"], fb["acceptors"]),
        (fb["halogens"], fa["acceptors"]),
    ):
        for x, cbonded in hals:
            for acc in accs:
                if _dist(x.coord, acc.coord) > c["halogen_dist"]:
                    continue
                if _angle_at(x.coord, cbonded.coord, acc.coord) < c["halogen_angle"]:
                    continue
                out.append(
                    _mk(
                        "halogen",
                        "",
                        x.label(),
                        acc.label(),
                        x.coord,
                        acc.coord,
                        x,
                        acc,
                    )
                )
    return out


def detect_metal(fa, fb):
    cut = CUTOFFS["metal_dist"]
    out = []
    for metals, accs in (
        (fa["metals"], fb["acceptors"]),
        (fb["metals"], fa["acceptors"]),
    ):
        for m in metals:
            for acc in accs:
                if _dist(m.coord, acc.coord) <= cut:
                    out.append(
                        _mk(
                            "metal",
                            "",
                            m.label(),
                            acc.label(),
                            m.coord,
                            acc.coord,
                            m,
                            acc,
                        )
                    )
    return out


def detect_pi_sulfur(fa, fb):
    cut = CUTOFFS["pi_sulfur_dist"]
    out = []
    for rings, sulfs in ((fa["rings"], fb["sulfurs"]), (fb["rings"], fa["sulfurs"])):
        for r in rings:
            for s in sulfs:
                if _dist(r.centroid, s.coord) <= cut:
                    out.append(
                        _mk(
                            "pi_sulfur", "", r.tag, s.label(), r.centroid, s.coord, r, s
                        )
                    )
    return out


def detect_pi_anion(fa, fb):
    c = CUTOFFS
    out = []
    for rings, anis in ((fa["rings"], fb["anions"]), (fb["rings"], fa["anions"])):
        for r in rings:
            for apt, albl, aatom in anis:
                if _dist(r.centroid, apt) > c["pi_anion_dist"]:
                    continue
                if _proj_offset(apt, r.centroid, r.normal) > c["pi_anion_offset"]:
                    continue
                out.append(_mk("pi_anion", "", r.tag, albl, r.centroid, apt, r, aatom))
    return out


def detect_water_bridge(fa, fb, waters):
    """Water-mediated H-bond. Ported. Emits two legs (partner--water) per bridge."""
    c = CUTOFFS

    def _partners(feat):
        seen, out = set(), []
        for donor, _hs in feat["donors"]:
            if donor.idx not in seen:
                seen.add(donor.idx)
                out.append(donor)
        for acc in feat["acceptors"]:
            if acc.idx not in seen:
                seen.add(acc.idx)
                out.append(acc)
        return out

    pa, pb = _partners(fa), _partners(fb)
    lo, hi = c["water_bridge_min"], c["water_bridge_max"]
    amin, amax = c["water_bridge_angle_min"], c["water_bridge_angle_max"]
    out = []
    for w in waters:
        near_a = [p for p in pa if lo <= _dist(w.coord, p.coord) <= hi]
        near_b = [p for p in pb if lo <= _dist(w.coord, p.coord) <= hi]
        for pai in near_a:
            for pbi in near_b:
                ang = _angle_at(w.coord, pai.coord, pbi.coord)
                if not (amin <= ang <= amax):
                    continue
                wlbl = w.label()
                for partner in (pai, pbi):
                    out.append(
                        _mk(
                            "water_bridge",
                            "",
                            partner.label(),
                            wlbl,
                            partner.coord,
                            w.coord,
                            partner,
                            w,
                        )
                    )
    return out


DETECTORS = {
    "hbond": lambda fa, fb, h: detect_hbond(fa, fb, h),
    "carbon_hbond": lambda fa, fb, h: detect_carbon_hbond(fa, fb, h),
    "saltbridge": lambda fa, fb, h: detect_saltbridge(fa, fb),
    "pipi": lambda fa, fb, h: detect_pipi(fa, fb),
    "pication": lambda fa, fb, h: detect_pication(fa, fb),
    "pialkyl": lambda fa, fb, h: detect_pialkyl(fa, fb),
    "alkyl": lambda fa, fb, h: detect_alkyl(fa, fb),
    "halogen": lambda fa, fb, h: detect_halogen(fa, fb),
    "metal": lambda fa, fb, h: detect_metal(fa, fb),
    "pi_sulfur": lambda fa, fb, h: detect_pi_sulfur(fa, fb),
    "pi_anion": lambda fa, fb, h: detect_pi_anion(fa, fb),
    # water_bridge handled separately (needs the water list)
}


# ===========================================================================
# High-level driver + endpoint metadata helpers
# ===========================================================================


def endpoint_side(obj):
    """'receptor' | 'ligand' | 'water' for an interaction endpoint (Atom/Ring)."""
    if isinstance(obj, Ring):
        return obj.atoms[0].side
    return obj.side


def endpoint_resid(obj):
    """Residue tag (resn+resi+chain) of an endpoint."""
    a = obj.atoms[0] if isinstance(obj, Ring) else obj
    return a.res_tag()


def endpoint_name(obj):
    """Atom name, or 'ring'/'ring2' for a ring endpoint."""
    if isinstance(obj, Ring):
        return obj.tag.rsplit("_", 1)[-1]
    return obj.name


def compute_interactions(receptor_atoms, ligand_atoms, waters=None, types=None):
    """Detect all requested interactions between receptor and ligand.

    Atoms must already have `.side` set ('receptor'/'ligand'/'water'). Returns a
    list of interaction dicts, each with an added 'dist' (endpoint separation).
    """
    waters = waters or []
    req = list(types) if types else list(VALID_TYPES)
    has_h = any(a.elem == "H" for a in receptor_atoms) or any(
        a.elem == "H" for a in ligand_atoms
    )

    feat_r = classify(receptor_atoms, _build_rings(receptor_atoms), has_h)
    feat_l = classify(ligand_atoms, _build_rings(ligand_atoms), has_h)

    inters = []
    for itype in req:
        if itype == "water_bridge":
            inters.extend(detect_water_bridge(feat_r, feat_l, waters))
        elif itype in DETECTORS:
            inters.extend(DETECTORS[itype](feat_r, feat_l, has_h))
    for it in inters:
        it["dist"] = _dist(it["a_point"], it["b_point"])
    return inters
