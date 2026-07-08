"""
End-to-end tests: mol2 acceptance cases + synthetic PDB/PDBQT + export.

Runnable directly (``python tests/test_pipeline.py``) or via pytest.
"""

import os
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from docklens import batch_runner as br  # noqa: E402
from docklens import export  # noqa: E402
from docklens.parser_pdb import parse_pdb  # noqa: E402
from docklens.parser_pdbqt import parse_pdbqt  # noqa: E402
from docklens.entity_resolver import resolve, set_sides  # noqa: E402
from docklens.interaction_core import compute_interactions  # noqa: E402

MOL2_A = r"C:\Users\adria\Downloads\1_peptide_ctxm_FDSVH_sol16.mol2"
MOL2_B = (
    r"C:\Users\adria\OneDrive\Dell_antigo_backup\Pós-doutorado\Pós doc IQ-UNESP"
    r"\AndreaChimarro\Proyecto\Docking\5l2fBestSolutions\B1_sol47_5l2f.mol2"
)


def test_mol2_acceptance():
    result = br.run([MOL2_A, MOL2_B])
    assert not result.pending, "CCDC files should not need confirmation"
    by_lig = {s.ligand_id: s for s in result.summaries}
    assert "1_FDSVH" in by_lig, list(by_lig)
    assert any("ZINC" in k for k in by_lig), list(by_lig)
    fdsvh = by_lig["1_FDSVH"]
    assert fdsvh.sol == 16, fdsvh.sol
    assert fdsvh.n_total_interactions > 0
    for d in result.details:
        assert d.ligand_atom, d
        assert d.receptor_residue, d
    return result


def test_key_residue_recompute():
    result = br.run([MOL2_A])
    some_res = result.details[0].receptor_residue
    br.recompute_key(result, [some_res])
    n_key = sum(1 for d in result.details if d.is_key_residue)
    assert n_key >= 1
    assert result.summaries[0].n_key_residue_interactions == n_key
    br.recompute_key(result, [])
    assert all(not d.is_key_residue for d in result.details)
    assert result.summaries[0].n_key_residue_interactions == 0


def test_export_roundtrip():
    result = br.run([MOL2_A])
    tmp = tempfile.mkdtemp()
    csvs = export.export_csv(result, os.path.join(tmp, "out"))
    assert all(os.path.exists(p) and os.path.getsize(p) > 0 for p in csvs), csvs
    xlsx = export.export_xlsx(result, os.path.join(tmp, "out.xlsx"))
    assert os.path.exists(xlsx) and os.path.getsize(xlsx) > 0

    import pandas as pd

    xls = pd.ExcelFile(xlsx)
    assert set(xls.sheet_names) == {"Summary", "Detail"}, xls.sheet_names


# --- synthetic PDB: SER backbone O (acceptor) + ligand N-H donor -> hbond ---
_PDB = """\
ATOM      1  N   SER A   1      -0.966   0.000   0.000  1.00  0.00           N
ATOM      2  CA  SER A   1       0.500   0.100   0.000  1.00  0.00           C
ATOM      3  C   SER A   1       1.100   1.400   0.000  1.00  0.00           C
ATOM      4  O   SER A   1       2.900   1.400   0.000  1.00  0.00           O
HETATM    5  N1  LIG B   1       5.700   1.400   0.000  1.00  0.00           N
HETATM    6  H1  LIG B   1       4.750   1.400   0.000  1.00  0.00           H
HETATM    7  C1  LIG B   1       6.300   2.700   0.000  1.00  0.00           C
END
"""


def test_pdb_synthetic():
    tmp = tempfile.mkdtemp()
    path = os.path.join(tmp, "synthetic.pdb")
    with open(path, "w") as fh:
        fh.write(_PDB)
    poses = parse_pdb(path)
    assert len(poses) == 1
    pose = poses[0]
    assert pose.is_hetatm[5] is True and pose.is_hetatm[1] is False
    res = resolve(pose)[0]
    assert res.method == "hetatm"
    assert {a.name for a in res.ligand_atoms} == {"N1", "H1", "C1"}
    set_sides(res)
    inters = compute_interactions(res.receptor_atoms, res.ligand_atoms, res.waters)
    types = {it["type"] for it in inters}
    assert "hbond" in types, types


# --- synthetic PDBQT: MODEL with Vina score + two atoms ---
_PDBQT = """\
MODEL 1
REMARK VINA RESULT:    -7.5      0.000      0.000
ROOT
HETATM    1  N1  LIG d   1       5.700   1.400   0.000  1.00  0.00    -0.350 N
HETATM    2  C1  LIG d   1       6.300   2.700   0.000  1.00  0.00     0.100 C
ENDROOT
TORSDOF 0
ENDMDL
"""


def test_pdbqt_synthetic_score():
    tmp = tempfile.mkdtemp()
    path = os.path.join(tmp, "pose_sol3.pdbqt")
    with open(path, "w") as fh:
        fh.write(_PDBQT)
    poses = parse_pdbqt(path)
    assert len(poses) == 1
    assert poses[0].score == -7.5, poses[0].score
    assert poses[0].sol == 3
    assert poses[0].atoms[0].elem == "N"


if __name__ == "__main__":
    test_mol2_acceptance()
    print("mol2 acceptance OK")
    test_key_residue_recompute()
    print("key-residue recompute OK")
    test_export_roundtrip()
    print("export roundtrip OK")
    test_pdb_synthetic()
    print("pdb synthetic OK")
    test_pdbqt_synthetic_score()
    print("pdbqt synthetic OK")
    print("\nALL PIPELINE TESTS PASSED")
