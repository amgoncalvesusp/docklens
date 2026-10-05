# Tutorial: analysing a virtual-screening campaign with DockLens

> **Version pin.** Every number and control name below belongs to DockLens
> **1.0.0** (tag `v1.0.0`, commit `b76fa11`), the version evaluated in the
> paper. Later releases changed the detection core and the interface, so the
> counts differ. To reproduce the numbers, install the 1.0.0 release. For the
> current release, see [TUTORIAL.md](TUTORIAL.md).

This walkthrough reproduces Section 4.3 of the BSB 2026 paper: 150 docking
solutions of 110 ZINC compounds docked into the metallo-β-lactamase BcII
(PDB 2M5D), analysed with the **DS-calibrated** chemistry profile and the
**Discovery Studio-like** view. Every count below was produced by DockLens
1.0.0 on that set; with another dataset the numbers differ, the sequence of
actions does not.

The `tests/` directory of this repository is **not** example data. It holds
small synthetic PDB/PDBQT/MOL2 fixtures used only by the automated test suite.

## What you need

* DockLens (the single-file executable from the Releases page, or
  `python -m docklens.app` from a source checkout).
* A folder of structure files. Each file must contain the receptor **and**
  one ligand pose. Accepted formats: `.mol2`, `.pdb`, `.pdbqt`. In the paper's
  set the ligand is identified by the MOL2 `SUBSTRUCTURE` `GROUP` record, so
  no manual receptor/ligand split is needed.
* Alternatively, one receptor file plus one multi-model ligand file (for
  example AutoDock Vina output) through **Open files**.

The 2M5D screening set used in the paper is available from the authors on
request. It is the same corpus that calibrated the Discovery Studio-like view
(Section 4.1 of the paper), so it is a demonstration set, not a blind test.

## Step by step

### 1. Load the campaign

**Open folder** → select the directory that holds the 150 `.mol2` files.
DockLens scans it recursively; the status bar reports the number of files
selected. Nothing is computed yet.

### 2. Choose the chemistry

In the context bar:

| Control | Set to | Why |
|---|---|---|
| **Criteria** | *DS-calibrated beta (explicit-H geometry)* | Explicit-hydrogen H-bond geometry and conservative donor/acceptor rules. The default, *PLIP*, is the permissive legacy geometry kept for cross-checking with the PyMOL plug-in. |
| **Profile** | *Discovery Studio-like* | Reversible post-detection view: keeps the DS interaction families and rejects salt bridges longer than 4.0 Å. The *Complete* result is never discarded. |
| **Evidence** | *Docking poses* | Poses are a population, not a time series. Switch to *Molecular dynamics* only for saved MD frames. |
| **Chart labels** | *Ligand name* | Presentation only; stable pose and source identifiers are unchanged. |

### 3. Run detection

**Run detection**. On a 16-thread desktop the 150 files (3,616 atoms each)
take about 16 s and peak at ~94 MB. Expected status line:

```
150 pose(s), 1211 interaction(s), 0 error(s) [DS-calibrated beta (explicit-H geometry)].
```

The complete result holds 1,220 rows; the Discovery Studio-like view removes
nine salt bridges longer than 4.0 Å, hence 1,211. If a file fails to parse,
the error count is non-zero and the **Input QC** sheet of the XLSX export
says why.

### 4. Read the residue profile (Residues workspace)

The header strip should read `21 residues · 150 poses · 61 residue/type
channels`. Each bar segment counts **once per pose × residue × interaction
type**, so several atom pairs to the same residue do not inflate it. Stacked
segments can exceed 100 % because one residue may be contacted through more
than one family in the same pose.

Expected top of the list (prevalence = fraction of poses with at least one
contact):

| Residue | Prevalence |
|---|---|
| HIS88 | 86.7 % |
| ALA89 | 80.0 % |
| ASP90 | 69.3 % |
| HIS149 | 66.7 % |
| HIS210 | 63.3 % |
| TRP59 | 55.3 % |
| ZN301 | 52.7 % |
| PHE34 | 46.0 % |
| ZN302 | 36.0 % |

Both catalytic zinc ions appear because receptor metals stay on the receptor
side. 53 poses coordinate both ions in the same pose.

Pick a residue in the **DockLens Lens** selector; the inspector on the right
shows its frequency by channel (for example `ALA89A · hbond 73.3 % (110/150)`)
and stays synchronised with the tables.

### 5. Triage the population (Fingerprint workspace)

* **Interaction comparison heatmap** — *Compare rows: Ligands / uploaded
  files*, *Columns: Residue × interaction type*, *Top features: 40*. One row
  per ligand file with its own denominator; zero-contact poses are kept. The
  status strip states row, feature and observation counts and whether the
  feature limit was applied.
* **Similarity and pose families** — Jaccard/Tanimoto similarity of binary
  fingerprints and complete-link clustering at threshold 0.65. Expected:
  **91 pose families, 60 of them singletons; largest family 7 poses (4.7 %)**;
  90 poses (60 %) fall in families with more than one member. With fewer than
  300 observations the whole set is used for training and no pose is labelled
  `OUTLIER`.
* The representative of the largest family is a solution of compound D22; its
  defining contacts are hydrogen bonds to ALA89, ASP90, HIS149 and HIS210.
  **Open representative structure** hands that file to your viewer.

### 6. Export the evidence

* **Export figure** (Fingerprint header, *Export view* selector) writes a
  bundle: `<name>.png` (and/or `.svg`, `.pdf`), `<name>_data.csv` with the
  exact rows drawn, and `<name>_manifest.json` recording profile, view,
  counting unit, denominator and version. Figures 1 and 2 of the paper are
  such bundles.
* **XLSX** writes the six-sheet workbook (Summary, Residue Matrix, Key Residue
  Coverage, Detail, Parameters, Input QC). *Detail* is the canonical
  atom-level evidence; *Parameters* records every cutoff and filter in force.
* **Save project** stores the complete pre-filter result, settings and SHA-256
  digests of the inputs as a `.docklens` file that reopens without
  re-detection.

### 7. (Optional) Draw the same interactions in PyMOL

Install the companion plug-in (`pymol_interactions_plugin.zip`, Plugin
Manager) and, for one pose:

```
load D22_sol18_2m5d.mol2, pose
detect_interactions not resn ZINC1, resn ZINC1, show_residues=1
interactions_parity_status
```

`resn ZINC1` is how PyMOL names the ligand substructure of these MOL2 files;
check with `iterate` if your files differ. The plug-in draws nine objects for
this pose — four hydrogen bonds, one carbon hydrogen bond, two metal contacts
and two π-alkyl contacts — and reports the active parity contract
(`docklens-dsv-2026.07`), meaning the dashes come from the same filtered set
that DockLens plotted.

## Reading the agreement with Discovery Studio

Section 4.1 of the paper compares these 150 poses with a Discovery Studio
annotation matrix on one unit, the (pose, residue, family) triple. The script
that reproduces Table 3 of the paper needs the annotation workbook and the pose
folder, neither of which is distributed with this repository; both, with the
script, are available from the authors on request.

## Where the numbers come from

| Quantity | Value | Where to see it |
|---|---|---|
| Poses / files | 150 | status bar after Run detection |
| Rows, complete → DS-like | 1,220 → 1,211 | Profile switch; Summary sheet |
| Residues / channels | 21 / 61 | Residues header strip |
| Pose families | 91 (60 singletons) | Fingerprint → families table |
| Metal rows | 133 | Detail sheet, filter `metal` |
