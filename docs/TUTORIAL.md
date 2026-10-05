# DockLens tutorial

This tutorial walks through DockLens 1.4.1 from installation to exported
evidence. It covers docking poses and saved molecular-dynamics (MD) frames,
explains every control of the main window, and ends with troubleshooting. For a
worked virtual-screening campaign with published numbers, see
[tutorial-2m5d-screening.md](tutorial-2m5d-screening.md), which is pinned to
release 1.0.0.

Contents

1. [What DockLens does](#1-what-docklens-does)
2. [Install](#2-install)
3. [Core concepts](#3-core-concepts)
4. [Quick start with the bundled fixture](#4-quick-start-with-the-bundled-fixture)
5. [Loading structures](#5-loading-structures)
6. [Choosing the chemistry](#6-choosing-the-chemistry)
7. [Running detection and reading the workspaces](#7-running-detection-and-reading-the-workspaces)
8. [Key residues and recalculation](#8-key-residues-and-recalculation)
9. [Docking campaigns](#9-docking-campaigns)
10. [Molecular-dynamics frames](#10-molecular-dynamics-frames)
11. [Comparing two systems](#11-comparing-two-systems)
12. [Exporting evidence](#12-exporting-evidence)
13. [Projects](#13-projects)
14. [Command-line and scripted use](#14-command-line-and-scripted-use)
15. [Troubleshooting](#15-troubleshooting)

## 1. What DockLens does

DockLens detects non-covalent intermolecular interactions between a receptor
and a ligand in `.mol2`, `.pdb` and `.pdbqt` files. It keeps the atom-level
evidence for every contact and summarizes it as residue profiles, interaction
fingerprints, pose families, comparisons and auditable tables. It does not
depend on PyMOL, RDKit or Open Babel.

It detects 18 interaction types: `hbond`, `carbon_hbond`, `saltbridge`,
`attractive_charge`, `charge_repulsion`, `pipi`, `pication`, `pialkyl`,
`pi_sigma`, `alkyl`, `halogen`, `metal`, `water_bridge`, `pi_sulfur`,
`pi_anion`, `pi_donor_hbond`, `pi_lone_pair` and `chalcogen`.

## 2. Install

**Prebuilt executable (no Python).** Download the file for your system from the
[Releases](https://github.com/amgoncalvesusp/docklens/releases) page:

| System | File | Run |
|---|---|---|
| Windows, portable | `DockLens-windows-x86_64.exe` | double-click |
| Windows, installer | `DockLens-setup-windows-x86_64.exe` | run the installer |
| Linux | `DockLens-linux-x86_64` | `chmod +x DockLens-linux-x86_64 && ./DockLens-linux-x86_64` |

**From source.** Python 3.9 or newer is required.

```bash
git clone https://github.com/amgoncalvesusp/docklens.git
cd docklens
pip install -r requirements.txt
python -m docklens.app
```

**Check the installation.** The self-check writes and reads back a small XLSX
workbook and exits with status 0 when it succeeds:

```bash
python -m docklens.app --self-check
```

Memory use, timing and hardware guidance are in
[SYSTEM_REQUIREMENTS.md](SYSTEM_REQUIREMENTS.md).

## 3. Core concepts

**Observation.** One pose (docking) or one saved frame (MD). Every observation
has a stable `pose_id`, and every uploaded file has a stable `source_id`. Chart
labels can be changed freely; these identifiers never change.

**Receptor and ligand.** Each input must contain a receptor and a ligand.
DockLens separates them with a fixed priority:

1. MOL2 `CCDC_LIGAND` and `CCDC_AMINOACID` tags.
2. The MOL2 `SUBSTRUCTURE` record whose type is `GROUP`, taken as the ligand.
3. PDB or PDBQT `HETATM` records (excluding water and metals) as ligand, and
   `ATOM` records as receptor.
4. A fallback to the smallest connected component or chain. DockLens always
   asks for confirmation before applying it.
5. A manual override, which is always available.

**Scientific profile versus analysis view.** These are two different controls.

* The *Scientific profile* chooses the detector geometry and is applied when
  you press **Run detection**.
* The *Profile* selector (Complete or Discovery Studio-like) is a reversible
  filter applied after detection. Complete keeps every detected interaction;
  Discovery Studio-like keeps the Discovery Studio families and rejects salt
  bridges longer than 4.0 Å. Switching it never discards data.

For new comparisons, choose the detector profile directly and normally keep the
Complete view.

**Counting unit.** In the Residues workspace, a channel counts at most once per
observation, receptor residue and interaction type. Several atom pairs to the
same residue do not inflate prevalence. Stacked bars can exceed 100% because one
residue can be contacted through more than one interaction type in the same
pose.

## 4. Quick start with the bundled fixture

The repository ships a three-atom fixture, `tests/fixtures/two_pose_complex_sol7.pdb`,
that is small enough to check by hand. It is a test fixture and not scientific
data. It has two `MODEL` records. In model 1 the ligand nitrogen lies 2.9 Å from
the oxygen of SER70; in model 2 it lies 9.0 Å away.

1. Start DockLens.
2. Click **Open complex files** and select `two_pose_complex_sol7.pdb`.
3. Leave **Scientific profile** at *PLIP (legacy default)*, **Evidence** at
   *Docking poses* and **Profile** at *Complete*.
4. Type `SER70` in **Key residues** and press Enter.
5. Click **Run detection**.

Expected result (counts verified by running the detection core on this fixture;
the interface calls the same code):

| Scientific profile | Poses | Interaction rows | Why |
|---|---|---|---|
| PLIP (legacy default) | 2 | 2 | Pose 1 contains one N···O pair. The legacy profile reports it twice, once with each atom as donor and as acceptor. Pose 2 has no contact. |
| DSV-like | 2 | 1 | The DSV-like profile deduplicates semantically equivalent atom pairs, so the same pair is counted once. |

Open **Tables**, then **Detail**. Each row shows the ligand and receptor
endpoint, atom names and serials, the distance (2.9 Å), the `pose_id` and
whether the receptor residue is a key residue. Open **Summary** to see that pose
2 is present with zero interactions: observations without contacts are never
dropped from denominators.

Change **Scientific profile** to *DSV-like* and press **Run detection** again
to see the row count fall from 2 to 1. This is the simplest demonstration of why
the profile is recorded with every export.

## 5. Loading structures

The main window has three loading buttons.

| Button | Use it when |
|---|---|
| **Protein + ligands** | The receptor and the ligand poses are separate files, for example one receptor and the multi-model output of a docking run. |
| **Open complex files** | Each file already contains receptor and ligand. |
| **Open folder** | A whole campaign is in one directory of complex files. The folder is scanned recursively. |

**Protein + ligands** opens a dialog with two tabs, *Protein + ligands* and
*Complex files*. Choose one receptor (`PDB`, `PDBQT` or `MOL2`), add one or more
ligand files or folders, and use **Add another receptor group** when different
ligand sets belong to different receptors. **Remove selected** removes an entry.
Reopening the dialog appends files within each group.

A file larger than 256 MB is refused and reported in the *Input QC* sheet.
Invalid or unparsable files are recorded there with their error and do not stop
the other files.

## 6. Choosing the chemistry

Set these controls before **Run detection**.

| Control | Options | Notes |
|---|---|---|
| **Chart labels** | Ligand name (recommended); Uploaded file name; Pose / frame index | Presentation only. It can be changed after detection without rerunning. Repeated names receive a pose or frame suffix. |
| **Scientific profile** | PLIP (legacy default); LUNA 0.14 defaults; DSV-like (Discovery Studio 2024 defaults); LUNA × DSV conservative | Chooses the detector geometry. See the table below. |
| **Evidence** | Docking poses; Molecular dynamics frames | Changes labels and enables the saved-frame statistics. It does not reinterpret the raw rows. |
| **Profile** | Complete; Discovery Studio-like | Reversible post-detection view (Section 3). |

The four scientific profiles in short (the full parameter table is in the
README):

| Profile | When to choose it |
|---|---|
| PLIP (legacy default) | You need results that match earlier DockLens versions or the PyMOL plug-in's PLIP geometry. It is permissive and can report the same atom pair in both donor and acceptor roles. |
| LUNA 0.14 defaults | You want LUNA's published defaults. It requires explicit donor hydrogens. |
| DSV-like | You want parity with the criteria exposed by Discovery Studio Visualizer 2024, with semantic deduplication. |
| LUNA × DSV conservative | You want the union of the represented families using the stricter limit wherever both systems define one. |

The LUNA and DSV profiles claim parameter compatibility, not identical output,
with those programs.

## 7. Running detection and reading the workspaces

Click **Run detection**. Detection runs on the loaded files; the status area
reports the number of poses, interactions and errors. The navigation rail on the
left opens four workspaces.

### Residues

A stacked bar chart of residue by interaction type and a pose or frame barcode.
Use the **DockLens Lens** selector to pick a residue. The inspector shows its
frequency (docking) or occupancy (MD) by channel and stays synchronized with the
tables.

### Fingerprint

* **Interaction comparison heatmap.** Choose **Compare rows**
  (*Ligands / uploaded files* or *Individual poses*), **Columns**
  (*Residue × interaction type* or *Residue (any interaction)*) and **Top
  features** (20, 40, 80 or All). Each row has its own denominator, and
  zero-contact observations are included. The status line states how many rows
  and features are shown and whether a limit was applied.
* **Pose families.** Binary fingerprints are compared with Jaccard/Tanimoto
  similarity and grouped by complete-link clustering. Set **Similarity
  threshold** to control the grouping. With more than 300 observations, the model
  is trained on a disclosed, evenly spaced sample, and patterns outside the
  trained threshold are labeled `OUTLIER`.
* **Representatives and defining contacts.** The medoid observation of each
  family and the contacts that define it.
* **Export view.** Selects which Fingerprint figure **Export figure** writes
  (Section 12).

### Compare

Differences between System B and System A with independent denominators. See
Section 11.

### Tables

The auditable tables behind every chart: **Summary**, **Key Residue Coverage**
and **Detail**. Click a column header to sort. Use the filters described next.

**Filters.** *Search* matches residue, file or ligand text. *Key residues only*
restricts the Detail table to the key residues. The *Interaction types* check
boxes filter the Detail table by type. Filters change what is displayed and
what a filtered export contains; they never change the stored result.

**Chart scope — ligand/file A** restricts every analytical chart to one uploaded
ligand source (selected by its internal `source_id`, so two files that both
contain a ligand labeled `LIG` stay distinct). *All ligands / uploaded files*
restores the pooled view, which is weighted by the number of observations per
ligand.

## 8. Key residues and recalculation

Key residues mark positions you care about, for example catalytic residues.

1. Type them in **Key residues** (`SER70; LYS73; GLU166`; spaces, commas,
   semicolons and line breaks all work), or tick them in **Pick key residues
   from detected protein**, which lists the residues found in the receptor. Use
   the filter box above the list to find a residue.
2. DockLens reports identifiers that are invalid, unmatched or ambiguous
   between chains. It never discards them silently.
3. Press **Recalculate** to apply the changes to both systems without repeating
   detection, or **Discard changes** to revert. Pending changes must be
   recalculated or discarded before you save or export.

Key-residue counts in the Summary are raw interaction-pair counts. For rankings
based on conserved positions, use the `distinct_key_residue_count` and
`key_residue_coverage` fields of the Key Residue Coverage table instead.

## 9. Docking campaigns

A typical campaign, from a folder of poses to a figure:

1. **Open folder** (complex files) or **Protein + ligands** (one receptor plus a
   directory of ligand files).
2. Set **Chart labels** to *Ligand name*, **Scientific profile** to the profile
   you will report, **Evidence** to *Docking poses*.
3. **Run detection** and check the error count. A non-zero count means files
   failed to parse; the *Input QC* sheet of an XLSX export says why.
4. In **Residues**, read the prevalence of each residue. In **Fingerprint**,
   inspect the heatmap and the pose families, and open the representative of the
   largest family in your molecular viewer.
5. Export the figure bundle and the workbook (Section 12) and save the project
   (Section 13).

Docking poses are a population, not a time series. DockLens never applies
transition or bootstrap terminology to them.

## 10. Molecular-dynamics frames

Switch **Evidence** to *Molecular dynamics frames* before detection when the
inputs are saved MD frames. Then:

1. Set **Saved-frame step (ns)** to the time between consecutive saved frames
   before interpreting any duration.
2. Frames are treated as one contiguous series by default. For several replicas
   or for non-consecutive frames, press **Load trajectory map** and choose a CSV
   with the columns `observation_id`, `replica_id` and `frame_index`, plus an
   optional `time_ns`. The map must list exactly the analyzed observations. With
   it, episodes, transitions and resampling respect replica boundaries and
   frame gaps.
3. In **Fingerprint**, the four tabs are **Population** (state populations),
   **Timeline** (the interaction state of each frame in order), **Transitions**
   (observed counts and probabilities) and **Confidence**. For docking data the
   second tab is called **Pose order**.
4. Press **Compute confidence intervals** for a circular moving-block bootstrap
   of occupancy. Choose the number of resamples (**Bootstrap**) and the **Block**
   length. Frames next to each other are resampled together because they are not
   independent. The block length, seed, confidence level, number of resamples and
   any short-trajectory warning are stored with the plotted rows.

Transitions are descriptive and are not a validated Markov model. DockLens reads
saved frames as separate structure files and has no native trajectory reader;
convert a trajectory to per-frame PDB or MOL2 files first. The MD path has been
exercised on synthetic frame series by the automated tests and has not yet been
validated on a real trajectory.

## 11. Comparing two systems

Press **Load system B** to load a second set of inputs. In **Compare**, set the
role of each system (*Docking* or *MD frames*) with **A role** and **B role**.

* The difference view shows B minus A, each with its own denominator.
* A docking A to MD B comparison adds a retention analysis that labels each
  contact `retained`, `intermittent`, `lost` or `gained`.
* For MD, the block bootstrap can add intervals for the B minus A occupancy
  difference.

System A and System B keep independent chart scopes.

## 12. Exporting evidence

**Tables.** Choose all interactions or only the current filtered interactions,
then use **CSV** or **XLSX**.

* CSV writes `<prefix>_summary.csv`, `<prefix>_key_residue_coverage.csv` and
  `<prefix>_detail.csv`.
* The XLSX workbook holds `Summary`, `Ligand Index`, `Residue Matrix`,
  `Key Residue Coverage`, `Detail`, `Parameters` and `Input QC`. `Detail` is the
  canonical atom-level source from which the matrix and the filtered summaries
  are derived. `Parameters` records the DockLens version, profile, cutoffs,
  export filters and the matched, unmatched and ambiguous key residues.
  The matrix can store interaction counts or binary presence values.
* Every analyzed pose stays in the Summary and the Matrix, with zero counts when
  no interaction survives a filter. "All interactions" always exports the
  complete result.
* Text taken from input files is neutralized before writing so that it cannot
  become a spreadsheet formula.

**Figures.** Pick the figure in **Export view** (Interaction fingerprint,
Interaction comparison heatmap, Similarity matrix, State populations, State
timeline / pose order, State transitions or Confidence intervals) and press
**Export figure**. It writes a bundle: the figure (PNG, SVG, PDF or TIFF), a CSV
with the exact rows drawn, and a JSON manifest with the profile, view, counting
unit, denominator and version. **Export all TIFFs** writes every chart available
in the current state at 600 dpi with a per-chart manifest and a collection
manifest. Quote the manifest values in your methods section.

## 13. Projects

**Save project** stores the complete result before the Complete or Discovery
Studio-like view is applied, the active settings, SHA-256 digests of the source
files and a methods record, as a `.docklens` file. **Open project** resumes
without repeating detection. If an external structure file has changed or is
missing, DockLens reports it and disables rerunning, while the verified cached
evidence stays inspectable. **Reset** clears the session so that you can start a
new analysis.

## 14. Command-line and scripted use

The detection core runs without a display and can be called from Python. This
example reproduces the quick start of Section 4:

```python
from docklens import batch_runner

result = batch_runner.run(
    ["tests/fixtures/two_pose_complex_sol7.pdb"],
    key_residues=["SER70"],
    hbond_preset="dsv",   # "plip", "luna", "dsv" or "luna_dsv"
)
print(len(result.summaries), "poses,", len(result.details), "interaction rows")
for row in result.details:
    print(row.pose_id, row.interaction_type, row.receptor.resname,
          row.receptor.resseq, row.distance_A)
```

With the DSV-like profile this prints two poses and one interaction row.
Use `run_paired(receptor, poses_path, ...)` for one receptor plus a pose file or
folder. Set `QT_QPA_PLATFORM=offscreen` to run the test suite or other
Qt-dependent code on a machine without a display:

```bash
pip install -r requirements-dev.txt
QT_QPA_PLATFORM=offscreen pytest -q
```

## 15. Troubleshooting

| Symptom | Likely cause and fix |
|---|---|
| The status area reports errors. | One or more files failed to parse. Export the workbook and read the `Input QC` sheet for the file, the resolution method and the error text. |
| DockLens asks you to confirm a ligand. | No tag, `GROUP` record or `HETATM` block identified the ligand, so the smallest component was proposed. Confirm it only if it is correct, or use the manual override. |
| A pose has no interactions. | It is kept on purpose. Zero-contact observations stay in every denominator. Stricter profiles (DSV-like, LUNA 0.14) report fewer contacts than PLIP, and LUNA 0.14 needs explicit donor hydrogens, so a pose can lose all its contacts when you change the profile. |
| Two profiles give different counts. | Expected. The profiles use different cutoffs and deduplication. Report the profile name with every number. |
| A residue you typed is not applied. | Read the key-residue status line. The identifier may be invalid, absent from the receptor, or ambiguous between chains (write the chain, for example `SER70A`). Press **Recalculate** afterward. |
| Save or export is blocked. | Key-residue changes are pending. Press **Recalculate** or **Discard changes**. |
| The heatmap shows fewer rows or columns than expected. | A top-feature limit or the hard cell budget was applied. The status line and the manifest state the reduction; choose **All** in **Top features** for the complete set. |
| Reopening a project disables rerunning. | A source structure moved or changed, so its SHA-256 digest no longer matches. The cached evidence remains available; restore the files to rerun. |
| Counts differ from another program on the same poses. | The two programs apply different geometric rules and perceive chemistry differently. In the BSB 2026 case study, DockLens and Discovery Studio reported 1,077 and 1,079 (pose, residue, family) triples but disagreed on 464 of them. Compare on a stated unit, not on totals. |

Questions and bug reports: open an issue at
<https://github.com/amgoncalvesusp/docklens/issues>.
