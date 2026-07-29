# DockLens

![DockLens](docklens/assets/docklens_logo.png)

*Visual Intermolecular Interaction Analytics*

Standalone desktop tool (PyQt5) that detects and audits non-covalent
intermolecular interactions in docking poses or saved molecular-dynamics
frames (`.mol2`, `.pdb`, `.pdbqt`). DockLens separates receptor from ligand,
preserves atom-level evidence and presents residue profiles, interaction
fingerprints, comparisons and sortable/filterable tables.

The legacy PLIP profile in `interaction_core.py` preserves the geometry and
cutoffs ported from the PyMOL plugin `interactions_plugin.py`. A separate DSV
profile adds chemistry-aware, corpus-calibrated rules without changing legacy
results. DockLens has **no PyMOL dependency**.

**Inventor:** Adriano Marques Gonçalves — Universidade de Araraquara (UNIARA).

## Interaction types (15)

`hbond`, `carbon_hbond`, `saltbridge`, `pipi` (sandwich/T-shaped), `pication`,
`pialkyl`, `pi_sigma`, `alkyl`, `halogen`, `metal`, `water_bridge`, `pi_sulfur`,
`pi_anion`, `pi_donor_hbond`, `pi_lone_pair`.
Colours use the same colour-blind-safe Okabe-Ito palette as the PyMOL plugin.

## Download (standalone, no Python needed)

Grab the executable for your OS from the [Releases](../../releases) page:

| OS | File | How to run |
|----|------|------------|
| Windows | `DockLens-windows-x86_64.exe` | double-click |
| Linux | `DockLens-linux-x86_64` | `chmod +x DockLens-linux-x86_64 && ./DockLens-linux-x86_64` |

Both are single-file bundles built with PyInstaller — no Python or dependencies
to install. Builds are produced by GitHub Actions (`.github/workflows/build.yml`)
on Windows and Ubuntu runners.

## Run from source

Requirements: Python 3.9+, `numpy`, `pandas`, `openpyxl`, `matplotlib`, `PyQt5`
(no RDKit/OpenBabel).

```
pip install -r requirements.txt
python -m docklens.app
```

Build your own executable:

```
pip install pyinstaller
pyinstaller --noconfirm --onefile --windowed --name DockLens run_docklens.py
```

## DockingHub integration

DockingHub can open a completed docking directly in DockLens. The applications
remain in separate processes so their Qt runtimes do not conflict. DockingHub
creates a `docklens-launch-v1` JSON manifest that explicitly identifies one
receptor and one multipose docking file, confines both paths to the DockingHub
project and records a SHA-256 digest for each input. DockLens validates the
contract before parsing any structure.

When the paired analysis finishes, DockLens atomically writes the
`docklens-result-v1` path declared by DockingHub. The response includes manifest
and input hashes, parameters, per-pose and per-type counts, interaction detail
and input QC. DockingHub validates that chain before adding the analysis to its
project and consolidated report.

The desktop route is:

```
DockLens.exe --manifest <project>/reports/docklens_launch_<run-id>.json
```

For packaging and diagnostic checks without opening the GUI:

```
DockLens.exe --check-manifest <manifest.json>
```

This prints a compact JSON summary and returns a nonzero exit code when the
manifest, hashes or paired analysis are invalid.

## DockLens 1.0 analytical atlas

The v1 interface has four peer workspaces:

- **Residues:** stacked residue/type bars and a pose/frame barcode. Every
  channel counts at most once per observation × receptor residue × interaction
  type, preventing multiple atom pairs from inflating prevalence or occupancy.
- **Fingerprint:** binary pose/frame fingerprints, Jaccard/Tanimoto similarity,
  an interaction comparison heatmap, deterministic complete-link pose
  families or interaction states, medoid observations, defining contacts,
  population charts and an ordered ribbon. The comparison heatmap can use
  ligand/uploaded-file rows with independently normalized frequency/occupancy,
  or individual pose/frame rows with binary presence.
  Above 300 observations, model training uses a disclosed evenly spaced
  sample; patterns outside the trained threshold remain labeled `OUTLIER`.
- **Compare:** System B minus System A differences with independent
  denominators, plus an explicit docking A → MD B retention analysis
  (`retained`, `intermittent`, `lost`, `gained`). MD comparisons can add
  independent block-bootstrap intervals for the B − A occupancy difference.
- **Tables:** the complete Summary, Key Residue Coverage and atom-level Detail
  views from earlier DockLens versions.

The **DockLens Lens** keeps the selected residue synchronized with its
frequency/occupancy evidence. Docking poses are never described as a temporal
series; switching to molecular dynamics changes labels and enables saved-frame
episode statistics without reinterpreting the raw interaction rows. In MD mode,
Fingerprint also provides a state timeline, observed transition
counts/probabilities and representative frames. By default, frames are treated
as one contiguous series. **Load trajectory map** accepts a CSV with
`observation_id`, `replica_id`, `frame_index` and optional `time_ns`; with that
map, transitions, episodes and resampling preserve declared replica boundaries
and frame gaps. Transitions are descriptive, not a validated Markov model.

MD occupancy intervals use a circular moving-block bootstrap so consecutive
saved frames are resampled together instead of being treated as independent.
The block length, seed, confidence level, number of resamples and any
short-trajectory warning are retained with the plotted rows. Docking data never
uses temporal transition or block-bootstrap terminology.

The **Chart scope — ligand/file** selectors distinguish uploaded ligand
sources by their internal `source_id`, even when different files contain the
same generic ligand label such as `LIG` or `RES1`. Selecting a source
recalculates every analytical chart from all its poses or frames, including
observations with zero surviving interactions. System A and System B have
independent scopes. **All ligands / uploaded files** retains the pooled view;
that view is weighted by observation count, so a ligand with more poses or
frames contributes more to the aggregate.

**Chart labels** is available before detection and selects the human-readable
identity used on analytical axes: detected ligand name, uploaded filename or
pose/frame index. Stable `pose_id` and `source_id` values remain unchanged and
are retained in exported source rows. Repeated names receive a pose/frame
suffix so multipose inputs stay unambiguous.

The **Interaction comparison heatmap** uses all loaded ligand/file groups in
its aggregate mode, even when another chart is focused on a single source.
Each row has its own pose/frame denominator and includes zero-contact
observations. Its individual mode follows the active chart scope. Columns can
represent residue × interaction type or binary presence of any interaction
with a residue; a visible top-feature limit keeps large analyses readable.
A hard cell budget protects the desktop from oversized projects. When it is
reached, rows are selected deterministically across the complete order and the
reduction is stated in the status and reproducibility manifest.

Analyses can be saved as versioned `.docklens` projects. A project stores the
complete immutable result before the Complete/Discovery Studio-like view is
applied, active settings, source SHA-256 digests and a methods record. Cached
result members and an explicit trajectory map are integrity checked when
reopened. Changed or missing external structures are reported and rerunning is
disabled while the verified cached evidence remains inspectable.

Every analytical chart can export a publication bundle containing PNG, SVG or
PDF output, the exact tidy CSV rows used to draw it and a JSON reproducibility
manifest with profile, criteria, denominator and counting unit. The established
CSV/XLSX export remains backward compatible.

The **Complete** and **Discovery Studio-like** views remain visible in the
global profile selector and are applied consistently to charts and tables.

## Using the app

1. **Open file(s)** or **Open folder** (recursive scan of `.mol2/.pdb/.pdbqt`).
2. (optional) set key residues — type them (`SER70; LYS73; GLU166`) **or** tick
   them from the checkbox list of the detected protein residues. Spaces,
   commas, semicolons and line breaks are accepted. DockLens reports invalid,
   unmatched and chain-ambiguous identifiers instead of silently discarding
   them.
3. Pick **Chart labels**: ligand name, uploaded filename or pose/frame index.
   This controls presentation only and can be changed later without rerunning.
4. Pick the **H-bond criteria** preset (see below).
5. Pick the **Analysis view**: Complete preserves every detected interaction;
   Discovery Studio-like is a reversible post-detection view that keeps the
   DSV interaction families and rejects salt bridges longer than 4.0 Å.
6. **Run detection**.
7. Explore **Residues**, **Fingerprint** and **Compare**. Choose whether the
   observations represent docking poses or saved molecular-dynamics frames.
   Load a separate System B when a differential or retention analysis is needed.
   Use **Chart scope — ligand/file A** and **System B** to inspect or compare
   individual uploaded ligands without removing other results from the tables.
   Fingerprint contains population, timeline, transition and confidence tabs.
   Set the saved-frame interval before interpreting MD durations. For multiple
   replicas or nonconsecutive saved frames, load the trajectory-map CSV for
   System A and, when comparing MD systems, for System B.
8. Sort by clicking a column; filter by interaction type, search text, or
   "key residues only". Edit key residues any time — counts recompute without
   re-running detection.
9. Choose the desired **Export view** in Fingerprint, then use **Export
   figure** to write that figure, its source rows and reproducibility manifest.
10. **Export CSV** or **Export XLSX**. Choose all interactions or the current
   filtered interactions; every analyzed pose remains in Summary/Matrix, with
   zero counts when no interaction survives the filter. XLSX matrices can use
   interaction counts or binary presence values. "All interactions" always
   exports the complete result; a filtered export records the selected analysis
   view in the Parameters sheet.
11. Use **Save project** to preserve cached evidence, settings and methods;
    use **Open project** to resume without redetection.
11. **Reset** clears everything to start a new analysis.

## H-bond criteria — DockLens and the strict profile

DockLens ships two H-bond profiles:

| Preset | Explicit-H geometry | Chemistry policy |
|--------|---------------------|------------------|
| **PLIP** (default) | D···A ≤ 4.1 Å; legacy angle ≥ 100° | Legacy behavior, retained for reproducibility. |
| **DS-calibrated beta** | Conventional: H···A ≤ 3.1 Å, D–H···A ≥ 90°, H–A–Y ≥ 90°. Carbon H-bond: H···A ≤ 3.0 Å, D–H···A ≥ 90°, H–A–Y ≥ 90°. | Uses retained MOL2 atom/bond types, explicit hydrogens and conservative donor/acceptor rules. |

The default (PLIP) is **deliberately permissive** — it matches the companion
PyMOL plugin so historical numbers remain cross-checkable. The strict profile
also excludes known non-acceptor nitrogens, protonated hydroxyl acceptors and
non-donor carbonyl/carboxylate oxygens. Carbon H-bonds require a polarized
carbon. Its H-bond rows export the participating hydrogen, H···A distance,
D–H···A angle, H–A–Y angle, inference basis and confidence.

The beta profile was recalibrated with 150 matched 2m5d docking poses annotated
in Discovery Studio. The visible annotation matrix contains 1,150 interactions
in eight families; DockLens v0.6.0 reports 1,151 in those same mapped families
(a 0.1% aggregate difference). At the residue-and-type level, the transparent
approximation has 0.658 precision/recall against that matrix. The most important
structural fix is preserving unassigned receptor metals: the calibrated run
recovers all 133 annotated metal contacts. The aromatic/hydrophobic classes
remain less exact and should still be reviewed in the auditable Detail sheet.

This calibration adds `pi_sigma`, `pi_donor_hbond` and `pi_lone_pair`; applies
semantic deduplication to `pialkyl`; and retains generic ligand identifiers from
the source filename instead of exporting every CCDC ligand as `RES1`. The
**Discovery Studio-like** analysis view does not alter or discard the complete
result. It is a transparent screening aid, not a claim that DockLens reproduces
every proprietary or user-customized Discovery Studio rule.

## Ligand vs. receptor resolution (fixed priority)

1. mol2 CCDC tags (`CCDC_LIGAND` / `CCDC_AMINOACID`) — source of truth.
2. mol2 `SUBSTRUCTURE` `GROUP` record → that group is the ligand.
3. pdb/pdbqt `HETATM` (excluding water & metals) = ligand; `ATOM` = receptor.
4. Fallback: smallest connected component / smallest chain — **asks for
   confirmation** (never applied silently).
5. Manual override always available.

## XLSX export schema

DockLens exports a versioned, six-sheet workbook:

| Sheet | Contents |
|-------|----------|
| `Summary` | One row per unambiguous pose, with pair counts, distinct key-residue coverage and per-type counts. |
| `Residue Matrix` | One row per pose and grouped columns for residue × interaction type. |
| `Key Residue Coverage` | Ranking-audit view separating raw atomic pairs, distinct conserved residues, conventional H-bond coverage and interaction diversity. |
| `Detail` | Auditable interaction endpoints, pose IDs, atom serials, participating hydrogen, H···A/DHA/HAY geometry and angular descriptors for π contacts. |
| `Parameters` | DockLens/schema version, preset, cutoffs, export filters and matched/unmatched/ambiguous key-residue identifiers. |
| `Input QC` | Status, resolution method, atom counts, warnings and parse errors. |

The matrix uses full residue identifiers such as `SER70A`, so chains do not
collide. **Count** mode stores the number of interactions; **Presence** mode
stores only `0` or `1`. `Detail` is the canonical source of truth from which the
matrix and filtered summaries are derived.

Each source receives a deterministic `source_id`; every pose/resolution receives
a unique `pose_id`. Water bridges are represented as one semantic
receptor-water-ligand interaction and therefore count once.

`n_key_residue_interactions` is a raw interaction-pair count. It must not be
used as a synonym for residue coverage. `distinct_key_residue_count` and
`key_residue_coverage` are the appropriate fields for rankings based on
conserved positions.

CSV writes `<prefix>_summary.csv`, `<prefix>_key_residue_coverage.csv` and
`<prefix>_detail.csv`. Text controlled by input files is neutralized before
CSV/XLSX writing so it cannot become an Excel formula.

## Modules

| File | Role |
|------|------|
| `interaction_core.py` | Ported detection core (Atom, Ring, classify, detect_*, CUTOFFS). No PyMOL. |
| `structures.py` | `ParsedPose` schema + covalent-radii bond inference. |
| `parser_mol2.py` / `parser_pdb.py` / `parser_pdbqt.py` | Format readers. |
| `entity_resolver.py` | Ligand/receptor split (priority above). |
| `batch_runner.py` | Input modes, detection driver, Summary/Detail rows, key-residue recompute. |
| `integration_manifest.py` | Validated, hashed and project-confined DockingHub launch contract. |
| `integration_result.py` | Atomic `docklens-result-v1` round-trip contract for DockingHub. |
| `results.py` | Immutable schema-v3 result, endpoint, QC and parameter contracts. |
| `result_analysis.py` | Immutable conserved-residue coverage and ranking-audit metrics. |
| `residue_keys.py` | Canonical residue-list parsing, validation and chain-aware matching. |
| `export_views.py` | Pure Summary/Detail/Coverage/Residue Matrix/Parameters/QC transformations. |
| `export.py` | Atomic CSV / XLSX writers with filtering and Okabe-Ito shading. |
| `analytics.py` | Consolidated prevalence/occupancy, fingerprints, similarity, clustering, episodes, comparison and retention. |
| `observation_series.py` | Explicit saved-frame order, time, replica boundaries and gap-aware transitions. |
| `dynamic_states.py` / `dynamic_plotting.py` | Complete-link pose families, MD states, representatives, timelines, observed transitions and publication artifacts. |
| `uncertainty.py` | Circular moving-block bootstrap for MD occupancy and independent B − A differences. |
| `project_session.py` / `project_controller.py` | Integrity-checked `.docklens` projects, cached results, provenance, methods and desktop restoration. |
| `analysis_tasks.py` | Background bootstrap execution without blocking the desktop UI. |
| `plotting.py` / `figure_export.py` | Publication figure builders and atomic figure/data/manifest bundles. |
| `analytics_widgets.py` / `main_window_ui.py` | Responsive analytical workspaces and desktop shell construction. |
| `main_window.py` / `app.py` | Desktop behavior and entry point. |

## Tests

```
pip install -r requirements-dev.txt
pytest --cov=docklens --cov-branch --cov-fail-under=80
```

Tests use repository-owned synthetic PDB/PDBQT/MOL2 fixtures and cover parsing,
pose identity, immutable key-residue recomputation, water bridges, filtered
exports, formula neutralization, the six-sheet XLSX schema and offscreen UI
flows. No test depends on files outside the repository.
