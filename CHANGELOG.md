# Changelog

## 1.4.0 — 2026-09-07

### Added

- Explicit Recalculate and Discard changes controls for staged key-residue edits.
- Chart ranking limited to the top 100 ligands, with selectable mean interaction
  counts or mean distinct fingerprint features per pose and deterministic ties.
- Visible chart scope and ranking disclosures, also included in figure exports.

### Changed

- Analytical figures are prepared on demand and redundant refreshes are avoided.
- Complete table exports retain every analyzed ligand and pose regardless of
  chart ranking. Limited figure data remain explicitly distinguished.

### Fixed

- Apply key-residue changes consistently to both comparison systems.
- Preserve chain-unspecified residue selections when editing checkboxes.
- Prevent saving/exporting pending parameters with previously calculated results.
- Keep the host Qt event loop active when running an embedded smoke check.

### Validation

- 316 local tests passed with 83.63% branch-aware coverage; Ruff checks passed.
- Synthetic GUI checks with 1,000 ligands retained all 3,000 poses while limiting
  chart rendering to the selected top 100 ligands.

## 1.3.0 — 2026-08-26

### Added

- Preserve source-native ligand and pose identities throughout analysis.
- Parse molecule names and GOLD `Gold.*.Fitness` docking scores from MOL2.
- Add source pose metadata to results, QC, interface labels and XLSX exports.
- Add the `Ligand Index` worksheet and scalable chunked `.docklens` schema 5.

### Changed

- Paired workflows now prefer source ligand identities over generated pose names.
- Analytical exports repeat direct ligand-to-pose traceability.
- Large project results are written and verified as incremental NDJSON chunks.

### Compatibility

- Project schemas 1–4 remain readable.

## 1.2.0 — 2026-08-18

### Input architecture and scientific correctness

- Added the immutable `InputPlan` / `InputJob` contract and one generic
  runner for combined complexes and explicit protein-plus-ligand jobs.
- Added multi-group receptor loading, shared-receptor caching and append-only
  ligand-file selection. System A and System B use the same plan and runner.
- Every uploaded source receives a stable `source_id`; every resolved pose is
  identified as `Sxxxxxx:Pxxxx:Rxxx`, with receptor provenance retained in
  Input QC, exports and reproducible projects.
- External receptors must contain exactly one structural model. A malformed
  source only affects its dependent jobs; a bad receptor does not silently
  select a model or cancel unrelated groups.
- Corrected MOL2 multi-group resolution so ligand groups are never treated as
  receptor atoms for one another. PDB/PDBQT HETATM candidates are split by
  connectivity, and PDBQT ATOM/HETATM classification is preserved.
- `.docklens` projects now use schema 4 while retaining migrations for
  schemas 1–4. Explicit input plans, external receptor paths and input modes
  are included in the methods record.

### Desktop release

- Added a dedicated **Protein + ligands** workflow with multiple receptor
  groups, multiple ligand/pose files per group and separate complex-file
  loading.
- Windows distribution now includes both the portable executable and an
  Inno Setup installer with Start Menu and Desktop shortcuts.
- Publication TIFF export, PLIP/LUNA/DSV-like profiles, the conservative
  LUNA × DSV profile and the shared Okabe–Ito color contract remain available
  in this release.

### Validation

- Added regression coverage for input plans, receptor caching and isolation,
  multi-model receptor rejection, PDB/PDBQT/MOL2 resolution, project schema
  round trips and the multi-group Qt dialog.
- The release CI requires lint, full tests with 80% branch coverage,
  dependency audit, executable smoke tests and installer asset validation.

## 1.1.0 — 2026-08-13

### Added

- Publication-quality TIFF export for every analytical chart.
- LZW compression and 600 DPI metadata for individual TIFF bundles.
- Batch **Export all TIFFs** command with stable chart names, per-chart CSV and
  JSON sidecars, and a collection manifest.
- Atomic staging and rollback for individual bundles and complete collections.
- Pillow as the explicit TIFF runtime dependency.

### Scientific profiles

- PLIP legacy profile remains available for reproducibility.
- LUNA 0.14 defaults and corrected Discovery Studio Visualizer-like criteria
  are available in DockLens and the PyMOL plugin.
- Conservative LUNA x DSV combines representable interaction families using
  stricter shared geometry and semantic deduplication.
- Added `chalcogen`, `attractive_charge` and `charge_repulsion` families while
  preventing salt-bridge duplication and generic contact inflation.
- Unified the 18-type Okabe-Ito-derived color contract across charts, tables,
  exports and PyMOL.

### Validation

- 247 tests passing with 83% total branch coverage.
- PyInstaller executable build and packaged `--self-check` passing.
- Dependency audit reports no known vulnerabilities.
