# Product

<!-- impeccable:product-schema 1 -->

## Platform

adaptive

## Users

Primary users are computational chemistry, structural biology, medicinal
chemistry and molecular-modeling researchers working on docking campaigns and
molecular-dynamics trajectories.

**Inferred from the current product and brief:** users commonly work on
Windows or Linux desktop systems, inspect dense scientific datasets and need
publication-ready evidence rather than simplified consumer-style summaries.

## Product Purpose

DockLens detects, audits, compares and communicates non-covalent
protein-ligand interactions. Success means that a researcher can move from
structures or frames to residue-level conclusions while retaining the exact
poses, atoms, geometry, parameters and filters behind every result.

## Positioning

DockLens joins transparent interaction detection, a calibrated
Discovery Studio-like option, docking-pose analysis and molecular-dynamics
interaction persistence in one reproducible desktop workflow.

## Operating Context

- Inputs include docking poses and structure files or extracted frames.
- Results are reviewed as pose/frame summaries, residue coverage, interaction
  details and matrices.
- Researchers compare ligands, poses, time ranges, trajectories and key
  residues.
- Uploaded ligand files can be charted individually with independent System A
  and System B denominators; the pooled view discloses observation weighting.
- Docking fingerprints can form pose families; ordered MD fingerprints can
  form interaction states with descriptive transitions and representatives.
- Tables and figures are exported for audit, collaboration and publication.
- The companion PyMOL plug-in remains the molecular-viewing surface.

## Capabilities and Constraints

- The complete and Discovery Studio-like analysis views must remain available.
- Discovery Studio-like is an empirical, transparent approximation; the
  product must not claim proprietary algorithm parity.
- Docking frequency and molecular-dynamics occupancy are distinct scientific
  concepts and must be labeled separately.
- Raw atom-pair events remain available even when charts use consolidated
  residue/type/pose-or-frame counting units.
- Existing CSV/XLSX results and v0.6.0 behavior must remain reproducible.
- State clustering must disclose method, threshold, training sample and
  outliers.
- MD intervals must preserve temporal dependence through block resampling,
  disclose every parameter and never be applied to docking poses.
- MD analyses must visibly disclose the single-series fallback. An explicit
  trajectory map defines observation order, replica boundaries, frame gaps and
  time values and must survive project save/open.
- Saved projects retain the unprofiled result so Complete and Discovery
  Studio-like remain reversible without redetection.
- Chart scopes use the uploaded source identity rather than generic ligand
  labels, preserve zero-contact observations and survive project save/open.
- Version 1.0.0 is developed on an isolated branch and does not replace the
  current stable release until explicitly approved.

## Brand Commitments

- Product name: DockLens.
- Scientific, direct and restrained voice.
- Okabe-Ito interaction colors remain the semantic color authority.
- Existing DockLens icon and authorship attribution are preserved.

## Evidence on Hand

- Existing detector, result contracts, residue matrix and export pipeline in
  `docklens/`.
- Discovery Studio-calibrated 2m5d corpus analysis documented in `README.md`.
- Automated tests covering detection, export, profiles and desktop behavior.
- No user study, accessibility audit or external visualization benchmark is
  currently recorded; future work must not fabricate these.

## Product Principles

1. Every visual statement must be traceable to auditable rows and parameters.
2. Preserve scientific meaning across docking and molecular dynamics.
3. Use interaction fingerprints to reveal patterns, not to hide assumptions.
4. Keep dense workflows scannable and selections synchronized.
5. Export the data behind every publication figure.
6. Express uncertainty and temporal assumptions instead of implying that
   saved frames are independent.

## Accessibility & Inclusion

Use color-blind-safe interaction colors, never encode meaning by color alone,
support keyboard navigation and readable high-DPI scaling, and retain usable
layouts on compact Windows and Linux screens.
