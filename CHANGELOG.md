# Changelog

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
