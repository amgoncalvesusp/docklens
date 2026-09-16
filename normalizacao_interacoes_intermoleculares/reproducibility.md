# DockLens 1.4.1 / PyMOL plugin 0.7.2 compatibility

The two applications ship the same reviewed interaction core. Match chemistry
profile (`plip`, `luna`, `dsv`, `luna_dsv`), cutoffs, enabled families, analysis
view (`complete` / `ds_like`), frame, receptor/ligand selections and export filters
before comparing results. Application defaults are not a substitute for recorded
parameters. Reuse the original molecular files and preserve their hashes.

## Corrected causes

1. Ring labels and members previously followed sets of temporary atom indices.
   They now follow residue identity, atom names and geometry, with the serial
   only breaking ties. Reindexing, reordered inputs and renumbered serials do
   not swap the physical rings in the regression fixtures.
2. A formal charge on a protein group atom was also counted alongside the
   group's centre. Each recognized charged protein group now contributes one
   centre, including when PyMOL has inferred a charge on an atom in that group.
3. Standard PDB protein rings lacked SYBYL aromatic evidence. Complete named
   PHE/TYR/TRP/HIS rings now supply this evidence after the existing topology and
   planarity checks. Explicit non-aromatic typing overrides this fallback.
   This restores the missing tyrosine–sulfur contact in the audited case.
4. PyMOL's `??` means unavailable atom type. The adapter now treats it as empty,
   keeps PDB connectivity and formal charges, and does not interpret inferred
   PDB bond orders as source MOL2 evidence. Typed MOL2 bonds retain their orders.
5. Chalcogen donor geometry previously used only the first bonded neighbor.
   The detector now checks all bonded directions and uses the best qualifying
   donor angle, so bond enumeration cannot hide a valid contact.

The standard aromatic templates follow the
[wwPDB Chemical Component Dictionary](https://www.wwpdb.org/data/ccd), including
[TYR](https://www.rcsb.org/ligand/TYR). PDB
[CONECT records](https://www.wwpdb.org/documentation/file-format-content/format33/sect10.html)
describe connectivity; external chemistry perception must not be silently
confused with retained source typing.

## Verification and limits

`tests/test_reproducibility.py` covers the corrected contracts with synthetic,
portable fixtures. Existing chemistry, parser, profile-isolation and project
save/load tests remain part of the release gate.

An additional local audit compares five supplied ligand structures against the
same receptor in native PyMOL 3.1.4.1 and DockLens, in all four profiles. Comparison
uses interaction family/subtype, endpoint atom membership, roles, hydrogen
identity and distances at the CSV precision of 0.01 angstrom. This is semantic
equality, not binary equality of floating-point coordinates: PyMOL stores
coordinates at a different precision. Audit inputs and detailed evidence remain
local because they contain unpublished research data.

Final local results: 330 DockLens tests passed (83.84% branch-aware coverage),
100 plugin tests passed, and all 365 contacts across the 20 native cases matched.
Maximum native distance difference was 0.000002893 angstrom, below the exported
0.01 angstrom precision. Two drawing/export passes and native PSE save/reload
were identical in every case. DockLens fresh runs and `.docklens` project
reload/recalculation produced byte-identical CSVs for all four profiles.

Different source types, added/removed hydrogens, different states, changed cutoffs
or re-exported molecular chemistry do not meet the identical-input contract.
In particular, a PDB cannot preserve all MOL2 chemistry. This release does not
claim universal equivalence after arbitrary format conversions, or reconstruct
missing original inputs from exported interaction tables.
