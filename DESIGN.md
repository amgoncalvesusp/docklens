---
name: DockLens
description: A reciprocal-space atlas for auditable molecular interaction analysis.
colors:
  instrument-frame: "#071A2E"
  navigation-accent: "#0072B2"
  analytical-paper: "#F7F9FA"
  panel-paper: "#FFFFFF"
  secondary-surface: "#EEF2F4"
  structural-rule: "#D8E0E5"
  primary-ink: "#142330"
  secondary-ink: "#5B6976"
  strong-heading-ink: "#101C26"
  metric-ink: "#294151"
  control-border: "#AEBBC4"
  navigation-text: "#DCE7EE"
  navigation-hover: "#102E49"
  table-heading-ink: "#203644"
typography:
  title:
    fontFamily: "Aptos, Segoe UI Variable, Segoe UI, Noto Sans, sans-serif"
    fontSize: "20px"
    fontWeight: 700
  section:
    fontFamily: "Aptos, Segoe UI Variable, Segoe UI, Noto Sans, sans-serif"
    fontSize: "14px"
    fontWeight: 700
  body:
    fontFamily: "Aptos, Segoe UI Variable, Segoe UI, Noto Sans, sans-serif"
    fontSize: "12px"
    fontWeight: 400
  label:
    fontFamily: "Aptos, Segoe UI Variable, Segoe UI, Noto Sans, sans-serif"
    fontSize: "11px"
    fontWeight: 600
rounded:
  analytical: "3px"
  control: "4px"
spacing:
  control-y: "7px"
  control-x: "11px"
  workspace: "24px"
  section: "16px"
components:
  button-primary:
    backgroundColor: "{colors.instrument-frame}"
    textColor: "{colors.panel-paper}"
    rounded: "{rounded.control}"
    padding: "{spacing.control-y} {spacing.control-x}"
  field:
    backgroundColor: "{colors.panel-paper}"
    textColor: "{colors.primary-ink}"
    rounded: "{rounded.control}"
    padding: "5px 7px"
---

# Design System: DockLens

<!--
THESIS: A reciprocal-space atlas: dense molecular evidence becomes calm,
traceable residue, fingerprint and comparison maps.
OWN-WORLD: Deep-navy instrument frame, paper-white analytical fields,
Okabe-Ito evidence channels, precise one-pixel rules and the DockLens Lens as
the synchronized selection signature.
STORY: Start with residue prevalence, reveal pose/frame fingerprints, compare
systems or intervals, inspect the underlying geometry, then export the figure,
rows and parameters together.
FIRST VIEWPORT: Compact rail and context bar frame a dominant residue chart and
barcode; the Lens inspector explains the current residue without obscuring
evidence. Fingerprint and Compare are peer workspaces, not secondary dialogs.
FORM: Desktop research instrument with squared plot fields, shallow controls,
independent scrolling and collapsible inspector. Compact windows collapse labels
before charts; interaction colors never decorate navigation. The approved
compositions are structural references, not pixel-literal specifications.
-->

## Overview

**Creative North Star: "The Reciprocal-Space Atlas"**

DockLens should feel like a contemporary structural-biology instrument: precise
enough for dense scientific evidence, calm enough for long analysis sessions,
and unmistakably organized around residues, interaction channels and temporal
states. Contact maps, reciprocal-space plots and publication figures inform the
visual grammar without turning the application into a literal laboratory prop.

The reusable signature is the **DockLens Lens**: one selected residue, pose,
frame or time interval becomes the shared focus across every chart, table and
inspector. Selection is expressed through position, labels and shape as well as
color.

**Key Characteristics:**

- Restrained analytical surfaces with one strong navigational accent.
- Okabe-Ito interaction colors used only as scientific category semantics.
- Wide data viewports, compact controls and a persistent contextual inspector.
- Crisp plot fields, restrained contours and sparse depth.
- Familiar desktop behavior with high-DPI and compact-screen adaptation.

## Colors

The palette uses cool paper-like surfaces, deep blue-black structure and a
controlled teal-blue navigational accent. Interaction colors inherit the
Okabe-Ito authority already established by the detector.

**The Semantic Color Rule.** Interaction colors never decorate containers,
navigation or headings; they identify interaction evidence only.

**The Single Focus Rule.** Primary accent marks the current task, selection and
keyboard focus, never multiple competing calls to action.

## Typography

One workhorse system sans-serif family carries navigation, controls and prose.
Tabular numerals and a compact monospaced fallback are reserved for residue
identifiers, frame indices, distances, angles and fingerprints.

Hierarchy is tight and operational: titles establish location, labels support
scanning and data remains the dominant visual voice. No display typefaces or
ornamental scientific lettering appear in controls.

## Layout

The spatial grammar is an atlas workspace: compact navigation, one dominant
analysis field and an optional contextual inspector. Toolbars remain shallow;
filters use progressive disclosure rather than occupying permanent vertical
space.

The main analysis field supports a responsive two-column dashboard that becomes
one column on compact windows. Tables and heatmaps retain independent scrolling.
The inspector collapses before the analysis field is reduced below a useful
scientific width.

**The Data Owns the Viewport Rule.** Controls frame the evidence; they do not
consume the area required to compare residues, poses or frames.

## Elevation & Depth

Depth is tonal and structural. Most surfaces are flat and separated by subtle
changes in background or one-pixel rules. Shadows appear only for temporary
overlays, floating inspectors or drag feedback.

## Shapes

Containers use gently softened corners; plots, matrices and tables remain more
rectilinear. Pills are reserved for removable filters and compact status, not
used as a default component silhouette. Focus rings are visible and never
replaced by color-only state.

## Components

### Buttons

Controls are compact and operational. Primary detection uses the instrument
frame; secondary actions use white analytical paper with a structural border.
Hover is tonal and keyboard focus uses a two-pixel navigation-accent outline.
Disabled controls retain readable text and do not rely on opacity alone.

### Inputs and Selectors

Fields use panel paper, a structural border and the control corner. Profile and
evidence selectors remain in the shallow context bar because they change the
meaning of every workspace.

### Navigation

The deep-navy rail is the application frame. The current peer workspace receives
the navigation accent; hover uses a darker tonal layer. Labels remain visible
until compact width makes horizontal scrolling preferable.

### Analytical Fields

Plots use paper-white backgrounds and one-pixel structural rules. Interaction
colors belong to evidence marks only. Every chart owns a tidy dataframe and
reproducibility metadata, exported together through the figure bundle.

### Pose Families & Dynamic States

Advanced fingerprint analysis remains inside the Fingerprint workspace.
Population, timeline, observed-transition and confidence views use compact
internal tabs so they do not compete with the four peer workspaces. State
colors form a separate categorical channel and never reuse interaction-type
colors. `OUTLIER` is always labeled and rendered neutrally.

Docking uses “pose family” and “frequency”; MD uses “interaction state”,
“occupancy”, “episode” and “observed transition”. Transition plots must state
that they are descriptive rather than validated Markov models. Long bootstrap
calculations run in background workers while existing evidence remains visible.
The MD controls expose a trajectory-map action beside the saved-frame step.
Status text distinguishes the implicit contiguous fallback from an explicit
map and summarizes replica and gap counts. The Fingerprint header contains an
explicit export-view selector so export never depends on whichever internal
panel happened to be active.

### Reproducible Projects

Open and Save project actions live in the navigation rail as persistent
workflow actions. Stale external sources produce a visible warning while
verified cached evidence remains readable. Integrity verification must never
be described as cryptographic authorship or authenticity.

### DockLens Lens

The Lens is the signature inspector. It names the selected residue, reports the
current frequency or occupancy by channel and states the active analytical
profile and counting unit. It collapses on compact windows without changing the
selection.

### Ligand / Uploaded-file Scope

A persistent two-row context bar separates detection settings from chart
scope. System A and System B use independent selectors. Each option shows the
uploaded filename, detected ligand label, observation count and stable source
identifier; generic labels such as `LIG` never merge different uploads. The
metric text states whether charts show one source or a pooled,
observation-weighted view. Scope changes update analytical charts and the Lens
while the complete tables remain available.

### Observation Identity & Interaction Heatmap

The first action row exposes **Chart labels** before detection. Ligand name is
the recommended default; uploaded filename and pose/frame index remain
available. These values are presentation labels only. Stable observation and
source identifiers remain the basis of filtering, joins, state analysis and
export, and repeated display names receive a concise pose/frame suffix.

Fingerprint contains a dedicated interaction comparison heatmap in addition
to the Jaccard/Tanimoto similarity matrix. Aggregate rows compare all loaded
ligand/uploaded-file groups, normalized by each row's own observations.
Individual rows follow the active chart scope. Columns switch between
residue × interaction type and any-interaction residue presence. The status
strip always states row count, feature count, observation count, scope,
zero-contact convention and whether the visible top-feature or hard cell limit
was applied. Safety reductions are deterministic and recorded in exported
metadata.

Header labels on the navy top bar always render on a transparent background.
White title text and pale-blue context text are reserved for this dark surface
to maintain readable contrast.

## Do's and Don'ts

### Do:

- **Do** keep filters, charts and tables synchronized through the DockLens Lens.
- **Do** label docking frequency and molecular-dynamics occupancy distinctly.
- **Do** expose the counting unit and denominator beside every aggregate chart.
- **Do** preserve the Discovery Studio-like profile as a visible analytical
  choice with its parameters available for audit.
- **Do** export the data and parameters behind every figure.

### Don't:

- **Don't** use interaction colors as general UI decoration.
- **Don't** hide raw evidence behind a composite score.
- **Don't** place every analysis inside identical cards.
- **Don't** rely on hover for essential scientific information.
- **Don't** imply that docking poses are temporal observations.
