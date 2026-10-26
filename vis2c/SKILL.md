---
name: vis2c
description: Perform end-to-end wavefunction analysis, quantitative comparison, and
  scientific visualization using Multiwfn, TRESC/vis2c (cub1c, real3damp,
  real3d, cuboper, alignmol, atomic, bond), and Blender. Trigger when a task
  calls for generating, extracting, validating, analyzing, comparing, or
  visualizing wavefunction-derived properties or 3D real-space fields,
  preparing the numerical data for these tasks, or refining existing
  scientific figures. Common tasks include plotting canonical, natural,
  natural-transition, natural-bond, singly occupied natural, or localized
  orbitals (HOMO/LUMO/SOMO, NO, NTO, NBO, SONO) from suitable Molden files
  with configurable orbital grids ($TRESC/gridsettings.ini); analyzing or
  rendering scalar fields and mapped surfaces such as electron/spin
  density, ESP, IGMH, ELF, LOL, and Fukui/dual descriptors based on common
  wavefunction file formats (e.g., .fch, .molden.input, etc.); investigating
  substituent effects or conformational property differences through
  aligned field comparisons or atomic/bond-property tables; aligning
  molecular structures or fields; and improving figure views, colors,
  legends, lighting, and layout. These examples are not exhaustive,
  route other wavefunction-analysis quantities by their scientific
  definition and required representation. Reuse valid existing results,
  complete supported figure requests through rendering and actual visual
  review, and stop at validated numerical data when requested.
---

# vis2c: From Electronic-Structure Results to Scientific Figures

## 1. One workflow, selectable stopping point

This skill connects numerical preparation, comparison, and visualization for wavefunction-derived analyses. Its scientific examples are open-ended, but execution and completion are limited by verified numerical generators and representation capabilities. Orbital images, IGMH, electrostatic potential (ESP), density/spin-density surfaces, atomic/bond properties and substituent comparisons below are illustrative scientific questions, not a closed list of supported applications. Apply the same reasoning to other wavefunction-derived quantities: establish their definitions and source requirements, compute the required data, choose a scientifically appropriate representation, and complete the requested visualization where the tools support it.

Distinguish broad scientific scope from concrete exporter capabilities. The three supplied model exporters implement signed real orbitals, one scalar contour, and a scalar contour colored by a second field. These are reusable representation primitives for many analyses, not three permitted scientific topics. A new descriptor fitting one of these primitives needs no new rendering adapter. A vector/tensor display, topology graph, complex phase representation, section, or other representation outside those primitives may require a task-local extension or a documented data handoff; never silently replace it with an unrelated scalar image.

Treat numerical preparation, field analysis and rendering as connected stages, not unrelated skills.

```text
Define the scientific question and representation
-> reuse suitable results, or run an authorized calculation/export
-> obtain wavefunction/Molden data and the required numerical properties
-> generate and validate cubes or atom/bond-property tables
-> align and combine fields when the question requires it
-> export the correct model representation
-> render, actually inspect, and refine bounded previews
-> deliver final candidate figures, numerical data and reproduction records
```

Select the stopping point from the request:

- A figure request means complete the supported chain through rendering; do not stop at a cube and tell the user to draw it.
- A data-only request ends at validated numerical files. Do not render against the user's stated preference.
- A geometry-alignment request may end at aligned XYZ and transformation/fit records.
- Atomic/bond-property requests use the data-export paths documented below; their GUI scripts are not automatically compatible with the OBJ/Blender surface pipeline.
- Unsupported representations must be disclosed rather than silently replaced with a superficially similar picture.

Reuse any already valid stage. Do not rerun SCF when an appropriate wavefunction exists, generate a duplicate cube before cub1c generates it internally, or re-export a model just to change lighting.

### Trigger boundaries and route-specific reading

Trigger for a concrete request to generate, validate, compare, export, or visualize supported numerical data. Do not start numerical jobs solely because a conceptual question mentions wavefunction analysis. Do not turn every substituent or conformer comparison into cube subtraction. A data-only request is a first-class endpoint, not an incomplete figure task.

Read sections 1-3 and 14 for task scope, deployment, routing, authority, and honest completion. Then read the applicable sections completely before execution; retrieve further passages if a tool truncates them. Do not load unrelated API examples merely to execute a simple route.

| Selected route | Additional sections to read |
|---|---|
| Numerical generation or validation | Relevant parts of section 4 and any selected downstream route |
| Molden real orbital | Sections 4.3-4.5, 7, and shared bundle rules in 9.4 |
| Geometry-only fitting | Section 5 |
| Aligned scalar-field arithmetic | Sections 4.4-4.5, 5-6; include property-specific definitions from section 4 |
| Single-contour surface | Sections 4.4-4.5, 8, and 9.4 |
| Quantitatively mapped surface | Sections 4.4-4.6 and 9; include 4.2 for IGMH |
| Any supported rendered figure | Section 10 and the applicable exporter section; section 13 when its playbook is relevant |
| Atomic or bond-property handoff | Relevant numerical preparation plus section 11 or 12 |

Configured deployment paths do not require repeated discovery. However, this document does not override observed runtime behavior: verify relevant help, signatures, or implementation when compatibility is uncertain or an actual mismatch appears. Record the discrepancy and do not invent flags or weaken validation to continue. Reuse version-appropriate verified recipes within a session when their identity and applicability remain established.

## 2. Deployment and task environment

Use the guaranteed deployment facts directly, without repeatedly checking script existence, locating alternate copies, enumerating Python environments or inspecting shebangs:

- `$TRESC` is the installation root.
- `$TRESC/gridsettings.ini` is the documented orbital-grid configuration for `cub1c.py` / `tshell.sh -cub1c`, with `rgrid` and `edge` settings. Read and record applicable values when generating orbital cubes. This is a shared installation resource, not unconditional permission to edit it; follow section 7.4. Verify parameter semantics against the installed implementation before changing values.
- The plotting/export tools are `$TRESC/vis2c/cub1c.py`, `real3damp.py`, `real3d.py` and `render_obj.py`.
- `$VIS2C_PYTHON` is the executable for project numerical/export scripts, with their dependencies configured.
- `$BLENDER_BIN` is the user-provided Blender executable path, with its matching bundled resources and color management.
- `$MULTIWFN_EXE` is the user-provided executable entrypoint for Multiwfn. It may be the Multiwfn binary OR an executable user-written shell wrapper. Invoke this entrypoint as supplied; do not replace it with a discovered binary or bypass its initialization.
- Numerical helpers and user-facing atomic/bond tools are addressed by their own names under `$TRESC/vis2c/`: `cuboper.py`, `alignmol.py`, `atomic.py`, `bond.py`.

These four environment variables are guaranteed user inputs. Use them directly; do not ask the user to supply additional configuration or repeatedly probe these guaranteed paths. `$MULTIWFN_EXE` replaces the former Multiwfn variable name throughout this workflow. Quote it as one executable path, not a shell command string; do not use `eval`, split it into words, or force `bash`/`sh` onto a possible binary.

Use `"$VIS2C_PYTHON"` for the vis2c tools. Blender executes rendering code with its own Python, not VIS2C_PYTHON. Call Multiwfn with `"$MULTIWFN_EXE"` for both terminal and batch workflows.

### Agent-selected additional configuration

The agent determines ALL remaining operational settings from the task and the available runtime. The user is not expected to define extra environment variables, choose a display runner, tune threads, or provide output-directory defaults. Use ordinary task-local variables/configuration and record the choices. Respect explicit scientific or presentation preferences already given by the user.

| Decision | Agent responsibility |
|---|---|
| Output workspace | Choose a writable absolute task parent, normally `figures/` under the project; use a fresh task directory. |
| CPU and memory | Choose a positive thread limit within the allocation, normally start with at most 4 rendering threads; set preview/final costs and numerical grid sizes appropriately. |
| Mayavi display | Use an available working display; if necessary and available, choose a working `xvfb-run` wrapper. Diagnose actual backend/authentication failures without requesting a new user-defined configuration variable. |
| Numerical and display settings | Choose property-appropriate grids, thresholds, mapping limits, LUT, view, lighting, layout, resolution and a bounded refinement budget. Explain choices that affect scientific interpretation. |
| Task reproducibility | Save actual commands, local variable/configuration values, tool identities, logs and review records. |

Invoke only the tools needed by the chosen route. Existing cubes do not require new SCF or orbital generation; data-only work does not require Blender. If a guaranteed entrypoint actually fails, diagnose that execution failure and retain its log rather than silently choosing a different installation. Ask for help only for a real blocker the agent cannot resolve within its authority.

The Mayavi export branches use offscreen=False and skip interactive controllers. They may briefly create windows and still require usable display backends. Blender background mode is a separate stage and does not provide that display context. Do not change shared backend code or install drivers/packages to solve an aesthetic issue.

Portable deployment means configured paths, native dependencies and resources, not a promise that Markdown alone runs on an unconfigured host. Do not embed usernames, hosts, Conda paths or fixed Blender installation paths.

Create an isolated task. The following is a Bash pattern; choose the work parent, thread count and runner for the actual task before execution:

```bash
WORK_ROOT="$PWD/figures"
mkdir -p "$WORK_ROOT"
TASK="$(mktemp -d "$WORK_ROOT/vis2c-XXXXXXXX")"
mkdir -p "$TASK"/{inputs,calculations,data,work,models,scripts,previews,final,logs}
THREADS=4  # Lower if the allocation is smaller; select the actual limit first.
EXPORT_RUNNER=()
# If the agent determines a virtual display is needed and available, use:
# EXPORT_RUNNER=(xvfb-run -a -s "-screen 0 1280x1024x24")
```

Carry absolute task paths and settings explicitly between independent shell/tool calls. Use separate subdirectories for each molecule and each Multiwfn run. Add `workflow.json`, `render_review.json` and `reproduction.md`. Do not dump the entire user environment, which may contain credentials, into records.

All command examples below are patterns to instantiate for the actual task, not commands to paste with unset variables. Bind input paths, property definitions and numerical/display settings first; create any named per-stage work directory before `cd`, and copy/link the intended inputs into it where the example uses local filenames. Use a Bash process for Bash array examples even if the login shell differs. Execute only the stages needed for the requested endpoint.

## 3. Plan from the scientific question

| Goal | Numerical/model chain | Endpoint |
|---|---|---|
| Selected real orbital (canonical, natural, natural-transition, natural-bond, singly occupied natural, or localized, when faithfully exported) | Suitable result -> validated Molden and index mapping -> verify orbital-grid settings (section 7.4) -> cub1c -> render_obj only if requested | Signed orbital cube and provenance; rendered image when requested |
| One density/spin-density/other field surface | Wavefunction -> numerical generator -> cube -> real3damp -> render_obj | One-contour image |
| IGMH surface with sign(lambda_2)*rho colors | Wavefunction -> Multiwfn paired cubes -> real3d -> render_obj | Quantitative mapped image and legend |
| ESP mapped onto a density-defined molecular surface | Wavefunction -> Multiwfn density + ESP on matching grids -> real3d -> render_obj | ESP-colored surface, actual density contour and ESP legend/unit |
| Another wavefunction-analysis quantity | Establish definition/source -> Multiwfn or a verified suitable generator -> appropriate data/model representation | Requested analysis visualization, or a precise representation-extension requirement |
| Aligned pointwise scalar-field comparison for a defined substituent or conformer question | Comparable fields -> justified shared-scaffold fit -> cuboper arithmetic -> validate difference, overlap, and interpolation -> render only if requested | Signed difference cube, domain-qualified statistics, comparison record, and optional sign-complete views |
| Common-fragment geometry fit | XYZ pair -> alignmol numerical API | Aligned XYZ, transformation and fit record |
| Atomic scalar-property visualization | Calculate/extract properties -> XYZQ | Data plus user-run atomic.py command |
| Bond/pair-property visualization | Calculate/extract properties -> XYZ + square matrix | Data plus user-run bond.py command |
| Two-component complex orbital | Appropriate wavefunction -> TRESC component cubes | Numerical data unless a separate spinor-rendering route is authorized |
| 3D-field cross-section | Generate full cube(s) and plane definition | Data unless an explicit section-rendering extension is scoped |

The table is a set of examples. Route an unlisted quantity by its mathematical character and intended representation rather than reject it for not appearing here. Before running anything, specify relevant structure/conformer/state, charge/multiplicity, computational method, orbital/spin, field definition, fragments, units, comparison direction and intended region. Use project-approved choices when unambiguous; ask only for consequential missing scientific information. Operational additional configuration is the agent's responsibility.

### Calculation stage

If suitable results exist, validate and reuse them. If the task explicitly requires a new calculation and the project grants execution authority, prepare/run the necessary calculation with the installed quantum-chemistry workflow. Use an existing validated program recipe rather than inventing executable flags.

Record geometry, charge/spin, method/basis, solvent/model and relevant convergence/state checks. Export a wavefunction/Molden file adequate for the downstream property. Do not claim that a successful file conversion proves SCF convergence, correct spin state or the physical meaning of the orbitals.

If a calculation is expensive, undefined or outside permissions, present the plan/resource need and obtain the required decision before execution. Once inputs and authority are sufficient, do not unnecessarily delegate executable calculation/postprocessing steps to the user.

## 4. Wavefunction and cube generation

### 4.1 Multiwfn: run the numerical workflow, do not guess it

Multiwfn is the core analysis engine for wavefunction analyses beyond orbital visualization in this workflow. Actively use it to compute the required quantities, including but not limited to IGMH, ESP, density/spin density, and applicable atomic/bond properties. Reuse valid results or a verified more suitable generator when justified; do not force every property through the orbital-generation route or treat the examples as the extent of Multiwfn's use.

Consult the [online Multiwfn manual](http://sobereva.com/multiwfn/manual/) for the requested analysis, its theoretical definition, supported inputs, menu workflow, numerical settings and exports. This is an available reference for any Multiwfn analysis, not only the examples here. Match online instructions to the actual runtime version/prompts; local documentation and verified project recipes are also useful. If the website is inaccessible, use a trustworthy local/version-appropriate source and observed prompts; do not claim to have read inaccessible material or invent menu sequences.

A filename extension alone does not guarantee that the source contains the basis, coefficients, occupations, density or spin information needed by the requested quantity. Confirm the input's adequacy and perform a documented numerical conversion if needed. Identify the Multiwfn version/build from actual startup output, logs or available documentation; do not invent version flags. A wrapper is part of the configured runtime and may set libraries, resources, threads or working-directory behavior: retain it and check where the run actually writes its outputs.

1. Read the relevant section of the online manual above, applicable installed manual or a verified version-appropriate project recipe. Inspect actual prompts with a controlled process/PTY when necessary.
2. Select the exact quantity, orbital/spin, fragment definition and numerical output. Interactive terminal menus are acceptable; do not launch graphical viewers as a data-generation workaround.
3. Choose documented numerical settings. For volumetric output, specify grid extent, origin/axes and resolution; for paired mapped fields, request identical physical grids. For atomic/bond data or other non-grid analyses, request the appropriate numerical tables/results instead of inventing a cube requirement.
4. Run in an isolated per-molecule directory. Save the actual input recipe/transcript, command, log and property definitions. Use bounded runtime/resources.
5. Check process/log status and output provenance, then cube headers/counts/finite values/grids or the appropriate table dimensions/atom order/units before copying to descriptive deliverable names. Reject stale default-named outputs after a failed run. A wrapper returning success alone does not establish successful analysis.

Usage example AFTER the recipe has been verified and created:

```bash
mkdir -p "$TASK/work/multiwfn_A"
cd "$TASK/work/multiwfn_A"
"$MULTIWFN_EXE" "$WAVEFUNCTION_A" \
  < "$TASK/scripts/multiwfn_A_verified.in" \
  > "$TASK/logs/multiwfn_A.log" 2>&1
```

The directory, input file and recipe must actually exist. This is an invocation pattern, not a universal Multiwfn menu sequence. Do not hard-code menu numbers from memory or silently reuse another property/version's recipe. If prompts differ, inspect rather than feed repeated guesses.

For an observed terminal-menu workflow, start the same configured entrypoint in a PTY:

```bash
cd "$TASK/work/multiwfn_A"
"$MULTIWFN_EXE" "$WAVEFUNCTION_A"
```

Respond to actual prompts and save the transcript and final recipe. Set the numerical extent/resolution/output through the documented analysis workflow, not guessed command-line flags. The agent performs the numerical work when needed; simply telling the user to run Multiwfn does not complete an authorized data-generation or figure task.

### 4.2 IGMH definitions must be explicit

Resolve what `dg` means: total delta_g, delta_g_intra, or delta_g_inter with a particular atom-to-fragment assignment. Intramolecular interactions can be studied between defined fragments inside a molecule; do not equate interfragment with necessarily intermolecular.

The geometry-fitting atom selection and the IGMH fragment partition serve DIFFERENT purposes. Record both. Across substituted molecules, preserve chemically corresponding fragment definitions, not identical raw index strings when numbering/composition changed.

For a conventional mapped IGMH figure, export the intended delta_g field for geometry and sign(lambda_2)*rho (sl2r) for coloring from the appropriate actual-density analysis. Preserve the correct quantity, units and grid. Do not substitute RDG, promolecular IGM or a different fragment component without a deliberate scientific decision.

For methyl-substitution comparisons, choose the question: total change after substitution, change in a common region, or a more controlled comparison. Optimized-geometry differences mix structural and electronic effects. A difference of IGMH descriptors is a descriptor change, not directly an interaction energy or proof of an isolated causal contribution of the methyl group.

### 4.3 TRESC cube-only routes

For a data-only request, example commands are:

```bash
cd "$TASK/work"
# Signed real scalar orbital data only.
tshell.sh -cub1c "$MOLDEN" 35A

# Two-component complex orbital data only.
tshell.sh -cub2c "$SPINOR_MOLDEN_DIR" 15
```

Replace example orbital indices using actual metadata. Expected component names for the examples are `35A-real.cub` and `15-rar.cub`, `15-rai.cub`, `15-rbr.cub`, `15-rbi.cub`; verify actual generated files.

For a scalar-orbital FIGURE, call cub1c.py --export-obj directly instead of generating the same cube once more beforehand. Its normal no-export CLI is interactive, not a cube-only command. Complex spinor components must not be reduced silently to one modulus/density cube.

### 4.4 Numerical validation before analysis or export

Check geometry/state/source, orbital numbering, units, origin, grid vectors/counts, scalar-data counts, finite values, range and spatial coverage. Preserve sign. Do not edit unit labels without actual unit conversion or smooth/rescale scientific inputs for aesthetics.

Paired fields need matching physical coordinates, not merely matching shape. Record deliberate resampling and its errors/limits. Full cube volume plus a requested plane definition is the default cross-section handoff; -slice in the plotting CLIs is not unattended image export.

### 4.5 Orbital identity, components and surface definitions

- Determine HOMO/LUMO/SOMO from the actual energies, occupations, spin and orbital type in the source. Neither cub1c nor the renderer infers them. Do not use electron_count/2 or unverified file order as a universal selection rule.
- For open-shell systems, identify alpha/beta channels and state which channel or cross-channel convention defines the requested orbital. Resolve degeneracy, fractional occupations, natural/localized orbital labels and noncanonical ordering explicitly when relevant.
- Map the selected orbital to the numerical generator's index convention. A label in a caption does not verify that the correct coefficients were exported.
- Real signed amplitude, absolute amplitude, squared modulus, electron density and spin density are different quantities. Preserve all required real/imaginary and spin components of complex/two-component data. Define any derived scalar before computing it.
- Orbital positive/negative colors indicate wavefunction sign under a chosen phase convention, not positive/negative charge. An orbital's global sign is arbitrary; inter-system orbital differences need a defined correspondence and phase convention before subtraction has the intended meaning.
- Specify which field defines a surface and which field, if any, defines its color. ESP on a density isosurface is different from an ESP isosurface. Record both definitions, contour and units separately.
- Choose property-specific contours/ranges. The exporters' numerical default of 0.05 is a software default, not a universal scientifically appropriate threshold for orbitals, IGMH, density, ESP or differences.

### 4.6 ESP and other analyses: generalize the numerical plan

For an ESP-colored molecular surface, compute the surface-defining electron-density cube and the ESP mapping cube with Multiwfn using the same grid-quality/resolution settings and the same physical grid: origin, grid vectors, point counts and spatial extent must match. Selecting the same nominal quality preset alone is not sufficient; validate the exported grids as described in section 9. Verify that both fields refer to the intended structure/electronic state and that the selected ESP definition/unit is the one used in labels.

For this ESP-on-density workflow, explicitly set the density-surface `--isovalue=0.001` as the recommended starting value in the density cube's standard atomic-unit convention; do not leave the generic software default of 0.05 in effect. Verify the density unit and record any scientifically justified alternative. This threshold belongs to electron density, not ESP; choose the ESP color range separately. Under the user-tested installed workflow, ionic systems use the same preparation and `real3d` rendering route as neutral systems: the script already accounts for them, so no ion-specific preprocessing, neutralization, or ad hoc ESP offset is required. Preserve the intended charged-state source data and normal validation.

An ESP isosurface instead uses ESP as the contour field, and a negative ESP contour must follow a route that actually supports signed levels, such as real3damp for a fixed-color single contour. The recommended density threshold above does not apply to an ESP contour.

For an unlisted analysis, consult the relevant Multiwfn manual section and determine the mathematical output and meaningful representation first. A scalar field can use the single-contour or mapped-field routes; atom-resolved scalar values can use XYZQ; pair-resolved values can use a matrix plus XYZ. For vector/tensor/complex outputs, retain the required components and explicitly define any scalar projection, component or invariant being visualized. For topology or other non-field outputs, retain the actual analysis data and choose a suitable task-local representation instead of forcing it into a misleading cube. The examples guide this decision; they do not restrict the analyses the skill can undertake.

## 5. alignmol.py: fit geometry through the numerical API

Use for a whole-molecule or common-scaffold rigid fit. Reference is fixed; probe is rotated/translated into its frame. XYZ coordinates and transformation translation are in angstrom.

Manual rules: 1-based atom numbers excluding XYZ headers; inclusive ascending ranges such as `1-9,16-23,25`; preserve supplied pairing order; equal lists of at least three element-matched, unique, valid, non-collinear pairs. Do not sort selections or permit reflection/nonuniform scaling.

Manual input needs two strings, REF then PROBE, for example `['1,2,4', '8,9,11']` only after verifying those atom correspondences. Entering `mcs` at the first prompt selects both mappings automatically; entering it at the second prompt replaces BOTH manual selections. The non-collinearity test rejects a first centered-coordinate singular value below 1e-8 or a second/first ratio below 1e-3.

MCS uses element/connectivity matching with explicit H and connectivity inferred by distance <= 1.15 times the covalent-radius sum. It does not certify actual bond orders, valence or stereochemical equivalence. Inspect matches for metal coordination, symmetry and reaction geometries. Timeout is 60 seconds; embeddings are capped at 1000 per structure and 100000 tested pairs. Preserve warnings, SMARTS, actual pairs and fitting-fragment RMSD; capped searches do not guarantee a global symmetry optimum.

The CLI `alignmol.py REF_XYZ PROBE_XYZ` writes an aligned XYZ then opens a GUI. Do not use it for unattended numerical work. `align_xyz_fragment()` returns a 4 x 4 transform without writing or plotting, but still prompts for selections.

The exact numerical API is `align_xyz_fragment(ref_xyz: str, probe_xyz: str) -> numpy.ndarray`; it does not accept matching selections as keyword arguments. Its CLI takes only the two positional XYZ paths; there is no `--mcs`, `--ref-atoms`, `--probe-atoms` or `--no-gui`. The CLI output is `<probe-stem>_aligned.xyz` in CWD and may overwrite an unrelated existing result, so use the checked numerical wrapper below for agent work. The intended installed `alignmol.py` must be importable; adding its directory to `sys.path` in the wrapper avoids relying on a hard-coded shebang or a coincidentally installed module.

Task-local usage example; save as `$TASK/scripts/align_task.py` and run in an isolated work directory containing the two input XYZ files:

```python
import os
from pathlib import Path
import sys
from unittest.mock import patch

sys.path.insert(0, str(Path(os.environ['TRESC']) / 'vis2c'))
import numpy as np
from rdkit import Chem
from alignmol import align_xyz_fragment

ref = Path('reference.xyz').resolve()
probe = Path('probe.xyz').resolve()
out = Path('probe_aligned.xyz').resolve()
transform_file = Path('probe_to_reference.txt').resolve()
if not ref.is_file() or not probe.is_file():
  raise FileNotFoundError('Missing input XYZ')
if any(p in (ref, probe) or p.exists() for p in (out, transform_file)):
  raise FileExistsError('Use new output paths')

# Choose MCS only after deciding it is suitable for the scientific comparison.
# Manual alternative: two verified paired-selection strings in the same order.
answers = ['mcs']
with patch('builtins.input', side_effect=answers):
  T = align_xyz_fragment(str(ref), str(probe))

mol = Chem.MolFromXYZFile(str(probe))
if mol is None:
  raise ValueError('Cannot reload probe XYZ')
positions = np.asarray(mol.GetConformer().GetPositions())
aligned = positions @ T[:3, :3].T + T[:3, 3]
if not np.all(np.isfinite(T)) or not np.all(np.isfinite(aligned)):
  raise ValueError('Nonfinite transformation or coordinates')
for i, xyz in enumerate(aligned):
  mol.GetConformer().SetAtomPosition(i, xyz)
Chem.MolToXYZFile(mol, str(out))
np.savetxt(transform_file, T, fmt='%.12g',
           header='probe -> reference; translation in angstrom')
```

```bash
cd "$TASK/work/alignment"
"$VIS2C_PYTHON" "$TASK/scripts/align_task.py" \
  > "$TASK/logs/alignment.log" 2>&1
```

Re-read XYZ, check atom count/order, finite coordinates, rigid-transform orthogonality/determinant, preserved internal geometry and fit RMSD. RMSD is for fitting atoms, not automatically every atom. A transform alone is not the complete fit record.

Do not apply an angstrom translation directly to Bohr cube coordinates. Rotating an XYZ while reusing unchanged field coordinates does not align the field. For cube arithmetic use cuboper's built-in field transform/interpolation; do not automatically run a second separate geometry alignment first.

## 6. cuboper.py: align, interpolate and combine scalar fields

Use `cuboper.py` when the scientific question explicitly requires pointwise arithmetic on real scalar fields after a justified rigid alignment. Substituent and conformer studies are possible applications, not automatic triggers. Compare atom/bond tables, existing numerical summaries, or other observables directly when they answer the question without volumetric subtraction.

Before subtraction, establish compatible field definitions, units, computational conventions, state, and a meaningful reference frame. Record why the fitting region is suitable and how registration affects interpretation. For ESP comparisons, establish compatible potential-reference conventions; for orbital amplitudes, establish orbital correspondence and phase conventions. When these conditions cannot be established, do not present subtraction as a physically resolved property change.

Different composition, electron count, optimized geometries, or fragment definitions can change the meaning of a difference. Explicitly distinguish the requested total change from an attempted isolated substituent contribution. Range and truncated-domain integral statistics describe the processed field and sampling domain; they do not automatically establish charge transfer, interaction energy, or a causal substituent effect.

First input = reference R; second input = probe P. Probe is aligned to reference before sampling and arithmetic. Differences mean aligned probe MINUS reference.

| CLI | Python operation | Formula |
|---|---|---|
| --add | sum | P + R |
| --diff | difference | P - R |
| --avg | average | (P + R)/2 |
| --posdiff | positive_difference | max(P - R, 0) |
| --negdiff | negative_difference | max(R - P, 0) |
| --absdiff | absolute_difference | abs(P - R) |

--negdiff is a nonnegative magnitude, not signed negative data. Choose one operation. Compare fields with compatible physical definitions, methods, fragment conventions and units.

Full CLI shape:

```text
cuboper.py (--add | --diff | --avg | --posdiff | --negdiff | --absdiff)
           REF_CUBE PROBE_CUBE [-o OUTPUT | --output OUTPUT]
           [--chunk-size POSITIVE_INT]
```

Inputs must be distinct existing `.cub`/`.cube` paths. Output must use a cube suffix, have an existing parent, be different from both inputs and not already exist. Without `-o/--output`, the CLI writes `<probe-stem>_<label>_<ref-stem>.cub` in CWD, where label is `add`, `diff`, `avg`, `posdiff`, `negdiff` or `absdiff`. Chunk size defaults to 200000 points. Use descriptive explicit output names in reproducible tasks.

PTY CLI usage example:

```bash
# Run through an execution tool with a real terminal/PTY on stdin.
cd "$TASK/work/cube_ops"
"$VIS2C_PYTHON" "$TRESC/vis2c/cuboper.py" \
  --diff reference_dg.cub probe_dg.cub \
  -o probe_minus_reference_dg.cub --chunk-size 200000
```

Answer prompts with verified paired indices or mcs; save the terminal transcript. The CLI checks stdin.isatty(), so an ordinary pipe or redirected input file supplying mcs is rejected. There are no --mcs/--ref-atoms/--probe-atoms flags.

Without a PTY, use a local single-process API wrapper. First copy inputs into a new writable cube_ops directory under distinct stems, protecting their source-side XYZ sidecars. Example `$TASK/scripts/cube_operation.py`:

```python
import os
from pathlib import Path
import sys
from unittest.mock import patch

sys.path.insert(0, str(Path(os.environ['TRESC']) / 'vis2c'))
from cuboper import operate_cubes

ref = Path('reference_dg.cub').resolve()
probe = Path('probe_dg.cub').resolve()
out = Path('probe_minus_reference_dg.cub').resolve()
if not ref.is_file() or not probe.is_file():
  raise FileNotFoundError('Missing input cube')
if ref == probe or out in (ref, probe) or out.exists():
  raise ValueError('Inputs/output must be distinct and output must be new')
sidecars = [ref.with_suffix('.xyz'), probe.with_suffix('.xyz')]
if sidecars[0] == sidecars[1] or any(p.exists() for p in sidecars):
  raise FileExistsError('Cube-to-XYZ sidecars would collide or overwrite files')

# For a scaffold-specific question, replace with two verified manual selections.
answers = ['mcs']
with patch('builtins.input', side_effect=answers):
  result = operate_cubes(
    str(ref), str(probe), 'difference', str(out), chunk_size=200000,
  )
if not Path(result).is_file() or Path(result).stat().st_size == 0:
  raise RuntimeError('No nonempty output cube')
```

```bash
cd "$TASK/work/cube_ops"
"$VIS2C_PYTHON" "$TASK/scripts/cube_operation.py" \
  > "$TASK/logs/cube_difference.log" 2>&1
```

Numerical and safety limits:

- cuboper already calls alignmol's numerical function; no GUI is opened by that API. Geometry and field are transformed consistently inside this route.
- The API signature is `operate_cubes(ref_cube: str, probe_cube: str, operation: str, output_cube: str, chunk_size: int = 200000) -> str`. It returns the written cube path, not a manifest or transform object. Save its alignment/operation log and independently record the actual pairs/transform if the task needs a reusable transform file.
- It writes same-stem XYZ beside EACH INPUT cube, possibly overwriting existing sidecars. Changing CWD alone is insufficient if original inputs are elsewhere. Use copies as in the example.
- CLI refuses existing output and input overwrite; direct API can replace existing unrelated output, so the wrapper must check first. Do not run concurrent tasks sharing paths or patched input.
- The metadata parser expects standard single-scalar Bohr cubes, positive atom/grid counts, at least two points/axis, and exactly four fields in the origin/grid header lines. Validate external cube compatibility; do not assume all Multiwfn cube variants satisfy it.
- Inputs may have different atom counts/grids, but must admit a meaningful non-collinear fragment fit and full-dimensional sampling overlap.
- Output uses reference geometry/frame and an approximate minimum-volume bounding box of the domain intersection, generally with non-axis-aligned vectors. It is not necessarily the original reference grid.
- Target spacing is the largest input step length; actual axis steps do not exceed it. Both fields undergo first-order interpolation via `map_coordinates(order=1, mode='nearest', prefilter=False)` at shared-domain points. This is not a finest-resolution or conservative remapping guarantee.
- Outside-overlap points inside the box are zero padding for visualization, not physically evaluated zeros. Keep overlap fraction, shape/spacing, range, boundary-value and truncated-domain integral logs.
- --chunk-size limits temporary interpolation batches, not full input/output resident memory.
- Re-read/validate output. Do not interpret the printed truncated-domain integral automatically as full-space electron transfer or interaction energy.

### Signed-difference rendering

A difference cube is not fully represented by one positive real3damp contour. Inspect its range and requested scientific meaning.

Default complete handoff: retain the signed difference cube, and visualize both positive and negative regions with a common justified magnitude c. If one sign is genuinely absent, report it instead of fabricating a surface. Example model exports:

```bash
"${EXPORT_RUNNER[@]}" "$VIS2C_PYTHON" "$TRESC/vis2c/real3damp.py" \
  "$DIFFERENCE_CUBE" --isovalue="$POS_LEVEL" \
  --export-obj "$TASK/models/change_positive.obj" \
  > "$TASK/logs/export_change_positive.log" 2>&1

"${EXPORT_RUNNER[@]}" "$VIS2C_PYTHON" "$TRESC/vis2c/real3damp.py" \
  "$DIFFERENCE_CUBE" --isovalue="$NEG_LEVEL" \
  --export-obj "$TASK/models/change_negative.obj" \
  > "$TASK/logs/export_change_negative.log" 2>&1
```

POS_LEVEL is +c and NEG_LEVEL is -c, not two arbitrary defaults. Check BOTH effective contours after plotter adjustments. If magnitudes differ, choose a valid shared magnitude and re-export or disclose the non-comparable thresholds.

The existing render CLI accepts ONE bundle, not an undocumented --input2. Default to two clearly labeled sign panels. Use a common display frame/scale/camera derived from the union of the models/reference; independent auto-fit must not create misleading size comparisons. For a combined 3D overlay, adapt only the task renderer to validate both bundles, namespace their materials and normalize them together once. Do not merge unverified source files or align them again.

Label signs as descriptor increase/decrease for P-R, not as orbital phase. Color choices and captions must preserve that distinction. A user-requested one-sign plot is allowed if explicitly labeled as partial.

## 7. cub1c.py: Molden -> signed orbital model

### 7.1 Purpose, arguments and execution modes

Use this route for a selected real scalar orbital faithfully represented in a suitable Molden source (`.molden.input` or `.molden`). Possible source orbital types include canonical orbitals (including HOMO/LUMO/SOMO, core, and virtual orbitals), natural orbitals (NO), natural transition orbitals (NTO), natural bond orbitals (NBO), singly occupied natural orbitals (SONO), and localized molecular orbitals. These are examples of source identities, not a claim that the exporter computes or recognizes each analysis automatically.

Verify that the source actually contains the requested coefficients and that basis order, normalization, angular-function conventions, spin labels, and 1-based generator index are compatible. Preserve analysis-specific metadata, including NTO hole/particle pairing and weights when relevant; a rendered orbital alone does not supply those relationships. An unsupported or unverified conversion is a blocker, not permission to relabel another orbital.

`cub1c.py` generates the selected orbital cube and constructs signed contours plus the molecular skeleton. It does not accept a pre-existing cube as input, select HOMO/LUMO automatically, generate the requested orbital family from canonical coefficients, produce electron density by relabeling amplitude, or visualize a complex spinor.

Interface and direct usage example:

```text
cub1c.py MOLDEN INDEX [-slice] [--isovalue FLOAT] [--export-obj FILE.obj]
```

| Argument/option | Meaning, default and constraints |
|---|---|
| `MOLDEN` | Existing `.molden.input` or `.molden` path; a basename is also accepted and tried with `.molden.input` then `.molden`. Use an absolute source path. |
| `INDEX` | Explicit 1-based positive integer followed by A/B, e.g. `35A` or `35b`; case becomes uppercase. A bare integer, HOMO string or zero is invalid. |
| `--isovalue FLOAT` | Positive finite amplitude magnitude; default 0.05. Requested contours are +c and -c; the plotting function may adjust the effective magnitude. |
| `--export-obj FILE.obj` | Export OBJ/MTL/scene JSON and exit without a GUI controller. No overwrite; parent must already exist. Blender is not run. |
| `-slice` | Add interactive sections in the ordinary viewer; incompatible with `--export-obj`. This is one hyphen, not an invented `--slice`. |
| `-h`, `--help` | Print CLI usage and exit without generating scientific data. |

### 7.2 Unattended model export

```bash
cd "$TASK/work/orbital"
"${EXPORT_RUNNER[@]}" "$VIS2C_PYTHON" "$TRESC/vis2c/cub1c.py" \
  "$MOLDEN" "$ORBITAL_INDEX" --isovalue="$ISOVALUE" \
  --export-obj "$TASK/models/orbital.obj" \
  > "$TASK/logs/export_orbital.log" 2>&1
```

Resolve ORBITAL_INDEX and ISOVALUE scientifically. Index is positive and includes A/B, case-insensitive. Input is a full .molden.input/.molden path or a resolvable basename. The script internally performs the TRESC cube stage described above; do not duplicate it.

The generated `<INDEX>-real.cub` is written/updated in the process CWD, independently of the model destination. Use a separate writable directory for each molecule/orbital/run. A pre-existing cube does not cause cube reuse: the generator runs again. The script checks that the result is nonempty and changed rather than accepting an unchanged stale cube, validates 3D finite coordinate/value arrays with at least two points per axis, and checks source hashes during export. Retain the generated cube; a model bundle alone is not the full numerical result.

The scene is constructed with `plot3d.real_orb_plot_cub`. The reviewed pipeline applies Gaussian smoothing with sigma=1 in grid-index units and inherited contour adjustment. This is existing construction behavior, not an invitation to smooth input data for beauty; record it with the actual plot3d identity/hash. For a changed installation, inspect the relevant plot3d implementation if its behavior matters rather than assume an old review describes it exactly.

The model contains explicitly mapped orbital_positive/orbital_negative/molecular_geometry materials. Effective opacity includes the original constant LUT alpha and actor opacity; do not apply alpha twice. Keep requested/effective contours, signed roles and cube provenance. Check nonempty_surface_roles. The plotter can reject unavailable contours; do not manufacture missing lobes.

The exporter identifies signs by the actual IsoSurface contour and VTK actor identity, not actor order or color. It requires both sign modules to be identifiable, but a sign may contribute no polygon geometry; at least one nonempty orbital surface is required. Report actual exported signs and investigate a missing sign rather than assume that a successful command guarantees a complete two-sign picture. Textured actors and unclassified scalar-colored geometry are rejected by this fixed-material route.

-slice is interactive-only and incompatible with --export-obj. Always use export mode for unattended work; the ordinary CLI opens a controller. This route is not a two-component spinor renderer.

### 7.3 Python API and interactive handoff

Public APIs:

```text
cub1c(molden: str, index: str, isovalue: float = 0.05,
      slice: bool = False) -> None
export_obj(molden: str, index: str, isovalue: float = 0.05, *, filename)
```

`cub1c()` is the interactive API. `export_obj()` returns the manifest dictionary after successful publication. Use a task-local wrapper with the common import setup in section 9.4, for example:

```python
from cub1c import export_obj

manifest = export_obj(
  molden_path, orbital_index, isovalue, filename=output_obj,
)
```

Bind those variables to actual validated paths/values and run the wrapper with `"$VIS2C_PYTHON"` in the isolated orbital CWD, with the selected display runner if required. API export generates the cube too and has the same display, path and no-overwrite requirements as the CLI. Python exceptions propagate to the caller; do not suppress them and consume partial files.

Only for a requested user-operated interactive session:

```bash
"$VIS2C_PYTHON" "$TRESC/vis2c/cub1c.py" "$MOLDEN" "$ORBITAL_INDEX" \
  --isovalue="$ISOVALUE"
# Add -slice only when interactive cross-sections are wanted.
```

Changing the controller redraws the viewer; it is not an automated figure-export or cube-only mode. For cube-only work use section 4.3 instead.

### 7.4 Grid accuracy and box-margin control ($TRESC/gridsettings.ini)

The documented orbital-cube generator uses `rgrid` for grid resolution and `edge` for bounding-box margin. Verify the installed parser's syntax, units, allowed values, and actual behavior before adjustment; do not assume that an undocumented value or numeric direction gives the intended spacing. Record both configuration values and the resulting cube's physical spacing and extent.

Preserve current settings unless the task requires a change. When choosing a new margin, normally prefer documented `edge = auto`, but do not reset an existing explicit margin merely because it is not auto. Respect user/project settings and check actual spatial coverage, especially for diffuse orbitals. Finer image pixels or more Blender samples do not improve the numerical orbital grid.

Use the following configuration precedence:

1. Prefer an installed, documented task-local configuration or override mechanism, if one exists. Verify that the generator reads it. Do not invent flags, environment variables, alternate-root behavior, or local-file discovery rules.
2. If only the shared installation file is supported, modify it only within explicit user/project authorization and a reliable exclusive-use arrangement covering all relevant readers and writers, including other sessions. A lock used only by this agent does not protect jobs that ignore it. Do not hold a shared edit across unattended or asynchronously launched work without controlled ownership and recovery.
3. If authorization or exclusivity cannot be established, do not edit the shared file. Use unchanged settings only if scientifically adequate; otherwise report the blocker and request a narrowly scoped decision. Existing validated cubes need no configuration edits.

For an authorized shared edit:

- Save the original complete file bytes and a hash in the task records, preserve permissions, and record the exact intended key changes. Inspect the INI section structure and change only required `rgrid` and/or `edge` entries; retain unrelated values and formatting. Do not replace the file with a generic template.
- Establish cleanup before editing. Use controlled synchronous execution with restoration in a `finally` block or equivalent shell traps for normal completion and catchable errors/signals. Such mechanisms cannot guarantee recovery after SIGKILL, host failure, or power loss; retain the backup, ownership record, and explicit recovery instructions.
- Hold the exclusive-use arrangement until all affected cube-generation processes have terminated and restoration has been verified. Do not run conflicting grid jobs concurrently.
- Before restoration, verify the live file still equals the agent's expected edited content. If another actor changed it, do not overwrite their changes. Preserve both versions, report the conflict immediately, and resolve it within authority.
- With no conflict, restore the original bytes and verify the original hash. Record restoration success or failure explicitly. Never claim cleanup succeeded merely because the analysis finished.
- Record original, requested, effective, and restored configuration states, generator identity, commands, and actual cube-grid metadata in `workflow.json` and `reproduction.md`. Archive an exact copy of the effective configuration for reproducibility.

For quantitative comparisons, keep numerical settings comparable where justified and disclose resampling or extent differences. For surface figures, inspect coverage and contour stability where they affect the conclusion. Do not change a scientific contour or smooth the numerical input merely to hide an inadequate grid.

## 8. real3damp.py: one scalar-field contour -> model

### 8.1 Purpose and full CLI

Use for any scientifically defined real scalar field when one contour is the intended representation: the field need not be an orbital or one of the examples. This route supplies a fixed-color surface, not quantitative colors from a second property.

Interface and usage example:

```text
real3damp.py AMPCUBE [--isovalue FLOAT] [--export-obj FILE.obj]
```

| Argument/option | Meaning, default and constraints |
|---|---|
| `AMPCUBE` | Existing real scalar `.cub`/`.cube` file; suffix check is case-insensitive. Use an absolute path. |
| `--isovalue FLOAT` | Finite signed contour level; default 0.05. Positive, negative and zero values pass the wrapper's finite-value check; actual useful contours and adjustments depend on plot3d and the field. |
| `--export-obj FILE.obj` | Export the single surface and skeleton as an OBJ/MTL/scene JSON bundle; no GUI controller and no Blender invocation. |
| `-h`, `--help` | Print CLI help. |

There is no `-slice`, map-range, colormap, output-image or automatic two-sign option in this script. Without `--export-obj`, it opens an interactive contour controller.

### 8.2 Model export and effective scientific representation

```bash
"${EXPORT_RUNNER[@]}" "$VIS2C_PYTHON" "$TRESC/vis2c/real3damp.py" \
  "$CUBE" --isovalue="$ISOVALUE" \
  --export-obj "$TASK/models/scalar.obj" \
  > "$TASK/logs/export_scalar.log" 2>&1
```

The finite level may be signed; handling/adjustment follows the existing plotting function. This is one contour, not an automatic signed pair or second-field map. Do not universally reuse orbital/density thresholds for IGMH or difference data.

For an electron localization function (ELF) isosurface, use the validated ELF cube directly through this normal `real3damp` route; no ELF-specific preprocessing is needed. Based on the user-tested workflow, explicitly choose an ELF `--isovalue` in the recommended range 0.6-0.8, for example `--isovalue=0.7`, rather than inheriting the generic 0.05 default. Treat this range as a practical starting recommendation, not a mandatory cutoff; respect a specified scientific threshold and record the requested/effective value. Continue with the same model export, Blender rendering and visual review steps.

Exactly one nonempty amplitude surface is identified, with role isosurface; skeleton is molecular_geometry. Explicit obj_group/mtl_material, RGB/opacity and contour are recorded. Negative field contours stay generic surfaces, not orbital_negative. The wrapper passes the full physical x/y/z coordinate arrays to plot3d. For non-axis-aligned cuboper outputs, verify that the installed viskit/plot3d retain those coordinates correctly; the wrapper alone does not certify every external grid convention. Scientific quantity/unit must come from provenance when JSON fields are null.

No scalar-texture mapping is added by this route; ambiguous, textured or unhandled scalar-colored actors are rejected. Export mode skips its controller and does not invoke Blender.

The input cube is loaded through viskit; x/y/z/value arrays must be finite, have identical 3D shapes, and contain at least two points per axis. The scene uses `plot3d.real_func_ampplot_cub`, whose reviewed construction smooths the field with Gaussian sigma=1.0 in grid-index units and retains its own threshold-adjustment rules. Inspect the effective signed contour returned in the manifest, especially for small or negative requested levels. Record inherited processing, and never relabel an adjusted positive contour as the requested negative surface.

The surface is identified by its `amplitude` point-data array, not its color or position in the renderer's actor list. Exactly one nonempty polygonal surface is required; skeleton-only/ambiguous output is refused. Actor RGB and actual opacity become a fixed named material. JSON `quantity` and `field_unit` are null because the script cannot infer scientific identity from the filename: supply these in task provenance and captions without modifying the source bundle to forge metadata.

This route does not infer two signs. To represent a signed difference completely, use the two-exports/common-frame procedure in section 6. If the source is an existing signed orbital cube and no Molden is available, that procedure can retain both signs, but label them as orbital phase and document the paired representation; a single contour is not a complete two-sign orbital picture.

### 8.3 Python API and interactive handoff

```text
real3damp(ampcube: str, isovalue: float = 0.05) -> None
export_obj(ampcube: str, isovalue: float = 0.05, *, filename)
```

The first API opens the controller; the second writes the bundle and returns its manifest. With the import setup in section 9.4, a numerical-to-model wrapper may use:

```python
from real3damp import export_obj

manifest = export_obj(cube_path, isovalue, filename=output_obj)
```

Use `"$VIS2C_PYTHON"`, an appropriate display context, new output paths and normal exception propagation. CLI export progress goes to stderr; stdout is not a mesh/JSON/image stream. Read files and check process success rather than redirect stdout to a supposed model.

For a requested user-operated viewer, omit the export option:

```bash
"$VIS2C_PYTHON" "$TRESC/vis2c/real3damp.py" "$CUBE" --isovalue="$ISOVALUE"
```

## 9. real3d.py: paired fields -> quantitative mapped model

### 9.1 Ordered inputs, arguments and defaults

This generic mapping route applies to any suitable pair of real scalar fields. IGMH with sign(lambda_2)*rho and density surfaces with ESP colors are examples; select definitions from the scientific task.

Interface:

```text
real3d.py AMPCUBE MAPCUBE [-slice] [--isovalue FLOAT]
          [--export-obj FILE.obj]
          [--map-min FLOAT --map-max FLOAT]
          [--colormap NAME] [--map-label TEXT] [--map-unit TEXT]
```

| Argument/option | Meaning, default and constraints |
|---|---|
| `AMPCUBE` | First cube: defines the isosurface and supplies the scene's atom geometry. |
| `MAPCUBE` | Second cube: supplies the scalar used for surface colors; it does not define the contour. |
| `--isovalue FLOAT` | Positive finite surface-defining level; default 0.05. This wrapper does not accept negative/zero contours. |
| `--map-min FLOAT`, `--map-max FLOAT` | Optional paired finite mapping bounds, min < max. If neither is supplied, retain the plotting function's actual range. |
| `--colormap NAME` | Optional Mayavi LUT name; omission preserves the plotter's colormap. A Matplotlib-only name is not necessarily accepted. |
| `--map-label TEXT` | Meaning of the mapped quantity; export default is `Mapping value`, which should normally be replaced with an informative physical label. |
| `--map-unit TEXT` | Unit of the mapped quantity; omitted by default. Record the contour field's unit separately. |
| `--export-obj FILE.obj` | Write complete OBJ/MTL/JSON/LUT-texture bundle and exit without a controller. |
| `-slice` | Add interactive sections of the FIRST field; incompatible with `--export-obj`. It does not export a quantitatively mapped section image. |
| `-h`, `--help` | Print CLI help. |

Both cubes must be existing `.cub`/`.cube` files with finite, matching-shaped x/y/z/value arrays, at least two samples per axis, and matching physical sample positions. The script checks physical grids but does not establish that the two analyses use the same structure/state: the agent must verify that separately.

### 9.2 Export, validation and quantitative transfer

Usage example for an IGMH figure, after validating the actual fields and definitions:

```bash
"${EXPORT_RUNNER[@]}" "$VIS2C_PYTHON" "$TRESC/vis2c/real3d.py" \
  "$DG_CUBE" "$SL2R_CUBE" --isovalue="$DG_LEVEL" \
  --map-min="$MAP_MIN" --map-max="$MAP_MAX" \
  --colormap jet --map-label 'sign(lambda_2)*rho' --map-unit "$MAP_UNIT" \
  --export-obj "$TASK/models/igmh.obj" \
  > "$TASK/logs/export_igmh.log" 2>&1
```

Colormap is illustrative. Omit the unit option when unknown/inapplicable and disclose this; never invent it. Label the actual selected delta_g variant in provenance/caption.

- First cube determines contour geometry, second determines colors. Contour must be positive and finite.
- --map-min/--map-max must be supplied together, finite, min < max. If omitted, read the actual range retained by the plotter; it is not guaranteed to be raw min/max or a suitable comparison range.
- Array shapes AND physical x/y/z points must match within 1e-5 in the loader's Bohr convention. No automatic resampling is performed. Do not reshape or relabel a grid to force compatibility.
- Existing amplitude smoothing is preserved; the mapped array is not smoothed by this export path. Record processing and effective contour.
- The mapped role is mapped_isosurface, not a fixed-color isosurface. Array name phase does not make the property quantum phase.
- -slice is interactive-only and incompatible with model export.

Complete bundle, with exact texture basenames obtained from the manifest:

```text
scene.obj
scene.mtl
scene.scene.json
scene.<mapped-material>.color.png
scene.<mapped-material>.alpha.png   # present only when LUT alpha is nonopaque
```

Keep all manifest-referenced files together. Do not rename a single member or discard apparently auxiliary textures. Scene schema is version 3, with scene_type `real3d_mapped_isosurface`; the scalar-to-color definition is recorded in `color_mapping` and the relevant actor record. Fixed-material cub1c and real3damp manifests have schema versions 1 and 2 respectively; do not require all exporters to have one identical schema version.

Mapping semantics preserved by exporter/renderer:

- Vertex scalar is encoded into OBJ UV as `u = (0.5 + (N-1)*(value-min)/(max-min))/N`, with v=0.5.
- Do not clip vertex UVs; interpolate first and clamp at sampling through Linear/EXTEND, not REPEAT.
- color_texture has LUT RGB and opaque PNG alpha. alpha_texture, if present, stores LUT alpha in RGB with opaque PNG alpha; use Non-Color grayscale, not its Alpha output.
- Surface opacity is actor opacity times sampled LUT alpha once. Do not discard or duplicate either factor.
- A mapped actor's white color placeholder is not the intended single-color appearance.
- Preserve both field identities, actual range, LUT, label/unit, and below/above-range vertex counts where relevant. A descriptor color interpretation is not a quantitative interaction energy.

LUT RGB/alpha is sampled into lossless 8-bit PNG textures. Transfer retains the numerical range and scalar interpolation convention with that finite LUT precision; it does not promise pixel-identical VTK and Blender lighting or optical transparency. Arbitrary pre-textured actors, direct-RGB scalar arrays, nonlinear LUTs, special above/below-range colors and unclassified scalar-colored skeleton actors are not silently translated by this exporter. Diagnose the actual incompatibility rather than remove quantitative coloring.

For ESP mapping, use the validated electron-density cube as AMPCUBE (FIRST input) and the ESP cube as MAPCUBE (SECOND input). Generate both cubes at matching grid quality/resolution and verify matching physical grids, not merely identical nominal quality settings. Explicitly pass `--isovalue=0.001` for the density-defined surface as the recommended starting value in standard density atomic units; do not use the generic 0.05 default or copy the IGMH contour. Set the ESP mapping range, label and unit independently. As established for the user-tested workflow in section 4.6, ions follow the same route as neutral molecules; do not introduce special ion preprocessing or change the source charge. A series should use comparable surface definitions, mapping ranges, LUT and view/scale conventions; disclose any intentional differences.

### 9.3 Python API and interactive handoff

```text
real3d(ampcube: str, mapcube: str, isovalue: float = 0.05,
       slice: bool = False, *, map_min=None, map_max=None,
       colormap=None, map_label=None, map_unit=None) -> None

export_obj(ampcube, mapcube, isovalue=0.05, *, filename,
           map_min=None, map_max=None, colormap=None,
           map_label=None, map_unit=None)
```

The first is an interactive API; the export API returns the manifest. Example wrapper body after the common setup below:

```python
from real3d import export_obj

manifest = export_obj(
  amplitude_cube, mapping_cube, isovalue, filename=output_obj,
  map_min=mapping_min, map_max=mapping_max,
  colormap=lut_name, map_label=quantity_label, map_unit=quantity_unit,
)
```

Bind variables to actual task settings; optional values may be None. API export shares CLI validation, no-overwrite and display requirements and does not accept `slice`. Export progress is not a serialized result on stdout. A requested user-operated viewer uses the same positional cubes/options with `--export-obj` omitted; add `-slice` only for the intended interactive first-field sections.

### 9.4 Shared exporter API setup, bundle handling and limitations

CLI export is normally sufficient. Use a task-local Python wrapper when programmatic orchestration or direct manifest access is useful. Prepend this to the chosen API example above, save it under `$TASK/scripts/`, and execute it with `"$VIS2C_PYTHON"` and the same display runner as the corresponding CLI:

```python
import os
from pathlib import Path
import sys

sys.path.insert(0, str(Path(os.environ['TRESC']) / 'vis2c'))
# Bind the actual absolute input/output paths and property settings here.
# Then import and call only the exporter needed by this task.
```

All three modules expose `main(argv=None)` for normal CLI argument dispatch, in addition to the public viewer/export functions above. Passing `['--help']` to main prints help and exits; it is not model generation. Prefer public export APIs over private `_...` helpers. The identical function name `export_obj` belongs to different modules; use qualified imports or aliases if a wrapper needs more than one exporter.

Shared operational rules:

- Model output suffix must be lowercase `.obj`; basename/stem must contain no whitespace or `#`. Parent directories must already exist. Use distinct new prefixes for separate contours, molecules or mapping definitions.
- Each fixed-color export writes `PREFIX.obj`, `PREFIX.mtl`, `PREFIX.scene.json`. The mapped export also writes its LUT PNGs. Existing bundle destinations are refused; a collision involving any texture at publication also fails. Do not remove a previous result to reuse its prefix during refinement.
- Mayavi export explicitly uses `offscreen=False`, skips GUI controllers, and closes its own scene/restores the previous offscreen setting afterward. A functioning display backend remains necessary. No export call produces a final Blender image.
- The OBJ serializes transformed polygon geometry, triangles and per-point normals. Existing normals are retained; missing normals are generated for shading. Source scientific volumes are not replaced by mesh files. Axes widgets and interactive controls are not exported.
- Fixed materials use standard MTL `d` opacity (0 transparent, 1 opaque), not an assumed inverse `Tr`. Per-actor `obj_group`/`mtl_material` roles are explicit, unique within the bundle and unrelated to actor order or unverified atom/bond IDs.
- Preserve OBJ, MTL, manifest and referenced textures intact. Manifest inventories record file sizes/hashes; exporter and renderer input hashes protect provenance, not scientific correctness. Files listed in `files` need not include the manifest itself; record its hash separately for reproduction.
- CLI scientific/execution errors exit nonzero (the exporter CLIs use exit code 2). API errors raise exceptions. Check process completion, logs, expected nonempty files and manifest semantics before rendering; do not consume partial or stale outputs just because one file exists.
- Quantities, units, selected orbital/state and field pairing need independent task provenance. A successful parser or finite mesh does not establish them. Field units can remain null in generic manifests and must not be inferred from a filename.
- Reuse a valid bundle for lighting/camera/material iteration. New contours or scalar mapping ranges require a correct new export or an explicitly validated numerical remapping; changing a caption/manifest alone does not change the represented field.
- No `--no-gui`, `--output-image`, `--input2` or other undocumented flag is implied. Choose the documented export route, numerical-only generator or a scoped task extension according to the requested endpoint.

## 10. render_obj.py: native Blender reconstruction and adaptive rendering

### 10.1 Invocation, inputs and output files

All three exporter routes are supported natively in manifest mode, including mapped textures and the quantitative colorbar. No adapter-development phase is needed for an ordinary task.

Prepare the task copy once:

```bash
cp "$TRESC/vis2c/render_obj.py" "$TASK/scripts/render_task.py"
```

Execute that copy through `"$BLENDER_BIN" --background --python ... -- ...`. Blender command-line flags go BEFORE the `--`; renderer script options go AFTER it. Do not run this script with `"$VIS2C_PYTHON"`: actual rendering uses Blender's `bpy` runtime and bundled Python. `--python-exit-code 1` makes Python failures visible as a failing Blender process. Use `--factory-startup` as in the examples so task results do not depend on a saved personal startup scene.

The template supports Blender 3.0+ through version-specific OBJ import code and Principled shader socket lookup. Use the provided executable and its resources; if a task-local edit introduces an API error, fix that edit against the actual runtime rather than upgrade software speculatively.

Baseline usage example; MODEL is one validated complete bundle:

```bash
"$BLENDER_BIN" --background --factory-startup \
  --python-exit-code 1 --python "$TASK/scripts/render_task.py" -- \
  --input "$MODEL" --output "$TASK/previews/baseline.png" \
  --material-mode manifest --source-color-space srgb \
  --width 800 --height 640 --samples 24 --threads "$THREADS" \
  --azimuth 45 --elevation 25 \
  > "$TASK/logs/render_baseline.log" 2>&1
```

Native roles: orbital_positive, orbital_negative, isosurface, mapped_isosurface and molecular_geometry. The renderer verifies names, OBJ/MTL and referenced texture hashes; checks mapped UV presence; loads/packs LUT textures with explicit color convention; handles alpha and recorded mapped backfaces; records the resulting materials and image checks.

File and mode options:

| Option | Default | Meaning and constraints |
|---|---|---|
| `--input PATH` | Required | Existing OBJ with its referenced MTL and, for vis2c exports, complete matching manifest/textures. Input suffix check is case-insensitive. Accepts ONE bundle. |
| `--output PATH` | Required | PNG image path. Also writes the JSON record at the full image name plus `.json`, e.g. `figure.png.json`. Parent directories must exist; existing image or sidecar is refused. |
| `--scene-json PATH` | Input basename + `.scene.json` | Optional explicit manifest path; if supplied it must exist. It must still match the actual OBJ/MTL/textures and is not a way to borrow another model's metadata. |
| `--save-blend PATH` | No editable scene saved | Also save an editable `.blend`; parent must exist and destination must be new. Loaded LUT textures are packed. |
| `--material-mode auto\|manifest\|imported` | `auto` | Prefer `manifest` for vis2c exports. Detailed validation/fallback behavior is below. |
| `-h`, `--help` | — | Print renderer options; invoke through Blender with the script's argument separator. |

View, framing and image size:

| Option | Default | Meaning and constraints |
|---|---|---|
| `--azimuth FLOAT` | 45 | Degrees around source +Z, measured from source +X; finite. |
| `--elevation FLOAT` | 25 | Degrees above source XY; strictly -89 < elevation < 89. |
| `--margin FLOAT` | 1.15 | Finite framing factor > 1; larger values give more space around the model. |
| `--width INT` | 1200 | Positive PNG width in pixels. |
| `--height INT` | 1000 | Positive PNG height in pixels. |
| `--transparent` | Off | Transparent camera background instead of white; does not remove environment illumination or change surface opacity. |

Surface appearance, lighting and quantitative legend:

| Option | Default | Meaning and constraints |
|---|---|---|
| `--source-color-space srgb\|linear` | `srgb` | Source convention for reconstructed manifest RGB/LUT colors. Imported materials are not reconverted. Select by data convention, not perceived brightness. |
| `--surface-roughness FLOAT` | 0.36 | Finite [0,1], applied to reconstructed surface materials. Geometry/skeleton preset is 0.40. |
| `--surface-specular FLOAT` | 0.30 | Finite [0,1], applied to reconstructed surfaces using the appropriate Principled specular socket. Skeleton preset is 0.25. |
| `--surface-shadows transparent\|physical` | `transparent` | Surface shadow-ray policy: bypass surface shadow attenuation or retain shader shadows. Camera-facing opacity remains independently defined. |
| `--world-strength FLOAT` | 0.20 | Nonnegative finite environment illumination. White camera background strength remains 1. |
| `--key-energy FLOAT` | 350 | Nonnegative finite key area-light energy; zero omits that light. |
| `--fill-energy FLOAT` | 70 | Nonnegative finite fill area-light energy; zero omits it. |
| `--rim-energy FLOAT` | 0 | Nonnegative finite rim area-light energy; disabled by default. |
| `--colorbar-ticks INT` | 5 | Quantitative mapped colorbar ticks; integer 2–11, including endpoints. |
| `--no-colorbar` | Off | Omit the mapped legend only when justified by the task/shared-legend layout. |

Quality and resources:

| Option | Default | Meaning and constraints |
|---|---|---|
| `--samples INT` | 64 | Positive Cycles sample count; start lower for previews and raise for the selected final image. |
| `--threads INT` | 4 | Positive fixed CPU rendering thread count. Agent chooses within actual resource allocation. |
| `--no-denoise` | Off | Disable denoising; denoising is otherwise enabled. Compare real details if denoising changes thin surfaces. |

All numeric rendering values must satisfy the parser's finite/range checks. These are software defaults, not a scientific or aesthetic prescription. The reference uses CPU Cycles; there is no GPU/device CLI option. It uses orthographic projection; perspective, camera roll, an arbitrary pose, target/zoom controls, added annotations and multiple inputs require explicit task-copy code rather than guessed flags.

manifest is preferred for vis2c-generated models. auto is not a way around failed hashes or unsupported roles. imported deliberately skips reconstruction; mapped controls and the native quantitative colorbar are not restored automatically there. --scene-json otherwise defaults to the OBJ basename with .scene.json.

In detail, `manifest` requires a manifest with explicit `mtl_material` mappings and valid source-file hashes. `auto` can retain imported materials when there is no manifest or no explicit material map, but once mapped entries are present it validates them and raises on mismatches/unsupported roles rather than silently falling back. `imported` explicitly skips reconstruction. It is useful only when retaining an importer result is the intended representation; do not choose it to hide an invalid vis2c bundle.

The renderer checks every imported material against the manifest's supported roles, loads verified local LUT PNGs, requires active UVs on mapped meshes, and rebuilds shaders. A single-color MTL does not replace a mapped LUT. Fixed RGB is converted once according to the source convention; mapped color textures use that convention, and LUT alpha textures use Non-Color grayscale. Actor opacity and LUT alpha are multiplied once, then mixed between white Transparent BSDF and opaque Principled BSDF; no physical refraction or volume absorption is implied.

`build_scene` centers the combined imported bounds and uniformly scales the largest extent to 4 display units; it records the source-to-display matrix. A camera is placed 12 display units along the requested direction, aimed at the origin, and framed orthographically using the aspect ratio/margin. Lights are placed relative to camera right/up/view directions. These display transforms do not alter source OBJ coordinates, but independently auto-fitting two models can change apparent size, so a comparison must deliberately share scale and camera conventions.

The template uses an RGBA 8-bit PNG, fixed seed 0, Standard view transform, neutral look where available, exposure 0 and gamma 1. It records the active color settings rather than assume all installations behave identically. Changes to these display conventions belong in the task copy and reproduction record.

### 10.2 Quantitative colorbar and mapping fidelity

The mapped colorbar is ON by default; default tick count 5, allowed 2-11. It uses the same LUT/range unlit, with low at bottom, high at top, label and unit. It reserves right-side frame space by camera adjustment, not molecular displacement. Example optional control is `--colorbar-ticks 7`.

Use --no-colorbar only for an explicit request or a documented shared-legend layout. Do not remove the only quantitative legend for aesthetics. One mapped colorbar is supported per figure; multi-mapping layouts are an explicit task-local extension.

Non-mapped scenes do not acquire a quantitative legend merely by setting tick options. Mapped legend colors are the unlit LUT without surface illumination or actor opacity: this is intentional. Numeric ticks are generated across the actual mapping bounds with distinct labels. Check long labels, units, endpoint direction, tick spacing and available frame area in the actual image. If a label/layout needs improvement, adapt only the task renderer while retaining the same LUT/range and truthful labels. Mapping-range changes generally call for a new export; a different title or gradient alone does not remap surface UVs.

White camera background and environment illumination are independent. Source sRGB/linear choice is a technical convention, not a saturation knob; avoid double conversion. Alpha is applied once in transparent/opaque shader mixing. Shadow bypass is a scientific-illustration choice, not a physical optical prediction. Uniform model normalization is for display; the source coordinates remain in the OBJ.

The renderer does not recover a Mayavi camera or automatically find the best chemical view. New palette/opacity controls, arbitrary camera poses or combined-bundle rendering require justified task-local code, not imaginary CLI flags.

### 10.3 Task-local visual autonomy and extension points

Copy the template; never edit shared scripts. Freely adapt camera, projection, lights, roughness/highlights, background, framing, layout and quality to the scientific brief. Existing defaults are starting points, not guaranteed attractive settings.

Preserve palette/effective opacity by default; when the user grants redesign, record source and final values and keep scientific role/phase conventions consistent. Quantitative maps must retain a valid scalar range/LUT/legend. Changing transfer settings is a representation decision, not hidden image beautification.

Validate source bundles BEFORE approved rendering overrides. Apply overrides after validation by verified roles/materials in the task copy. Do not edit source MTL/JSON/textures, forge hashes or weaken safety checks. Update output notes if colors/alpha/shaders change. Keep comments/docstrings English and Python block indentation two spaces; inspect diffs and syntax before execution.

Useful existing implementation locations in the task copy:

| Function | Role in an intentional task-local adaptation |
|---|---|
| `parse_args` | Actual supported CLI/defaults; add a real validated option here if a new control is needed. |
| `load_manifest`, `material_dependencies` | Source/mapping integrity; retain checks before overrides. |
| `import_obj`, `rebuild_materials` | Verified role-based shader construction and LUT/alpha handling; keep quantitative semantics. |
| `build_scene`, `make_world` | Camera, uniform display transform, lighting, CPU Cycles and color-management presets. |
| `add_colorbar`, `_colorbar_values` | Legend placement, text and numeric ticks; preserve matching LUT/range. |
| `inspect_png`, `publish`, `main` | Pixel/file validation, no-overwrite publication and orchestration; preserve error propagation and provenance. |

The renderer's `main()` reads script arguments; it is not an ordinary-Python image-rendering API with `input=`/`output=` keyword parameters. Use Blender to execute the copied script. When adding an explicit new display feature, update its task-local parser/configuration and metadata as needed so the executed settings remain reproducible. Do not duplicate native mapped-material/legend support before checking the existing functions.

For a combined signed overlay or multi-bundle panel layout, validate every complete bundle against its own manifest, namespace object/material names, retain their common physical frame, and compute display normalization/camera once from the combined bounds. Keep signs, contour magnitudes, field identities and any shared legend explicit. The native CLI handles one input; such an extension is code in the task copy, not a second `--input` or an imaginary `--input2`.

Do not change scientific geometry, source normals/fields, atom positions or contour labels for beauty. No lobe deletion, nonuniform stretching or scientific-data smoothing. Distinguish approved detail crops, backface/shadow changes and other illustration choices from physical conclusions.

### 10.4 Preview and final QA

Actually open each image claimed as reviewed. Logs, JSON and pixel extrema are not vision. Without image-reading tools, deliver limited candidates as unreviewed rather than claim autonomous optimization.

Suggested budget: one baseline, at most two refinement rounds with up to three candidates each, then a selected final render. Use low-cost previews, one process at a time, configured threads, bounded runtime and an explicit total budget for multi-panel/multi-molecule tasks. Do not endlessly rerun or repeatedly regenerate numerical data.

Search view and material changes separately when practical, preserve the best candidate, and record concrete reasons. Check scientific identity and labels first, then complete geometry/visible ROI, correct phase/scalar encoding, skeleton contrast, glare/haze/noise, framing and legend legibility. A nonblank colorbar can make pixel checks pass with a missing model: inspect the MODEL region too.

Final usage example with the actually selected settings replacing these illustrative values:

```bash
"$BLENDER_BIN" --background --factory-startup \
  --python-exit-code 1 --python "$TASK/scripts/render_task.py" -- \
  --input "$MODEL" --output "$TASK/final/figure.png" \
  --save-blend "$TASK/final/figure.blend" \
  --material-mode manifest --source-color-space srgb \
  --width 2400 --height 1920 --samples 128 --threads "$THREADS" \
  --azimuth 45 --elevation 25 \
  > "$TASK/logs/render_final.log" 2>&1
```

Reinspect the final image; denoising, resolution and sampling can change appearance. Keep PNG plus renderer JSON, task code/config, logs, selected settings and optional .blend. Native LUT textures are packed; ensure any custom resources also remain portable. Do not equate file/pixel validation with visual/scientific approval.

The renderer stages outputs beside their intended destinations and refuses any existing PNG, `PNG.json` or requested `.blend`. On success, the JSON records Blender version, input/model/MTL/texture/script hashes, rebuilt roles/materials, source/display transform, camera, lighting, mapping/colorbar, image settings, decoded pixel checks and output hash. It rechecks source OBJ/MTL/texture hashes before publication. It does not save the original cube data into the image or .blend; preserve the numerical sources and source bundle separately.

File checks cover completed rendering/saving, decoded dimensions/channels, finite pixels and entirely black/white/transparent failures. They do not prove that a lobe, weak interaction region, legend or model is visible/correct. Keep template metadata's `visual_review_performed=false` as an execution fact; save the agent's actual image review separately in `render_review.json`, tied to the exact opened image hash. A modified final render needs its own review.

## 11. atomic.py: atomic scalar-property data handoff

This script's current role in this workflow is USER-OPERATED visualization. Agent prepares validated XYZQ; it must not assume --export-obj, --column or a noninteractive atomic model exporter exists.

Its CLI takes one positional input path and then prompts for the 1-based property column. Accepted suffixes in the supplied script are `.xyzq`, `.txt`, `.chg`, `.dat`, `.xyz`, `.pdb` and `.pqr` (case-sensitive checks). A `.xyz` suffix does not make a standard XYZ with no property columns adequate: the actual data still need property values. Prefer an explicit clean `.xyzq` to avoid ambiguous extension-based interpretation.

XYZQ format, one row per atom:

```text
ELEMENT X Y Z PROPERTY_1 [PROPERTY_2 ...]
```

Coordinates are angstrom, ELEMENT is a symbol/atomic number, not a serial index. Use the exact source atom order and geometry/state. No XYZ header, comments, blank lines, row numbers or embedded labels. Require at least one property, fixed column count and finite values. Use ordinary decimal/E notation and never fill missing scientific data with silent zero/NaN.

Format example ONLY, not user calculation data:

```text
C  0.000000  0.000000  0.000000  -0.120000  0.015000
H  0.000000  0.000000  1.080000   0.120000 -0.015000
```

User-run usage example, not an unattended agent command:

```bash
"$VIS2C_PYTHON" "$TRESC/vis2c/atomic.py" "$TASK/data/atomic_properties.xyzq"
```

Tell the user which column to select. Column 1 means the first PROPERTY after XYZ, not the first text column. Supply column names, units, population method, source/state and atom count separately. Distinguish Mulliken/Lowdin/Hirshfeld/NPA and vector component/magnitude definitions.

Validate every row; the parser can stop at an invalid/changed row and leave only part of the molecule. Prefer clean XYZQ over assuming PDB/PQR occupancy/temperature fields are charges. The CLI prompts and calls mlab.show(); do not execute it unattended. Do not label this data-only subtask as a completed rendered image.

The PDB/PQR branch uses fixed PDB columns and writes a same-stem `.xyzq`, potentially replacing an existing file. It requires sufficiently long ATOM/HETATM records, obtains the element from columns 77–78, and copies columns 31–72 into coordinates/property data. This is not a general semantic PQR charge/radius parser. If such a source is supplied, the agent should explicitly extract the intended atomic property and coordinates into validated XYZQ instead of relying on this convenience conversion.

The viewer uses a zero-centered range based on 1.3 times the largest absolute property value, applies its own diverging LUT, and draws connections by the covalent-distance criterion. It offers no CLI label/unit/range/colormap controls; provide the correct property definitions and verify the actual legend/sign convention rather than infer them from a code comment. An all-zero property set may give a degenerate range. Public `datainline` parses a row, but `pop_plot` is graphical and is not an unattended numerical-export API.

## 12. bond.py: pair-property matrix and XYZ handoff

Agent prepares data; the current GUI script is for the user's subsequent visualization, not a presumed OBJ export interface.

The exact CLI is `bond.py GEOMETRY.xyz MATRIX.{txt,dat,mat,bnd}` with two positional paths; the supplied suffix checks are case-sensitive. There are no column-selection, manual bond-list, colorbar-title or model-export options. `bond_plot(...)` is a graphical API ending in `mlab.show()`, not a hidden unattended exporter.

Deliver:

- geometry.xyz: standard XYZ with N, comment, N element/angstrom-coordinate records.
- bond_property.txt: N x N whitespace-separated numeric matrix without headers/indices/commas/comments.
- Property definition/unit/method/state, exact atom order, diagonal and missing/nonapplicable-entry conventions.

A list `i j value` must be converted to the required atom-ordered square matrix. The script reads the strict LOWER triangle, so upper-triangle-only output is insufficient. Symmetric properties should be symmetric; do not silently symmetrize directional data.

Format example ONLY:

```text
0.000000 0.910000 0.880000
0.910000 0.000000 0.020000
0.880000 0.020000 0.000000
```

User-run usage example:

```bash
"$VIS2C_PYTHON" "$TRESC/vis2c/bond.py" \
  "$TASK/data/geometry.xyz" "$TASK/data/bond_property.txt"
```

Check dimensions, finite values, units/order, expected symmetry and numerical readback with the installed viskit.load_matrix implementation, without invoking bond.py. Do not assume arbitrary paginated analysis output is directly readable.

Unknown is not zero. Preserve complete computed bond-order data; real pair distances may populate a distance matrix. For sparse energies, ensure all pairs the distance-based display would draw have valid data. Document placeholders only for nonapplicable pairs and never let them masquerade as measured zeros.

Connections are selected by distance <= 1.15 times the covalent-radius sum, not by any nonzero matrix value. Flag long coordination bonds/close contacts. No displayed edges or identical values can cause rendering/colorbar issues; do not add fake values. The supplied GUI's colorbar title is Bond Order; for other properties ensure the user's final labels/units are corrected. The CLI calls mlab.show(); do not run it as an unattended figure generator.

The viewer colors displayed connections with its `Greys` LUT and takes min/max from the values of those distance-selected connections, not necessarily all matrix entries. Preserve this distinction when preparing comparison data: correct matrices alone do not impose a shared visual range across separate viewer sessions. It reads the geometry through `viskit.load_geometry_xyz` and matrix through `viskit.load_matrix`; inspect actual numerical readback without invoking the viewer. An agent can prepare the data fully, but automatic atomic/bond rendering through Blender is a separately implemented representation extension.

## 13. Connected scientific playbooks: examples, not scope limits

These describe execution plans, not permission to invent missing numerical recipes. Reuse completed stages and bind variables to actual task files/settings. Every invoked script's concrete usage is given in its section above.

Keep these examples as guides to scientific reasoning. For a different wavefunction analysis, construct the corresponding connected plan from the requested quantity, the relevant Multiwfn workflow and the appropriate representation; do not require it to match one of the example titles.

### A. HOMO of one molecule

1. Establish structure, charge/spin and computational model. Reuse a converged suitable wavefunction or run the authorized calculation and export Molden using the program's validated workflow.
2. Determine HOMO identity/spin from energies/occupations and map it to the TRESC/cub1c index. Resolve open-shell ambiguity rather than guessing.
3. Run the cub1c export example in an isolated orbital directory. This performs cube generation internally; retain that cube and provenance.
4. Set MODEL to the exported orbital OBJ. Copy render_obj.py, run the baseline example in manifest mode, inspect both signs and the skeleton, refine the task copy if needed.
5. Final-render the selected settings. Deliver signed-phase labels, actual contour, source identity, figure and reproduction records.

```text
Geometry/results -> authorized QC + Molden if needed
-> explicit HOMO/spin selection -> cub1c --export-obj
-> OBJ/MTL/scene.json -> render_task.py previews -> final image
```

### B. Substituent or conformational effect on real-space descriptors (e.g., methyl-substitution IGMH comparison)

1. Define reference R (e.g., parent molecule or reference conformer) and probe P (e.g., substituted molecule or target conformer), unless the user specifies otherwise. State that the requested difference is P-R.
2. Obtain comparable wavefunctions using consistent computational choices. Define the exact IGMH dg variant and chemically corresponding fragment partitions for BOTH molecules. Define the shared-scaffold fitting region separately.
3. Run verified Multiwfn recipes in separate directories to generate the two comparable dg cubes. Match resolution/coverage as far as practical; cuboper can handle different physical grids but cannot make different definitions comparable.
4. Copy them to `$TASK/work/cube_ops/reference_dg.cub` and `probe_dg.cub`; protect XYZ sidecars. Run the cuboper CLI with a PTY or its local API wrapper. Inspect/save the actual atom mapping, RMSD and overlap/spacing/boundary logs. Do not run a redundant alignmol CLI first.
5. Validate the signed difference cube. Positive means increased selected descriptor in aligned P relative to R; negative means decreased. Do not call it an interaction-energy difference.
6. Export positive and negative contours of the SAME difference cube with real3damp, using matched actual magnitudes. No second alignment is needed. For a one-sign request, explicitly label partial coverage.
7. Render the two sign panels in a common frame/scale or implement a validated combined-bundle task renderer. Distinguish increase/decrease colors and captions; do not call them orbital phases. Inspect whether apparent features are boundary/registration artifacts or direct substituent-density contributions.
8. Deliver difference cube, alignment/definition records, sign-complete figure(s) and interpretation limits, including geometry/electronic effects and truncated sampling.

```text
Parent/P wavefunctions -> matching Multiwfn dg definitions
-> cuboper --diff R P (built-in fragment alignment and resampling)
-> signed dg difference cube -> real3damp at +c and -c
-> common-frame Blender sign panels or validated overlay -> final figure(s)
```

An integrand comparison or attractive-looking isosurface does not by itself isolate the methyl group's energy contribution. If a stronger causal claim is needed, propose a separately scoped controlled calculation/analysis instead of inferring it from colors.

### C. IGMH dg surface colored by sl2r

1. Obtain a suitable wavefunction and choose the intended dg variant/fragment partition.
2. Use a verified Multiwfn numerical export to generate dg and sign(lambda_2)*rho cubes on identical physical grids. Do not confuse a total, intra- or interfragment dg file.
3. Set DG_CUBE and SL2R_CUBE accordingly; select DG_LEVEL, a justified MAP_MIN/MAP_MAX, label and unit. For a series, use consistent thresholds and mapping conventions unless explicitly disclosed.
4. Run the real3d example: dg FIRST, sl2r SECOND. Retain all textures and effective range/contour metadata.
5. Set MODEL to the mapped OBJ; run the native renderer in manifest mode with its colorbar enabled. Inspect the model region, full mapped color variation, LUT alpha, legend endpoints/direction/label, and any out-of-range saturation.
6. Refine only task-local display choices while preserving quantitative mapping. Final-render and deliver source definitions, complete model bundle, image/legend and reproduction records.

```text
Wavefunction -> Multiwfn dg + sl2r on the same physical grid
-> real3d --export-obj (geometry + quantitative LUT/UV/textures)
-> native manifest rendering + matching unlit colorbar -> final image
```

### D. Electrostatic potential mapped onto a density surface

1. Establish the intended structure/electronic state and a wavefunction adequate for the requested density and ESP analysis. Reuse valid source data when available.
2. Consult the relevant online/local Multiwfn manual sections and execute verified numerical exports through `"$MULTIWFN_EXE"`. Generate density and ESP on the same physical grid and retain both fields' definitions and units.
3. Choose the density surface level, ESP color range, LUT and ESP label/unit from the scientific purpose. A density-contour unit is not the ESP unit. For comparisons, use consistent surface/mapping definitions and record any deliberate exceptions.
4. Set the following variables to the actual task data/settings and export:

```bash
"${EXPORT_RUNNER[@]}" "$VIS2C_PYTHON" "$TRESC/vis2c/real3d.py" \
  "$DENSITY_CUBE" "$ESP_CUBE" --isovalue="$DENSITY_LEVEL" \
  --map-min="$ESP_MIN" --map-max="$ESP_MAX" \
  --colormap "$ESP_COLORMAP" \
  --map-label 'Electrostatic potential' --map-unit "$ESP_UNIT" \
  --export-obj "$TASK/models/esp.obj" \
  > "$TASK/logs/export_esp.log" 2>&1
```

5. Validate the complete mapped bundle and set MODEL to `$TASK/models/esp.obj`. Use the section 10 Blender workflow in manifest mode, retaining the matching quantitative colorbar.
6. Inspect actual surface color variation, saturation, sign labels, visible chemical regions and legend legibility. Refine the task-local view/lighting/layout, final-render, and deliver both source fields, bundle, figure and reproduction records.

```text
Wavefunction -> Multiwfn density + ESP on matching physical grids
-> real3d (density defines geometry; ESP defines color)
-> native Blender material reconstruction + ESP legend -> final image
```

Use the same reasoning for any other descriptor-mapped surface: identify the two quantities independently and preserve their definitions through numerical generation, export, shading and labels. A shared representation primitive does not make their scientific interpretations interchangeable.

## 14. Safety, provenance and honest completion

Use unique paths, preserve source files and project authorization boundaries, and check each stage before passing output forward. Python API examples use local patched input in a single non-concurrent process; do not patch input in a multithreaded service.

Record actual executed commands and software identities, not proposed commands as completed work. Preserve numerical input/output hashes; calculation model/state/convergence; orbital identity; Multiwfn recipes; field/fragment definitions; grid units/shape/extent; alignment pairs/transform/RMSD; operation formula; interpolation/padding limits; effective contours; quantitative transfer/legend; task-render script/settings; and actual review status.

Do not interpret an OBJ or picture as a replacement for source data. For mapped transferable .blend files, verify texture packing or document remaining dependencies. Keep image plus renderer JSON together. Source metadata may initialize visual_review_performed=false; actual agent review belongs in render_review.json tied to the image genuinely opened.

The scripts publish staged outputs using no-overwrite hard links; chosen output storage must support them. Wait for successful process completion before consumption; multi-file publication is not atomic to concurrent readers.

Failure handling:

- Numerical result absent/invalid: stop downstream steps, inspect logs and explain missing evidence.
- cuboper alignment/overlap failure: fix selection/definitions/coverage, not by distorting structures or changing labels.
- Display/backend failure: use the agent-selected available working display/runner; Blender background does not fix the Mayavi stage.
- Hash/material/texture mismatch: locate the intact bundle or deliberately re-export; never forge metadata.
- Native mapped-rendering issue: inspect actual UV/LUT/alpha/legend and logs, not a single-color fallback.
- Color-management resource issue: resolve deployment configuration, not exposure/light compensation.
- Poor appearance: inspect actual previews, isolate meaningful changes and stay within budget.
- Unsupported atomic/bond/spinor/section rendering: deliver valid data and a precise limitation or separately scope an extension. Do not pretend existing flags implement it.

Final report should state:

```text
Scientific task and definition:
- Structure/state, field/orbital/spin, fragment selection or comparison formula.

Executed stages and validated data:
- Calculation/export, numerical preparation, alignment/arithmetic, model export.

Figures and visual decisions:
- Actual effective contours/ranges, views, palette/opacity/legend choices.
- Exact files, render script/configuration, editable scene if produced.

Verification and limitations:
- Separate numerical, file/pixel and actual visual checks.
- Missing stages, overlap/interpolation limits, partial signs or other caveats.

Reproduction:
- Commands, inputs, logs and records needed to repeat the chosen result.
```

For data-only or unsupported-rendering branches, explicitly state no image was generated and provide the user-run commands where applicable. For completed rendering, do not falsely retain a blanket "no visualization performed" statement. Call the output a publication-oriented candidate when appropriate, not guaranteed journal compliance or human approval.

## 15. Language maintenance

Keep this SKILL.md in English in future revisions, including frontmatter, headings, instructions, comments, and examples. Preserve literal filenames, paths, technical identifiers, and existing mathematical notation. This artifact-language requirement does not constrain conversation, captions, or report language.
