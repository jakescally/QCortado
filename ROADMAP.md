# QCortado Roadmap

## Direction

QCortado should deepen the workflow it already owns: structure import, SCF or relaxation, derived electronic or vibrational calculations, transport or property analysis, and portable exports. The near-term roadmap should favor features that reuse saved project artifacts, fit the dashboard-and-wizard model, and produce clear plots or tables inside the app.

This intentionally keeps QCortado focused on being a workflow tool, not a generic launcher for every Quantum ESPRESSO binary.

## Near-Term Priorities

### 1. Native transport analysis from Wannier results

**Summary**  
Add a real transport workflow on top of saved Wannier calculations. The first target should be bulk or semi-classical transport outputs such as conductivity, Seebeck coefficient, and electronic thermal conductivity across temperature and chemical-potential scans.

**Why this should be next**  
QCortado already has a solid Wannier foundation and already exports to Ludwig, so transport is the most natural way to turn that data into a higher-value in-app result. This would move QCortado from "prepare and export" toward "analyze and compare" without changing the core product shape.

**Basic implementation guide**  
Start by adding a transport action to saved Wannier calculations in the dashboard, then build a small wizard for the temperature range, chemical-potential scan, and output density. The backend can run a first transport engine against saved Wannier artifacts, parse tabulated results, and save them back as a new calculation type with reusable plots. The first version should stay tightly scoped: pick a saved Wannier run, configure a scan, run, view, and save.

### 2. Richer Fermi-surface analysis

**Summary**  
Expand the current Fermi-surface workflow from "generate FRMSF and open FermiSurfer" into a more informative analysis surface. The next step should be orbital projection support, better metadata, and stronger in-app summaries of what the generated surfaces mean.

**Why this should be next**  
QCortado already runs `fermi_velocity.x`, stores generated surface files, and has a dedicated wizard. That makes this a relatively low-friction upgrade with immediate value for transport-oriented users.

**Basic implementation guide**  
Extend the existing Fermi-surface pipeline to support `fermi_proj.x` and save projection-aware surface artifacts alongside the current outputs. In the frontend, add a viewer state that can switch between plain surface, velocity-oriented output, and orbital-projected output, with basic file and pocket summaries in the results panel. Keep FermiSurfer launch support as a secondary option rather than the main analysis path.

### 3. EPW-ready phonon pipeline, then full EPW workflow

**Summary**  
QCortado already exposes some EPW-prep controls in the phonon wizard, but it still needs a coherent bridge from phonon calculations to a true EPW workflow. The next step is to make artifact retention, saved-run labeling, and downstream prerequisites explicit and reliable.

**Why this should be next**  
This closes one of the biggest gaps between "phonons for plots" and "phonons for transport or superconductivity." It also builds directly on recent work around phonon recovery, compact saves, and EPW-related artifact handling.

**Basic implementation guide**  
First, formalize save policies such as `plot-only`, `EPW-ready`, and `full archive`, and show that status clearly in saved phonon calculations. Then tighten the backend so the exact required `dvscf`, `dyn*`, and `phsave` subsets are preserved consistently when EPW prep is enabled. Once those artifacts are stable, add an `epw.x` wizard that attaches to saved SCF and phonon prerequisites instead of rerunning setup loosely from scratch.

### 4. General post-processing viewers (`pp.x` and related outputs)

**Summary**  
Add a broader post-processing layer for projected DOS, charge density, spin density, electrostatic potential, and similar derived data. These are common endpoints in real QE workflows and are still mostly outside QCortado today.

**Why this should be next**  
This is one of the clearest "one-stop-shop" gaps in the current app. Many users finish the main calculation path but still need a separate tool for interpretation, which breaks the project-centered workflow QCortado is trying to own.

**Basic implementation guide**  
Start with projected DOS because it fits naturally beside the current DOS and band workflows and already aligns with existing `projwfc.x` usage. After that, add a preset-driven `pp.x` wizard with a small number of curated outputs rather than exposing the entire raw parameter space. Each preset should land in a dedicated viewer or export mode so the results feel like native QCortado features instead of loosely attached files.

### 5. `hp.x` for Hubbard parameter estimation

**Summary**  
Add a guided `hp.x` workflow so DFT+U values can be derived inside QCortado instead of entered entirely by hand. This would make the app more credible for correlated-material workflows.

**Why this should be next**  
QCortado already exposes DFT+U controls in the SCF flow, so `hp.x` is a clean way to close that loop. It is a focused addition that improves scientific defensibility without pulling the app into a completely different problem space.

**Basic implementation guide**  
Branch from a converged SCF calculation into a compact `hp.x` wizard that exposes only the core setup needed for a first version. Save computed Hubbard parameters as project artifacts with enough metadata to show provenance and source SCF. Then add a one-click path to apply the generated values into a new SCF or relaxation preset.

## Transport Direction

- **Near-term in-app path:** start with a Wannier-based transport workflow so QCortado can compute and visualize transport quantities from artifacts it already knows how to generate.
- **Advanced path:** build EPW after the phonon artifact pipeline is stable, so the app can support more serious electron-phonon transport and related workflows.
- **Existing external bridge:** keep Ludwig export as a useful external option, especially for 2D work, but do not treat it as QCortado's main transport story.
- **Not a near-term target:** ballistic transport tooling such as `pwcond.x` is valuable, but it does not match QCortado's current bulk-material and project-dashboard workflow as well as the items above.

## Explicitly Deferred From This Near-Term Roadmap

- NEB is out of scope for the current near-term plan.
- XSPECTRA and turboTDDFT are also out of scope for the current near-term plan.

