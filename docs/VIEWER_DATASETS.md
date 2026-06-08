# Viewer Datasets

## Purpose

Shared viewers should consume normalized datasets instead of engine-native input or output files. Engine modules should adapt their native results into these datasets.

This lets QE keep QE-specific parsing and provenance while preparing the UI for future engines such as remote-only Wien2k.

## Dataset Rule

Normalize viewer data, not engine inputs.

- QE inputs remain QE-specific.
- Future Wien2k inputs remain Wien2k-specific.
- Shared viewers receive stable, engine-neutral datasets.
- Engine-native result payloads can still be stored for provenance and specialized views.

## Current Viewer Payloads

### Bands

Current frontend type:

- `BandData` in `src/components/BandPlot.tsx`.

Current fields are already close to normalized:

- path coordinate distances.
- band energies.
- Fermi energy.
- high-symmetry markers.
- band gap.
- optional projection groups.

Migration target:

```ts
interface BandPathDataset {
  schema: "cortado.band_path.v1";
  quantity: "electronic_energy";
  x: number[];
  series: Array<{
    id: string;
    label: string;
    values: Array<number | null>;
  }>;
  xUnit: "path_distance";
  yUnit: "eV";
  referenceEnergyEv?: number | null;
  markers: Array<{
    x: number;
    label: string;
  }>;
  metadata?: Record<string, unknown>;
}
```

QE adapter sources:

- `src-tauri/src/engines/qe/bands.rs`
- `src/lib/engines/qe/bandDatasetAdapter.ts`
- current `BandData` result.
- optional `projwfc.x` projection groups.

Future Wien2k adapter sources:

- `spaghetti` output.
- Fermi energy from `case.scf`.
- k-path labels from the Wien2k workflow.

### Electronic DOS

Current frontend type:

- `ElectronicDOSData` in `src/components/ElectronicDOSPlot.tsx`.

Migration target:

```ts
interface DosDataset {
  schema: "cortado.dos.v1";
  x: number[];
  y: number[];
  xUnit: "eV";
  yUnit: string;
  referenceEnergyEv?: number | null;
  channels?: Array<{
    id: string;
    label: string;
    values: Array<number | null>;
  }>;
  metadata?: Record<string, unknown>;
}
```

QE adapter sources:

- `dos.x` output parsed by QE code.

Future Wien2k adapter sources:

- Wien2k DOS output from a future engine workflow.

### Phonon Dispersion

Current state:

- `src/components/PhononPlot.tsx` uses phonon-specific dispersion fields.
- `src/App.tsx` and `src/ViewerApp.tsx` adapt phonon dispersion into `BandData` for reuse in `BandPlot`.

Migration target:

```ts
interface LinePathDataset {
  schema: "cortado.line_path.v1";
  quantity: "phonon_frequency" | "electronic_energy";
  x: number[];
  series: Array<{
    id: string;
    label: string;
    values: Array<number | null>;
  }>;
  xUnit: "path_distance";
  yUnit: "cm-1" | "THz" | "eV";
  markers: Array<{
    x: number;
    label: string;
  }>;
  metadata?: Record<string, unknown>;
}
```

QE adapter sources:

- `ph.x`, `q2r.x`, and `matdyn.x` pipeline outputs.
- Current `PhononDispersion` payload.

Future Wien2k note:

- Wien2k phonon support is not part of the initial migration and should not drive current abstractions.

### Phonon DOS

Migration target:

```ts
interface PhononDosDataset {
  schema: "cortado.phonon_dos.v1";
  frequencies: number[];
  dos: number[];
  frequencyUnit: "cm-1" | "THz";
  metadata?: Record<string, unknown>;
}
```

QE adapter source:

- `matdyn.x` DOS output.

### Transport

Current frontend type:

- `TransportResult` in `src/lib/transport.ts`.

Current payload is BoltzWann-oriented but already includes an `engine` string.

Migration target:

```ts
interface TransportTensorDataset {
  schema: "cortado.transport_tensor.v1";
  engineSource: string;
  chemicalPotentialEv: number[];
  temperatureK: number[];
  tensors: Array<{
    id: string;
    label: string;
    unit: string;
    componentLabels: string[];
    values: Array<Array<Array<number | null>>>;
  }>;
  referenceEnergyEv?: number | null;
  metadata?: Record<string, unknown>;
}
```

QE adapter sources:

- Wannier90 and postw90 BoltzWann transport outputs.

Future Wien2k adapter sources:

- Future Wien2k transport workflow, if implemented.
- Do not assume it uses Wannier90 or BoltzWann.

### EPW

Current frontend type:

- `EpwViewerData` in `src/lib/engines/qe/epw.ts`.

EPW is QE-specific. Keep the EPW payload engine-native, but expose generic table and series views where useful.

Migration targets:

```ts
interface NumericTableDataset {
  schema: "cortado.numeric_table.v1";
  title: string;
  columns: string[];
  rows: Array<Array<number | string | null>>;
  metadata?: Record<string, unknown>;
}

interface SeriesDataset {
  schema: "cortado.series.v1";
  x: number[];
  y: Array<number | null>;
  xLabel: string;
  yLabel: string;
  metadata?: Record<string, unknown>;
}
```

QE adapter sources:

- EPW parsed tables.
- EPW mobility, superconductivity, and spectral summaries.

## Storage Strategy

During migration, saved calculations should support both native and normalized data.

Suggested intermediate shape:

```ts
interface SavedCalculationResultEnvelope {
  engine_id: "qe";
  native_result: unknown;
  viewer_datasets?: {
    bands?: BandPathDataset;
    dos?: DosDataset;
    phonons?: LinePathDataset;
    phononDos?: PhononDosDataset;
    transport?: TransportTensorDataset;
    tables?: NumericTableDataset[];
  };
}
```

Legacy saved calculations currently store `QEResult`. Existing records should continue to load. If a record has no `engine_id`, treat it as QE.

## Adapter Placement

Frontend adapters:

```text
src/lib/viewers/
  bands/
    types.ts
  dos/
    types.ts

src/lib/engines/qe/
  bandDatasetAdapter.ts
```

Future target, after more viewer datasets are normalized:

```text
src/lib/viewers/
  bands/
    types.ts
  dos/
    types.ts
  phonons/
    types.ts
  transport/
    types.ts
  tables/
    types.ts

src/lib/engines/qe/
  bandDatasetAdapter.ts
  dosDatasetAdapter.ts
  phononDatasetAdapter.ts
  transportDatasetAdapter.ts
```

Backend adapters:

```text
src-tauri/src/results/
  mod.rs
  bands.rs
  dos.rs
  phonons.rs
  transport.rs
  tables.rs

src-tauri/src/engines/qe/adapters.rs
```

Early PRs can keep adapters in existing files and re-export from the target path. The important boundary is that viewers consume normalized datasets through adapter functions.

## Migration Phases

### Phase 1: Type-only dataset definitions

- Add normalized TypeScript interfaces.
- Add no runtime conversion yet.
- Keep existing viewer props unchanged.

Validation:

```bash
npm run build
npm run test:unit
```

### Phase 2: QE adapters for existing viewer data

- Add functions that convert existing QE result payloads to normalized datasets.
- Use them in one low-risk viewer path first.
- Keep old payloads available.

Validation:

```bash
npm run build
npm run test:unit
```

### Phase 3: Shared viewer props

- Update `BandPlot` and DOS/phonon viewers to accept normalized datasets.
- Keep compatibility wrappers for existing call sites.
- Remove duplicated phonon-to-band adaptation from app shells.

Validation:

```bash
npm run build
npm run test:unit
```

### Phase 4: Backend result envelope

- Add optional normalized dataset fields to saved results.
- Preserve `QEResult` reads for legacy projects.
- Do not rewrite existing project files in place unless explicitly requested.

Validation:

```bash
npm run build
npm run test:unit
cargo check --manifest-path src-tauri/Cargo.toml
```

## Non-Goals

- Do not define a universal SCF input.
- Do not make EPW generic.
- Do not drop QE-native result payloads.
- Do not remove support for existing saved QE projects.
- Do not infer Wien2k workflows before the QE boundary is stable.
