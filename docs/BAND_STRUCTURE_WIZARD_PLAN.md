# Band Structure Wizard Implementation Plan

## Overview

Add band structure calculation and visualization to QCortado. Users can compute electronic band structures along high-symmetry k-paths and visualize them interactively.

## Prerequisites & Constraints

- **Requires completed SCF calculation** - Band structure is only accessible from the project dashboard after an SCF has been run for the selected structure
- **Two-step QE process**: SCF (already done) → NSCF along k-path → bands.x post-processing
- **Interactive visualization** required for v1
- **No 3D Brillouin zone** in v1 - use dropdown of standard paths + manual entry

---

## Architecture

### Data Flow

```
Dashboard (select SCF)
    → Band Structure Wizard
        → Step 1: Select source SCF
        → Step 2: Choose k-path (standard or manual)
        → Step 3: Configure parameters
        → Step 4: Run NSCF + bands.x
        → Step 5: Interactive band plot
    → Save to project (as new calculation type: "bands")
```

### New Files to Create

**Frontend:**
- `src/components/BandStructureWizard.tsx` - Main wizard component
- `src/components/BandPlot.tsx` - Interactive SVG band structure chart
- `src/lib/kpaths.ts` - K-path data for each Bravais lattice
- `src/lib/brillouinZone.ts` - Crystal system detection utilities

**Backend:**
- `src-tauri/src/qe/bands.rs` - Band structure input generation and output parsing

---

## Phase 1: K-Path Infrastructure

### 1.1 Crystal System Detection (`src/lib/brillouinZone.ts`)

Map space group number → Bravais lattice type:

```typescript
type BravaisLattice =
  | 'cubic-P'      // Primitive cubic (sc)
  | 'cubic-F'      // Face-centered cubic (fcc)
  | 'cubic-I'      // Body-centered cubic (bcc)
  | 'tetragonal-P' // Primitive tetragonal
  | 'tetragonal-I' // Body-centered tetragonal
  | 'orthorhombic-P' | 'orthorhombic-C' | 'orthorhombic-I' | 'orthorhombic-F'
  | 'hexagonal'    // Hexagonal
  | 'trigonal-R'   // Rhombohedral
  | 'monoclinic-P' | 'monoclinic-C'
  | 'triclinic';

function detectBravaisLattice(spaceGroupNumber: number): BravaisLattice;
```

Space group ranges:
- 1-2: Triclinic
- 3-15: Monoclinic
- 16-74: Orthorhombic
- 75-142: Tetragonal
- 143-167: Trigonal
- 168-194: Hexagonal
- 195-230: Cubic

Within cubic (195-230):
- P lattice: 195-199, 200-206 (some)
- F lattice: 196, 202, 203, 209, 210, 216, 219, 225-228
- I lattice: 197, 199, 204, 206, 211, 214, 217, 220, 229, 230

### 1.2 K-Path Database (`src/lib/kpaths.ts`)

```typescript
interface HighSymmetryPoint {
  label: string;           // e.g., "Γ", "X", "M"
  coords: [number, number, number];  // Crystal coordinates
  description?: string;    // e.g., "Zone center"
}

interface StandardKPath {
  lattice: BravaisLattice;
  points: HighSymmetryPoint[];
  defaultPath: string[];   // e.g., ["Γ", "X", "M", "Γ", "R"]
  description: string;
}

const KPATH_DATABASE: Record<BravaisLattice, StandardKPath>;
```

**Standard paths to implement:**

| Lattice | Points | Default Path |
|---------|--------|--------------|
| cubic-P | Γ(0,0,0), X(0.5,0,0), M(0.5,0.5,0), R(0.5,0.5,0.5) | Γ→X→M→Γ→R→X\|M→R |
| cubic-F | Γ, X, W, K, L, U | Γ→X→W→K→Γ→L→U→W\|L→K\|U→X |
| cubic-I | Γ, H, N, P | Γ→H→N→Γ→P→H\|P→N |
| hexagonal | Γ, M, K, A, L, H | Γ→M→K→Γ→A→L→H→A\|L→M\|K→H |
| tetragonal-P | Γ, X, M, Z, R, A | Γ→X→M→Γ→Z→R→A→Z\|X→R\|M→A |

Reference: Setyawan & Curtarolo, Comp. Mat. Sci. 49, 299 (2010)

### 1.3 K-Path Selection UI

```typescript
interface KPathConfig {
  mode: 'standard' | 'manual';
  // Standard mode
  selectedPath?: string[];  // List of point labels in order
  // Manual mode
  manualPoints?: Array<{
    label: string;
    coords: [number, number, number];
    npoints: number;  // Points to next segment
  }>;
  pointsPerSegment: number;  // Default: 20
}
```

UI Components:
1. **Lattice info display** - Show detected Bravais lattice
2. **Path dropdown** - "Recommended path" vs "Custom"
3. **Path editor** (for custom):
   - List of available high-symmetry points with checkboxes
   - Drag-and-drop reordering
   - Preview of path as text: "Γ → X → M → Γ"
4. **Manual entry** - Text input for direct coordinate entry

---

## Phase 2: Backend - Band Calculation

### 2.1 Input Generation (`src-tauri/src/qe/bands.rs`)

Generate NSCF input for bands:

```rust
pub struct BandsCalculation {
    // Inherit from SCF
    pub source_scf_dir: PathBuf,  // Where the SCF .save directory is
    pub prefix: String,
    pub pseudo_dir: String,

    // Band-specific
    pub k_path: Vec<KPathPoint>,
    pub nbnd: Option<u32>,  // Number of bands (auto if None)
}

pub struct KPathPoint {
    pub label: String,
    pub coords: [f64; 3],
    pub npoints: u32,  // 0 for last point
}

pub fn generate_bands_input(calc: &BandsCalculation, system: &SystemConfig) -> String;
```

Output format:
```
&CONTROL
  calculation = 'bands',
  prefix = 'qcortado_scf',
  outdir = './tmp',
  pseudo_dir = '...',
  verbosity = 'high',
/

&SYSTEM
  ... (same as SCF)
  nbnd = 20,
/

&ELECTRONS
/

ATOMIC_SPECIES
  ...

ATOMIC_POSITIONS {crystal}
  ...

K_POINTS {crystal_b}
5
0.000 0.000 0.000 30  Gamma
0.500 0.000 0.000 30  X
0.500 0.500 0.000 30  M
0.000 0.000 0.000 30  Gamma
0.000 0.000 0.500 0   Z

CELL_PARAMETERS {angstrom}
  ...
```

### 2.2 Bands.x Input Generation

After NSCF, run bands.x:

```
&BANDS
  prefix = 'qcortado_scf',
  outdir = './tmp',
  filband = 'bands.dat',
  lsym = .true.,
/
```

### 2.3 Output Parsing

Parse `bands.dat.gnu` (gnuplot format) output:

```rust
pub struct BandData {
    pub k_points: Vec<f64>,           // k-path distance
    pub energies: Vec<Vec<f64>>,      // [band_index][k_index]
    pub fermi_energy: f64,
    pub high_symmetry_points: Vec<(f64, String)>,  // (k_distance, label)
    pub n_bands: usize,
    pub n_kpoints: usize,
    pub band_gap: Option<BandGap>,
}

pub struct BandGap {
    pub value: f64,        // eV
    pub is_direct: bool,
    pub vbm_k: f64,        // k-point of valence band max
    pub cbm_k: f64,        // k-point of conduction band min
}
```

### 2.4 New Tauri Commands

```rust
#[tauri::command]
async fn run_bands_calculation(
    app: AppHandle,
    project_id: String,
    cif_id: String,
    source_calc_id: String,  // SCF calculation to use
    k_path: Vec<KPathPoint>,
    nbnd: Option<u32>,
    state: State<'_, AppState>,
) -> Result<BandData, String>;

#[tauri::command]
fn get_standard_kpath(
    lattice_type: String,
) -> Result<StandardKPath, String>;
```

---

## Phase 3: Band Structure Wizard UI

### 3.1 Wizard Flow

```
┌─────────────────────────────────────────────────────┐
│  Band Structure Wizard                    [← Back]  │
├─────────────────────────────────────────────────────┤
│  Step: ① Source  ② K-Path  ③ Parameters  ④ Run     │
├─────────────────────────────────────────────────────┤
│                                                     │
│  [Step-specific content]                            │
│                                                     │
└─────────────────────────────────────────────────────┘
```

### Step 1: Select Source SCF

- List available SCF calculations for current structure
- Show: date, total energy, convergence status
- Radio button selection
- "No SCF available" state with button to run SCF first

### Step 2: K-Path Selection

```
┌─────────────────────────────────────────────────────┐
│ Crystal System: Cubic (Face-centered)               │
│ Space Group: Fm-3m (#225)                           │
├─────────────────────────────────────────────────────┤
│ K-Path Mode:                                        │
│ ○ Recommended path for FCC                          │
│   Γ → X → W → K → Γ → L → U → W                    │
│                                                     │
│ ○ Custom path                                       │
│   [Γ] [X] [W] [K] [L] [U]  ← toggle points         │
│   Current: Γ → X → L → Γ   ← drag to reorder      │
│                                                     │
│ ○ Manual entry                                      │
│   ┌──────────────────────────────────────────┐     │
│   │ G 0.0 0.0 0.0 30                         │     │
│   │ X 0.5 0.0 0.0 30                         │     │
│   │ ...                                       │     │
│   └──────────────────────────────────────────┘     │
├─────────────────────────────────────────────────────┤
│ Points per segment: [20 ▼]                          │
└─────────────────────────────────────────────────────┘
```

### Step 3: Parameters

- **Number of bands**: Auto-calculated from electrons, with override
- **Points per segment**: 20 (default), affects resolution
- Show estimated k-points total

### Step 4: Run

- Show generated input
- Stream NSCF output (like SCF wizard)
- Then run bands.x
- Progress indicator

### Step 5: Results - Interactive Band Plot

See Phase 4 below.

---

## Phase 4: Interactive Band Plot Component

### 4.1 Component Props

```typescript
interface BandPlotProps {
  data: BandData;
  width?: number;
  height?: number;
  energyRange?: [number, number];  // Auto if not specified
  showFermiLevel?: boolean;
  colorScheme?: 'mono' | 'gradient';
}
```

### 4.2 Features

**Must have (v1):**
- SVG-based plot with proper scaling
- X-axis: k-path with labeled high-symmetry points (vertical lines)
- Y-axis: Energy (eV) relative to Fermi level
- Horizontal dashed line at E_F = 0
- Hover tooltip: band index, energy value, k-point
- Zoom: scroll wheel or pinch
- Pan: click and drag
- Band gap annotation (if semiconductor/insulator)

**Nice to have (v1.1+):**
- Click band to highlight it
- Color bands by orbital character (requires projwfc.x)
- Export as PNG/SVG
- Export raw data as CSV
- Comparison mode (overlay multiple calculations)

### 4.3 Implementation Approach

Use SVG with React for interactivity:

```typescript
function BandPlot({ data }: BandPlotProps) {
  const [zoom, setZoom] = useState(1);
  const [pan, setPan] = useState({ x: 0, y: 0 });
  const [hoveredPoint, setHoveredPoint] = useState<{band: number, k: number} | null>(null);

  // Transform data to SVG coordinates
  const xScale = (k: number) => ...;
  const yScale = (e: number) => ...;

  return (
    <svg viewBox="...">
      {/* Grid lines */}
      {/* High-symmetry point markers */}
      {/* Fermi level */}
      {/* Band lines */}
      {data.energies.map((band, bandIdx) => (
        <path
          key={bandIdx}
          d={bandToPath(band, data.k_points)}
          onMouseMove={(e) => handleHover(e, bandIdx)}
        />
      ))}
      {/* Tooltip */}
      {hoveredPoint && <Tooltip ... />}
    </svg>
  );
}
```

---

## Phase 5: Integration

### 5.1 Dashboard Integration

Add "Bands" button to calculation action grid:
- Enabled only if SCF calculations exist for current structure
- Opens BandStructureWizard as full-screen view (like SCF)

### 5.2 Saving Results

Save band calculations to project with type "bands":
- Store band data JSON
- Store input/output files
- Reference source SCF calculation

### 5.3 Calculation History

Band calculations appear in history with:
- Type: "BANDS"
- Show: band gap (if applicable), number of bands
- Expanded view: mini band plot preview + full details

---

## Implementation Order

1. **K-path infrastructure** (kpaths.ts, brillouinZone.ts)
2. **Backend band input generation** (bands.rs)
3. **Backend output parsing**
4. **Wizard UI** (steps 1-4)
5. **Basic band plot** (static SVG)
6. **Interactive features** (zoom, pan, hover)
7. **Dashboard integration**
8. **Save to project**

---

## Testing Plan

1. **K-path detection**: Test with various space groups
2. **Input generation**: Compare output with manual QE input
3. **Parsing**: Use known bands.dat.gnu files
4. **End-to-end**: Run on Si (diamond), Cu (fcc metal), NaCl (rocksalt insulator)

---

## References

- Setyawan & Curtarolo, "High-throughput electronic band structure calculations", Comp. Mat. Sci. 49, 299 (2010) - Standard k-paths
- QE documentation: https://www.quantum-espresso.org/Doc/INPUT_PW.html
- QE bands.x: https://www.quantum-espresso.org/Doc/INPUT_BANDS.html
