import { useId, useMemo, useState } from "react";
import { Canvas, ThreeEvent } from "@react-three/fiber";
import { GizmoHelper, GizmoViewport, OrbitControls } from "@react-three/drei";
import * as THREE from "three";
import type { SavedStructureData } from "../lib/types";
import type { AtomicMagneticMoment, MagneticMomentVector } from "../lib/magnetism";
import { getElementColor } from "../lib/elementColors";
import { InfoTooltip } from "./InfoTooltip";

type Vector3Tuple = [number, number, number];
const DEFAULT_ARROW_SCALE = 1;
const ARROW_BASE_SCALE = 2 / 3;
const WEAK_MOMENT_THRESHOLD_FRACTION = 0.05;

interface MagnetismViewerProps {
  structure: SavedStructureData;
  moments: AtomicMagneticMoment[];
  totalMagnetization?: number | null;
}

interface ViewerAtom {
  index: number;
  symbol: string;
  position: THREE.Vector3;
  moment: MagneticMomentVector;
  magnitude: number;
}

interface DisplayAtom {
  key: string;
  source: ViewerAtom;
  position: THREE.Vector3;
}

function toVector3(values: Vector3Tuple): THREE.Vector3 {
  return new THREE.Vector3(values[0], values[1], values[2]);
}

function magnitude(vector: THREE.Vector3): number {
  return Math.max(vector.length(), 1e-6);
}

function formatMoment(value: number): string {
  if (Math.abs(value) < 0.0005) return "0.000";
  return value.toFixed(3);
}

function shouldShowMomentArrow(
  atom: ViewerAtom,
  maxMoment: number,
  showWeakMoments: boolean,
): boolean {
  if (atom.magnitude <= 1e-8) return false;
  if (showWeakMoments) return true;
  return atom.magnitude / Math.max(maxMoment, 1e-6) >= WEAK_MOMENT_THRESHOLD_FRACTION;
}

function cartesianFromCrystal(
  position: Vector3Tuple,
  aVec: THREE.Vector3,
  bVec: THREE.Vector3,
  cVec: THREE.Vector3,
): THREE.Vector3 {
  return new THREE.Vector3()
    .addScaledVector(aVec, position[0])
    .addScaledVector(bVec, position[1])
    .addScaledVector(cVec, position[2]);
}

function fractionalFromCartesian(
  position: THREE.Vector3,
  aVec: THREE.Vector3,
  bVec: THREE.Vector3,
  cVec: THREE.Vector3,
): THREE.Vector3 {
  const bxc = bVec.clone().cross(cVec);
  const cxa = cVec.clone().cross(aVec);
  const axb = aVec.clone().cross(bVec);
  const cellVolume = aVec.dot(bxc);
  if (Math.abs(cellVolume) < 1e-10) return new THREE.Vector3(0, 0, 0);

  return new THREE.Vector3(
    position.dot(bxc) / cellVolume,
    position.dot(cxa) / cellVolume,
    position.dot(axb) / cellVolume,
  );
}

function wrapToUnitCell(value: number): number {
  const wrapped = value - Math.floor(value);
  return wrapped >= 1 - 1e-10 ? 0 : wrapped;
}

function expandAtomsOnCellFaces(
  atoms: ViewerAtom[],
  aVec: THREE.Vector3,
  bVec: THREE.Vector3,
  cVec: THREE.Vector3,
): DisplayAtom[] {
  const FACE_EPSILON = 1e-4;
  const expanded: DisplayAtom[] = [];

  atoms.forEach((atom) => {
    const fractional = fractionalFromCartesian(atom.position, aVec, bVec, cVec);
    const wrapped = new THREE.Vector3(
      wrapToUnitCell(fractional.x),
      wrapToUnitCell(fractional.y),
      wrapToUnitCell(fractional.z),
    );

    const onFaceA = wrapped.x <= FACE_EPSILON || wrapped.x >= 1 - FACE_EPSILON;
    const onFaceB = wrapped.y <= FACE_EPSILON || wrapped.y >= 1 - FACE_EPSILON;
    const onFaceC = wrapped.z <= FACE_EPSILON || wrapped.z >= 1 - FACE_EPSILON;
    const aShifts = onFaceA ? [0, 1] : [0];
    const bShifts = onFaceB ? [0, 1] : [0];
    const cShifts = onFaceC ? [0, 1] : [0];

    aShifts.forEach((aShift) => {
      bShifts.forEach((bShift) => {
        cShifts.forEach((cShift) => {
          const displayPosition = cartesianFromCrystal(
            [wrapped.x + aShift, wrapped.y + bShift, wrapped.z + cShift],
            aVec,
            bVec,
            cVec,
          );
          expanded.push({
            key: `${atom.index}:${aShift}${bShift}${cShift}`,
            source: atom,
            position: displayPosition,
          });
        });
      });
    });
  });

  return expanded;
}

function makeEdges(aVec: THREE.Vector3, bVec: THREE.Vector3, cVec: THREE.Vector3) {
  const o = new THREE.Vector3(0, 0, 0);
  const a = aVec.clone();
  const b = bVec.clone();
  const c = cVec.clone();
  const ab = a.clone().add(b);
  const ac = a.clone().add(c);
  const bc = b.clone().add(c);
  const abc = a.clone().add(b).add(c);

  return [
    [o, a],
    [a, ab],
    [ab, b],
    [b, o],
    [c, ac],
    [ac, abc],
    [abc, bc],
    [bc, c],
    [o, c],
    [a, ac],
    [b, bc],
    [ab, abc],
  ] as Array<[THREE.Vector3, THREE.Vector3]>;
}

function axisQuaternion(direction: THREE.Vector3): THREE.Quaternion {
  return new THREE.Quaternion().setFromUnitVectors(
    new THREE.Vector3(0, 1, 0),
    direction.clone().normalize(),
  );
}

function MomentArrow({
  atom,
  aDir,
  bDir,
  cDir,
  maxMoment,
  maxDimension,
  atomRadius,
  arrowScale,
}: {
  atom: ViewerAtom;
  aDir: THREE.Vector3;
  bDir: THREE.Vector3;
  cDir: THREE.Vector3;
  maxMoment: number;
  maxDimension: number;
  atomRadius: number;
  arrowScale: number;
}) {
  const direction = useMemo(() => {
    const vector = new THREE.Vector3()
      .addScaledVector(aDir, atom.moment[0])
      .addScaledVector(bDir, atom.moment[1])
      .addScaledVector(cDir, atom.moment[2]);
    return vector.lengthSq() > 1e-10 ? vector.normalize() : cDir.clone();
  }, [aDir, atom.moment, bDir, cDir]);

  const strength = atom.magnitude / Math.max(maxMoment, 1e-6);
  const scaleFactor = ARROW_BASE_SCALE * arrowScale;
  const totalLength = maxDimension * 0.44 * strength * scaleFactor;
  if (totalLength <= 1e-8) return null;

  const positiveLength = totalLength * 0.62;
  const tailLength = totalLength * 0.38;
  const coneLength = positiveLength * 0.36;
  const shaftLength = positiveLength - coneLength;
  const maxRadius = maxDimension * 0.022;
  const shaftRadius = maxRadius * strength * scaleFactor;
  const coneRadius = shaftRadius * 2.35;
  const sphereOverlap = Math.min(atomRadius * 0.28, shaftRadius * 1.25);
  const surfaceOffset = atomRadius - sphereOverlap;
  const orientation = axisQuaternion(direction);
  const materialColor = "#e24b35";

  const tailCenter = atom.position.clone().addScaledVector(direction, -(surfaceOffset + tailLength * 0.5));
  const shaftCenter = atom.position.clone().addScaledVector(direction, surfaceOffset + shaftLength * 0.5);
  const coneCenter = atom.position.clone().addScaledVector(direction, surfaceOffset + shaftLength + coneLength * 0.5);

  return (
    <group>
      {tailLength > 1e-8 && (
        <mesh position={tailCenter} quaternion={orientation}>
          <cylinderGeometry args={[shaftRadius, shaftRadius, tailLength, 20, 1]} />
          <meshStandardMaterial color={materialColor} roughness={0.42} metalness={0.08} />
        </mesh>
      )}
      {shaftLength > 1e-8 && (
        <mesh position={shaftCenter} quaternion={orientation}>
          <cylinderGeometry args={[shaftRadius, shaftRadius, shaftLength, 24, 1]} />
          <meshStandardMaterial color={materialColor} roughness={0.42} metalness={0.08} />
        </mesh>
      )}
      {coneLength > 1e-8 && (
        <mesh position={coneCenter} quaternion={orientation}>
          <coneGeometry args={[coneRadius, coneLength, 28, 1]} />
          <meshStandardMaterial color={materialColor} roughness={0.38} metalness={0.08} />
        </mesh>
      )}
    </group>
  );
}

function MagnetismScene({
  atoms,
  hiddenSpecies,
  selectedAtomIndex,
  onSelectAtom,
  aVec,
  bVec,
  cVec,
  center,
  maxDimension,
  arrowScale,
  showWeakMoments,
}: {
  atoms: ViewerAtom[];
  hiddenSpecies: Set<string>;
  selectedAtomIndex: number | null;
  onSelectAtom: (index: number | null) => void;
  aVec: THREE.Vector3;
  bVec: THREE.Vector3;
  cVec: THREE.Vector3;
  center: THREE.Vector3;
  maxDimension: number;
  arrowScale: number;
  showWeakMoments: boolean;
}) {
  const edges = useMemo(() => makeEdges(aVec, bVec, cVec), [aVec, bVec, cVec]);
  const visibleAtoms = atoms.filter((atom) => !hiddenSpecies.has(atom.symbol));
  const displayAtoms = useMemo(
    () => expandAtomsOnCellFaces(visibleAtoms, aVec, bVec, cVec),
    [visibleAtoms, aVec, bVec, cVec],
  );
  const maxMoment = Math.max(...visibleAtoms.map((atom) => atom.magnitude), 1e-6);
  const atomRadius = maxDimension * 0.045;
  const selectionRadius = atomRadius * 1.35;
  const aDir = aVec.clone().normalize();
  const bDir = bVec.clone().normalize();
  const cDir = cVec.clone().normalize();

  return (
    <>
      <OrbitControls
        makeDefault
        enablePan
        enableRotate
        enableZoom
        target={[0, 0, 0]}
        minDistance={maxDimension * 0.35}
        maxDistance={maxDimension * 8}
      />
      <GizmoHelper alignment="bottom-left" margin={[68, 68]}>
        <GizmoViewport
          axisColors={["#d94841", "#38a169", "#2563eb"]}
          labelColor="#f8fafc"
          labels={["a", "b", "c"]}
        />
      </GizmoHelper>
      <group position={[-center.x, -center.y, -center.z]}>
        {edges.map(([start, end], index) => (
          <line key={`cell-edge-${index}`}>
            <bufferGeometry>
              <bufferAttribute
                attach="attributes-position"
                args={[
                  new Float32Array([
                    start.x,
                    start.y,
                    start.z,
                    end.x,
                    end.y,
                    end.z,
                  ]),
                  3,
                ]}
              />
            </bufferGeometry>
            <lineBasicMaterial color="#9ca3af" transparent opacity={0.7} />
          </line>
        ))}

        {displayAtoms.map((displayAtom) => {
          const atom = displayAtom.source;
          const selected = selectedAtomIndex === atom.index;
          const showArrow = shouldShowMomentArrow(atom, maxMoment, showWeakMoments);
          return (
            <group key={`mag-atom-${displayAtom.key}`}>
              {showArrow && (
                <MomentArrow
                  atom={{ ...atom, position: displayAtom.position }}
                  aDir={aDir}
                  bDir={bDir}
                  cDir={cDir}
                  maxMoment={maxMoment}
                  maxDimension={maxDimension}
                  atomRadius={atomRadius}
                  arrowScale={arrowScale}
                />
              )}
              {selected && (
                <mesh position={displayAtom.position}>
                  <sphereGeometry args={[selectionRadius, 32, 32]} />
                  <meshBasicMaterial color="#f6c344" transparent opacity={0.25} depthWrite={false} />
                </mesh>
              )}
              <mesh
                position={displayAtom.position}
                onClick={(event: ThreeEvent<MouseEvent>) => {
                  event.stopPropagation();
                  onSelectAtom(atom.index);
                }}
              >
                <sphereGeometry args={[atomRadius, 40, 40]} />
                <meshStandardMaterial
                  color={getElementColor(atom.symbol)}
                  metalness={0.05}
                  roughness={0.48}
                  emissive={selected ? "#6b4a00" : "#000000"}
                  emissiveIntensity={selected ? 0.25 : 0}
                />
              </mesh>
            </group>
          );
        })}
      </group>
    </>
  );
}

export function MagnetismViewer({
  structure,
  moments,
  totalMagnetization,
}: MagnetismViewerProps) {
  const [hiddenSpecies, setHiddenSpecies] = useState<Set<string>>(() => new Set());
  const [selectedAtomIndex, setSelectedAtomIndex] = useState<number | null>(null);
  const [arrowScale, setArrowScale] = useState(DEFAULT_ARROW_SCALE);
  const [showWeakMoments, setShowWeakMoments] = useState(true);
  const arrowScaleInputId = useId();
  const showWeakMomentsInputId = useId();

  const geometry = useMemo(() => {
    const aVec = toVector3(structure.cell_parameters![0]);
    const bVec = toVector3(structure.cell_parameters![1]);
    const cVec = toVector3(structure.cell_parameters![2]);
    const center = aVec.clone().add(bVec).add(cVec).multiplyScalar(0.5);
    const maxDimension = Math.max(magnitude(aVec), magnitude(bVec), magnitude(cVec));
    const momentsByAtom = new Map(moments.map((moment) => [moment.atomIndex, moment]));
    const atoms: ViewerAtom[] = structure.atoms.map((atom, index) => {
      const moment = momentsByAtom.get(index);
      const position = structure.position_units === "crystal"
        ? cartesianFromCrystal(atom.position, aVec, bVec, cVec)
        : toVector3(atom.position);
      return {
        index,
        symbol: atom.symbol,
        position,
        moment: moment?.vector ?? [0, 0, 0],
        magnitude: moment?.magnitude ?? 0,
      };
    });

    return { atoms, aVec, bVec, cVec, center, maxDimension };
  }, [moments, structure]);

  const species = useMemo(
    () => Array.from(new Set(geometry.atoms.map((atom) => atom.symbol))).sort(),
    [geometry.atoms],
  );
  const selectedAtom = selectedAtomIndex == null
    ? null
    : geometry.atoms.find((atom) => atom.index === selectedAtomIndex) ?? null;
  const visibleAtoms = geometry.atoms.filter((atom) => !hiddenSpecies.has(atom.symbol));
  const maxVisibleMoment = Math.max(...visibleAtoms.map((atom) => atom.magnitude), 1e-6);
  const visibleMomentCount = visibleAtoms.filter((atom) =>
    shouldShowMomentArrow(atom, maxVisibleMoment, showWeakMoments)
  ).length;
  const cameraDistance = geometry.maxDimension * 2.6;

  function toggleSpecies(symbol: string) {
    setHiddenSpecies((current) => {
      const next = new Set(current);
      if (next.has(symbol)) {
        next.delete(symbol);
      } else {
        next.add(symbol);
        if (selectedAtom?.symbol === symbol) {
          setSelectedAtomIndex(null);
        }
      }
      return next;
    });
  }

  return (
    <div className="magnetism-viewer">
      <div className="magnetism-viewer-toolbar">
        <div className="magnetism-viewer-title">
          <span>Magnetism</span>
          <small>{visibleMomentCount} visible moments</small>
        </div>
        <div className="magnetism-species-toggles" aria-label="Visible atom species">
          {species.map((symbol) => (
            <label key={symbol} className="magnetism-species-toggle">
              <input
                type="checkbox"
                checked={!hiddenSpecies.has(symbol)}
                onChange={() => toggleSpecies(symbol)}
              />
              <span
                className="magnetism-species-swatch"
                style={{ backgroundColor: getElementColor(symbol) }}
                aria-hidden="true"
              />
              {symbol}
            </label>
          ))}
        </div>
      </div>
      <div className="magnetism-viewer-body">
        <div className="magnetism-canvas">
          <Canvas
            camera={{
              position: [0, 0, cameraDistance],
              up: [0, 1, 0],
              fov: 34,
              near: 0.01,
              far: cameraDistance * 12,
            }}
            onPointerMissed={() => setSelectedAtomIndex(null)}
          >
            <color attach="background" args={["#111827"]} />
            <ambientLight intensity={0.62} />
            <directionalLight position={[4, 6, 8]} intensity={1.1} />
            <directionalLight position={[-6, -3, -4]} intensity={0.35} />
            <MagnetismScene
              atoms={geometry.atoms}
              hiddenSpecies={hiddenSpecies}
              selectedAtomIndex={selectedAtomIndex}
              onSelectAtom={setSelectedAtomIndex}
              aVec={geometry.aVec}
              bVec={geometry.bVec}
              cVec={geometry.cVec}
              center={geometry.center}
              maxDimension={geometry.maxDimension}
              arrowScale={arrowScale}
              showWeakMoments={showWeakMoments}
            />
          </Canvas>
        </div>
        <aside className="magnetism-selection-panel">
          <div className="band-control-section magnetism-control-section">
            <div className="band-control-header magnetism-control-static-header">
              View Controls
            </div>
            <div className="band-control-grid">
              <div className="band-control-row">
                <label htmlFor={arrowScaleInputId}>
                  Arrow Scale
                  <InfoTooltip text="Multiplies the displayed arrow length and thickness. This does not change the saved magnetic moment values." />
                </label>
                <input
                  id={arrowScaleInputId}
                  type="range"
                  min={0.25}
                  max={3}
                  step={0.05}
                  value={arrowScale}
                  onChange={(event) => setArrowScale(Number.parseFloat(event.target.value))}
                />
                <span className="band-control-value">{arrowScale.toFixed(2)}x</span>
              </div>

              <div className="band-control-row">
                <label htmlFor={showWeakMomentsInputId}>
                  Show weak moments
                  <InfoTooltip text="When off, arrows below 5% of the strongest currently visible moment are hidden. Atom spheres remain selectable." />
                </label>
                <input
                  className="magnetism-weak-moments-toggle"
                  id={showWeakMomentsInputId}
                  type="checkbox"
                  checked={showWeakMoments}
                  onChange={(event) => setShowWeakMoments(event.target.checked)}
                />
                <span className="band-control-value">{showWeakMoments ? "On" : "Off"}</span>
              </div>
            </div>
          </div>

          {selectedAtom ? (
            <>
              <div className="magnetism-selection-heading">
                <span>{selectedAtom.symbol}{selectedAtom.index + 1}</span>
                <small>|m| = {formatMoment(selectedAtom.magnitude)} μB</small>
              </div>
              <dl className="magnetism-moment-grid">
                <div>
                  <dt>a</dt>
                  <dd>{formatMoment(selectedAtom.moment[0])}</dd>
                </div>
                <div>
                  <dt>b</dt>
                  <dd>{formatMoment(selectedAtom.moment[1])}</dd>
                </div>
                <div>
                  <dt>c</dt>
                  <dd>{formatMoment(selectedAtom.moment[2])}</dd>
                </div>
              </dl>
            </>
          ) : (
            <div className="magnetism-empty-selection">
              Select an atom to inspect its magnetic moment.
            </div>
          )}
          {Number.isFinite(Number(totalMagnetization)) && (
            <div className="magnetism-total">
              <span>Total</span>
              <strong>{formatMoment(Number(totalMagnetization))} μB/cell</strong>
            </div>
          )}
        </aside>
      </div>
    </div>
  );
}
