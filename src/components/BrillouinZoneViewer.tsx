/**
 * Interactive 3D Brillouin Zone Viewer with k-path selection.
 *
 * Displays the first Brillouin zone with high-symmetry points and allows
 * users to click points to build a k-path for band structure calculations.
 */

import { useRef, useState, useMemo, useCallback } from "react";
import { Canvas, useFrame } from "@react-three/fiber";
import { OrbitControls, Text, Line, OrthographicCamera, PerspectiveCamera, Cone, Billboard } from "@react-three/drei";
import * as THREE from "three";
import { CrystalData } from "../lib/types";
import {
  Vec3,
  realSpaceLatticeVectors,
  reciprocalLatticeVectors,
  fractionalToCartesian,
  calculateBrillouinZone,
  BrillouinZoneGeometry,
  conventionalToPrimitive,
  CenteringType,
} from "../lib/reciprocalLattice";
import {
  BrillouinZoneData,
  HighSymmetryPoint,
  getBrillouinZoneData,
  BravaisLatticeType,
} from "../lib/brillouinZoneData";
import { detectBravaisLattice, BravaisLattice } from "../lib/brillouinZone";

// ============================================================================
// Types
// ============================================================================

export interface KPathPoint {
  label: string;
  coords: Vec3;
  npoints: number;
}

interface BrillouinZoneViewerProps {
  crystalData: CrystalData;
  onPathChange: (path: KPathPoint[]) => void;
  initialPath?: KPathPoint[];
  pointsPerSegment?: number;
}

// ============================================================================
// Sub-components for Three.js scene
// ============================================================================

interface HighSymmetryPointMeshProps {
  point: HighSymmetryPoint;
  position: [number, number, number];
  isSelected: boolean;
  isInPath: boolean;
  pathIndex: number | null;
  onClick: () => void;
  onHover: (hovered: boolean) => void;
}

function HighSymmetryPointMesh({
  point,
  position,
  isSelected,
  isInPath,
  pathIndex,
  onClick,
  onHover,
}: HighSymmetryPointMeshProps) {
  const meshRef = useRef<THREE.Mesh>(null);
  const [hovered, setHovered] = useState(false);

  // Animate on hover
  useFrame(() => {
    if (meshRef.current) {
      const targetScale = hovered ? 1.3 : 1;
      meshRef.current.scale.lerp(
        new THREE.Vector3(targetScale, targetScale, targetScale),
        0.1
      );
    }
  });

  const color = isSelected
    ? "#ff4444"
    : isInPath
    ? "#44ff44"
    : hovered
    ? "#ffaa00"
    : "#4488ff";

  return (
    <group position={position}>
      {/* Point sphere */}
      <mesh
        ref={meshRef}
        onClick={(e) => {
          e.stopPropagation();
          onClick();
        }}
        onPointerOver={(e) => {
          e.stopPropagation();
          setHovered(true);
          onHover(true);
        }}
        onPointerOut={() => {
          setHovered(false);
          onHover(false);
        }}
      >
        <sphereGeometry args={[0.04, 16, 16]} />
        <meshStandardMaterial color={color} />
      </mesh>

      {/* Label - always faces camera */}
      <Billboard position={[0, 0.08, 0]}>
        <Text
          fontSize={0.06}
          color={isInPath ? "#44ff44" : "#ffffff"}
          anchorX="center"
          anchorY="bottom"
          outlineWidth={0.003}
          outlineColor="#000000"
        >
          {point.label}
          {pathIndex !== null ? ` (${pathIndex + 1})` : ""}
        </Text>
      </Billboard>
    </group>
  );
}

interface BZWireframeProps {
  vertices: Vec3[];
  edges: [number, number][];
  scale: number;  // Uniform scale factor for visualization
}

function BZWireframe({ vertices, edges, scale }: BZWireframeProps) {
  const lineSegments = useMemo(() => {
    const segments: [THREE.Vector3, THREE.Vector3][] = [];

    for (const [i, j] of edges) {
      const v1 = vertices[i];
      const v2 = vertices[j];
      if (v1 && v2) {
        // Vertices are already in Cartesian coordinates, just apply scale
        segments.push([
          new THREE.Vector3(v1[0] * scale, v1[1] * scale, v1[2] * scale),
          new THREE.Vector3(v2[0] * scale, v2[1] * scale, v2[2] * scale),
        ]);
      }
    }

    return segments;
  }, [vertices, edges, scale]);

  return (
    <group>
      {lineSegments.map((segment, i) => (
        <Line
          key={i}
          points={[segment[0], segment[1]]}
          color="#666666"
          lineWidth={1}
          transparent
          opacity={0.5}
        />
      ))}
    </group>
  );
}

interface KPathLinesProps {
  path: KPathPoint[];
  reciprocalBasis: [Vec3, Vec3, Vec3];
}

function KPathLines({ path, reciprocalBasis }: KPathLinesProps) {
  const points = useMemo(() => {
    return path.map((p) => {
      const cart = fractionalToCartesian(p.coords, reciprocalBasis);
      return new THREE.Vector3(...cart);
    });
  }, [path, reciprocalBasis]);

  if (points.length < 2) return null;

  return (
    <Line
      points={points}
      color="#ff6600"
      lineWidth={3}
    />
  );
}

/**
 * Helper component for a single axis with arrow
 */
interface AxisArrowProps {
  direction: Vec3;
  color: string;
  label: string;
  length: number;
}

function AxisArrow({ direction, color, label, length }: AxisArrowProps) {
  const mag = Math.sqrt(direction[0] ** 2 + direction[1] ** 2 + direction[2] ** 2);
  if (mag === 0) return null;

  // Normalize direction
  const norm: Vec3 = [direction[0] / mag, direction[1] / mag, direction[2] / mag];

  // End point of the line
  const endPoint = new THREE.Vector3(norm[0] * length, norm[1] * length, norm[2] * length);

  // Arrow cone position (at the end)
  const arrowPos = endPoint.clone();

  // Calculate rotation to point cone along axis direction
  const quaternion = new THREE.Quaternion();
  const up = new THREE.Vector3(0, 1, 0);
  const dir = new THREE.Vector3(norm[0], norm[1], norm[2]);
  quaternion.setFromUnitVectors(up, dir);

  // Label position (slightly past the arrow)
  const labelPos: [number, number, number] = [
    norm[0] * (length + 0.08),
    norm[1] * (length + 0.08),
    norm[2] * (length + 0.08),
  ];

  return (
    <group>
      <Line
        points={[new THREE.Vector3(0, 0, 0), endPoint]}
        color={color}
        lineWidth={1}
      />
      <Cone
        args={[0.02, 0.06, 8]}
        position={arrowPos}
        quaternion={quaternion}
      >
        <meshBasicMaterial color={color} />
      </Cone>
      <Billboard position={labelPos}>
        <Text
          fontSize={0.05}
          color={color}
          anchorX="center"
          anchorY="middle"
        >
          {label}
        </Text>
      </Billboard>
    </group>
  );
}

interface AxesHelperProps {
  reciprocalBasis: [Vec3, Vec3, Vec3];
}

function ReciprocallatticeAxes({ reciprocalBasis }: AxesHelperProps) {
  const [b1, b2, b3] = reciprocalBasis;

  // Axis length - extend well past the BZ (which is scaled to ~0.7)
  const axisLength = 1.0;

  return (
    <group>
      {/* Reciprocal lattice axes */}
      <AxisArrow direction={b1} color="#ff4444" label="b1" length={axisLength} />
      <AxisArrow direction={b2} color="#44ff44" label="b2" length={axisLength} />
      <AxisArrow direction={b3} color="#4444ff" label="b3" length={axisLength} />

      {/* Cartesian axes in dark gray */}
      <AxisArrow direction={[1, 0, 0]} color="#666666" label="x" length={axisLength} />
      <AxisArrow direction={[0, 1, 0]} color="#666666" label="y" length={axisLength} />
      <AxisArrow direction={[0, 0, 1]} color="#666666" label="z" length={axisLength} />
    </group>
  );
}

// ============================================================================
// Main Scene Component
// ============================================================================

interface BZSceneProps {
  bzData: BrillouinZoneData;
  bzGeometry: BrillouinZoneGeometry;
  primitiveRecipBasis: [Vec3, Vec3, Vec3];  // Primitive reciprocal lattice for k-points
  path: KPathPoint[];
  onPointClick: (point: HighSymmetryPoint) => void;
  selectedPoint: HighSymmetryPoint | null;
  useOrthographic: boolean;
}

function BZScene({
  bzData,
  bzGeometry,
  primitiveRecipBasis,
  path,
  onPointClick,
  selectedPoint,
  useOrthographic,
}: BZSceneProps) {
  const [, setHoveredPoint] = useState<string | null>(null);

  // Calculate a uniform scale factor for visualization
  // Use the BZ geometry extent to determine the scale
  const visualScale = useMemo(() => {
    if (bzGeometry.vertices.length === 0) {
      // Fallback: use reciprocal lattice vector magnitude
      const mags = primitiveRecipBasis.map((v) =>
        Math.sqrt(v[0] ** 2 + v[1] ** 2 + v[2] ** 2)
      );
      return 1 / Math.max(...mags);
    }

    // Find the maximum distance of any BZ vertex from origin
    let maxDist = 0;
    for (const v of bzGeometry.vertices) {
      const dist = Math.sqrt(v[0] ** 2 + v[1] ** 2 + v[2] ** 2);
      if (dist > maxDist) maxDist = dist;
    }

    // Scale so the BZ fits nicely in view (max vertex at distance ~0.7)
    return maxDist > 0 ? 0.7 / maxDist : 1;
  }, [bzGeometry.vertices, primitiveRecipBasis]);

  // Scale the primitive reciprocal basis for k-point placement
  // K-point coordinates in standard tables are in primitive reciprocal lattice coords
  const scaledBasis = useMemo((): [Vec3, Vec3, Vec3] => {
    return primitiveRecipBasis.map((v) => [
      v[0] * visualScale,
      v[1] * visualScale,
      v[2] * visualScale,
    ]) as [Vec3, Vec3, Vec3];
  }, [primitiveRecipBasis, visualScale]);

  return (
    <>
      {/* Lighting */}
      <ambientLight intensity={0.6} />
      <pointLight position={[10, 10, 10]} intensity={0.8} />
      <pointLight position={[-10, -10, -10]} intensity={0.4} />

      {/* BZ wireframe - calculated geometry */}
      {bzGeometry.vertices.length > 0 && bzGeometry.edges.length > 0 && (
        <BZWireframe
          vertices={bzGeometry.vertices}
          edges={bzGeometry.edges}
          scale={visualScale}
        />
      )}

      {/* Reciprocal lattice axes */}
      <ReciprocallatticeAxes reciprocalBasis={scaledBasis} />

      {/* High-symmetry points */}
      {bzData.points.map((point) => {
        const cartesian = fractionalToCartesian(point.coords, scaledBasis);
        const pathIndex = path.findIndex((p) => p.label === point.label);

        return (
          <HighSymmetryPointMesh
            key={point.label}
            point={point}
            position={cartesian as [number, number, number]}
            isSelected={selectedPoint?.label === point.label}
            isInPath={pathIndex >= 0}
            pathIndex={pathIndex >= 0 ? pathIndex : null}
            onClick={() => onPointClick(point)}
            onHover={(hovered) =>
              setHoveredPoint(hovered ? point.label : null)
            }
          />
        );
      })}

      {/* K-path lines */}
      <KPathLines path={path} reciprocalBasis={scaledBasis} />

      {/* Camera */}
      {useOrthographic ? (
        <OrthographicCamera makeDefault position={[2, 2, 2]} zoom={200} />
      ) : (
        <PerspectiveCamera makeDefault position={[1.5, 1.5, 1.5]} fov={50} />
      )}

      {/* Camera controls */}
      <OrbitControls
        enablePan={true}
        enableZoom={true}
        enableRotate={true}
        minDistance={0.5}
        maxDistance={5}
      />
    </>
  );
}

// ============================================================================
// Main Component
// ============================================================================

export function BrillouinZoneViewer({
  crystalData,
  onPathChange,
  initialPath = [],
  pointsPerSegment = 20,
}: BrillouinZoneViewerProps) {
  const [path, setPath] = useState<KPathPoint[]>(initialPath);
  const [selectedPoint, setSelectedPoint] = useState<HighSymmetryPoint | null>(null);
  const [useOrthographic, setUseOrthographic] = useState(true);

  // Determine Bravais lattice type first (needed for centering)
  const bravaisInfo = useMemo(() => {
    const spaceGroup = crystalData.space_group_IT_number;
    const bravaisType: BravaisLattice = spaceGroup
      ? detectBravaisLattice(spaceGroup)
      : "triclinic";

    // Map BravaisLattice type to BravaisLatticeType and centering
    const typeMap: Record<BravaisLattice, { latticeType: BravaisLatticeType; centering: CenteringType }> = {
      "cubic-P": { latticeType: "cP", centering: "P" },
      "cubic-F": { latticeType: "cF", centering: "F" },
      "cubic-I": { latticeType: "cI", centering: "I" },
      "tetragonal-P": { latticeType: "tP", centering: "P" },
      "tetragonal-I": { latticeType: "tI", centering: "I" },
      "orthorhombic-P": { latticeType: "oP", centering: "P" },
      "orthorhombic-C": { latticeType: "oC", centering: "C" },
      "orthorhombic-I": { latticeType: "oI", centering: "I" },
      "orthorhombic-F": { latticeType: "oF", centering: "F" },
      "hexagonal": { latticeType: "hP", centering: "P" },
      "trigonal-R": { latticeType: "hR", centering: "R" },
      "monoclinic-P": { latticeType: "mP", centering: "P" },
      "monoclinic-C": { latticeType: "mC", centering: "C" },
      "triclinic": { latticeType: "aP", centering: "P" },
    };

    return typeMap[bravaisType] || { latticeType: "aP" as BravaisLatticeType, centering: "P" as CenteringType };
  }, [crystalData.space_group_IT_number]);

  // Compute lattice data
  const latticeData = useMemo(() => {
    const a = crystalData.cell_length_a.value;
    const b = crystalData.cell_length_b.value;
    const c = crystalData.cell_length_c.value;
    const alpha = crystalData.cell_angle_alpha.value;
    const beta = crystalData.cell_angle_beta.value;
    const gamma = crystalData.cell_angle_gamma.value;

    // Conventional cell vectors
    const conventionalLattice = realSpaceLatticeVectors(a, b, c, alpha, beta, gamma);

    // Convert to primitive cell for BZ calculation
    const primitiveLattice = conventionalToPrimitive(conventionalLattice, bravaisInfo.centering);

    // Reciprocal lattice of the PRIMITIVE cell (for BZ calculation)
    const primitiveRecipLattice = reciprocalLatticeVectors(primitiveLattice);

    // Reciprocal lattice of the CONVENTIONAL cell (for k-point coordinates)
    // Standard k-point tables use conventional reciprocal lattice coordinates
    const conventionalRecipLattice = reciprocalLatticeVectors(conventionalLattice);

    return {
      conventionalLattice,
      primitiveLattice,
      primitiveRecipLattice,
      conventionalRecipLattice,
      a, b, c, alpha, beta, gamma,
    };
  }, [crystalData, bravaisInfo.centering]);

  // Get BZ high-symmetry point data (uses conventional coordinates)
  const bzData = useMemo(() => {
    return getBrillouinZoneData(bravaisInfo.latticeType, {
      a: latticeData.a,
      b: latticeData.b,
      c: latticeData.c,
      alpha: latticeData.alpha,
      beta: latticeData.beta,
      gamma: latticeData.gamma,
    });
  }, [bravaisInfo.latticeType, latticeData]);

  // Calculate BZ geometry (Wigner-Seitz cell) from PRIMITIVE reciprocal lattice
  const bzGeometry = useMemo(() => {
    return calculateBrillouinZone(latticeData.primitiveRecipLattice);
  }, [latticeData.primitiveRecipLattice]);

  // Handle point click - add to path
  const handlePointClick = useCallback(
    (point: HighSymmetryPoint) => {
      setSelectedPoint(point);

      // Check if point is already last in path
      if (path.length > 0 && path[path.length - 1].label === point.label) {
        return;
      }

      // Add point to path
      // npoints = number of k-points from THIS point to the NEXT point
      // Last point should have npoints = 0 (no segment after it)
      const newPoint: KPathPoint = {
        label: point.label,
        coords: point.coords,
        npoints: 0,  // New point is last, so no segment after it
      };

      // Update previous last point to have pointsPerSegment (now has a segment after it)
      const newPath = path.map((p, i) =>
        i === path.length - 1 ? { ...p, npoints: pointsPerSegment } : p
      );
      newPath.push(newPoint);

      setPath(newPath);
      onPathChange(newPath);
    },
    [path, pointsPerSegment, onPathChange]
  );

  // Remove last point from path
  const handleUndo = useCallback(() => {
    if (path.length > 0) {
      const newPath = path.slice(0, -1);
      setPath(newPath);
      onPathChange(newPath);
      setSelectedPoint(null);
    }
  }, [path, onPathChange]);

  // Clear entire path
  const handleClear = useCallback(() => {
    setPath([]);
    onPathChange([]);
    setSelectedPoint(null);
  }, [onPathChange]);

  // Use recommended path
  const handleUseRecommended = useCallback(() => {
    const recommendedPath: KPathPoint[] = [];

    for (let i = 0; i < bzData.recommendedPath.length; i++) {
      const [from, to] = bzData.recommendedPath[i];

      // Add 'from' point if it's the first or different from last added
      if (recommendedPath.length === 0 || recommendedPath[recommendedPath.length - 1].label !== from) {
        const fromPoint = bzData.points.find((p) => p.label === from);
        if (fromPoint) {
          recommendedPath.push({
            label: fromPoint.label,
            coords: fromPoint.coords,
            npoints: pointsPerSegment,
          });
        }
      }

      // Add 'to' point
      const toPoint = bzData.points.find((p) => p.label === to);
      if (toPoint) {
        recommendedPath.push({
          label: toPoint.label,
          coords: toPoint.coords,
          npoints: pointsPerSegment,
        });
      }
    }

    // Remove duplicate consecutive points and set last npoints to 0
    const cleanedPath: KPathPoint[] = [];
    for (const point of recommendedPath) {
      if (cleanedPath.length === 0 || cleanedPath[cleanedPath.length - 1].label !== point.label) {
        cleanedPath.push(point);
      }
    }
    if (cleanedPath.length > 0) {
      cleanedPath[cleanedPath.length - 1].npoints = 0;
    }

    setPath(cleanedPath);
    onPathChange(cleanedPath);
  }, [bzData, pointsPerSegment, onPathChange]);

  // Format path for display
  const pathString = path.map((p) => p.label).join(" → ");

  return (
    <div className="bz-viewer">
      <div className="bz-viewer-header">
        <h4>Brillouin Zone - {bzData.name}</h4>
        <div className="bz-viewer-info">
          <span>Click points to build k-path</span>
        </div>
      </div>

      <div className="bz-viewer-canvas">
        <Canvas style={{ background: "#1a1a2e" }}>
          <BZScene
            bzData={bzData}
            bzGeometry={bzGeometry}
            primitiveRecipBasis={latticeData.primitiveRecipLattice}
            path={path}
            onPointClick={handlePointClick}
            selectedPoint={selectedPoint}
            useOrthographic={useOrthographic}
          />
        </Canvas>
      </div>

      <div className="bz-viewer-controls">
        <div className="bz-path-display">
          <label>Current Path:</label>
          <span className="bz-path-string">
            {pathString || "(none - click points to build path)"}
          </span>
        </div>

        <div className="bz-viewer-buttons">
          <button onClick={handleUseRecommended} className="bz-btn recommended">
            Use Recommended Path
          </button>
          <button onClick={handleUndo} disabled={path.length === 0} className="bz-btn">
            Undo
          </button>
          <button onClick={handleClear} disabled={path.length === 0} className="bz-btn">
            Clear
          </button>
          <button
            onClick={() => setUseOrthographic(!useOrthographic)}
            className="bz-btn"
          >
            {useOrthographic ? "Perspective" : "Orthographic"}
          </button>
        </div>
      </div>

      <div className="bz-points-list">
        <h5>High-Symmetry Points</h5>
        <div className="bz-points-grid">
          {bzData.points.map((point) => (
            <div
              key={point.label}
              className={`bz-point-item ${
                path.some((p) => p.label === point.label) ? "in-path" : ""
              }`}
              onClick={() => handlePointClick(point)}
            >
              <span className="bz-point-label">{point.label}</span>
              <span className="bz-point-coords">
                ({point.coords.map((c) => c.toFixed(3)).join(", ")})
              </span>
              <span className="bz-point-desc">{point.description}</span>
            </div>
          ))}
        </div>
      </div>
    </div>
  );
}

export default BrillouinZoneViewer;
