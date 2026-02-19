/**
 * Interactive 3D Brillouin Zone Viewer with k-path selection.
 *
 * Displays the first Brillouin zone with high-symmetry points and allows
 * users to click points to build a k-path for band structure calculations.
 */

import { useRef, useState, useMemo, useCallback, useEffect } from "react";
import { Canvas, useFrame } from "@react-three/fiber";
import { OrbitControls, Text, Line, OrthographicCamera, PerspectiveCamera, Cone, Billboard } from "@react-three/drei";
import * as THREE from "three";
import { CrystalData } from "../lib/types";
import {
  Vec3,
  dot,
  magnitude,
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
  findHighSymmetryPoint,
} from "../lib/brillouinZoneData";
import { detectBravaisLattice, BravaisLattice } from "../lib/brillouinZone";
import { SymmetryTransformResult } from "../lib/symmetryTransform";

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
  symmetryTransform?: SymmetryTransformResult | null;
}

function coerceSpaceGroupNumber(value: unknown): number | null {
  if (value == null) return null;
  const parsed = typeof value === "number"
    ? value
    : Number.parseInt(String(value).trim(), 10);
  if (!Number.isInteger(parsed) || parsed < 1 || parsed > 230) {
    return null;
  }
  return parsed;
}

function extractSpaceGroupNumberFromHM(hm: string | undefined): number | null {
  if (!hm) return null;
  const match = hm.match(/#\s*(\d{1,3})/);
  if (!match) return null;
  return coerceSpaceGroupNumber(match[1]);
}

function normalizeSpaceGroupHM(value: string | undefined): string {
  if (!value) return "";
  return value
    .toLowerCase()
    .replace(/[−–—]/g, "-")
    .replace(/[\s_:'"]/g, "");
}

function detectBravaisFromHMSymbol(value: string | undefined): BravaisLattice | null {
  const normalized = normalizeSpaceGroupHM(value);
  if (!normalized) return null;

  // Rhombohedral setting (e.g., "R -3 m H", "R-3m", "R3c")
  if (normalized.startsWith("r")) {
    return "trigonal-R";
  }

  // Common hexagonal HM symbols in P setting.
  if (
    normalized.startsWith("p6") ||
    normalized.startsWith("p-6") ||
    normalized.startsWith("p63")
  ) {
    return "hexagonal";
  }

  return null;
}

const SUBSCRIPT_TO_ASCII: Record<string, string> = {
  "₀": "0",
  "₁": "1",
  "₂": "2",
  "₃": "3",
  "₄": "4",
  "₅": "5",
  "₆": "6",
  "₇": "7",
  "₈": "8",
  "₉": "9",
};

function formatLabelForDisplay(label: string): string {
  return label.replace(/[₀₁₂₃₄₅₆₇₈₉]/g, (char) => SUBSCRIPT_TO_ASCII[char] ?? char);
}

function crystalSystemRankFromSpaceGroup(spaceGroupNumber: number): number {
  if (spaceGroupNumber >= 1 && spaceGroupNumber <= 2) return 0; // triclinic
  if (spaceGroupNumber >= 3 && spaceGroupNumber <= 15) return 1; // monoclinic
  if (spaceGroupNumber >= 16 && spaceGroupNumber <= 74) return 2; // orthorhombic
  if (spaceGroupNumber >= 75 && spaceGroupNumber <= 142) return 3; // tetragonal
  if (spaceGroupNumber >= 143 && spaceGroupNumber <= 167) return 4; // trigonal
  if (spaceGroupNumber >= 168 && spaceGroupNumber <= 194) return 5; // hexagonal
  if (spaceGroupNumber >= 195 && spaceGroupNumber <= 230) return 6; // cubic
  return -1;
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

  const displayLabel = formatLabelForDisplay(point.label);
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
          {displayLabel}
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
  const polylineSegments = useMemo(() => {
    const segments: THREE.Vector3[][] = [];
    let currentSegment: THREE.Vector3[] = [];

    for (let i = 0; i < path.length; i++) {
      const point = path[i];
      const cart = fractionalToCartesian(point.coords, reciprocalBasis);
      currentSegment.push(new THREE.Vector3(...cart));

      // npoints=0 marks the end of a segment (path discontinuity)
      if (i < path.length - 1 && point.npoints === 0) {
        if (currentSegment.length >= 2) {
          segments.push(currentSegment);
        }
        currentSegment = [];
      }
    }

    if (currentSegment.length >= 2) {
      segments.push(currentSegment);
    }

    return segments;
  }, [path, reciprocalBasis]);

  if (polylineSegments.length === 0) return null;

  return (
    <group>
      {polylineSegments.map((segment, i) => (
        <Line
          key={i}
          points={segment}
          color="#ff6600"
          lineWidth={3}
        />
      ))}
    </group>
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
  symmetryTransform = null,
}: BrillouinZoneViewerProps) {
  const [path, setPath] = useState<KPathPoint[]>(initialPath);
  const [selectedPoint, setSelectedPoint] = useState<HighSymmetryPoint | null>(null);
  const [useOrthographic, setUseOrthographic] = useState(true);

  useEffect(() => {
    setPath(initialPath);
  }, [initialPath]);

  // Determine Bravais lattice type first (needed for centering)
  const bravaisInfo = useMemo(() => {
    const cifSpaceGroup =
      coerceSpaceGroupNumber(crystalData.space_group_IT_number) ??
      extractSpaceGroupNumberFromHM(crystalData.space_group_HM);
    const symmetrySpaceGroup = coerceSpaceGroupNumber(symmetryTransform?.spacegroupNumber);
    const preferCifSpaceGroup =
      cifSpaceGroup != null &&
      symmetrySpaceGroup != null &&
      crystalSystemRankFromSpaceGroup(symmetrySpaceGroup) <
        crystalSystemRankFromSpaceGroup(cifSpaceGroup);

    if (
      preferCifSpaceGroup &&
      symmetrySpaceGroup != null &&
      cifSpaceGroup != null &&
      symmetrySpaceGroup !== cifSpaceGroup
    ) {
      console.warn(
        `Ignoring downgraded spglib space group ${symmetrySpaceGroup} in favor of CIF space group ${cifSpaceGroup} for BZ labeling.`,
      );
    }

    const spaceGroup = preferCifSpaceGroup
      ? cifSpaceGroup
      : (symmetrySpaceGroup ?? cifSpaceGroup);

    const hmFallback = preferCifSpaceGroup
      ? (detectBravaisFromHMSymbol(crystalData.space_group_HM) ??
        detectBravaisFromHMSymbol(symmetryTransform?.internationalSymbol))
      : (detectBravaisFromHMSymbol(symmetryTransform?.internationalSymbol) ??
        detectBravaisFromHMSymbol(crystalData.space_group_HM));

    let bravaisType: BravaisLattice = hmFallback ?? "triclinic";
    if (spaceGroup != null) {
      try {
        bravaisType = detectBravaisLattice(spaceGroup);
      } catch (error) {
        console.warn("Failed to detect Bravais lattice from space group:", spaceGroup, error);
        if (hmFallback) {
          bravaisType = hmFallback;
        }
      }
    }

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

    const mapped =
      typeMap[bravaisType] || {
        latticeType: "aP" as BravaisLatticeType,
        centering: "P" as CenteringType,
      };

    const useSymmetryTransform = symmetryTransform != null && !preferCifSpaceGroup;
    return { ...mapped, useSymmetryTransform };
  }, [
    crystalData.space_group_HM,
    crystalData.space_group_IT_number,
    symmetryTransform?.internationalSymbol,
    symmetryTransform?.spacegroupNumber,
  ]);

  // Compute lattice data
  const latticeData = useMemo(() => {
    const angleDegrees = (u: Vec3, v: Vec3): number => {
      const denom = magnitude(u) * magnitude(v);
      if (denom <= 0) return 0;
      const cosine = Math.max(-1, Math.min(1, dot(u, v) / denom));
      return (Math.acos(cosine) * 180) / Math.PI;
    };

    const useSymmetryTransform = bravaisInfo.useSymmetryTransform && symmetryTransform != null;
    const conventionalLattice = useSymmetryTransform
      ? symmetryTransform.standardizedConventionalLattice
      : realSpaceLatticeVectors(
          crystalData.cell_length_a.value,
          crystalData.cell_length_b.value,
          crystalData.cell_length_c.value,
          crystalData.cell_angle_alpha.value,
          crystalData.cell_angle_beta.value,
          crystalData.cell_angle_gamma.value,
        );

    // Convert to primitive cell for BZ calculation unless backend provided one.
    const primitiveLattice = useSymmetryTransform
      ? symmetryTransform.standardizedPrimitiveLattice
      : conventionalToPrimitive(conventionalLattice, bravaisInfo.centering);
    const primitiveAlpha = angleDegrees(primitiveLattice[1], primitiveLattice[2]);

    // Reciprocal lattice of the PRIMITIVE cell (for BZ calculation)
    const primitiveRecipLattice = reciprocalLatticeVectors(primitiveLattice);

    // Reciprocal lattice of the CONVENTIONAL cell (for k-point coordinates)
    // Standard k-point tables use conventional reciprocal lattice coordinates
    const conventionalRecipLattice = reciprocalLatticeVectors(conventionalLattice);

    const a = magnitude(conventionalLattice[0]);
    const b = magnitude(conventionalLattice[1]);
    const c = magnitude(conventionalLattice[2]);
    const alpha = angleDegrees(conventionalLattice[1], conventionalLattice[2]);
    const beta = angleDegrees(conventionalLattice[0], conventionalLattice[2]);
    const gamma = angleDegrees(conventionalLattice[0], conventionalLattice[1]);

    return {
      conventionalLattice,
      primitiveLattice,
      primitiveRecipLattice,
      conventionalRecipLattice,
      primitiveAlpha,
      a, b, c, alpha, beta, gamma,
    };
  }, [crystalData, bravaisInfo.centering, bravaisInfo.useSymmetryTransform, symmetryTransform]);

  // Get BZ high-symmetry point data (uses conventional coordinates)
  const bzData = useMemo(() => {
    const alphaForPath =
      bravaisInfo.latticeType === "hR"
        ? latticeData.primitiveAlpha
        : latticeData.alpha;

    return getBrillouinZoneData(bravaisInfo.latticeType, {
      a: latticeData.a,
      b: latticeData.b,
      c: latticeData.c,
      alpha: alphaForPath,
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
      const newPath = path.slice(0, -1).map((point, i, points) =>
        i === points.length - 1 ? { ...point, npoints: 0 } : point
      );
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

    for (const [fromLabel, toLabel] of bzData.recommendedPath) {
      const fromPoint = findHighSymmetryPoint(bzData, fromLabel);
      const toPoint = findHighSymmetryPoint(bzData, toLabel);
      if (!fromPoint || !toPoint) {
        continue;
      }

      if (recommendedPath.length === 0) {
        recommendedPath.push({
          label: fromPoint.label,
          coords: fromPoint.coords,
          npoints: pointsPerSegment,
        });
      } else {
        const lastPoint = recommendedPath[recommendedPath.length - 1];
        if (lastPoint.label !== fromLabel) {
          // Disconnected segment (e.g. ... | U→X): end previous segment and restart.
          lastPoint.npoints = 0;
          recommendedPath.push({
            label: fromPoint.label,
            coords: fromPoint.coords,
            npoints: pointsPerSegment,
          });
        } else if (lastPoint.npoints === 0) {
          // Re-entering from a segment endpoint: reopen interpolation to next point.
          lastPoint.npoints = pointsPerSegment;
        }
      }

      // Avoid duplicate point entries for degenerate segments.
      if (recommendedPath[recommendedPath.length - 1].label !== toLabel) {
        recommendedPath.push({
          label: toPoint.label,
          coords: toPoint.coords,
          npoints: pointsPerSegment,
        });
      }
    }

    if (recommendedPath.length > 0) {
      recommendedPath[recommendedPath.length - 1].npoints = 0;
    }

    setPath(recommendedPath);
    onPathChange(recommendedPath);
  }, [bzData, pointsPerSegment, onPathChange]);

  // Format path for display
  const pathString = useMemo(() => {
    if (path.length === 0) return "";
    let result = formatLabelForDisplay(path[0].label);
    for (let i = 1; i < path.length; i++) {
      const separator = path[i - 1].npoints === 0 ? " | " : " → ";
      result += `${separator}${formatLabelForDisplay(path[i].label)}`;
    }
    return result;
  }, [path]);

  return (
    <div className="bz-viewer">
      <div className="bz-viewer-header">
        <h4>Brillouin Zone - {bzData.name}</h4>
        <div className="bz-viewer-info">
          <span>Click points to build k-path (Setyawan-Curtarolo convention)</span>
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
              <span className="bz-point-label">{formatLabelForDisplay(point.label)}</span>
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
