// 3D Unit Cell Viewer - ported from Reciprocal project

import { useRef, useMemo, useState } from "react";
import { Canvas, useFrame, useThree } from "@react-three/fiber";
import { OrbitControls } from "@react-three/drei";
import * as THREE from "three";
import { CrystalData } from "../lib/types";
import { expandAtomSites } from "../lib/symmetry";

// Element colors based on CPK coloring convention
const ELEMENT_COLORS: Record<string, string> = {
  H: "#FFFFFF",
  He: "#D9FFFF",
  Li: "#CC80FF",
  Be: "#C2FF00",
  B: "#FFB5B5",
  C: "#909090",
  N: "#3050F8",
  O: "#FF0D0D",
  F: "#90E050",
  Ne: "#B3E3F5",
  Na: "#AB5CF2",
  Mg: "#8AFF00",
  Al: "#BFA6A6",
  Si: "#F0C8A0",
  P: "#FF8000",
  S: "#FFFF30",
  Cl: "#1FF01F",
  Ar: "#80D1E3",
  K: "#8F40D4",
  Ca: "#3DFF00",
  Sc: "#E6E6E6",
  Ti: "#BFC2C7",
  V: "#A6A6AB",
  Cr: "#8A99C7",
  Mn: "#9C7AC7",
  Fe: "#E06633",
  Co: "#F090A0",
  Ni: "#50D050",
  Cu: "#C88033",
  Zn: "#7D80B0",
  Ga: "#C28F8F",
  Ge: "#668F8F",
  As: "#BD80E3",
  Se: "#FFA100",
  Br: "#A62929",
  Kr: "#5CB8D1",
  Rb: "#702EB0",
  Sr: "#00FF00",
  Y: "#94FFFF",
  Zr: "#94E0E0",
  Nb: "#73C2C9",
  Mo: "#54B5B5",
  Tc: "#3B9E9E",
  Ru: "#248F8F",
  Rh: "#0A7D8C",
  Pd: "#006985",
  Ag: "#C0C0C0",
  Cd: "#FFD98F",
  In: "#A67573",
  Sn: "#668080",
  Sb: "#9E63B5",
  Te: "#D47A00",
  I: "#940094",
  Xe: "#429EB0",
  Cs: "#57178F",
  Ba: "#00C900",
  La: "#70D4FF",
  Ce: "#FFFFC7",
  Pr: "#D9FFC7",
  Nd: "#C7FFC7",
  Pm: "#A3FFC7",
  Sm: "#8FFFC7",
  Eu: "#61FFC7",
  Gd: "#45FFC7",
  Tb: "#30FFC7",
  Dy: "#1FFFC7",
  Ho: "#00FF9C",
  Er: "#00E675",
  Tm: "#00D452",
  Yb: "#00BF38",
  Lu: "#00AB24",
  Hf: "#4DC2FF",
  Ta: "#4DA6FF",
  W: "#2194D6",
  Re: "#267DAB",
  Os: "#266696",
  Ir: "#175487",
  Pt: "#D0D0E0",
  Au: "#FFD123",
  Hg: "#B8B8D0",
  Tl: "#A6544D",
  Pb: "#575961",
  Bi: "#9E4FB5",
  Po: "#AB5C00",
  At: "#754F45",
  Rn: "#428296",
  Fr: "#420066",
  Ra: "#007D00",
  Ac: "#70ABFA",
  Th: "#00BAFF",
  Pa: "#00A1FF",
  U: "#008FFF",
  Np: "#0080FF",
  Pu: "#006BFF",
  Am: "#545CF2",
  Cm: "#785CE3",
  Bk: "#8A4FE3",
  Cf: "#A136D4",
  Es: "#B31FD4",
  Fm: "#B31FBA",
  Md: "#B30DA6",
  No: "#BD0D87",
  Lr: "#C70066",
};

const DEFAULT_COLOR = "#FF00FF";

function getElementColor(symbol: string): string {
  const element = symbol.replace(/\d+[+-]?$/, "").replace(/[+-]$/, "");
  return ELEMENT_COLORS[element] || DEFAULT_COLOR;
}

function latticeToCartesian(
  a: number,
  b: number,
  c: number,
  alpha: number,
  beta: number,
  gamma: number
): { aVec: THREE.Vector3; bVec: THREE.Vector3; cVec: THREE.Vector3 } {
  const alphaRad = (alpha * Math.PI) / 180;
  const betaRad = (beta * Math.PI) / 180;
  const gammaRad = (gamma * Math.PI) / 180;

  const aVec = new THREE.Vector3(a, 0, 0);
  const bVec = new THREE.Vector3(
    b * Math.cos(gammaRad),
    b * Math.sin(gammaRad),
    0
  );

  const cx = c * Math.cos(betaRad);
  const cy =
    (c * (Math.cos(alphaRad) - Math.cos(betaRad) * Math.cos(gammaRad))) /
    Math.sin(gammaRad);
  const cz = Math.sqrt(c * c - cx * cx - cy * cy);
  const cVec = new THREE.Vector3(cx, cy, cz);

  return { aVec, bVec, cVec };
}

function fractToCartesian(
  fx: number,
  fy: number,
  fz: number,
  aVec: THREE.Vector3,
  bVec: THREE.Vector3,
  cVec: THREE.Vector3
): THREE.Vector3 {
  return new THREE.Vector3()
    .addScaledVector(aVec, fx)
    .addScaledVector(bVec, fy)
    .addScaledVector(cVec, fz);
}

interface SceneProps {
  crystalData: CrystalData;
  autoRotate: boolean;
  onInteraction: () => void;
}

function Scene({ crystalData, autoRotate, onInteraction }: SceneProps) {
  const mainGroupRef = useRef<THREE.Group>(null);
  const axisGroupRef = useRef<THREE.Group>(null);
  const { camera } = useThree();

  const { aDir, bDir, cDir } = useMemo(() => {
    const a = crystalData.cell_length_a.value;
    const b = crystalData.cell_length_b.value;
    const c = crystalData.cell_length_c.value;
    const alpha = crystalData.cell_angle_alpha.value;
    const beta = crystalData.cell_angle_beta.value;
    const gamma = crystalData.cell_angle_gamma.value;

    const { aVec, bVec, cVec } = latticeToCartesian(a, b, c, alpha, beta, gamma);

    return {
      aDir: aVec.clone().normalize(),
      bDir: bVec.clone().normalize(),
      cDir: cVec.clone().normalize(),
    };
  }, [crystalData]);

  const { aVec, bVec, cVec, center, maxDim } = useMemo(() => {
    const a = crystalData.cell_length_a.value;
    const b = crystalData.cell_length_b.value;
    const c = crystalData.cell_length_c.value;
    const alpha = crystalData.cell_angle_alpha.value;
    const beta = crystalData.cell_angle_beta.value;
    const gamma = crystalData.cell_angle_gamma.value;

    const vectors = latticeToCartesian(a, b, c, alpha, beta, gamma);

    const center = new THREE.Vector3()
      .add(vectors.aVec)
      .add(vectors.bVec)
      .add(vectors.cVec)
      .multiplyScalar(0.5);

    const maxDim = Math.max(
      vectors.aVec.length(),
      vectors.bVec.length(),
      vectors.cVec.length()
    );

    return { ...vectors, center, maxDim };
  }, [crystalData]);

  const wireframePoints = useMemo(() => {
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
    ];
  }, [aVec, bVec, cVec]);

  const expandedAtoms = useMemo(() => {
    return expandAtomSites(
      crystalData.atom_sites,
      crystalData.symmetry_operations
    );
  }, [crystalData.atom_sites, crystalData.symmetry_operations]);

  // Calculate proper camera distance to fit the unit cell
  useMemo(() => {
    const fov = 35 * (Math.PI / 180);
    const fitDistance = (maxDim * 1.5) / Math.tan(fov / 2);
    camera.position.set(
      Math.sin((10 * Math.PI) / 180) * fitDistance,
      fitDistance * 0.2,
      Math.cos((10 * Math.PI) / 180) * fitDistance
    );
    camera.lookAt(0, 0, 0);
    camera.updateProjectionMatrix();
  }, [maxDim, camera]);

  useFrame((_, delta) => {
    if (autoRotate) {
      if (mainGroupRef.current) {
        mainGroupRef.current.rotation.y += delta * 0.3;
      }
      if (axisGroupRef.current) {
        axisGroupRef.current.rotation.y += delta * 0.3;
      }
    }
  });

  const scale = 4 / maxDim;
  const axisLength = 0.7;

  return (
    <>
      <OrbitControls
        enablePan={true}
        enableZoom={true}
        enableRotate={true}
        onStart={onInteraction}
        minDistance={1}
        maxDistance={50}
      />
      <group
        ref={mainGroupRef}
        scale={[scale, scale, scale]}
        position={[0, 0, 0]}
      >
        <group position={[-center.x, -center.y, -center.z]}>
          {wireframePoints.map((edge, i) => (
            <line key={`edge-${i}`}>
              <bufferGeometry>
                <bufferAttribute
                  attach="attributes-position"
                  args={[
                    new Float32Array([
                      edge[0].x,
                      edge[0].y,
                      edge[0].z,
                      edge[1].x,
                      edge[1].y,
                      edge[1].z,
                    ]),
                    3,
                  ]}
                />
              </bufferGeometry>
              <lineBasicMaterial color="#666666" linewidth={1} />
            </line>
          ))}

          {expandedAtoms.map((atom, i) => {
            const pos = fractToCartesian(
              atom.fract_x,
              atom.fract_y,
              atom.fract_z,
              aVec,
              bVec,
              cVec
            );
            const color = getElementColor(atom.type_symbol);
            const radius = maxDim * 0.04;

            return (
              <mesh key={`atom-${i}`} position={[pos.x, pos.y, pos.z]}>
                <sphereGeometry args={[radius, 32, 32]} />
                <meshStandardMaterial color={color} />
              </mesh>
            );
          })}
        </group>
      </group>

      <group position={[-2.5, -2.0, 0]}>
        <group ref={axisGroupRef}>
          <arrowHelper
            args={[
              aDir,
              new THREE.Vector3(0, 0, 0),
              axisLength,
              0xff4444,
              0.15,
              0.08,
            ]}
          />
          <arrowHelper
            args={[
              bDir,
              new THREE.Vector3(0, 0, 0),
              axisLength,
              0x44ff44,
              0.15,
              0.08,
            ]}
          />
          <arrowHelper
            args={[
              cDir,
              new THREE.Vector3(0, 0, 0),
              axisLength,
              0x4444ff,
              0.15,
              0.08,
            ]}
          />
        </group>
        <AxisLabels
          aDir={aDir}
          bDir={bDir}
          cDir={cDir}
          axisLength={axisLength}
          groupRef={axisGroupRef}
        />
      </group>
    </>
  );
}

function AxisLabels({
  aDir,
  bDir,
  cDir,
  axisLength,
  groupRef,
}: {
  aDir: THREE.Vector3;
  bDir: THREE.Vector3;
  cDir: THREE.Vector3;
  axisLength: number;
  groupRef: React.RefObject<THREE.Group | null>;
}) {
  const aLabelRef = useRef<THREE.Sprite>(null);
  const bLabelRef = useRef<THREE.Sprite>(null);
  const cLabelRef = useRef<THREE.Sprite>(null);

  const createTextTexture = (text: string, color: string) => {
    const canvas = document.createElement("canvas");
    canvas.width = 64;
    canvas.height = 64;
    const ctx = canvas.getContext("2d")!;
    ctx.fillStyle = color;
    ctx.font = "bold 48px Arial";
    ctx.textAlign = "center";
    ctx.textBaseline = "middle";
    ctx.fillText(text, 32, 32);
    const texture = new THREE.CanvasTexture(canvas);
    return texture;
  };

  const aTexture = useMemo(() => createTextTexture("a", "#ff4444"), []);
  const bTexture = useMemo(() => createTextTexture("b", "#44ff44"), []);
  const cTexture = useMemo(() => createTextTexture("c", "#4444ff"), []);

  useFrame(() => {
    if (groupRef.current) {
      const rotation = groupRef.current.rotation.clone();
      const labelOffset = axisLength + 0.18;

      if (aLabelRef.current) {
        const pos = aDir
          .clone()
          .multiplyScalar(labelOffset)
          .applyEuler(rotation);
        aLabelRef.current.position.copy(pos);
      }
      if (bLabelRef.current) {
        const pos = bDir
          .clone()
          .multiplyScalar(labelOffset)
          .applyEuler(rotation);
        bLabelRef.current.position.copy(pos);
      }
      if (cLabelRef.current) {
        const pos = cDir
          .clone()
          .multiplyScalar(labelOffset)
          .applyEuler(rotation);
        cLabelRef.current.position.copy(pos);
      }
    }
  });

  return (
    <>
      <sprite ref={aLabelRef} scale={[0.22, 0.22, 1]}>
        <spriteMaterial map={aTexture} transparent />
      </sprite>
      <sprite ref={bLabelRef} scale={[0.22, 0.22, 1]}>
        <spriteMaterial map={bTexture} transparent />
      </sprite>
      <sprite ref={cLabelRef} scale={[0.22, 0.22, 1]}>
        <spriteMaterial map={cTexture} transparent />
      </sprite>
    </>
  );
}

interface UnitCellViewerProps {
  crystalData: CrystalData;
}

export function UnitCellViewer({ crystalData }: UnitCellViewerProps) {
  const [autoRotate, setAutoRotate] = useState(true);

  const handleInteraction = () => {
    setAutoRotate(false);
  };

  return (
    <div style={{ width: "100%", height: "100%", background: "#1a1a2e", borderRadius: "12px", overflow: "hidden" }}>
      <Canvas
        camera={{
          position: [10, 2, 10],
          fov: 35,
          near: 0.1,
          far: 100,
        }}
      >
        <ambientLight intensity={0.6} />
        <directionalLight position={[5, 5, 5]} intensity={0.8} />
        <directionalLight position={[-5, -5, -5]} intensity={0.3} />
        <Scene
          crystalData={crystalData}
          autoRotate={autoRotate}
          onInteraction={handleInteraction}
        />
      </Canvas>
    </div>
  );
}
