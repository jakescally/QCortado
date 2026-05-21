import assert from "node:assert/strict";
import test from "node:test";
import {
  getMagneticSpeciesFields,
  getMagnetismViewerData,
  parseAtomicMagneticMomentsFromOutput,
} from "../../src/lib/magnetism";

test("parses final collinear site moments as c-axis vectors", () => {
  const moments = parseAtomicMagneticMomentsFromOutput(`
     Magnetic moment per site:
     atom    1 (R=0.108)  charge= 14.8805  magn=  2.1250
     atom:   2            charge: 12.1000  magn: -1.8750
  `);

  assert.deepEqual(moments.map((moment) => moment.vector), [
    [0, 0, 2.125],
    [0, 0, -1.875],
  ]);
});

test("uses stored magnetic moments and viewer structure when available", () => {
  const viewerData = getMagnetismViewerData({
    parameters: {
      magnetism_viewer_structure: {
        position_units: "crystal",
        atoms: [
          { symbol: "Eu", position: [0, 0, 0] },
          { symbol: "Eu", position: [0.5, 0.5, 0.5] },
        ],
        cell_parameters: [
          [4, 0, 0],
          [0, 5, 0],
          [0, 0, 6],
        ],
        cell_units: "angstrom",
      },
    },
    result: {
      atomic_magnetic_moments: [
        [1, 0, 0],
        [-1, 0, 0],
      ],
    },
  });

  assert.equal(viewerData?.structure.atoms.length, 2);
  assert.deepEqual(viewerData?.moments.map((moment) => moment.vector), [
    [1, 0, 0],
    [-1, 0, 0],
  ]);
});

test("parses atom-number magnetization blocks from noncollinear QE output", () => {
  const moments = parseAtomicMagneticMomentsFromOutput(`
 ==============================================================================
     atom number    1 relative position :    0.0000   0.0000   0.0000
     magnetization :          2.840143    0.099180    2.841836
 ==============================================================================
     atom number    2 relative position :   -0.5000   0.5000   0.5000
     magnetization :         -0.019422   -0.000678   -0.019434
  `);

  assert.deepEqual(moments.map((moment) => moment.vector), [
    [2.840143, 0.09918, 2.841836],
    [-0.019422, -0.000678, -0.019434],
  ]);
});

test("prefers source structure and remaps primitive moments onto conventional atoms", () => {
  const viewerData = getMagnetismViewerData({
    parameters: {
      source_structure: {
        position_units: "crystal",
        atoms: [
          { symbol: "A", position: [0, 0, 0] },
          { symbol: "A", position: [0, 0.5, 0.5] },
          { symbol: "A", position: [0.5, 0, 0.5] },
          { symbol: "A", position: [0.5, 0.5, 0] },
          { symbol: "B", position: [0.5, 0.5, 0.5] },
          { symbol: "B", position: [0.5, 0, 0] },
          { symbol: "B", position: [0, 0.5, 0] },
          { symbol: "B", position: [0, 0, 0.5] },
        ],
        cell_parameters: [
          [1, 0, 0],
          [0, 1, 0],
          [0, 0, 1],
        ],
        cell_units: "angstrom",
      },
      magnetism_viewer_structure: {
        position_units: "crystal",
        atoms: [
          { symbol: "A", position: [0, 0, 0] },
          { symbol: "B", position: [0.5, 0.5, 0.5] },
        ],
        cell_parameters: [
          [0, 0.5, 0.5],
          [0.5, 0, 0.5],
          [0.5, 0.5, 0],
        ],
        cell_units: "angstrom",
      },
    },
    result: {
      atomic_magnetic_moments: [
        [1, 2, 3],
        [0.1, 0.2, 0.3],
      ],
    },
  });

  assert.ok(viewerData);
  assert.equal(viewerData?.structure.atoms.length, 8);
  assert.deepEqual(viewerData?.moments.map((moment) => moment.vector), [
    [1, 2, 3],
    [1, 2, 3],
    [1, 2, 3],
    [1, 2, 3],
    [0.1, 0.2, 0.3],
    [0.1, 0.2, 0.3],
    [0.1, 0.2, 0.3],
    [0.1, 0.2, 0.3],
  ]);
});

test("carries noncollinear species theta and phi from saved SCF parameters", () => {
  assert.deepEqual(
    getMagneticSpeciesFields({
      starting_magnetization: { Fe: 0.7 },
      starting_magnetization_theta: { Fe: 90 },
      starting_magnetization_phi: { Fe: 45 },
    }, "Fe"),
    {
      starting_magnetization: 0.7,
      theta: 90,
      phi: 45,
    },
  );
});
