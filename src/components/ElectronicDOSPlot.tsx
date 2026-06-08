import { useCallback, useEffect, useLayoutEffect, useMemo, useRef, useState } from "react";

export interface ElectronicDOSData {
  energies: number[];
  dos: number[];
  fermi_energy: number | null;
  energy_range: [number, number];
  max_dos: number;
  points: number;
}

interface ElectronicDOSPlotProps {
  data: ElectronicDOSData;
  width?: number;
  height?: number;
}

function formatAxisInputValue(value: number): string {
  if (!Number.isFinite(value)) return "";
  return Number.parseFloat(value.toFixed(6)).toString();
}

function getTickDecimals(step: number): number {
  if (!Number.isFinite(step) || step <= 0) {
    return 1;
  }
  const normalized = formatAxisInputValue(step);
  const decimalPart = normalized.split(".")[1] ?? "";
  return Math.min(6, decimalPart.length);
}

function getTickValues(min: number, max: number, count: number): number[] {
  if (!Number.isFinite(min) || !Number.isFinite(max) || max < min) {
    return [0];
  }
  const steps = Math.max(1, count);
  if (Math.abs(max - min) < 1e-12) {
    return [min];
  }
  return Array.from({ length: steps + 1 }, (_, index) => min + (index / steps) * (max - min));
}

export function ElectronicDOSPlot({
  data,
  width = 900,
  height = 500,
}: ElectronicDOSPlotProps) {
  const canvasRef = useRef<HTMLDivElement>(null);
  const [containerSize, setContainerSize] = useState<{ width: number; height: number } | null>(null);

  const [xMin, setXMin] = useState<number | null>(null);
  const [xMax, setXMax] = useState<number | null>(null);
  const [yMin, setYMin] = useState<number | null>(null);
  const [yMax, setYMax] = useState<number | null>(null);
  const [manualXMinInput, setManualXMinInput] = useState("0");
  const [manualXMaxInput, setManualXMaxInput] = useState("1");
  const [manualYMinInput, setManualYMinInput] = useState("0");
  const [manualYMaxInput, setManualYMaxInput] = useState("1");
  const [xRangeError, setXRangeError] = useState<string | null>(null);
  const [yRangeError, setYRangeError] = useState<string | null>(null);
  const [axisExpanded, setAxisExpanded] = useState(true);

  useLayoutEffect(() => {
    const el = canvasRef.current;
    if (!el) return;

    const observer = new ResizeObserver((entries) => {
      const entry = entries[0];
      if (!entry) return;
      const nextWidth = Math.floor(entry.contentRect.width);
      const nextHeight = Math.floor(entry.contentRect.height);
      setContainerSize((prev) => {
        if (prev && prev.width === nextWidth && prev.height === nextHeight) {
          return prev;
        }
        return {
          width: nextWidth,
          height: nextHeight,
        };
      });
    });

    observer.observe(el);
    return () => observer.disconnect();
  }, []);

  const [autoXMin, autoXMax] = useMemo<[number, number]>(() => {
    const minFromPayload = Number.isFinite(data.energy_range?.[0]) ? data.energy_range[0] : Math.min(...data.energies);
    const maxFromPayload = Number.isFinite(data.energy_range?.[1]) ? data.energy_range[1] : Math.max(...data.energies);
    if (!Number.isFinite(minFromPayload) || !Number.isFinite(maxFromPayload)) {
      return [-10, 10];
    }
    if (Math.abs(maxFromPayload - minFromPayload) < 1e-9) {
      return [minFromPayload - 1, maxFromPayload + 1];
    }
    return [minFromPayload, maxFromPayload];
  }, [data.energy_range, data.energies]);

  const [autoYMin, autoYMax] = useMemo<[number, number]>(() => {
    const maxFromPayload = Number(data.max_dos);
    const maxFromSeries = data.dos.reduce((acc, value) => Math.max(acc, value), 0);
    const max = Math.max(maxFromPayload, maxFromSeries);
    if (!Number.isFinite(max) || max <= 0) {
      return [0, 1];
    }
    return [0, max];
  }, [data.dos, data.max_dos]);

  useEffect(() => {
    if (xMin === null && xMax === null) {
      setManualXMinInput(formatAxisInputValue(autoXMin));
      setManualXMaxInput(formatAxisInputValue(autoXMax));
    }
  }, [autoXMax, autoXMin, xMax, xMin]);

  useEffect(() => {
    if (yMin === null && yMax === null) {
      setManualYMinInput(formatAxisInputValue(autoYMin));
      setManualYMaxInput(formatAxisInputValue(autoYMax));
    }
  }, [autoYMax, autoYMin, yMax, yMin]);

  const resolvedXMin = xMin ?? autoXMin;
  const resolvedXMax = xMax ?? autoXMax;
  const resolvedYMin = yMin ?? autoYMin;
  const resolvedYMax = yMax ?? autoYMax;

  const resolvedWidth = containerSize ? Math.max(1, containerSize.width) : Math.max(320, width);
  const resolvedHeight = containerSize ? Math.max(1, containerSize.height) : Math.max(260, height);

  const margins = { top: 28, right: 26, bottom: 56, left: 64 };
  const innerWidth = Math.max(1, resolvedWidth - margins.left - margins.right);
  const innerHeight = Math.max(1, resolvedHeight - margins.top - margins.bottom);

  const xSpan = Math.max(1e-9, resolvedXMax - resolvedXMin);
  const ySpan = Math.max(1e-9, resolvedYMax - resolvedYMin);
  const toX = (energy: number) => ((energy - resolvedXMin) / xSpan) * innerWidth;
  const toY = (value: number) => innerHeight - ((value - resolvedYMin) / ySpan) * innerHeight;

  const pathD = useMemo(() => {
    const count = Math.min(data.energies.length, data.dos.length);
    if (count === 0) return "";

    let path = `M ${toX(data.energies[0]).toFixed(3)} ${toY(data.dos[0]).toFixed(3)}`;
    for (let i = 1; i < count; i += 1) {
      path += ` L ${toX(data.energies[i]).toFixed(3)} ${toY(data.dos[i]).toFixed(3)}`;
    }
    return path;
  }, [data.dos, data.energies, innerHeight, innerWidth, resolvedXMax, resolvedXMin, resolvedYMax, resolvedYMin]);

  const fermiX = useMemo(() => {
    if (data.fermi_energy == null || !Number.isFinite(data.fermi_energy)) return null;
    if (data.fermi_energy < resolvedXMin || data.fermi_energy > resolvedXMax) return null;
    return toX(data.fermi_energy);
  }, [data.fermi_energy, resolvedXMax, resolvedXMin, innerWidth]);

  const xTicks = useMemo(() => {
    return getTickValues(resolvedXMin, resolvedXMax, 5);
  }, [resolvedXMax, resolvedXMin]);

  const yTicks = useMemo(() => {
    return getTickValues(resolvedYMin, resolvedYMax, 4);
  }, [resolvedYMax, resolvedYMin]);

  const xTickStep = useMemo(() => {
    if (xTicks.length < 2) return 1;
    return Math.abs(xTicks[1] - xTicks[0]);
  }, [xTicks]);

  const yTickStep = useMemo(() => {
    if (yTicks.length < 2) return 1;
    return Math.abs(yTicks[1] - yTicks[0]);
  }, [yTicks]);

  const xTickDecimals = useMemo(() => getTickDecimals(xTickStep), [xTickStep]);
  const yTickDecimals = useMemo(() => getTickDecimals(yTickStep), [yTickStep]);

  const applyXRange = useCallback(() => {
    const parsedMin = Number.parseFloat(manualXMinInput.trim());
    const parsedMax = Number.parseFloat(manualXMaxInput.trim());
    if (!Number.isFinite(parsedMin) || !Number.isFinite(parsedMax)) {
      setXRangeError("Enter valid numeric values for X min and X max.");
      return;
    }
    if (parsedMax <= parsedMin) {
      setXRangeError("X max must be greater than X min.");
      return;
    }
    setXMin(parsedMin);
    setXMax(parsedMax);
    setXRangeError(null);
  }, [manualXMaxInput, manualXMinInput]);

  const applyYRange = useCallback(() => {
    const parsedMin = Number.parseFloat(manualYMinInput.trim());
    const parsedMax = Number.parseFloat(manualYMaxInput.trim());
    if (!Number.isFinite(parsedMin) || !Number.isFinite(parsedMax)) {
      setYRangeError("Enter valid numeric values for Y min and Y max.");
      return;
    }
    if (parsedMax <= parsedMin) {
      setYRangeError("Y max must be greater than Y min.");
      return;
    }
    setYMin(parsedMin);
    setYMax(parsedMax);
    setYRangeError(null);
  }, [manualYMaxInput, manualYMinInput]);

  const handleRangeKeyDown = useCallback(
    (event: React.KeyboardEvent<HTMLInputElement>, applyRange: () => void) => {
      if (event.key !== "Enter") return;
      event.preventDefault();
      applyRange();
    },
    [],
  );

  const resetView = useCallback(() => {
    setXMin(null);
    setXMax(null);
    setYMin(null);
    setYMax(null);
    setXRangeError(null);
    setYRangeError(null);
    setManualXMinInput(formatAxisInputValue(autoXMin));
    setManualXMaxInput(formatAxisInputValue(autoXMax));
    setManualYMinInput(formatAxisInputValue(autoYMin));
    setManualYMaxInput(formatAxisInputValue(autoYMax));
  }, [autoXMax, autoXMin, autoYMax, autoYMin]);

  return (
    <div className="electronic-dos-viewer-layout">
      <div className="electronic-dos-plot-panel">
        <div className="electronic-dos-plot">
          <div className="electronic-dos-plot-canvas" ref={canvasRef}>
            <svg width={resolvedWidth} height={resolvedHeight} viewBox={`0 0 ${resolvedWidth} ${resolvedHeight}`}>
              <g transform={`translate(${margins.left},${margins.top})`}>
                <rect x={0} y={0} width={innerWidth} height={innerHeight} fill="#ffffff" />

                {xTicks.map((tick) => (
                  <g key={`x-${tick.toFixed(4)}`} transform={`translate(${toX(tick)},0)`}>
                    <line y1={0} y2={innerHeight} stroke="#e2e8f0" strokeWidth={1} />
                    <text y={innerHeight + 22} textAnchor="middle" fontSize={12} fill="#4a5568">
                      {tick.toFixed(xTickDecimals)}
                    </text>
                  </g>
                ))}

                {yTicks.map((tick) => (
                  <g key={`y-${tick.toFixed(4)}`} transform={`translate(0,${toY(tick)})`}>
                    <line x1={0} x2={innerWidth} stroke="#edf2f7" strokeWidth={1} />
                    <text x={-10} y={4} textAnchor="end" fontSize={12} fill="#4a5568">
                      {tick.toFixed(yTickDecimals)}
                    </text>
                  </g>
                ))}

                <line x1={0} y1={innerHeight} x2={innerWidth} y2={innerHeight} stroke="#2d3748" strokeWidth={1.4} />
                <line x1={0} y1={0} x2={0} y2={innerHeight} stroke="#2d3748" strokeWidth={1.4} />

                {pathD && (
                  <path d={pathD} fill="none" stroke="#1f77b4" strokeWidth={1.8} />
                )}

                {fermiX !== null && (
                  <>
                    <line
                      x1={fermiX}
                      y1={0}
                      x2={fermiX}
                      y2={innerHeight}
                      stroke="#d53f8c"
                      strokeWidth={1.25}
                      strokeDasharray="4 3"
                    />
                    <text
                      x={Math.min(innerWidth - 6, fermiX + 4)}
                      y={14}
                      textAnchor="start"
                      fontSize={11}
                      fill="#97266d"
                    >
                      EF
                    </text>
                  </>
                )}

                <text
                  x={innerWidth / 2}
                  y={innerHeight + 44}
                  textAnchor="middle"
                  fontSize={13}
                  fill="#2d3748"
                >
                  Energy (eV)
                </text>
                <text
                  transform={`translate(${-46},${innerHeight / 2}) rotate(-90)`}
                  textAnchor="middle"
                  fontSize={13}
                  fill="#2d3748"
                >
                  DOS (states/eV)
                </text>
              </g>
            </svg>
          </div>
        </div>
      </div>

      <aside className="band-plot-sidebar electronic-dos-sidebar">
        <div className="band-plot-controls">
          <button onClick={resetView} className="band-plot-reset">
            Reset View
          </button>
          <span className="band-plot-hint">Enter ranges to rescale both axes.</span>
        </div>

        <div className="band-plot-control-panel">
          <div className="band-control-section">
            <button
              type="button"
              className="band-control-header"
              onClick={() => setAxisExpanded((prev) => !prev)}
            >
              <span className={`collapse-icon ${axisExpanded ? "expanded" : ""}`}>▶</span>
              Axis
            </button>
            {axisExpanded && (
              <div className="band-control-grid">
                <div className="band-control-row">
                  <label>Energy Range (eV)</label>
                  <div className="band-control-range-inputs">
                    <input
                      type="number"
                      step="any"
                      value={manualXMinInput}
                      onChange={(event) => {
                        setManualXMinInput(event.target.value);
                        setXRangeError(null);
                      }}
                      onKeyDown={(event) => handleRangeKeyDown(event, applyXRange)}
                      aria-label="Energy minimum"
                    />
                    <span className="band-control-range-separator">to</span>
                    <input
                      type="number"
                      step="any"
                      value={manualXMaxInput}
                      onChange={(event) => {
                        setManualXMaxInput(event.target.value);
                        setXRangeError(null);
                      }}
                      onKeyDown={(event) => handleRangeKeyDown(event, applyXRange)}
                      aria-label="Energy maximum"
                    />
                  </div>
                  <button
                    type="button"
                    className="band-control-apply"
                    onClick={applyXRange}
                  >
                    Apply X Range
                  </button>
                  {xRangeError && (
                    <span className="band-control-range-error">{xRangeError}</span>
                  )}
                </div>

                <div className="band-control-row">
                  <label>DOS Range (states/eV)</label>
                  <div className="band-control-range-inputs">
                    <input
                      type="number"
                      step="any"
                      value={manualYMinInput}
                      onChange={(event) => {
                        setManualYMinInput(event.target.value);
                        setYRangeError(null);
                      }}
                      onKeyDown={(event) => handleRangeKeyDown(event, applyYRange)}
                      aria-label="DOS minimum"
                    />
                    <span className="band-control-range-separator">to</span>
                    <input
                      type="number"
                      step="any"
                      value={manualYMaxInput}
                      onChange={(event) => {
                        setManualYMaxInput(event.target.value);
                        setYRangeError(null);
                      }}
                      onKeyDown={(event) => handleRangeKeyDown(event, applyYRange)}
                      aria-label="DOS maximum"
                    />
                  </div>
                  <button
                    type="button"
                    className="band-control-apply"
                    onClick={applyYRange}
                  >
                    Apply Y Range
                  </button>
                  {yRangeError && (
                    <span className="band-control-range-error">{yRangeError}</span>
                  )}
                </div>

                <div className="band-control-note">
                  X scales the energy axis, Y scales the DOS axis. Reset View restores automatic scaling.
                </div>
              </div>
            )}
          </div>
        </div>
      </aside>
    </div>
  );
}
