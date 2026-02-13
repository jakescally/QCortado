import { useLayoutEffect, useMemo, useRef, useState } from "react";

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

export function ElectronicDOSPlot({
  data,
  width = 900,
  height = 500,
}: ElectronicDOSPlotProps) {
  const canvasRef = useRef<HTMLDivElement>(null);
  const [containerSize, setContainerSize] = useState<{ width: number; height: number } | null>(null);

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

  const resolvedWidth = containerSize ? Math.max(1, containerSize.width) : Math.max(320, width);
  const resolvedHeight = containerSize ? Math.max(1, containerSize.height) : Math.max(260, height);

  const margins = { top: 28, right: 26, bottom: 56, left: 64 };
  const innerWidth = Math.max(1, resolvedWidth - margins.left - margins.right);
  const innerHeight = Math.max(1, resolvedHeight - margins.top - margins.bottom);

  const [energyMin, energyMax] = useMemo<[number, number]>(() => {
    const min = Number.isFinite(data.energy_range?.[0]) ? data.energy_range[0] : Math.min(...data.energies);
    const max = Number.isFinite(data.energy_range?.[1]) ? data.energy_range[1] : Math.max(...data.energies);
    if (!Number.isFinite(min) || !Number.isFinite(max)) return [-10, 10];
    if (Math.abs(max - min) < 1e-9) return [min - 1, max + 1];
    return [min, max];
  }, [data.energy_range, data.energies]);

  const dosMax = useMemo(() => {
    const maxFromPayload = Number(data.max_dos);
    const maxFromSeries = data.dos.reduce((acc, value) => Math.max(acc, value), 0);
    const max = Math.max(maxFromPayload, maxFromSeries);
    return max > 0 ? max : 1;
  }, [data.dos, data.max_dos]);

  const energySpan = Math.max(1e-9, energyMax - energyMin);
  const toX = (energy: number) => ((energy - energyMin) / energySpan) * innerWidth;
  const toY = (value: number) => innerHeight - (Math.max(0, value) / dosMax) * innerHeight;

  const pathD = useMemo(() => {
    const count = Math.min(data.energies.length, data.dos.length);
    if (count === 0) return "";

    let path = `M ${toX(data.energies[0]).toFixed(3)} ${toY(data.dos[0]).toFixed(3)}`;
    for (let i = 1; i < count; i += 1) {
      path += ` L ${toX(data.energies[i]).toFixed(3)} ${toY(data.dos[i]).toFixed(3)}`;
    }
    return path;
  }, [data.dos, data.energies, innerHeight, innerWidth]);

  const fermiX = useMemo(() => {
    if (data.fermi_energy == null || !Number.isFinite(data.fermi_energy)) return null;
    if (data.fermi_energy < energyMin || data.fermi_energy > energyMax) return null;
    return toX(data.fermi_energy);
  }, [data.fermi_energy, energyMin, energyMax, innerWidth]);

  const xTicks = useMemo(() => {
    const ticks: number[] = [];
    const count = 5;
    for (let i = 0; i <= count; i += 1) {
      ticks.push(energyMin + (i / count) * (energyMax - energyMin));
    }
    return ticks;
  }, [energyMin, energyMax]);

  const yTicks = useMemo(() => {
    const ticks: number[] = [];
    const count = 4;
    for (let i = 0; i <= count; i += 1) {
      ticks.push((i / count) * dosMax);
    }
    return ticks;
  }, [dosMax]);

  return (
    <div className="electronic-dos-plot">
      <div className="electronic-dos-plot-canvas" ref={canvasRef}>
        <svg width={resolvedWidth} height={resolvedHeight} viewBox={`0 0 ${resolvedWidth} ${resolvedHeight}`}>
          <g transform={`translate(${margins.left},${margins.top})`}>
            <rect x={0} y={0} width={innerWidth} height={innerHeight} fill="#ffffff" />

            {xTicks.map((tick) => (
              <g key={`x-${tick.toFixed(4)}`} transform={`translate(${toX(tick)},0)`}>
                <line y1={0} y2={innerHeight} stroke="#e2e8f0" strokeWidth={1} />
                <text y={innerHeight + 22} textAnchor="middle" fontSize={12} fill="#4a5568">
                  {tick.toFixed(1)}
                </text>
              </g>
            ))}

            {yTicks.map((tick) => (
              <g key={`y-${tick.toFixed(4)}`} transform={`translate(0,${toY(tick)})`}>
                <line x1={0} x2={innerWidth} stroke="#edf2f7" strokeWidth={1} />
                <text x={-10} y={4} textAnchor="end" fontSize={12} fill="#4a5568">
                  {tick.toFixed(1)}
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
  );
}
