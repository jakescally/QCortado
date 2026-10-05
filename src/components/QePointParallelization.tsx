import { InfoTooltip } from "./InfoTooltip";
import type { QeAutoPlan } from "../lib/qeAutoParallelization";
import { balancedPoolCounts, PointEstimateState, QePointEstimates } from "../lib/qePointEstimate";

export function PointEstimateSummary({ estimate, label }: { estimate: PointEstimateState; label: string }) {
  return <p className="param-hint qe-point-summary" aria-live="polite">
    {estimate.result ? `${estimate.explicit ? "Explicit" : estimate.conservative ? "Conservative upper bound" : "Expected irreducible"} ${label}: ${estimate.result.irreduciblePoints} of ${estimate.result.fullPoints}`
      : estimate.loading ? `Estimating irreducible ${label}…` : `Irreducible ${label}: unavailable`}
    <InfoTooltip text={estimate.error ?? (estimate.explicit ? "Band paths use explicitly sampled points, including endpoints and path breaks. They are not reduced as a uniform reciprocal mesh. These sampled points determine pool load balance." : "Computed from the calculation cell and mesh offsets using spglib. This is a preflight estimate; QE's symmetry tolerance, magnetic symmetries and mesh handling can change the final count. For noncollinear/SOC configurations the preview conservatively disables symmetry reduction. Collinear spin channels do not double the displayed number of distinct k-points.")} />
  </p>;
}

interface Props {
  ranks: number;
  pools: number;
  images?: number;
  phonon?: boolean;
  estimates?: QePointEstimates;
  autoPlan?: QeAutoPlan;
  disabled?: boolean;
  onChange: (pools: number, images: number) => void;
}

export function QePointParallelization({ ranks, pools, images = 1, phonon = false, estimates, autoPlan, disabled = false, onChange }: Props) {
  const valid = Number.isInteger(pools) && Number.isInteger(images) && pools > 0 && images > 0 && ranks % (pools * images) === 0;
  const balancedEstimate = phonon ? estimates?.q?.result : estimates?.k.result;
  const choices = balancedEstimate ? balancedPoolCounts(ranks, phonon ? pools : images, balancedEstimate.irreduciblePoints) : [];
  const selectedCount = phonon ? images : pools;
  return <section className="option-section qe-point-parallelization">
    <h4 className="qe-point-label">QE point parallelization <InfoTooltip text="Split MPI tasks into k-point pools and, for phonons, images that distribute q-points and irreducible representations. A count of 1 disables that level. The existing plane-wave and FFT parallelization continues inside each pool." /></h4>
    {estimates && <PointEstimateSummary estimate={estimates.k} label={estimates.kLabel ?? "k-points"} />}
    {estimates?.q && <PointEstimateSummary estimate={estimates.q} label="q-points" />}
    {autoPlan && <div className="qe-auto-breakdown" aria-live="polite">
      <span>{autoPlan.ranks} MPI ranks × {autoPlan.threads} CPU{autoPlan.threads === 1 ? "" : "s"} / rank</span>
      {phonon && <span>{images} phonon image{images === 1 ? "" : "s"} · {autoPlan.ranksPerImage} ranks / image</span>}
      <span>{pools} electronic pool{pools === 1 ? "" : "s"}{phonon ? " / image" : ""} · {autoPlan.ranksPerPool} ranks / pool</span>
      {autoPlan.notes.length > 0 && <span className="qe-point-label">Auto details <InfoTooltip text={autoPlan.notes.join(" ")} /></span>}
    </div>}
    {!autoPlan && <div className="hpc-grid qe-point-fields">
      {phonon && <label>
        <span className="qe-point-label">Phonon images <InfoTooltip text="QE calls these phonon images (-nimage). They distribute q-points and irreducible representations, so their useful count can exceed the q-point count. Images require shared scratch and automatically collect dynamical matrices before post-processing. Use a compatible source SCF pool/rank layout; keep the image and k-point pool layout unchanged for recovery." /></span>
        <input type="number" min={1} max={ranks} step={1} value={images} disabled={disabled}
          onChange={(event) => onChange(pools, Math.max(1, Number(event.target.value) || 1))} />
      </label>}
      <label>
        <span className="qe-point-label">{phonon ? "K-point pools / image" : "K-point pools"} <InfoTooltip text="QE's -nk setting distributes electronic k-points within each image. Use no more pools than irreducible k-points. Balanced choices divide both MPI tasks per image and the expected k-point count. Keep each pool within a node where practical, and benchmark the tradeoff with FFT performance. For phonons the k-point symmetry changes at each q, so the source SCF count is guidance rather than a strict cap." /></span>
        <input type="number" min={1} max={ranks} step={1} value={pools} disabled={disabled}
          onChange={(event) => onChange(Math.max(1, Number(event.target.value) || 1), images)} />
      </label>
    </div>}
    {!autoPlan && choices.length > 0 && <div className="qe-point-balanced">
      <span className="qe-point-label">{phonon ? "Balanced image counts" : "Balanced pool counts"} <InfoTooltip text={phonon
        ? "These image counts divide both the expected irreducible q-point count and the MPI tasks available after allowing for k-point pools. They balance q-points; representation costs can differ and images can also distribute individual representations."
        : "These counts divide both the expected distinct k-point count and the available MPI tasks per image. They describe load balance, not a guarantee of the fastest layout."} /></span>
      <div className="qe-point-counts">
        {choices.map((count) => <button key={count} type="button" className="btn btn-secondary qe-point-count" disabled={disabled}
          aria-pressed={selectedCount === count}
          onClick={() => phonon ? onChange(pools, count) : onChange(count, images)}>{count}</button>)}
      </div>
    </div>}
    {!autoPlan && <p>{valid ? `${ranks / (pools * images)} MPI tasks per pool; ${ranks / images} per image.` : "MPI tasks must be divisible by pools × images. Adjust the counts before running."}</p>}
  </section>;
}
