interface InfoTooltipProps {
  text: string;
  className?: string;
}

export function InfoTooltip({ text, className }: InfoTooltipProps) {
  return (
    <span className={className ? `tooltip-container ${className}` : "tooltip-container"}>
      <span className="tooltip-icon" aria-hidden="true">?</span>
      <span className="tooltip-text" role="tooltip">{text}</span>
    </span>
  );
}
