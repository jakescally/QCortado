import type { ReactNode } from "react";
import { InfoTooltip } from "../InfoTooltip";

interface Wien2kFieldLabelProps {
  children: ReactNode;
  tooltip: string;
}

export function Wien2kFieldLabel({ children, tooltip }: Wien2kFieldLabelProps) {
  return (
    <span className="wien2k-field-label-row">
      {children}
      <InfoTooltip text={tooltip} />
    </span>
  );
}
