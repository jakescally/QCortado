export const TOOLTIP_OPEN_DELAY_MS = 750;

export function getTooltipOpenDelay(hasCustomTrigger: boolean): number {
  return hasCustomTrigger ? TOOLTIP_OPEN_DELAY_MS : 0;
}
