import { useCallback, useEffect, useId, useLayoutEffect, useRef, useState } from "react";
import type { ReactNode } from "react";
import { createPortal } from "react-dom";
import { getTooltipOpenDelay } from "./tooltipTiming";

interface InfoTooltipProps {
  text: string;
  className?: string;
  children?: ReactNode;
}

type TooltipPlacement = "top" | "bottom";

interface TooltipPosition {
  top: number;
  left: number;
  arrowLeft: number;
  placement: TooltipPlacement;
}

const TOOLTIP_GAP = 8;
const VIEWPORT_PADDING = 10;
const CLOSE_DELAY_MS = 120;

export function InfoTooltip({ text, className, children }: InfoTooltipProps) {
  const triggerRef = useRef<HTMLSpanElement | null>(null);
  const bubbleRef = useRef<HTMLDivElement | null>(null);
  const openTimerRef = useRef<number | null>(null);
  const closeTimerRef = useRef<number | null>(null);
  const tooltipId = useId();

  const [open, setOpen] = useState(false);
  const [position, setPosition] = useState<TooltipPosition | null>(null);

  const clearCloseTimer = useCallback(() => {
    if (closeTimerRef.current != null) {
      window.clearTimeout(closeTimerRef.current);
      closeTimerRef.current = null;
    }
  }, []);

  const clearOpenTimer = useCallback(() => {
    if (openTimerRef.current != null) {
      window.clearTimeout(openTimerRef.current);
      openTimerRef.current = null;
    }
  }, []);

  const showTooltip = useCallback(() => {
    clearOpenTimer();
    clearCloseTimer();
    setOpen(true);
  }, [clearCloseTimer, clearOpenTimer]);

  const showTooltipOnHover = useCallback(() => {
    const openDelay = getTooltipOpenDelay(children != null);
    if (openDelay === 0) {
      showTooltip();
      return;
    }

    clearOpenTimer();
    clearCloseTimer();
    openTimerRef.current = window.setTimeout(() => {
      openTimerRef.current = null;
      setOpen(true);
    }, openDelay);
  }, [children, clearCloseTimer, clearOpenTimer, showTooltip]);

  const hideTooltip = useCallback(() => {
    clearOpenTimer();
    clearCloseTimer();
    closeTimerRef.current = window.setTimeout(() => {
      setOpen(false);
    }, CLOSE_DELAY_MS);
  }, [clearCloseTimer, clearOpenTimer]);

  const updatePosition = useCallback(() => {
    const trigger = triggerRef.current;
    const bubble = bubbleRef.current;
    if (!trigger || !bubble) return;

    const triggerRect = trigger.getBoundingClientRect();
    const bubbleRect = bubble.getBoundingClientRect();
    const viewportWidth = window.innerWidth;
    const viewportHeight = window.innerHeight;

    let left = triggerRect.left + (triggerRect.width / 2) - (bubbleRect.width / 2);
    const minLeft = VIEWPORT_PADDING;
    const maxLeft = viewportWidth - bubbleRect.width - VIEWPORT_PADDING;
    if (maxLeft < minLeft) {
      left = minLeft;
    } else {
      left = Math.min(maxLeft, Math.max(minLeft, left));
    }

    let placement: TooltipPlacement = "top";
    let top = triggerRect.top - bubbleRect.height - TOOLTIP_GAP;
    if (top < VIEWPORT_PADDING) {
      placement = "bottom";
      top = triggerRect.bottom + TOOLTIP_GAP;
    }
    const maxTop = viewportHeight - bubbleRect.height - VIEWPORT_PADDING;
    top = Math.min(maxTop, Math.max(VIEWPORT_PADDING, top));

    const triggerCenterX = triggerRect.left + (triggerRect.width / 2);
    const arrowLeft = Math.max(
      12,
      Math.min(bubbleRect.width - 12, triggerCenterX - left),
    );

    setPosition({ top, left, arrowLeft, placement });
  }, []);

  useLayoutEffect(() => {
    if (!open) return;
    updatePosition();
  }, [open, text, updatePosition]);

  useEffect(() => {
    if (!open) return;
    const handleViewportChange = () => updatePosition();
    window.addEventListener("resize", handleViewportChange);
    window.addEventListener("scroll", handleViewportChange, true);
    return () => {
      window.removeEventListener("resize", handleViewportChange);
      window.removeEventListener("scroll", handleViewportChange, true);
    };
  }, [open, updatePosition]);

  useEffect(() => {
    return () => {
      clearOpenTimer();
      clearCloseTimer();
    };
  }, [clearCloseTimer, clearOpenTimer]);

  const bubble = open ? createPortal(
    <div
      ref={bubbleRef}
      id={tooltipId}
      role="tooltip"
      className={className ? `tooltip-layer ${className}` : "tooltip-layer"}
      data-placement={position?.placement ?? "top"}
      style={position
        ? { top: `${position.top}px`, left: `${position.left}px` }
        : { top: "-9999px", left: "-9999px", visibility: "hidden" }}
      onMouseEnter={showTooltip}
      onMouseLeave={hideTooltip}
    >
      <div className="tooltip-layer-content">{text}</div>
      <span
        className="tooltip-arrow"
        style={position ? { left: `${position.arrowLeft}px` } : undefined}
        aria-hidden="true"
      />
    </div>,
    document.body,
  ) : null;

  return (
    <>
      <span
        ref={triggerRef}
        className={children ? "tooltip-trigger-wrap" : "tooltip-container"}
        onMouseEnter={showTooltipOnHover}
        onMouseLeave={hideTooltip}
        onFocus={showTooltip}
        onBlur={hideTooltip}
        aria-describedby={open ? tooltipId : undefined}
      >
        {children || <span className="tooltip-icon" aria-hidden="true">?</span>}
      </span>
      {bubble}
    </>
  );
}
