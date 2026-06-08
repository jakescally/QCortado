import { useEffect, useRef } from "react";

const NARROW_DRAWER_BREAKPOINT = 900;

export function useOverlayDrawer<T extends HTMLElement>(
  isOpen: boolean,
  onClose: () => void,
) {
  const drawerRef = useRef<T | null>(null);

  useEffect(() => {
    if (!isOpen) return;

    const drawer = drawerRef.current;
    const firstFocusable = drawer?.querySelector<HTMLElement>(
      "button:not(:disabled), a[href], input:not(:disabled), select:not(:disabled), textarea:not(:disabled), [tabindex]:not([tabindex='-1'])",
    );
    firstFocusable?.focus();

    function handlePointerDown(event: MouseEvent) {
      if (drawerRef.current && !drawerRef.current.contains(event.target as Node)) {
        onClose();
      }
    }

    function handleKeyDown(event: KeyboardEvent) {
      if (event.key === "Escape") {
        event.preventDefault();
        onClose();
        return;
      }
      if (event.key !== "Tab" || window.innerWidth > NARROW_DRAWER_BREAKPOINT) return;

      const focusable = Array.from(drawerRef.current?.querySelectorAll<HTMLElement>(
        "button:not(:disabled), a[href], input:not(:disabled), select:not(:disabled), textarea:not(:disabled), [tabindex]:not([tabindex='-1'])",
      ) ?? []);
      if (focusable.length === 0) return;
      const first = focusable[0];
      const last = focusable[focusable.length - 1];
      if (event.shiftKey && document.activeElement === first) {
        event.preventDefault();
        last.focus();
      } else if (!event.shiftKey && document.activeElement === last) {
        event.preventDefault();
        first.focus();
      }
    }

    document.addEventListener("mousedown", handlePointerDown);
    document.addEventListener("keydown", handleKeyDown);
    return () => {
      document.removeEventListener("mousedown", handlePointerDown);
      document.removeEventListener("keydown", handleKeyDown);
    };
  }, [isOpen, onClose]);

  return drawerRef;
}
