import { useLayoutEffect, useState, type HTMLAttributes, type ReactNode } from "react";
import { createPortal } from "react-dom";

const APP_HEADER_SLOT_ID = "app-dynamic-header-slot";

interface AppHeaderPortalProps extends HTMLAttributes<HTMLDivElement> {
  children: ReactNode;
}

export function AppHeaderPortal({ children, className = "", ...props }: AppHeaderPortalProps) {
  const [target, setTarget] = useState<HTMLElement | null>(null);

  useLayoutEffect(() => {
    setTarget(document.getElementById(APP_HEADER_SLOT_ID));
  }, []);

  if (!target) return null;

  return createPortal(
    <div className={`app-dynamic-header ${className}`.trim()} {...props}>
      {children}
    </div>,
    target,
  );
}
