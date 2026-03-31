import { useEffect } from "react";

export function useViewportScrollLock(locked: boolean) {
  useEffect(() => {
    const root = document.getElementById("root");
    if (!root) {
      return;
    }

    const previousOverflowY = root.style.overflowY;
    const previousOverscrollBehavior = root.style.overscrollBehavior;

    if (locked) {
      root.style.overflowY = "hidden";
      root.style.overscrollBehavior = "contain";
    }

    return () => {
      root.style.overflowY = previousOverflowY;
      root.style.overscrollBehavior = previousOverscrollBehavior;
    };
  }, [locked]);
}
