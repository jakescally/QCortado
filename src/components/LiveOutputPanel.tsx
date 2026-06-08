import { useCallback, useEffect, useLayoutEffect, useRef, useState } from "react";

interface LiveOutputPanelProps {
  title: string;
  output: string;
  placeholder?: string;
  totalLineCount?: number;
  visibleLineCount?: number;
  panelClassName?: string;
  outputClassName?: string;
  showTools?: boolean;
}

export function LiveOutputPanel({
  title,
  output,
  placeholder = "Waiting for output...",
  totalLineCount = 0,
  visibleLineCount = 0,
  panelClassName = "output-panel",
  outputClassName = "output-text",
  showTools = true,
}: LiveOutputPanelProps) {
  const outputRef = useRef<HTMLPreElement>(null);
  const scrollFrameRef = useRef<number | null>(null);
  const [isFollowing, setIsFollowing] = useState(true);
  const hiddenLineCount = Math.max(0, totalLineCount - visibleLineCount);

  const scrollToBottom = useCallback(() => {
    const el = outputRef.current;
    if (!el) return;
    el.scrollTop = el.scrollHeight;
  }, []);

  const scheduleScrollToBottom = useCallback(() => {
    if (scrollFrameRef.current !== null) {
      window.cancelAnimationFrame(scrollFrameRef.current);
    }
    scrollFrameRef.current = window.requestAnimationFrame(() => {
      scrollFrameRef.current = null;
      scrollToBottom();
    });
  }, [scrollToBottom]);

  const handleScroll = useCallback(() => {
    const el = outputRef.current;
    if (!el) return;
    const distanceToBottom = el.scrollHeight - el.scrollTop - el.clientHeight;
    setIsFollowing(distanceToBottom <= 16);
  }, []);

  const handleFollowModeChange = useCallback((nextFollow: boolean) => {
    setIsFollowing(nextFollow);
    if (nextFollow) {
      scheduleScrollToBottom();
    }
  }, [scheduleScrollToBottom]);

  useLayoutEffect(() => {
    if (!isFollowing) return;
    scheduleScrollToBottom();
  }, [isFollowing, output, totalLineCount, scheduleScrollToBottom]);

  useEffect(() => {
    if (!output && totalLineCount === 0) {
      setIsFollowing(true);
    }
  }, [output, totalLineCount]);

  useEffect(() => {
    return () => {
      if (scrollFrameRef.current !== null) {
        window.cancelAnimationFrame(scrollFrameRef.current);
      }
    };
  }, []);

  return (
    <div className={panelClassName}>
      <div className="output-panel-header">
        <h3>{title}</h3>
        {showTools && (
          <div className="output-panel-tools">
            {hiddenLineCount > 0 && (
              <span
                className="output-panel-note"
                title={`${hiddenLineCount.toLocaleString()} earlier lines are hidden to keep the live view responsive.`}
              >
                Newest {visibleLineCount.toLocaleString()} of {totalLineCount.toLocaleString()} lines
              </span>
            )}
            <div className="output-follow-toggle" role="group" aria-label="Output scrolling mode">
              <button
                type="button"
                className={isFollowing ? "active" : ""}
                aria-pressed={isFollowing}
                onClick={() => handleFollowModeChange(true)}
              >
                Follow
              </button>
              <button
                type="button"
                className={!isFollowing ? "active" : ""}
                aria-pressed={!isFollowing}
                onClick={() => handleFollowModeChange(false)}
              >
                Pause
              </button>
            </div>
          </div>
        )}
      </div>
      <pre ref={outputRef} className={outputClassName} onScroll={handleScroll}>
        {output || placeholder}
      </pre>
    </div>
  );
}
