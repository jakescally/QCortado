import { useCallback, useEffect, useRef, useState } from "react";

interface LiveOutputPanelProps {
  title: string;
  output: string;
  placeholder?: string;
  totalLineCount?: number;
  visibleLineCount?: number;
  panelClassName?: string;
  outputClassName?: string;
}

export function LiveOutputPanel({
  title,
  output,
  placeholder = "Waiting for output...",
  totalLineCount = 0,
  visibleLineCount = 0,
  panelClassName = "output-panel",
  outputClassName = "output-text",
}: LiveOutputPanelProps) {
  const outputRef = useRef<HTMLPreElement>(null);
  const [isFollowing, setIsFollowing] = useState(true);
  const hiddenLineCount = Math.max(0, totalLineCount - visibleLineCount);

  const scrollToBottom = useCallback(() => {
    const el = outputRef.current;
    if (!el) return;
    el.scrollTop = el.scrollHeight;
  }, []);

  const handleScroll = useCallback(() => {
    if (!isFollowing) return;
    const el = outputRef.current;
    if (!el) return;
    const distanceToBottom = el.scrollHeight - el.scrollTop - el.clientHeight;
    if (distanceToBottom > 16) {
      setIsFollowing(false);
    }
  }, [isFollowing]);

  const handleFollowModeChange = useCallback((nextFollow: boolean) => {
    setIsFollowing(nextFollow);
    if (nextFollow) {
      requestAnimationFrame(() => {
        scrollToBottom();
      });
    }
  }, [scrollToBottom]);

  useEffect(() => {
    if (!isFollowing) return;
    scrollToBottom();
  }, [isFollowing, output, scrollToBottom]);

  useEffect(() => {
    if (!output && totalLineCount === 0) {
      setIsFollowing(true);
    }
  }, [output, totalLineCount]);

  return (
    <div className={panelClassName}>
      <div className="output-panel-header">
        <h3>{title}</h3>
        <div className="output-panel-tools">
          {hiddenLineCount > 0 && (
            <span
              className="output-panel-note"
              title={`${hiddenLineCount.toLocaleString()} earlier lines are hidden to keep the live view responsive.`}
            >
              Showing last {visibleLineCount.toLocaleString()} lines
            </span>
          )}
          <div className="output-follow-toggle" role="group" aria-label="Output scrolling mode">
            <button
              type="button"
              className={isFollowing ? "active" : ""}
              aria-pressed={isFollowing}
              onClick={() => handleFollowModeChange(true)}
            >
              Lock
            </button>
            <button
              type="button"
              className={!isFollowing ? "active" : ""}
              aria-pressed={!isFollowing}
              onClick={() => handleFollowModeChange(false)}
            >
              Free
            </button>
          </div>
        </div>
      </div>
      <pre ref={outputRef} className={outputClassName} onScroll={handleScroll}>
        {output || placeholder}
      </pre>
    </div>
  );
}
