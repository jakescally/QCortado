export const MAX_LIVE_OUTPUT_LINES = 1600;
export const MAX_LIVE_OUTPUT_CHARS = 200_000;

export interface VisibleOutputWindow {
  lines: string[];
  text: string;
  totalLineCount: number;
}

export function buildVisibleOutputWindow(
  lines: string[],
  maxLines = MAX_LIVE_OUTPUT_LINES,
  maxChars = MAX_LIVE_OUTPUT_CHARS,
): VisibleOutputWindow {
  const totalLineCount = lines.length;
  if (totalLineCount === 0) {
    return {
      lines: [],
      text: "",
      totalLineCount: 0,
    };
  }

  let visibleLines = lines.slice(Math.max(0, totalLineCount - maxLines));
  let visibleText = visibleLines.join("\n");
  while (visibleLines.length > 1 && visibleText.length > maxChars) {
    visibleLines = visibleLines.slice(1);
    visibleText = visibleLines.join("\n");
  }

  return {
    lines: visibleLines,
    text: visibleText,
    totalLineCount,
  };
}

export function countVisibleOutputLines(text: string): number {
  if (!text) {
    return 0;
  }
  return text.split("\n").length - Number(text.endsWith("\n"));
}
