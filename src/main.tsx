import React from "react";
import ReactDOM from "react-dom/client";
import App from "./App";
import { validateBrillouinZoneFixtures } from "./lib/brillouinZoneFixtures";

if (import.meta.env.DEV) {
  validateBrillouinZoneFixtures();
}

const userAgent = navigator.userAgent.toLowerCase();
const platform = userAgent.includes("linux")
  ? "linux"
  : userAgent.includes("mac")
    ? "macos"
    : userAgent.includes("win")
      ? "windows"
      : "other";
document.documentElement.setAttribute("data-platform", platform);

ReactDOM.createRoot(document.getElementById("root") as HTMLElement).render(
  <React.StrictMode>
    <App />
  </React.StrictMode>,
);
