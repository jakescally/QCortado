import { useEffect, useState } from "react";
import { invoke } from "@tauri-apps/api/core";
import "./App.css";

function App() {
  const [qePath, setQePath] = useState<string | null>(null);
  const [availableExecutables, setAvailableExecutables] = useState<string[]>([]);
  const [status, setStatus] = useState<string>("Not configured");
  const [error, setError] = useState<string | null>(null);

  // Check for existing QE configuration on startup
  useEffect(() => {
    checkQEPath();
  }, []);

  async function checkQEPath() {
    try {
      const path = await invoke<string | null>("get_qe_path");
      if (path) {
        setQePath(path);
        await loadExecutables();
        setStatus("Ready");
      }
    } catch (e) {
      console.log("No QE path configured yet");
    }
  }

  async function loadExecutables() {
    try {
      const exes = await invoke<string[]>("check_qe_executables");
      setAvailableExecutables(exes);
      setError(null);
    } catch (e) {
      setError(String(e));
    }
  }

  async function selectQEPath() {
    try {
      // For now, use a prompt - in full version we'd use Tauri's dialog
      const path = prompt("Enter path to QE bin directory:", "/path/to/qe-7.5/bin");
      if (path) {
        await invoke("set_qe_path", { path });
        setQePath(path);
        await loadExecutables();
        setStatus("Ready");
        setError(null);
      }
    } catch (e) {
      setError(String(e));
    }
  }

  return (
    <main className="container">
      <header className="header">
        <h1>QCortado</h1>
        <p className="subtitle">A Modern Interface for Quantum ESPRESSO</p>
      </header>

      <section className="config-section">
        <h2>Configuration</h2>
        <div className="config-row">
          <label>QE Installation:</label>
          {qePath ? (
            <span className="path">{qePath}</span>
          ) : (
            <span className="not-set">Not configured</span>
          )}
          <button onClick={selectQEPath}>
            {qePath ? "Change" : "Configure"}
          </button>
        </div>

        <div className="status-row">
          <label>Status:</label>
          <span className={`status ${status === "Ready" ? "ready" : "pending"}`}>
            {status}
          </span>
        </div>

        {error && <div className="error">{error}</div>}

        {availableExecutables.length > 0 && (
          <div className="executables">
            <label>Available programs:</label>
            <div className="exe-list">
              {availableExecutables.map((exe) => (
                <span key={exe} className="exe-tag">{exe}</span>
              ))}
            </div>
          </div>
        )}
      </section>

      {qePath && (
        <section className="actions-section">
          <h2>Quick Actions</h2>
          <div className="action-grid">
            <button className="action-btn" disabled>
              <span className="action-icon">+</span>
              <span className="action-label">New Project</span>
              <span className="action-hint">Coming soon</span>
            </button>
            <button className="action-btn" disabled>
              <span className="action-icon">SCF</span>
              <span className="action-label">SCF Calculation</span>
              <span className="action-hint">Coming soon</span>
            </button>
            <button className="action-btn" disabled>
              <span className="action-icon">Band</span>
              <span className="action-label">Band Structure</span>
              <span className="action-hint">Coming soon</span>
            </button>
            <button className="action-btn" disabled>
              <span className="action-icon">DOS</span>
              <span className="action-label">Density of States</span>
              <span className="action-hint">Coming soon</span>
            </button>
          </div>
        </section>
      )}

      <footer className="footer">
        <p>QCortado v0.1.0 - Quantum ESPRESSO 7.5 Interface</p>
      </footer>
    </main>
  );
}

export default App;
