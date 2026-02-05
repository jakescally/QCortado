import { useEffect, useState } from "react";
import { invoke } from "@tauri-apps/api/core";
import { open } from "@tauri-apps/plugin-dialog";
import "./App.css";
import { SCFWizard } from "./components/SCFWizard";
import { ProjectBrowser } from "./components/ProjectBrowser";
import { CreateProjectDialog } from "./components/CreateProjectDialog";

interface ProjectSummary {
  id: string;
  name: string;
  description: string | null;
  created_at: string;
  formula: string | null;
  calculation_count: number;
  last_activity: string;
}

type AppView = "home" | "scf-wizard" | "project-browser";

function App() {
  const [qePath, setQePath] = useState<string | null>(null);
  const [availableExecutables, setAvailableExecutables] = useState<string[]>([]);
  const [status, setStatus] = useState<string>("Not configured");
  const [error, setError] = useState<string | null>(null);
  const [currentView, setCurrentView] = useState<AppView>("home");
  const [projectCount, setProjectCount] = useState<number>(0);
  const [showCreateDialog, setShowCreateDialog] = useState(false);

  // Check for existing QE configuration on startup
  useEffect(() => {
    checkQEPath();
    loadProjectCount();
  }, []);

  async function loadProjectCount() {
    try {
      const projects = await invoke<ProjectSummary[]>("list_projects");
      setProjectCount(projects.length);
    } catch (e) {
      console.log("Failed to load project count:", e);
    }
  }

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
      const selected = await open({
        directory: true,
        multiple: false,
        title: "Select Quantum ESPRESSO bin directory",
      });

      if (selected && typeof selected === "string") {
        await invoke("set_qe_path", { path: selected });
        setQePath(selected);
        await loadExecutables();
        setStatus("Ready");
        setError(null);
      }
    } catch (e) {
      setError(String(e));
    }
  }

  if (currentView === "scf-wizard" && qePath) {
    return <SCFWizard qePath={qePath} onBack={() => setCurrentView("home")} />;
  }

  if (currentView === "project-browser") {
    return (
      <ProjectBrowser
        onBack={() => {
          setCurrentView("home");
          loadProjectCount(); // Refresh count in case projects were added/deleted
        }}
        onSelectProject={(projectId) => {
          // Placeholder: In the future, navigate to project dashboard
          console.log("Selected project:", projectId);
        }}
      />
    );
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
        <>
          <section className="actions-section">
            <h2>Projects</h2>
            <div className="action-grid">
              <button className="action-btn" onClick={() => setShowCreateDialog(true)}>
                <span className="action-icon">+</span>
                <span className="action-label">New Project</span>
                <span className="action-hint">Create a project</span>
              </button>
              <button
                className="action-btn"
                onClick={() => setCurrentView("project-browser")}
                disabled={projectCount === 0}
              >
                <span className="action-icon">{projectCount}</span>
                <span className="action-label">
                  {projectCount === 1 ? "View Project" : "View Projects"}
                </span>
                <span className="action-hint">
                  {projectCount === 0 ? "No projects yet" : "Browse & manage"}
                </span>
              </button>
            </div>
          </section>

          <section className="actions-section">
            <h2>Quick Calculations</h2>
            <div className="action-grid">
              <button className="action-btn" onClick={() => setCurrentView("scf-wizard")}>
                <span className="action-icon">SCF</span>
                <span className="action-label">SCF Calculation</span>
                <span className="action-hint">Import CIF & run</span>
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
        </>
      )}

      <footer className="footer">
        <p>QCortado v0.1.0 - Quantum ESPRESSO 7.5 Interface</p>
      </footer>

      <CreateProjectDialog
        isOpen={showCreateDialog}
        onClose={() => setShowCreateDialog(false)}
        onCreated={() => {
          setShowCreateDialog(false);
          loadProjectCount();
        }}
      />
    </main>
  );
}

export default App;
