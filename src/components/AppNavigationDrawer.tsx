import {
  Activity,
  Database,
  FolderOpen,
  HardDrive,
  RotateCcw,
  Settings,
  X,
} from "lucide-react";
import type { ExecutionMode, HpcProfile } from "../lib/types";
import type { TaskManagerSummary } from "../lib/taskManager";
import { useOverlayDrawer } from "../lib/useOverlayDrawer";

interface AppNavigationDrawerProps {
  isOpen: boolean;
  executionMode: ExecutionMode;
  activeHpcProfile: HpcProfile | null;
  taskSummary: TaskManagerSummary;
  onClose: () => void;
  onProjects: () => void;
  onStorage: () => void;
  onNodeActivity: () => void;
  onRecoverHeadless: () => void;
  onSettings: () => void;
  onHpcTasks: () => void;
}

export function AppNavigationDrawer({
  isOpen,
  executionMode,
  activeHpcProfile,
  taskSummary,
  onClose,
  onProjects,
  onStorage,
  onNodeActivity,
  onRecoverHeadless,
  onSettings,
  onHpcTasks,
}: AppNavigationDrawerProps) {
  const drawerRef = useOverlayDrawer<HTMLElement>(isOpen, onClose);

  return (
    <div className={`app-drawer-layer left ${isOpen ? "open" : ""}`} aria-hidden={!isOpen} inert={!isOpen}>
      <aside ref={drawerRef} className="app-drawer app-navigation-drawer" role="dialog" aria-modal="false" aria-label="Application menu">
        <div className="app-drawer-header">
          <div>
            <span className="app-drawer-eyebrow">QCortado</span>
            <h2>Workspace</h2>
          </div>
          <button type="button" className="chrome-icon-btn" onClick={onClose} aria-label="Close application menu">
            <X size={19} />
          </button>
        </div>

        <div className="app-drawer-scroll">
          <section className="drawer-status-panel">
            <div className="drawer-status-row">
              <span>Execution</span>
              <strong>{executionMode === "hpc" ? "HPC" : "Local"}</strong>
            </div>
            {executionMode === "hpc" && (
              <div className="drawer-status-row">
                <span>Profile</span>
                <strong>{activeHpcProfile?.name ?? "Not configured"}</strong>
              </div>
            )}
            <button type="button" className="drawer-hpc-summary" onClick={onHpcTasks}>
              <Activity size={17} />
              <span>
                <strong>{taskSummary.hpc} HPC tasks</strong>
                <small>{taskSummary.hpcRunning} running, {taskSummary.hpcQueued} queued</small>
              </span>
            </button>
          </section>

          <nav className="drawer-nav" aria-label="Application navigation">
            <button type="button" onClick={onProjects}><FolderOpen size={18} /><span>Projects</span></button>
            <button type="button" onClick={onStorage}><HardDrive size={18} /><span>Storage Manager</span></button>
            {executionMode === "hpc" && (
              <>
                <button type="button" onClick={onNodeActivity}><Database size={18} /><span>Node Activity</span></button>
                <button type="button" onClick={onRecoverHeadless}><RotateCcw size={18} /><span>Recover Headless Job</span></button>
              </>
            )}
            <button type="button" onClick={onSettings}><Settings size={18} /><span>Settings</span></button>
          </nav>
        </div>
      </aside>
    </div>
  );
}
