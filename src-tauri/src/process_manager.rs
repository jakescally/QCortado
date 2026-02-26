use std::collections::HashMap;
use std::sync::Arc;
use tokio::sync::Mutex;

#[derive(Debug, Clone, serde::Serialize, serde::Deserialize, PartialEq)]
#[serde(rename_all = "lowercase")]
pub enum TaskStatus {
    Running,
    Completed,
    Failed,
    Cancelled,
}

#[derive(Debug, Clone, serde::Serialize)]
pub struct TaskInfo {
    pub task_id: String,
    pub task_type: String,
    pub label: String,
    pub started_at: String,
    pub status: TaskStatus,
    pub result: Option<serde_json::Value>,
    pub error: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub backend: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub remote_job_id: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub scheduler_state: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub remote_node: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub remote_workdir: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub remote_project_path: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub remote_storage_bytes: Option<u64>,
}

#[derive(Debug, Clone, serde::Serialize)]
pub struct TaskSummary {
    pub task_id: String,
    pub task_type: String,
    pub label: String,
    pub started_at: String,
    pub status: TaskStatus,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub backend: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub remote_job_id: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub scheduler_state: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub remote_node: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub remote_workdir: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub remote_project_path: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub remote_storage_bytes: Option<u64>,
}

pub struct RunningTask {
    pub task_id: String,
    pub task_type: String,
    pub label: String,
    pub started_at: String,
    pub output_buffer: Vec<String>,
    pub status: TaskStatus,
    pub result: Option<serde_json::Value>,
    pub error: Option<String>,
    pub child_id: Option<u32>,
    pub cancel_flag: Arc<std::sync::atomic::AtomicBool>,
    pub backend: Option<String>,
    pub remote_job_id: Option<String>,
    pub scheduler_state: Option<String>,
    pub remote_node: Option<String>,
    pub remote_workdir: Option<String>,
    pub remote_project_path: Option<String>,
    pub remote_storage_bytes: Option<u64>,
    pub hpc_profile_id: Option<String>,
    pub local_sync_dir: Option<String>,
}

#[derive(Debug, Clone)]
pub struct HpcTransferContext {
    pub task_id: String,
    pub task_type: String,
    pub status: TaskStatus,
    pub backend: Option<String>,
    pub remote_workdir: Option<String>,
    pub hpc_profile_id: Option<String>,
    pub local_sync_dir: Option<String>,
}

/// Thread-safe process manager. Clone is cheap (shared Arc).
#[derive(Clone)]
pub struct ProcessManager {
    tasks: Arc<Mutex<HashMap<String, RunningTask>>>,
}

impl ProcessManager {
    pub fn new() -> Self {
        Self {
            tasks: Arc::new(Mutex::new(HashMap::new())),
        }
    }

    pub async fn register(
        &self,
        task_type: String,
        label: String,
    ) -> (String, Arc<std::sync::atomic::AtomicBool>) {
        let task_id = uuid::Uuid::new_v4().to_string();
        let cancel_flag = Arc::new(std::sync::atomic::AtomicBool::new(false));
        let task = RunningTask {
            task_id: task_id.clone(),
            task_type,
            label,
            started_at: chrono::Utc::now().to_rfc3339(),
            output_buffer: Vec::new(),
            status: TaskStatus::Running,
            result: None,
            error: None,
            child_id: None,
            cancel_flag: cancel_flag.clone(),
            backend: Some("local".to_string()),
            remote_job_id: None,
            scheduler_state: None,
            remote_node: None,
            remote_workdir: None,
            remote_project_path: None,
            remote_storage_bytes: None,
            hpc_profile_id: None,
            local_sync_dir: None,
        };
        let mut tasks = self.tasks.lock().await;
        tasks.insert(task_id.clone(), task);
        (task_id, cancel_flag)
    }

    pub async fn has_running_tasks(&self) -> bool {
        let tasks = self.tasks.lock().await;
        tasks.values().any(|t| t.status == TaskStatus::Running)
    }

    pub async fn append_output(&self, task_id: &str, line: String) {
        let mut tasks = self.tasks.lock().await;
        if let Some(task) = tasks.get_mut(task_id) {
            task.output_buffer.push(line);
        }
    }

    pub async fn set_child_id(&self, task_id: &str, pid: u32) {
        let mut tasks = self.tasks.lock().await;
        if let Some(task) = tasks.get_mut(task_id) {
            task.child_id = Some(pid);
        }
    }

    pub async fn complete(&self, task_id: &str, result: serde_json::Value) {
        let mut tasks = self.tasks.lock().await;
        if let Some(task) = tasks.get_mut(task_id) {
            task.status = TaskStatus::Completed;
            task.result = Some(result);
        }
    }

    pub async fn fail(&self, task_id: &str, error: String) {
        let mut tasks = self.tasks.lock().await;
        if let Some(task) = tasks.get_mut(task_id) {
            task.status = TaskStatus::Failed;
            task.error = Some(error);
        }
    }

    pub async fn get_task(&self, task_id: &str) -> Option<TaskInfo> {
        let tasks = self.tasks.lock().await;
        tasks.get(task_id).map(|t| TaskInfo {
            task_id: t.task_id.clone(),
            task_type: t.task_type.clone(),
            label: t.label.clone(),
            started_at: t.started_at.clone(),
            status: t.status.clone(),
            result: t.result.clone(),
            error: t.error.clone(),
            backend: t.backend.clone(),
            remote_job_id: t.remote_job_id.clone(),
            scheduler_state: t.scheduler_state.clone(),
            remote_node: t.remote_node.clone(),
            remote_workdir: t.remote_workdir.clone(),
            remote_project_path: t.remote_project_path.clone(),
            remote_storage_bytes: t.remote_storage_bytes,
        })
    }

    pub async fn list_tasks(&self) -> Vec<TaskSummary> {
        let tasks = self.tasks.lock().await;
        tasks
            .values()
            .map(|t| TaskSummary {
                task_id: t.task_id.clone(),
                task_type: t.task_type.clone(),
                label: t.label.clone(),
                started_at: t.started_at.clone(),
                status: t.status.clone(),
                backend: t.backend.clone(),
                remote_job_id: t.remote_job_id.clone(),
                scheduler_state: t.scheduler_state.clone(),
                remote_node: t.remote_node.clone(),
                remote_workdir: t.remote_workdir.clone(),
                remote_project_path: t.remote_project_path.clone(),
                remote_storage_bytes: t.remote_storage_bytes,
            })
            .collect()
    }

    pub async fn get_output(&self, task_id: &str) -> Option<Vec<String>> {
        let tasks = self.tasks.lock().await;
        tasks.get(task_id).map(|t| t.output_buffer.clone())
    }

    pub async fn cancel(&self, task_id: &str) -> Result<(), String> {
        let mut tasks = self.tasks.lock().await;
        let task = tasks
            .get_mut(task_id)
            .ok_or_else(|| format!("Task not found: {}", task_id))?;

        if task.status != TaskStatus::Running {
            return Err("Task is not running".to_string());
        }

        task.cancel_flag
            .store(true, std::sync::atomic::Ordering::SeqCst);

        if let Some(pid) = task.child_id {
            kill_process(pid);
        }

        task.status = TaskStatus::Cancelled;
        task.error = Some("Cancelled by user".to_string());
        Ok(())
    }

    pub async fn remove(&self, task_id: &str) {
        let mut tasks = self.tasks.lock().await;
        tasks.remove(task_id);
    }

    pub async fn kill_all(&self) {
        let mut tasks = self.tasks.lock().await;
        for task in tasks.values_mut() {
            if task.status == TaskStatus::Running {
                task.cancel_flag
                    .store(true, std::sync::atomic::Ordering::SeqCst);
                if let Some(pid) = task.child_id {
                    kill_process(pid);
                }
                task.status = TaskStatus::Cancelled;
            }
        }
    }

    pub async fn set_task_backend(&self, task_id: &str, backend: Option<String>) {
        let mut tasks = self.tasks.lock().await;
        if let Some(task) = tasks.get_mut(task_id) {
            task.backend = backend;
        }
    }

    pub async fn set_remote_job_id(&self, task_id: &str, remote_job_id: Option<String>) {
        let mut tasks = self.tasks.lock().await;
        if let Some(task) = tasks.get_mut(task_id) {
            task.remote_job_id = remote_job_id;
        }
    }

    pub async fn set_scheduler_state(&self, task_id: &str, scheduler_state: Option<String>) {
        let mut tasks = self.tasks.lock().await;
        if let Some(task) = tasks.get_mut(task_id) {
            task.scheduler_state = scheduler_state;
        }
    }

    pub async fn set_remote_node(&self, task_id: &str, remote_node: Option<String>) {
        let mut tasks = self.tasks.lock().await;
        if let Some(task) = tasks.get_mut(task_id) {
            task.remote_node = remote_node;
        }
    }

    pub async fn set_remote_workdir(&self, task_id: &str, remote_workdir: Option<String>) {
        let mut tasks = self.tasks.lock().await;
        if let Some(task) = tasks.get_mut(task_id) {
            task.remote_workdir = remote_workdir;
        }
    }

    pub async fn set_remote_project_path(
        &self,
        task_id: &str,
        remote_project_path: Option<String>,
    ) {
        let mut tasks = self.tasks.lock().await;
        if let Some(task) = tasks.get_mut(task_id) {
            task.remote_project_path = remote_project_path;
        }
    }

    pub async fn set_remote_storage_bytes(&self, task_id: &str, remote_storage_bytes: Option<u64>) {
        let mut tasks = self.tasks.lock().await;
        if let Some(task) = tasks.get_mut(task_id) {
            task.remote_storage_bytes = remote_storage_bytes;
        }
    }

    pub async fn set_hpc_profile_id(&self, task_id: &str, hpc_profile_id: Option<String>) {
        let mut tasks = self.tasks.lock().await;
        if let Some(task) = tasks.get_mut(task_id) {
            task.hpc_profile_id = hpc_profile_id;
        }
    }

    pub async fn set_local_sync_dir(&self, task_id: &str, local_sync_dir: Option<String>) {
        let mut tasks = self.tasks.lock().await;
        if let Some(task) = tasks.get_mut(task_id) {
            task.local_sync_dir = local_sync_dir;
        }
    }

    pub async fn get_hpc_transfer_context(&self, task_id: &str) -> Option<HpcTransferContext> {
        let tasks = self.tasks.lock().await;
        tasks.get(task_id).map(|task| HpcTransferContext {
            task_id: task.task_id.clone(),
            task_type: task.task_type.clone(),
            status: task.status.clone(),
            backend: task.backend.clone(),
            remote_workdir: task.remote_workdir.clone(),
            hpc_profile_id: task.hpc_profile_id.clone(),
            local_sync_dir: task.local_sync_dir.clone(),
        })
    }
}

fn kill_process(pid: u32) {
    #[cfg(unix)]
    {
        unsafe {
            // Kill the process group to ensure MPI children are also killed
            libc::kill(-(pid as i32), libc::SIGTERM);
            // Also kill the process directly as fallback
            libc::kill(pid as i32, libc::SIGTERM);
        }
    }
    #[cfg(not(unix))]
    {
        let _ = pid;
    }
}
