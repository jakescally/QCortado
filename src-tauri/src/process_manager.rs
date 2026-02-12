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
}

#[derive(Debug, Clone, serde::Serialize)]
pub struct TaskSummary {
    pub task_id: String,
    pub task_type: String,
    pub label: String,
    pub started_at: String,
    pub status: TaskStatus,
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
