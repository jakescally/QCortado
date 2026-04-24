use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, Default)]
#[serde(rename_all = "lowercase")]
pub enum ExecutionMode {
    #[default]
    Local,
    Hpc,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum HpcAuthMethod {
    SshKey,
    Password,
}

impl Default for HpcAuthMethod {
    fn default() -> Self {
        Self::SshKey
    }
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum ResourceType {
    Cpu,
    Gpu,
}

impl Default for ResourceType {
    fn default() -> Self {
        Self::Cpu
    }
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, Default)]
#[serde(rename_all = "snake_case")]
pub enum HpcResourceMode {
    CpuOnly,
    GpuOnly,
    #[default]
    Both,
}

impl HpcResourceMode {
    pub fn supports(self, resource_type: ResourceType) -> bool {
        match (self, resource_type) {
            (Self::Both, _) => true,
            (Self::CpuOnly, ResourceType::Cpu) => true,
            (Self::GpuOnly, ResourceType::Gpu) => true,
            _ => false,
        }
    }

    pub fn preferred_resource_type(self) -> ResourceType {
        match self {
            Self::GpuOnly => ResourceType::Gpu,
            _ => ResourceType::Cpu,
        }
    }
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, Default)]
#[serde(rename_all = "snake_case")]
pub enum HpcLauncher {
    #[default]
    Srun,
    Mpirun,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct SlurmResourceRequest {
    #[serde(default)]
    pub resource_type: ResourceType,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub partition: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub walltime: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub nodes: Option<u32>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub ntasks: Option<u32>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub cpus_per_task: Option<u32>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub memory_gb: Option<u32>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub gpus: Option<u32>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub qos: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub account: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub constraint: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub module_preamble: Option<String>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub additional_sbatch: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct HpcProfile {
    pub id: String,
    pub name: String,
    #[serde(default = "default_cluster")]
    pub cluster: String,
    pub host: String,
    #[serde(default = "default_port")]
    pub port: u16,
    pub username: String,
    #[serde(default)]
    pub auth_method: HpcAuthMethod,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub ssh_key_path: Option<String>,
    pub remote_qe_bin_dir: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub remote_qe_cpu_bin_dir: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub remote_qe_gpu_bin_dir: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub remote_epw_path: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub remote_wannier90_path: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub remote_postw90_path: Option<String>,
    pub remote_pseudo_dir: String,
    pub remote_workspace_root: String,
    pub remote_project_root: String,
    #[serde(default)]
    pub resource_mode: HpcResourceMode,
    #[serde(default)]
    pub launcher: HpcLauncher,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub launcher_cpu_extra_args: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub launcher_gpu_extra_args: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub launcher_extra_args: Option<String>,
    #[serde(default = "default_cpu_resources")]
    pub default_cpu_resources: SlurmResourceRequest,
    #[serde(default = "default_gpu_resources")]
    pub default_gpu_resources: SlurmResourceRequest,
    #[serde(default)]
    pub credential_persisted: bool,
    #[serde(default)]
    pub warning_acknowledged: bool,
    #[serde(default)]
    pub created_at: String,
    #[serde(default)]
    pub updated_at: String,
}

impl HpcProfile {
    pub fn supports_resource_type(&self, resource_type: ResourceType) -> bool {
        self.resource_mode.supports(resource_type)
    }

    pub fn preferred_resource_type(&self) -> ResourceType {
        self.resource_mode.preferred_resource_type()
    }

    pub fn remote_qe_bin_dir_for_resource(&self, resource_type: ResourceType) -> &str {
        let fallback = self.remote_qe_bin_dir.trim();
        match resource_type {
            ResourceType::Cpu => self
                .remote_qe_cpu_bin_dir
                .as_deref()
                .map(str::trim)
                .filter(|value| !value.is_empty())
                .unwrap_or(fallback),
            ResourceType::Gpu => self
                .remote_qe_gpu_bin_dir
                .as_deref()
                .map(str::trim)
                .filter(|value| !value.is_empty())
                .unwrap_or(fallback),
        }
    }

    pub fn launcher_extra_args_for_resource(&self, resource_type: ResourceType) -> Option<&str> {
        let specific = match resource_type {
            ResourceType::Cpu => self.launcher_cpu_extra_args.as_deref(),
            ResourceType::Gpu => self.launcher_gpu_extra_args.as_deref(),
        };
        specific
            .or(self.launcher_extra_args.as_deref())
            .map(str::trim)
            .filter(|value| !value.is_empty())
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct HpcExecutionTarget {
    #[serde(skip_serializing_if = "Option::is_none")]
    pub profile_id: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub resources: Option<SlurmResourceRequest>,
    #[serde(default)]
    pub interactive_debug: bool,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub recovery_save: Option<HpcRecoverySaveSpec>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct ExecutionTarget {
    #[serde(default)]
    pub mode: ExecutionMode,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub hpc: Option<HpcExecutionTarget>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct HpcRecoverySaveSpec {
    pub project_id: String,
    pub cif_id: String,
    pub calc_type: String,
    #[serde(default)]
    pub parameters: serde_json::Value,
    #[serde(default)]
    pub tags: Vec<String>,
    #[serde(default)]
    pub input_content: String,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct ResourceValidation {
    #[serde(default)]
    pub errors: Vec<String>,
    #[serde(default)]
    pub warnings: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct HpcEnvironmentValidation {
    pub reachable: bool,
    pub sbatch_available: bool,
    pub squeue_available: bool,
    pub sacct_available: bool,
    pub qe_pw_available: bool,
    pub qe_epw_available: bool,
    pub qe_pw2wannier_available: bool,
    pub wannier90_available: bool,
    pub postw90_available: bool,
    pub workspace_writable: bool,
    #[serde(default)]
    pub messages: Vec<String>,
}

pub fn default_cluster() -> String {
    "andromeda".to_string()
}

fn default_port() -> u16 {
    22
}

pub fn default_cpu_resources() -> SlurmResourceRequest {
    SlurmResourceRequest {
        resource_type: ResourceType::Cpu,
        partition: Some("short".to_string()),
        walltime: Some("02:00:00".to_string()),
        nodes: Some(1),
        ntasks: Some(4),
        cpus_per_task: Some(1),
        memory_gb: Some(16),
        gpus: Some(0),
        qos: None,
        account: None,
        constraint: None,
        module_preamble: None,
        additional_sbatch: Vec::new(),
    }
}

pub fn default_gpu_resources() -> SlurmResourceRequest {
    SlurmResourceRequest {
        resource_type: ResourceType::Gpu,
        partition: Some("short".to_string()),
        walltime: Some("02:00:00".to_string()),
        nodes: Some(1),
        ntasks: Some(1),
        cpus_per_task: Some(8),
        memory_gb: Some(32),
        gpus: Some(1),
        qos: None,
        account: None,
        constraint: None,
        module_preamble: None,
        additional_sbatch: Vec::new(),
    }
}

pub fn normalize_walltime(input: &str) -> Option<String> {
    let trimmed = input.trim();
    if trimmed.is_empty() {
        return None;
    }

    let parts: Vec<&str> = trimmed.split(':').collect();
    if parts.len() != 3 {
        return None;
    }
    let hours = parts[0].parse::<u32>().ok()?;
    let minutes = parts[1].parse::<u32>().ok()?;
    let seconds = parts[2].parse::<u32>().ok()?;
    if minutes > 59 || seconds > 59 {
        return None;
    }
    Some(format!("{:02}:{:02}:{:02}", hours, minutes, seconds))
}

pub fn walltime_seconds(input: &str) -> Option<u64> {
    let normalized = normalize_walltime(input)?;
    let parts: Vec<&str> = normalized.split(':').collect();
    let hours = parts.first()?.parse::<u64>().ok()?;
    let minutes = parts.get(1)?.parse::<u64>().ok()?;
    let seconds = parts.get(2)?.parse::<u64>().ok()?;
    Some((hours * 3600) + (minutes * 60) + seconds)
}

pub fn andromeda_partition_limit_seconds(partition: &str) -> Option<u64> {
    match partition.trim().to_ascii_lowercase().as_str() {
        "interactive" => Some(12 * 3600),
        "short" => Some(12 * 3600),
        "medium" => Some(48 * 3600),
        "long" => Some(5 * 24 * 3600),
        _ => None,
    }
}

pub fn validate_andromeda_resources(resources: &SlurmResourceRequest) -> ResourceValidation {
    let mut validation = ResourceValidation::default();
    let ntasks = resources.ntasks.unwrap_or(1);
    let cpus_per_task = resources.cpus_per_task.unwrap_or(1);

    let partition = resources
        .partition
        .as_deref()
        .map(|value| value.trim().to_ascii_lowercase())
        .unwrap_or_else(|| "short".to_string());

    if let Some(walltime) = resources.walltime.as_deref() {
        if let Some(seconds) = walltime_seconds(walltime) {
            if let Some(limit) = andromeda_partition_limit_seconds(&partition) {
                if seconds > limit {
                    validation.errors.push(format!(
                        "Requested walltime {} exceeds '{}' partition limit.",
                        walltime, partition
                    ));
                }
            }
        } else {
            validation
                .errors
                .push("Walltime must be in HH:MM:SS format.".to_string());
        }
    }

    if matches!(resources.resource_type, ResourceType::Gpu) {
        let gpus = resources.gpus.unwrap_or(0);
        if gpus == 0 {
            validation
                .errors
                .push("GPU jobs must request at least one GPU (--gres=gpu:N).".to_string());
        } else if ntasks > gpus {
            validation.warnings.push(format!(
                "GPU oversubscription risk: tasks ({}) exceed requested GPUs ({}). Prefer one MPI rank per GPU.",
                ntasks, gpus
            ));
        }
        if resources.nodes.unwrap_or(1) != 1 {
            validation.warnings.push(
                "Andromeda GPU examples default to one node; verify multi-node GPU usage."
                    .to_string(),
            );
        }

        if cpus_per_task < 4 {
            validation.warnings.push(
                "GPU runs with CPUs/task below 4 can bottleneck host-side work; start with 8-16 CPUs/task."
                    .to_string(),
            );
        }
    }

    let total_cpu = ntasks.saturating_mul(cpus_per_task);
    let mem = resources.memory_gb.unwrap_or(0);
    if total_cpu > 16 {
        validation.warnings.push(
            "Time-partition usage guidance is up to 16 aggregate cores per user.".to_string(),
        );
    }
    if mem > 64 {
        validation.warnings.push(
            "Time-partition usage guidance is up to 64GB aggregate memory per user.".to_string(),
        );
    }

    validation.warnings.push(
        "Andromeda time partitions typically allow up to 2 jobs per user concurrently.".to_string(),
    );

    validation
}
