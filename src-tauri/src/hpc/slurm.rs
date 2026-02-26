use super::profile::{
    validate_andromeda_resources, HpcProfile, ResourceType, ResourceValidation,
    SlurmResourceRequest,
};

#[derive(Debug, Clone)]
pub struct SlurmScript {
    pub script: String,
    pub sbatch_preview: String,
    pub effective_resources: SlurmResourceRequest,
    pub validation: ResourceValidation,
}

#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct SchedulerSnapshot {
    pub state: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub node: Option<String>,
    pub source: String,
}

pub fn merge_resources(
    profile: &HpcProfile,
    override_resources: Option<SlurmResourceRequest>,
) -> SlurmResourceRequest {
    let mut base = if let Some(ref request) = override_resources {
        match request.resource_type {
            ResourceType::Gpu => profile.default_gpu_resources.clone(),
            ResourceType::Cpu => profile.default_cpu_resources.clone(),
        }
    } else {
        match profile.preferred_resource_type() {
            ResourceType::Gpu => profile.default_gpu_resources.clone(),
            ResourceType::Cpu => profile.default_cpu_resources.clone(),
        }
    };

    if let Some(request) = override_resources {
        base.resource_type = request.resource_type;
        if request.partition.is_some() {
            base.partition = request.partition;
        }
        if request.walltime.is_some() {
            base.walltime = request.walltime;
        }
        if request.nodes.is_some() {
            base.nodes = request.nodes;
        }
        if request.ntasks.is_some() {
            base.ntasks = request.ntasks;
        }
        if request.cpus_per_task.is_some() {
            base.cpus_per_task = request.cpus_per_task;
        }
        if request.memory_gb.is_some() {
            base.memory_gb = request.memory_gb;
        }
        if request.gpus.is_some() {
            base.gpus = request.gpus;
        }
        if request.qos.is_some() {
            base.qos = request.qos;
        }
        if request.account.is_some() {
            base.account = request.account;
        }
        if request.constraint.is_some() {
            base.constraint = request.constraint;
        }
        if request.module_preamble.is_some() {
            base.module_preamble = request.module_preamble;
        }
        if !request.additional_sbatch.is_empty() {
            base.additional_sbatch = request.additional_sbatch;
        }
    }

    if matches!(base.resource_type, ResourceType::Gpu) {
        if base.gpus.unwrap_or(0) == 0 {
            base.gpus = Some(1);
        }
        if base.nodes.is_none() {
            base.nodes = Some(1);
        }
    }

    if base.nodes.is_none() {
        base.nodes = Some(1);
    }
    if base.ntasks.is_none() {
        base.ntasks = Some(1);
    }
    if base.cpus_per_task.is_none() {
        base.cpus_per_task = Some(1);
    }
    if base.partition.is_none() {
        base.partition = Some("short".to_string());
    }
    if base.walltime.is_none() {
        base.walltime = Some("02:00:00".to_string());
    }

    base
}

pub fn build_slurm_script(
    profile: &HpcProfile,
    job_name: &str,
    script_body_commands: &[String],
    override_resources: Option<SlurmResourceRequest>,
) -> SlurmScript {
    let resources = merge_resources(profile, override_resources);
    let mut validation = validate_andromeda_resources(&resources);
    if !profile.supports_resource_type(resources.resource_type) {
        let requested = match resources.resource_type {
            ResourceType::Cpu => "CPU",
            ResourceType::Gpu => "GPU",
        };
        let allowed = match profile.resource_mode {
            super::profile::HpcResourceMode::CpuOnly => "CPU-only",
            super::profile::HpcResourceMode::GpuOnly => "GPU-only",
            super::profile::HpcResourceMode::Both => "CPU+GPU",
        };
        validation.errors.push(format!(
            "Profile '{}' is {} and does not support {} runs.",
            profile.name, allowed, requested
        ));
    }

    let mut header = vec![
        "#!/bin/bash".to_string(),
        format!("#SBATCH --job-name={}", job_name),
        "#SBATCH --output=slurm.out".to_string(),
        "#SBATCH --error=slurm.err".to_string(),
    ];

    if let Some(partition) = resources.partition.as_deref() {
        header.push(format!("#SBATCH --partition={}", partition.trim()));
    }
    if let Some(walltime) = resources.walltime.as_deref() {
        header.push(format!("#SBATCH --time={}", walltime.trim()));
    }
    if let Some(nodes) = resources.nodes {
        header.push(format!("#SBATCH --nodes={}", nodes.max(1)));
    }
    if let Some(ntasks) = resources.ntasks {
        header.push(format!("#SBATCH --ntasks={}", ntasks.max(1)));
    }
    if let Some(cpus_per_task) = resources.cpus_per_task {
        header.push(format!("#SBATCH --cpus-per-task={}", cpus_per_task.max(1)));
    }
    if let Some(memory_gb) = resources.memory_gb {
        header.push(format!("#SBATCH --mem={}G", memory_gb.max(1)));
    }
    if matches!(resources.resource_type, ResourceType::Gpu) {
        header.push(format!(
            "#SBATCH --gres=gpu:{}",
            resources.gpus.unwrap_or(1).max(1)
        ));
    }
    if let Some(qos) = resources.qos.as_deref() {
        if !qos.trim().is_empty() {
            header.push(format!("#SBATCH --qos={}", qos.trim()));
        }
    }
    if let Some(account) = resources.account.as_deref() {
        if !account.trim().is_empty() {
            header.push(format!("#SBATCH --account={}", account.trim()));
        }
    }
    if let Some(constraint) = resources.constraint.as_deref() {
        if !constraint.trim().is_empty() {
            header.push(format!("#SBATCH --constraint={}", constraint.trim()));
        }
    }
    for additional in &resources.additional_sbatch {
        let trimmed = additional.trim();
        if !trimmed.is_empty() {
            if trimmed.starts_with("#SBATCH") {
                header.push(trimmed.to_string());
            } else {
                header.push(format!("#SBATCH {}", trimmed));
            }
        }
    }

    let mut script_lines = header.clone();
    script_lines.push("set -euo pipefail".to_string());
    script_lines.push("echo \"[QCortado] Job started on $(date -u)\"".to_string());
    script_lines.push("echo \"[QCortado] Host: $(hostname)\"".to_string());
    script_lines.push("echo \"[QCortado] Working directory: $(pwd)\"".to_string());

    if let Some(module_preamble) = resources.module_preamble.as_deref() {
        let trimmed = module_preamble.trim();
        if !trimmed.is_empty() {
            script_lines.push(String::new());
            script_lines.push("# User module/setup preamble".to_string());
            script_lines.extend(trimmed.lines().map(|line| line.to_string()));
        }
    }

    script_lines.push(String::new());
    script_lines.push("# Calculation commands".to_string());
    for cmd in script_body_commands {
        script_lines.push(cmd.to_string());
    }

    script_lines.push(String::new());
    script_lines.push("echo \"[QCortado] Job finished on $(date -u)\"".to_string());

    let preview = format!(
        "sbatch --partition={} --time={} --nodes={} --ntasks={} --cpus-per-task={}{} run.sbatch",
        resources.partition.as_deref().unwrap_or("short"),
        resources.walltime.as_deref().unwrap_or("02:00:00"),
        resources.nodes.unwrap_or(1),
        resources.ntasks.unwrap_or(1),
        resources.cpus_per_task.unwrap_or(1),
        if matches!(resources.resource_type, ResourceType::Gpu) {
            format!(" --gres=gpu:{}", resources.gpus.unwrap_or(1))
        } else {
            String::new()
        }
    );

    SlurmScript {
        script: script_lines.join("\n"),
        sbatch_preview: preview,
        effective_resources: resources,
        validation,
    }
}

pub fn parse_sbatch_job_id(output: &str) -> Option<String> {
    output
        .split_whitespace()
        .rev()
        .find(|segment| segment.chars().all(|ch| ch.is_ascii_digit()))
        .map(|value| value.to_string())
}

pub fn parse_squeue_snapshot(output: &str) -> Option<SchedulerSnapshot> {
    let line = output
        .lines()
        .find(|line| !line.trim().is_empty())
        .map(|line| line.trim().to_string())?;

    let mut parts = line.split('|');
    let state = parts.next()?.trim().to_string();
    let node = parts
        .next()
        .map(|value| value.trim().to_string())
        .filter(|value| !value.is_empty() && value != "(null)");

    Some(SchedulerSnapshot {
        state,
        node,
        source: "squeue".to_string(),
    })
}

pub fn parse_sacct_snapshot(output: &str) -> Option<SchedulerSnapshot> {
    let line = output
        .lines()
        .map(|line| line.trim())
        .find(|line| !line.is_empty() && !line.contains("batch") && !line.contains("extern"))?;
    let mut parts = line.split('|');
    let state = parts.next()?.trim().to_string();
    let node = parts
        .next()
        .map(|value| value.trim().to_string())
        .filter(|value| !value.is_empty() && value != "None");

    Some(SchedulerSnapshot {
        state,
        node,
        source: "sacct".to_string(),
    })
}

pub fn normalize_scheduler_state(state: &str) -> String {
    let trimmed = state.trim();
    if trimmed.is_empty() {
        return "UNKNOWN".to_string();
    }
    let upper = trimmed.to_ascii_uppercase();
    upper
        .split_whitespace()
        .next()
        .unwrap_or("UNKNOWN")
        .split('+')
        .next()
        .unwrap_or("UNKNOWN")
        .to_string()
}

pub fn is_terminal_state(state: &str) -> bool {
    matches!(
        normalize_scheduler_state(state).as_str(),
        "COMPLETED"
            | "FAILED"
            | "CANCELLED"
            | "TIMEOUT"
            | "OUT_OF_MEMORY"
            | "NODE_FAIL"
            | "PREEMPTED"
            | "BOOT_FAIL"
            | "DEADLINE"
    )
}
