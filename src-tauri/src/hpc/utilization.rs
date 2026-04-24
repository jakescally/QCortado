use serde::Serialize;

use super::profile::ResourceType;

const SQUEUE_BEGIN: &str = "__QCORTADO_SQUEUE_BEGIN__";
const SQUEUE_END: &str = "__QCORTADO_SQUEUE_END__";
const SCONTROL_BEGIN: &str = "__QCORTADO_SCONTROL_BEGIN__";
const SCONTROL_END: &str = "__QCORTADO_SCONTROL_END__";
const SSTAT_BEGIN: &str = "__QCORTADO_SSTAT_BEGIN__";
const SSTAT_END: &str = "__QCORTADO_SSTAT_END__";
const GPU_BEGIN: &str = "__QCORTADO_GPU_BEGIN__";
const GPU_END: &str = "__QCORTADO_GPU_END__";
const ERROR_PREFIX: &str = "__QCORTADO_ERROR__:";

#[derive(Debug, Clone, Serialize)]
pub struct HpcUtilizationSample {
    pub captured_at: String,
    pub resource_type: ResourceType,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub job_id: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub node: Option<String>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub sources: Vec<String>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub warnings: Vec<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub scheduler: Option<HpcSchedulerTelemetry>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub cpu: Option<HpcCpuTelemetry>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub memory: Option<HpcMemoryTelemetry>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub gpu: Option<HpcGpuTelemetry>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub raw: Option<String>,
}

#[derive(Debug, Clone, Serialize, Default)]
pub struct HpcSchedulerTelemetry {
    #[serde(skip_serializing_if = "Option::is_none")]
    pub state: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub node: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub allocated_cpus: Option<u32>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub requested_memory_bytes: Option<u64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub nodes: Option<u32>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub elapsed: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub time_limit: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub reason: Option<String>,
}

#[derive(Debug, Clone, Serialize)]
pub struct HpcCpuTelemetry {
    #[serde(skip_serializing_if = "Option::is_none")]
    pub allocated_cpus: Option<u32>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub total_cpu_seconds: Option<f64>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub steps: Vec<HpcCpuStepTelemetry>,
}

#[derive(Debug, Clone, Serialize)]
pub struct HpcCpuStepTelemetry {
    pub job_step: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub n_tasks: Option<u32>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub average_cpu_seconds: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub average_cpu_display: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub average_rss_bytes: Option<u64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub peak_rss_bytes: Option<u64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub average_vm_bytes: Option<u64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub peak_vm_bytes: Option<u64>,
}

#[derive(Debug, Clone, Serialize)]
pub struct HpcMemoryTelemetry {
    pub source: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub current_rss_bytes: Option<u64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub peak_rss_bytes: Option<u64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub average_vm_bytes: Option<u64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub peak_vm_bytes: Option<u64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub requested_bytes: Option<u64>,
}

#[derive(Debug, Clone, Serialize)]
pub struct HpcGpuTelemetry {
    pub source: String,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub devices: Vec<HpcGpuDeviceTelemetry>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub average_utilization_percent: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub memory_used_bytes: Option<u64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub memory_total_bytes: Option<u64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub memory_used_percent: Option<f64>,
}

#[derive(Debug, Clone, Serialize)]
pub struct HpcGpuDeviceTelemetry {
    pub index: u32,
    pub name: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub utilization_percent: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub memory_utilization_percent: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub memory_used_bytes: Option<u64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub memory_total_bytes: Option<u64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub temperature_c: Option<f64>,
}

#[derive(Debug, Clone)]
pub struct ParsedUtilization {
    pub sources: Vec<String>,
    pub warnings: Vec<String>,
    pub scheduler: Option<HpcSchedulerTelemetry>,
    pub cpu: Option<HpcCpuTelemetry>,
    pub memory: Option<HpcMemoryTelemetry>,
    pub gpu: Option<HpcGpuTelemetry>,
}

fn shell_single_quote(value: &str) -> String {
    if value.is_empty() {
        return "''".to_string();
    }
    let escaped = value.replace('\'', "'\"'\"'");
    format!("'{}'", escaped)
}

pub fn build_utilization_command(job_id: &str, resource_type: ResourceType) -> String {
    let quoted_job_id = shell_single_quote(job_id);
    let gpu_probe = r#"if command -v nvidia-smi >/dev/null 2>&1; then
  nvidia-smi --query-gpu=index,name,utilization.gpu,utilization.memory,memory.used,memory.total,temperature.gpu --format=csv,noheader,nounits 2>&1;
else
  printf '%s\n' '__QCORTADO_ERROR__:nvidia-smi not found inside allocation';
fi"#;
    let quoted_gpu_probe = shell_single_quote(gpu_probe);
    let gpu_section = if matches!(resource_type, ResourceType::Gpu) {
        format!(
            "echo \"{gpu_begin}\"; \
if command -v srun >/dev/null 2>&1; then \
  if command -v timeout >/dev/null 2>&1; then \
    timeout 8s srun --jobid={job_id} --overlap --nodes=1 --ntasks=1 --cpu-bind=none --quiet bash -lc {gpu_probe} 2>&1 || printf '%s\\n' '__QCORTADO_ERROR__:srun GPU probe failed or timed out'; \
  else \
    srun --jobid={job_id} --overlap --nodes=1 --ntasks=1 --cpu-bind=none --quiet bash -lc {gpu_probe} 2>&1 || printf '%s\\n' '__QCORTADO_ERROR__:srun GPU probe failed'; \
  fi; \
else \
  printf '%s\\n' '__QCORTADO_ERROR__:srun not available for allocation-scoped GPU probe'; \
fi; \
echo \"{gpu_end}\"; ",
            gpu_begin = GPU_BEGIN,
            gpu_end = GPU_END,
            job_id = quoted_job_id,
            gpu_probe = quoted_gpu_probe,
        )
    } else {
        String::new()
    };

    format!(
        "echo \"{squeue_begin}\"; \
if command -v squeue >/dev/null 2>&1; then \
  squeue -h -j {job_id} -o \"%T|%N|%C|%m|%D|%M|%l|%R\" 2>&1 || printf '%s\\n' '__QCORTADO_ERROR__:squeue failed'; \
else \
  printf '%s\\n' '__QCORTADO_ERROR__:squeue not found'; \
fi; \
echo \"{squeue_end}\"; \
echo \"{scontrol_begin}\"; \
if command -v scontrol >/dev/null 2>&1; then \
  scontrol show job -o {job_id} 2>&1 || printf '%s\\n' '__QCORTADO_ERROR__:scontrol failed'; \
else \
  printf '%s\\n' '__QCORTADO_ERROR__:scontrol not found'; \
fi; \
echo \"{scontrol_end}\"; \
echo \"{sstat_begin}\"; \
if command -v sstat >/dev/null 2>&1; then \
  sstat --allsteps --parsable2 --noheader --format=JobID,NTasks,AveCPU,AveRSS,MaxRSS,AveVMSize,MaxVMSize -j {job_id} 2>&1 || printf '%s\\n' '__QCORTADO_ERROR__:sstat failed'; \
else \
  printf '%s\\n' '__QCORTADO_ERROR__:sstat not found'; \
fi; \
echo \"{sstat_end}\"; \
{gpu_section}",
        squeue_begin = SQUEUE_BEGIN,
        squeue_end = SQUEUE_END,
        scontrol_begin = SCONTROL_BEGIN,
        scontrol_end = SCONTROL_END,
        sstat_begin = SSTAT_BEGIN,
        sstat_end = SSTAT_END,
        job_id = quoted_job_id,
        gpu_section = gpu_section,
    )
}

pub fn parse_utilization_output(raw: &str, resource_type: ResourceType) -> ParsedUtilization {
    let mut warnings = Vec::new();
    let mut sources = Vec::new();

    let squeue = section_lines(raw, SQUEUE_BEGIN, SQUEUE_END);
    collect_section_errors("squeue", &squeue, &mut warnings);
    let scontrol = section_lines(raw, SCONTROL_BEGIN, SCONTROL_END);
    collect_section_errors("scontrol", &scontrol, &mut warnings);
    let sstat = section_lines(raw, SSTAT_BEGIN, SSTAT_END);
    collect_section_errors("sstat", &sstat, &mut warnings);

    let mut scheduler = parse_squeue_section(&squeue);
    merge_scontrol_section(&mut scheduler, &scontrol);
    if scheduler.is_some() {
        sources.push("squeue/scontrol".to_string());
    }

    let steps = parse_sstat_steps(&sstat);
    if !steps.is_empty() {
        sources.push("sstat".to_string());
    } else if !sstat.iter().any(|line| line.starts_with(ERROR_PREFIX)) {
        warnings.push("Slurm accounting returned no running job-step rows.".to_string());
    }

    let allocated_cpus = scheduler.as_ref().and_then(|value| value.allocated_cpus);
    let cpu = if matches!(resource_type, ResourceType::Cpu) {
        let total_cpu_seconds = steps
            .iter()
            .filter_map(|step| {
                let cpu_seconds = step.average_cpu_seconds?;
                let tasks = step.n_tasks.unwrap_or(1).max(1) as f64;
                Some(cpu_seconds * tasks)
            })
            .sum::<f64>();
        Some(HpcCpuTelemetry {
            allocated_cpus,
            total_cpu_seconds: if total_cpu_seconds > 0.0 {
                Some(total_cpu_seconds)
            } else {
                None
            },
            steps: steps.clone(),
        })
    } else {
        None
    };

    let requested_bytes = scheduler
        .as_ref()
        .and_then(|value| value.requested_memory_bytes);
    let memory = build_memory_telemetry(&steps, requested_bytes);

    let gpu = if matches!(resource_type, ResourceType::Gpu) {
        let gpu_lines = section_lines(raw, GPU_BEGIN, GPU_END);
        collect_section_errors("gpu", &gpu_lines, &mut warnings);
        let parsed = parse_gpu_section(&gpu_lines);
        if parsed
            .as_ref()
            .is_some_and(|value| !value.devices.is_empty())
        {
            sources.push("srun/nvidia-smi".to_string());
        }
        parsed
    } else {
        None
    };

    warnings.sort();
    warnings.dedup();
    sources.sort();
    sources.dedup();

    ParsedUtilization {
        sources,
        warnings,
        scheduler,
        cpu,
        memory,
        gpu,
    }
}

fn section_lines(raw: &str, begin: &str, end: &str) -> Vec<String> {
    let mut in_section = false;
    let mut lines = Vec::new();
    for line in raw.lines() {
        let trimmed = line.trim();
        if trimmed == begin {
            in_section = true;
            continue;
        }
        if trimmed == end {
            break;
        }
        if in_section && !trimmed.is_empty() {
            lines.push(trimmed.to_string());
        }
    }
    lines
}

fn collect_section_errors(section: &str, lines: &[String], warnings: &mut Vec<String>) {
    for line in lines {
        if let Some(message) = line.strip_prefix(ERROR_PREFIX) {
            warnings.push(format!("{}: {}", section, message.trim()));
        }
    }
}

fn parse_squeue_section(lines: &[String]) -> Option<HpcSchedulerTelemetry> {
    let line = lines
        .iter()
        .find(|line| !line.starts_with(ERROR_PREFIX) && line.contains('|'))?;
    let mut parts = line.split('|');
    let state = clean_optional(parts.next());
    let node = clean_optional(parts.next());
    let allocated_cpus = parts.next().and_then(parse_u32);
    let requested_memory_bytes = parts.next().and_then(parse_slurm_memory_to_bytes);
    let nodes = parts.next().and_then(parse_u32);
    let elapsed = clean_optional(parts.next());
    let time_limit = clean_optional(parts.next());
    let reason = clean_optional(parts.next());

    Some(HpcSchedulerTelemetry {
        state,
        node,
        allocated_cpus,
        requested_memory_bytes,
        nodes,
        elapsed,
        time_limit,
        reason,
    })
}

fn merge_scontrol_section(scheduler: &mut Option<HpcSchedulerTelemetry>, lines: &[String]) {
    let Some(line) = lines.iter().find(|line| !line.starts_with(ERROR_PREFIX)) else {
        return;
    };
    let target = scheduler.get_or_insert_with(HpcSchedulerTelemetry::default);
    for token in line.split_whitespace() {
        let Some((key, value)) = token.split_once('=') else {
            continue;
        };
        match key {
            "JobState" if target.state.is_none() => target.state = clean_optional(Some(value)),
            "NodeList" if target.node.is_none() => target.node = clean_optional(Some(value)),
            "NumCPUs" if target.allocated_cpus.is_none() => {
                target.allocated_cpus = parse_u32(value);
            }
            "NumNodes" if target.nodes.is_none() => target.nodes = parse_u32(value),
            "RunTime" if target.elapsed.is_none() => target.elapsed = clean_optional(Some(value)),
            "TimeLimit" if target.time_limit.is_none() => {
                target.time_limit = clean_optional(Some(value));
            }
            "Reason" if target.reason.is_none() => target.reason = clean_optional(Some(value)),
            "MinMemoryNode" | "MinMemoryCPU" if target.requested_memory_bytes.is_none() => {
                target.requested_memory_bytes = parse_slurm_memory_to_bytes(value);
            }
            _ => {}
        }
    }
}

fn parse_sstat_steps(lines: &[String]) -> Vec<HpcCpuStepTelemetry> {
    lines
        .iter()
        .filter(|line| !line.starts_with(ERROR_PREFIX) && line.contains('|'))
        .filter_map(|line| {
            let parts: Vec<&str> = line.split('|').collect();
            if parts.len() < 7 {
                return None;
            }
            let ave_cpu_display = clean_optional(Some(parts[2]));
            let average_cpu_seconds = ave_cpu_display
                .as_deref()
                .and_then(parse_slurm_cpu_time_seconds);
            Some(HpcCpuStepTelemetry {
                job_step: parts[0].trim().to_string(),
                n_tasks: parse_u32(parts[1]),
                average_cpu_seconds,
                average_cpu_display: ave_cpu_display,
                average_rss_bytes: parse_slurm_memory_to_bytes(parts[3]),
                peak_rss_bytes: parse_slurm_memory_to_bytes(parts[4]),
                average_vm_bytes: parse_slurm_memory_to_bytes(parts[5]),
                peak_vm_bytes: parse_slurm_memory_to_bytes(parts[6]),
            })
        })
        .collect()
}

fn build_memory_telemetry(
    steps: &[HpcCpuStepTelemetry],
    requested_bytes: Option<u64>,
) -> Option<HpcMemoryTelemetry> {
    let current_rss_bytes = steps.iter().filter_map(|step| step.average_rss_bytes).max();
    let peak_rss_bytes = steps.iter().filter_map(|step| step.peak_rss_bytes).max();
    let average_vm_bytes = steps.iter().filter_map(|step| step.average_vm_bytes).max();
    let peak_vm_bytes = steps.iter().filter_map(|step| step.peak_vm_bytes).max();
    if current_rss_bytes.is_none()
        && peak_rss_bytes.is_none()
        && average_vm_bytes.is_none()
        && peak_vm_bytes.is_none()
        && requested_bytes.is_none()
    {
        return None;
    }
    Some(HpcMemoryTelemetry {
        source: if current_rss_bytes.is_some() || peak_rss_bytes.is_some() {
            "sstat".to_string()
        } else {
            "scheduler".to_string()
        },
        current_rss_bytes,
        peak_rss_bytes,
        average_vm_bytes,
        peak_vm_bytes,
        requested_bytes,
    })
}

fn parse_gpu_section(lines: &[String]) -> Option<HpcGpuTelemetry> {
    let devices: Vec<HpcGpuDeviceTelemetry> = lines
        .iter()
        .filter(|line| !line.starts_with(ERROR_PREFIX))
        .filter_map(|line| parse_gpu_csv_line(line))
        .collect();
    if devices.is_empty() {
        return None;
    }

    let util_values: Vec<f64> = devices
        .iter()
        .filter_map(|device| device.utilization_percent)
        .collect();
    let average_utilization_percent = if util_values.is_empty() {
        None
    } else {
        Some(util_values.iter().sum::<f64>() / util_values.len() as f64)
    };
    let memory_used_bytes = sum_optional_u64(devices.iter().map(|device| device.memory_used_bytes));
    let memory_total_bytes =
        sum_optional_u64(devices.iter().map(|device| device.memory_total_bytes));
    let memory_used_percent = match (memory_used_bytes, memory_total_bytes) {
        (Some(used), Some(total)) if total > 0 => Some((used as f64 / total as f64) * 100.0),
        _ => None,
    };

    Some(HpcGpuTelemetry {
        source: "srun/nvidia-smi".to_string(),
        devices,
        average_utilization_percent,
        memory_used_bytes,
        memory_total_bytes,
        memory_used_percent,
    })
}

fn parse_gpu_csv_line(line: &str) -> Option<HpcGpuDeviceTelemetry> {
    let parts: Vec<&str> = line.split(',').map(|part| part.trim()).collect();
    if parts.len() < 7 {
        return None;
    }
    let index = parts[0].parse::<u32>().ok()?;
    Some(HpcGpuDeviceTelemetry {
        index,
        name: parts[1].to_string(),
        utilization_percent: parse_percent(parts[2]),
        memory_utilization_percent: parse_percent(parts[3]),
        memory_used_bytes: parse_mib_to_bytes(parts[4]),
        memory_total_bytes: parse_mib_to_bytes(parts[5]),
        temperature_c: parse_f64(parts[6]),
    })
}

fn clean_optional(value: Option<&str>) -> Option<String> {
    let trimmed = value?.trim();
    if trimmed.is_empty()
        || trimmed == "(null)"
        || trimmed == "None"
        || trimmed.eq_ignore_ascii_case("n/a")
        || trimmed.eq_ignore_ascii_case("unknown")
    {
        None
    } else {
        Some(trimmed.to_string())
    }
}

fn parse_u32(value: &str) -> Option<u32> {
    value.trim().parse::<u32>().ok()
}

fn parse_f64(value: &str) -> Option<f64> {
    let trimmed = normalize_numeric(value)?;
    trimmed.parse::<f64>().ok()
}

fn parse_percent(value: &str) -> Option<f64> {
    parse_f64(value.trim().trim_end_matches('%'))
}

fn parse_mib_to_bytes(value: &str) -> Option<u64> {
    let parsed = parse_f64(value)?;
    Some((parsed * 1024.0 * 1024.0).round() as u64)
}

fn sum_optional_u64(values: impl Iterator<Item = Option<u64>>) -> Option<u64> {
    let mut found = false;
    let mut total = 0u64;
    for value in values.flatten() {
        found = true;
        total = total.saturating_add(value);
    }
    found.then_some(total)
}

pub fn parse_slurm_memory_to_bytes(value: &str) -> Option<u64> {
    let mut trimmed = value.trim();
    if trimmed.is_empty()
        || trimmed == "0"
        || trimmed.eq_ignore_ascii_case("n/a")
        || trimmed.eq_ignore_ascii_case("none")
        || trimmed == "(null)"
    {
        return None;
    }
    trimmed = trimmed.trim_end_matches(['n', 'N', 'c', 'C']);
    let unit = trimmed.chars().last().filter(|ch| ch.is_ascii_alphabetic());
    let numeric = if unit.is_some() {
        &trimmed[..trimmed.len().saturating_sub(1)]
    } else {
        trimmed
    };
    let number = normalize_numeric(numeric)?.parse::<f64>().ok()?;
    let multiplier = match unit.map(|ch| ch.to_ascii_uppercase()) {
        Some('T') => 1024f64.powi(4),
        Some('G') => 1024f64.powi(3),
        Some('M') => 1024f64.powi(2),
        Some('K') | None => 1024f64,
        _ => return None,
    };
    Some((number * multiplier).round() as u64)
}

pub fn parse_slurm_cpu_time_seconds(value: &str) -> Option<f64> {
    let trimmed = value.trim();
    if trimmed.is_empty() || trimmed.eq_ignore_ascii_case("n/a") {
        return None;
    }

    let (days, time_part) = if let Some((day_part, rest)) = trimmed.split_once('-') {
        (day_part.parse::<f64>().ok()?, rest)
    } else {
        (0.0, trimmed)
    };

    let segments: Vec<&str> = time_part.split(':').collect();
    let seconds = match segments.as_slice() {
        [seconds] => seconds.parse::<f64>().ok()?,
        [minutes, seconds] => minutes.parse::<f64>().ok()? * 60.0 + seconds.parse::<f64>().ok()?,
        [hours, minutes, seconds] => {
            hours.parse::<f64>().ok()? * 3600.0
                + minutes.parse::<f64>().ok()? * 60.0
                + seconds.parse::<f64>().ok()?
        }
        _ => return None,
    };
    Some(days * 86_400.0 + seconds)
}

fn normalize_numeric(value: &str) -> Option<&str> {
    let trimmed = value.trim();
    if trimmed.is_empty()
        || trimmed.eq_ignore_ascii_case("n/a")
        || trimmed.eq_ignore_ascii_case("nan")
        || trimmed == "[Not Supported]"
    {
        None
    } else {
        Some(trimmed)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_nvidia_smi_csv_for_multiple_gpus() {
        let parsed = parse_gpu_section(&[
            "0, NVIDIA A100-SXM4-40GB, 87, 45, 20480, 40960, 62".to_string(),
            "1, NVIDIA A100-SXM4-40GB, 13, 15, 1024, 40960, 41".to_string(),
        ])
        .expect("gpu telemetry");

        assert_eq!(parsed.devices.len(), 2);
        assert_eq!(parsed.devices[0].index, 0);
        assert_eq!(parsed.devices[0].memory_used_bytes, Some(21_474_836_480));
        assert_eq!(parsed.memory_total_bytes, Some(85_899_345_920));
        assert_eq!(parsed.average_utilization_percent, Some(50.0));
    }

    #[test]
    fn skips_gpu_missing_values_without_failing_device() {
        let parsed = parse_gpu_csv_line("0, NVIDIA A40, N/A, [Not Supported], 128, 46068, N/A")
            .expect("device");
        assert_eq!(parsed.utilization_percent, None);
        assert_eq!(parsed.memory_utilization_percent, None);
        assert_eq!(parsed.memory_used_bytes, Some(134_217_728));
        assert_eq!(parsed.temperature_c, None);
    }

    #[test]
    fn parses_sstat_steps_and_memory_units() {
        let steps = parse_sstat_steps(&[
            "123.batch|1|00:00:12|1024K|2048K|10M|12M".to_string(),
            "123.0|4|00:02:30|1G|2G|3G|4G".to_string(),
        ]);

        assert_eq!(steps.len(), 2);
        assert_eq!(steps[1].n_tasks, Some(4));
        assert_eq!(steps[1].average_cpu_seconds, Some(150.0));
        assert_eq!(steps[1].peak_rss_bytes, Some(2_147_483_648));
    }

    #[test]
    fn parses_slurm_memory_units() {
        assert_eq!(parse_slurm_memory_to_bytes("1K"), Some(1024));
        assert_eq!(parse_slurm_memory_to_bytes("1M"), Some(1_048_576));
        assert_eq!(parse_slurm_memory_to_bytes("1G"), Some(1_073_741_824));
        assert_eq!(parse_slurm_memory_to_bytes("1T"), Some(1_099_511_627_776));
        assert_eq!(parse_slurm_memory_to_bytes("1024"), Some(1_048_576));
        assert_eq!(parse_slurm_memory_to_bytes(""), None);
        assert_eq!(parse_slurm_memory_to_bytes("0"), None);
    }

    #[test]
    fn parses_cpu_time_variants() {
        assert_eq!(parse_slurm_cpu_time_seconds("12"), Some(12.0));
        assert_eq!(parse_slurm_cpu_time_seconds("02:03"), Some(123.0));
        assert_eq!(parse_slurm_cpu_time_seconds("01:02:03"), Some(3723.0));
        assert_eq!(parse_slurm_cpu_time_seconds("1-01:00:00"), Some(90_000.0));
    }
}
