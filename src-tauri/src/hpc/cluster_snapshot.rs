use std::collections::{BTreeMap, BTreeSet};

use serde::{Deserialize, Serialize};

use super::{profile::HpcProfile, slurm::normalize_scheduler_state};

const NODE_SECTION_BEGIN: &str = "__QCORTADO_CLUSTER_NODES_BEGIN__";
const NODE_SECTION_END: &str = "__QCORTADO_CLUSTER_NODES_END__";
const QUEUE_SECTION_BEGIN: &str = "__QCORTADO_CLUSTER_QUEUE_BEGIN__";
const QUEUE_SECTION_END: &str = "__QCORTADO_CLUSTER_QUEUE_END__";
const NODE_ERROR_PREFIX: &str = "__QCORTADO_NODES_ERROR__:";
const QUEUE_ERROR_PREFIX: &str = "__QCORTADO_QUEUE_ERROR__:";

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum QueueScope {
    All,
    Mine,
}

impl QueueScope {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::All => "all",
            Self::Mine => "mine",
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct HpcNodeSnapshot {
    pub node_name: String,
    #[serde(default)]
    pub partitions: Vec<String>,
    pub state: String,
    pub raw_state: String,
    pub cpus: u32,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub cpu_name: Option<String>,
    pub memory_mb: u64,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub free_memory_mb: Option<u64>,
    pub gpus: u32,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub gpu_name: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub gres: Option<String>,
    #[serde(default)]
    pub features: Vec<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub reason: Option<String>,
    pub node_type: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct HpcQueueSnapshot {
    pub job_id: String,
    pub user: String,
    pub state: String,
    pub raw_state: String,
    pub partition: String,
    pub nodes: u32,
    pub elapsed: String,
    pub time_limit: String,
    pub reason: String,
    pub nodelist: String,
    pub name: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct HpcClusterSnapshot {
    pub captured_at: String,
    pub cluster: String,
    pub host: String,
    pub queue_scope: String,
    pub queue_included: bool,
    pub queue_limit: usize,
    #[serde(default)]
    pub nodes: Vec<HpcNodeSnapshot>,
    #[serde(default)]
    pub queue: Vec<HpcQueueSnapshot>,
    #[serde(default)]
    pub warnings: Vec<String>,
}

#[derive(Debug, Clone)]
pub struct ParsedClusterSnapshot {
    pub nodes: Vec<HpcNodeSnapshot>,
    pub queue: Vec<HpcQueueSnapshot>,
    pub warnings: Vec<String>,
}

pub fn normalize_queue_scope(queue_scope: Option<String>) -> QueueScope {
    match queue_scope
        .as_deref()
        .unwrap_or("")
        .trim()
        .to_ascii_lowercase()
        .as_str()
    {
        "mine" | "my" | "user" => QueueScope::Mine,
        _ => QueueScope::All,
    }
}

pub fn is_andromeda_profile(profile: &HpcProfile) -> bool {
    let cluster = profile.cluster.trim().to_ascii_lowercase();
    let host = profile.host.trim().to_ascii_lowercase();
    cluster == "andromeda" || host == "andromeda.bc.edu"
}

pub fn build_cluster_snapshot_command(
    queue_scope: QueueScope,
    include_queue: bool,
    queue_limit: usize,
) -> String {
    let queue_user_clause = match queue_scope {
        QueueScope::Mine => "-u \"$USER\"",
        QueueScope::All => "",
    };
    let include_queue_flag = if include_queue { 1 } else { 0 };

    format!(
        "echo \"{node_begin}\"; \
nodes_output=$(sinfo -N -h -o \"%N|%P|%T|%c|%m|%e|%G|%f|%b|%E\" 2>&1); \
nodes_status=$?; \
if [ $nodes_status -eq 0 ]; then \
  printf '%s\\n' \"$nodes_output\"; \
else \
  nodes_error=$(printf '%s' \"$nodes_output\" | tr '\\n' ';'); \
  printf '%s\\n' \"{node_err}$nodes_error\"; \
fi; \
echo \"{node_end}\"; \
if [ {include_queue_flag} -eq 1 ]; then \
  echo \"{queue_begin}\"; \
  queue_output=$(squeue -h {queue_user_clause} -o \"%i|%u|%T|%P|%D|%M|%l|%R|%N|%j\" 2>&1); \
  queue_status=$?; \
  if [ $queue_status -eq 0 ]; then \
    printf '%s\\n' \"$queue_output\" | head -n {queue_limit}; \
  else \
    queue_error=$(printf '%s' \"$queue_output\" | tr '\\n' ';'); \
    printf '%s\\n' \"{queue_err}$queue_error\"; \
  fi; \
  echo \"{queue_end}\"; \
fi",
        node_begin = NODE_SECTION_BEGIN,
        node_end = NODE_SECTION_END,
        node_err = NODE_ERROR_PREFIX,
        queue_begin = QUEUE_SECTION_BEGIN,
        queue_end = QUEUE_SECTION_END,
        queue_err = QUEUE_ERROR_PREFIX,
        include_queue_flag = include_queue_flag,
        queue_user_clause = queue_user_clause,
        queue_limit = queue_limit.max(1),
    )
}

pub fn parse_cluster_snapshot_output(
    output: &str,
    include_queue: bool,
) -> Result<ParsedClusterSnapshot, String> {
    let mut node_lines: Vec<String> = Vec::new();
    let mut queue_lines: Vec<String> = Vec::new();
    let mut warnings: Vec<String> = Vec::new();

    let mut in_nodes = false;
    let mut in_queue = false;
    let mut saw_nodes_begin = false;
    let mut saw_nodes_end = false;
    let mut saw_queue_begin = false;
    let mut saw_queue_end = false;

    for line in output.lines() {
        let trimmed = line.trim();
        if trimmed == NODE_SECTION_BEGIN {
            in_nodes = true;
            in_queue = false;
            saw_nodes_begin = true;
            continue;
        }
        if trimmed == NODE_SECTION_END {
            in_nodes = false;
            saw_nodes_end = true;
            continue;
        }
        if trimmed == QUEUE_SECTION_BEGIN {
            in_queue = true;
            in_nodes = false;
            saw_queue_begin = true;
            continue;
        }
        if trimmed == QUEUE_SECTION_END {
            in_queue = false;
            saw_queue_end = true;
            continue;
        }

        if in_nodes {
            if let Some(msg) = trimmed.strip_prefix(NODE_ERROR_PREFIX) {
                let message = msg.trim();
                warnings.push(if message.is_empty() {
                    "Node query failed.".to_string()
                } else {
                    format!("Node query failed: {}", message)
                });
            } else if !trimmed.is_empty() {
                node_lines.push(trimmed.to_string());
            }
            continue;
        }

        if in_queue {
            if let Some(msg) = trimmed.strip_prefix(QUEUE_ERROR_PREFIX) {
                let message = msg.trim();
                warnings.push(if message.is_empty() {
                    "Queue query failed.".to_string()
                } else {
                    format!("Queue query failed: {}", message)
                });
            } else if !trimmed.is_empty() {
                queue_lines.push(trimmed.to_string());
            }
        }
    }

    if !saw_nodes_begin || !saw_nodes_end {
        return Err(
            "Failed to parse cluster snapshot output (missing node section markers).".to_string(),
        );
    }
    if include_queue && (!saw_queue_begin || !saw_queue_end) {
        warnings.push("Queue section was missing from the snapshot response.".to_string());
    }

    Ok(ParsedClusterSnapshot {
        nodes: parse_node_lines(&node_lines),
        queue: parse_queue_lines(&queue_lines),
        warnings,
    })
}

fn parse_node_lines(lines: &[String]) -> Vec<HpcNodeSnapshot> {
    let mut merged: BTreeMap<String, HpcNodeSnapshot> = BTreeMap::new();

    for line in lines {
        let Some(parsed) = parse_node_line(line) else {
            continue;
        };
        let key = parsed.node_name.clone();
        if let Some(existing) = merged.get_mut(&key) {
            let merged_partitions = merge_sorted_unique(&existing.partitions, &parsed.partitions);
            existing.partitions = merged_partitions;
            existing.features = merge_sorted_unique(&existing.features, &parsed.features);
            existing.cpus = existing.cpus.max(parsed.cpus);
            if existing.cpu_name.is_none() {
                existing.cpu_name = parsed.cpu_name;
            }
            existing.memory_mb = existing.memory_mb.max(parsed.memory_mb);
            existing.free_memory_mb = match (existing.free_memory_mb, parsed.free_memory_mb) {
                (Some(a), Some(b)) => Some(a.max(b)),
                (Some(a), None) => Some(a),
                (None, Some(b)) => Some(b),
                (None, None) => None,
            };
            existing.gpus = existing.gpus.max(parsed.gpus);
            if existing.gpu_name.is_none() {
                existing.gpu_name = parsed.gpu_name;
            }
            if existing.gres.is_none() {
                existing.gres = parsed.gres;
            }
            if existing.reason.is_none() {
                existing.reason = parsed.reason;
            }
            if existing.state == "IDLE" && parsed.state != "IDLE" {
                existing.state = parsed.state;
                existing.raw_state = parsed.raw_state;
            }
            if existing.node_type == "cpu" && parsed.node_type == "gpu" {
                existing.node_type = "gpu".to_string();
            }
        } else {
            merged.insert(key, parsed);
        }
    }

    merged.into_values().collect()
}

fn parse_queue_lines(lines: &[String]) -> Vec<HpcQueueSnapshot> {
    let mut queue: Vec<HpcQueueSnapshot> = lines
        .iter()
        .filter_map(|line| parse_queue_line(line))
        .collect();
    queue.sort_by(|a, b| {
        let a_id = a.job_id.parse::<u64>().unwrap_or(0);
        let b_id = b.job_id.parse::<u64>().unwrap_or(0);
        b_id.cmp(&a_id).then_with(|| a.job_id.cmp(&b.job_id))
    });
    queue
}

fn parse_node_line(line: &str) -> Option<HpcNodeSnapshot> {
    let mut parts = line.splitn(10, '|');
    let node_name = sanitize_required(parts.next()?)?;
    let partition_raw = sanitize_optional(parts.next());
    let raw_state = sanitize_optional(parts.next()).unwrap_or_else(|| "UNKNOWN".to_string());
    let cpus = parse_u32(parts.next()).unwrap_or(0);
    let memory_mb = parse_u64(parts.next()).unwrap_or(0);
    let free_memory_mb = parse_u64_optional(parts.next());
    let gres_raw = sanitize_optional(parts.next());
    let features_raw = sanitize_optional(parts.next());
    let active_features_raw = sanitize_optional(parts.next());
    let reason = sanitize_optional(parts.next());

    let partitions = parse_csv_no_empty(&partition_raw.unwrap_or_default())
        .into_iter()
        .map(|entry| entry.trim_end_matches('*').trim().to_string())
        .filter(|entry| !entry.is_empty())
        .collect::<Vec<String>>();
    let gpus = parse_gpu_count(gres_raw.as_deref().unwrap_or(""));
    let combined_features = merge_sorted_unique(
        &parse_csv_no_empty(&features_raw.clone().unwrap_or_default()),
        &parse_csv_no_empty(&active_features_raw.unwrap_or_default()),
    );
    let mut cpu_name = derive_cpu_name(&combined_features);
    let gpu_name = derive_gpu_name(gres_raw.as_deref(), &combined_features)
        .or_else(|| derive_gpu_name(None, &combined_features));
    if cpu_name.is_none() {
        cpu_name = infer_andromeda_cpu_name(cpus, memory_mb, gpus, gpu_name.as_deref());
    }
    let node_type = if gpus > 0 { "gpu" } else { "cpu" }.to_string();

    Some(HpcNodeSnapshot {
        node_name,
        partitions,
        state: normalize_scheduler_state(&raw_state),
        raw_state,
        cpus,
        cpu_name,
        memory_mb,
        free_memory_mb,
        gpus,
        gpu_name,
        gres: gres_raw,
        features: combined_features,
        reason,
        node_type,
    })
}

fn parse_queue_line(line: &str) -> Option<HpcQueueSnapshot> {
    let mut parts = line.splitn(10, '|');
    let job_id = sanitize_required(parts.next()?)?;
    let user = sanitize_optional(parts.next()).unwrap_or_else(|| "unknown".to_string());
    let raw_state = sanitize_optional(parts.next()).unwrap_or_else(|| "UNKNOWN".to_string());
    let partition = sanitize_optional(parts.next()).unwrap_or_else(|| "unknown".to_string());
    let nodes = parse_u32(parts.next()).unwrap_or(0);
    let elapsed = sanitize_optional(parts.next()).unwrap_or_else(|| "00:00".to_string());
    let time_limit = sanitize_optional(parts.next()).unwrap_or_else(|| "unknown".to_string());
    let reason = sanitize_optional(parts.next()).unwrap_or_else(|| "unknown".to_string());
    let nodelist = sanitize_optional(parts.next()).unwrap_or_else(|| "unknown".to_string());
    let name = sanitize_optional(parts.next()).unwrap_or_else(|| "(unnamed)".to_string());

    Some(HpcQueueSnapshot {
        job_id,
        user,
        state: normalize_scheduler_state(&raw_state),
        raw_state,
        partition,
        nodes,
        elapsed,
        time_limit,
        reason,
        nodelist,
        name,
    })
}

fn sanitize_required(value: &str) -> Option<String> {
    let cleaned = sanitize_optional(Some(value))?;
    if cleaned.is_empty() {
        None
    } else {
        Some(cleaned)
    }
}

fn sanitize_optional(value: Option<&str>) -> Option<String> {
    let trimmed = value.unwrap_or("").trim();
    if trimmed.is_empty() || trimmed.eq_ignore_ascii_case("(null)") || trimmed == "None" {
        None
    } else {
        Some(trimmed.to_string())
    }
}

fn parse_u32(value: Option<&str>) -> Option<u32> {
    sanitize_optional(value)?.parse::<u32>().ok()
}

fn parse_u64(value: Option<&str>) -> Option<u64> {
    sanitize_optional(value)?.parse::<u64>().ok()
}

fn parse_u64_optional(value: Option<&str>) -> Option<u64> {
    parse_u64(value)
}

fn parse_csv_no_empty(value: &str) -> Vec<String> {
    value
        .split(',')
        .map(|part| part.trim())
        .filter(|part| !part.is_empty())
        .map(|part| part.to_string())
        .collect()
}

fn merge_sorted_unique(left: &[String], right: &[String]) -> Vec<String> {
    let mut set: BTreeSet<String> = left.iter().cloned().collect();
    for value in right {
        set.insert(value.clone());
    }
    set.into_iter().collect()
}

fn parse_gpu_count(gres: &str) -> u32 {
    let trimmed = gres.trim();
    if trimmed.is_empty() || trimmed.eq_ignore_ascii_case("(null)") {
        return 0;
    }

    let mut total: u32 = 0;
    for segment in trimmed.split(',') {
        let base = segment.split('(').next().unwrap_or("").trim();
        if base.is_empty() {
            continue;
        }
        let lower = base.to_ascii_lowercase();
        if !lower.contains("gpu") {
            continue;
        }
        let pieces: Vec<&str> = base.split(':').collect();
        let count = pieces
            .last()
            .and_then(|value| value.trim().parse::<u32>().ok())
            .unwrap_or(1);
        total = total.saturating_add(count);
    }

    if total == 0 && trimmed.to_ascii_lowercase().contains("gpu") {
        1
    } else {
        total
    }
}

fn derive_gpu_name(gres: Option<&str>, features: &[String]) -> Option<String> {
    if let Some(raw_gres) = gres {
        for segment in raw_gres.split(',') {
            let base = segment.split('(').next().unwrap_or("").trim();
            if base.is_empty() {
                continue;
            }
            let parts: Vec<&str> = base.split(':').collect();
            if parts.is_empty() {
                continue;
            }
            let first = parts[0].trim().to_ascii_lowercase();
            if first != "gpu" {
                continue;
            }
            if parts.len() >= 3 {
                return Some(normalize_model_label(parts[1]));
            }
            if parts.len() == 2 {
                let second = parts[1].trim();
                if second.parse::<u32>().is_err() {
                    return Some(normalize_model_label(second));
                }
            }
        }
    }

    for feature in features {
        let lower = feature.to_ascii_lowercase();
        if looks_like_gpu_model(&lower) {
            return Some(normalize_model_label(feature));
        }
    }
    None
}

fn derive_cpu_name(features: &[String]) -> Option<String> {
    for feature in features {
        let lower = feature.to_ascii_lowercase();
        if looks_like_gpu_model(&lower) {
            continue;
        }
        if looks_like_cpu_model(&lower) {
            return Some(normalize_model_label(feature));
        }
    }
    None
}

fn infer_andromeda_cpu_name(
    cpus: u32,
    memory_mb: u64,
    gpus: u32,
    gpu_name: Option<&str>,
) -> Option<String> {
    let memory_gb = memory_mb / 1024;
    let gpu_label = gpu_name.unwrap_or("").to_ascii_lowercase();

    if gpus > 0 {
        if gpu_label.contains("v100") {
            return Some("Intel Xeon Platinum 8260".to_string());
        }
        if gpu_label.contains("a100") {
            return Some("Intel Xeon Platinum 8358/8362".to_string());
        }
        if gpu_label.contains("a10") {
            return Some("Intel Xeon Platinum 8362".to_string());
        }
        if gpu_label.contains("l40s")
            || gpu_label.contains("l4")
            || gpu_label.contains("h100")
            || gpu_label.contains("h200")
        {
            return Some("Intel Xeon Platinum 8562Y".to_string());
        }
        if cpus == 48 {
            return Some("Intel Xeon Platinum 8260".to_string());
        }
        if cpus == 64 && memory_gb >= 1800 {
            return Some("Intel Xeon Platinum 8562Y".to_string());
        }
        if cpus == 64 {
            return Some("Intel Xeon Platinum (GPU host)".to_string());
        }
    } else {
        if cpus == 48 && (170..=220).contains(&memory_gb) {
            return Some("Intel Xeon Platinum 8260".to_string());
        }
        if cpus == 64 && (220..=300).contains(&memory_gb) {
            return Some("Intel Xeon Platinum 8352Y".to_string());
        }
        if cpus == 96 && memory_gb >= 1500 {
            return Some("Intel Xeon Platinum 8568Y+".to_string());
        }
        if cpus == 72 && (400..=700).contains(&memory_gb) {
            return Some("Intel Xeon Platinum 8452Y".to_string());
        }
        if cpus == 96 && (400..=700).contains(&memory_gb) {
            return Some("Intel Xeon Platinum 8568Y+".to_string());
        }
    }

    None
}

fn looks_like_gpu_model(lower: &str) -> bool {
    lower.contains("a100")
        || lower.contains("h100")
        || lower.contains("v100")
        || lower.contains("rtx")
        || lower.contains("tesla")
        || lower.contains("l40")
        || lower.contains("l4")
        || lower.contains("a40")
        || lower.contains("a30")
        || lower.contains("a10")
        || lower.contains("mi250")
        || lower.contains("mi300")
        || lower.contains("gpu")
}

fn looks_like_cpu_model(lower: &str) -> bool {
    lower.contains("xeon")
        || lower.contains("epyc")
        || lower.contains("ryzen")
        || lower.contains("icelake")
        || lower.contains("skylake")
        || lower.contains("cascadelake")
        || lower.contains("sapphirerapids")
        || lower.contains("intel")
        || lower.contains("amd")
        || lower.contains("w12")
        || lower.contains("w-")
        || lower.contains("gold")
        || lower.contains("silver")
        || lower.contains("platinum")
}

fn normalize_model_label(raw: &str) -> String {
    let cleaned = raw
        .trim()
        .trim_matches(|ch: char| !ch.is_ascii_alphanumeric() && ch != '_' && ch != '-');
    if cleaned.is_empty() {
        return raw.trim().to_string();
    }

    cleaned
        .replace('_', " ")
        .split_whitespace()
        .map(pretty_token)
        .collect::<Vec<String>>()
        .join(" ")
}

fn pretty_token(token: &str) -> String {
    let lower = token.to_ascii_lowercase();
    match lower.as_str() {
        "intel" => "Intel".to_string(),
        "amd" => "AMD".to_string(),
        "xeon" => "Xeon".to_string(),
        "epyc" => "EPYC".to_string(),
        "ryzen" => "Ryzen".to_string(),
        "gold" => "Gold".to_string(),
        "silver" => "Silver".to_string(),
        "platinum" => "Platinum".to_string(),
        _ => {
            if token
                .chars()
                .all(|ch| ch.is_ascii_lowercase() || ch.is_ascii_digit() || ch == '-')
            {
                token.to_ascii_uppercase()
            } else {
                token.to_string()
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::{
        is_andromeda_profile, normalize_queue_scope, parse_cluster_snapshot_output,
        parse_gpu_count, QueueScope,
    };
    use crate::hpc::profile::{
        default_cpu_resources, default_gpu_resources, HpcAuthMethod, HpcLauncher, HpcProfile,
        HpcResourceMode,
    };

    fn test_profile(cluster: &str, host: &str) -> HpcProfile {
        HpcProfile {
            id: "profile-1".to_string(),
            name: "Test".to_string(),
            cluster: cluster.to_string(),
            host: host.to_string(),
            port: 22,
            username: "user".to_string(),
            auth_method: HpcAuthMethod::SshKey,
            ssh_key_path: Some("~/.ssh/id_rsa".to_string()),
            remote_qe_bin_dir: "~/qe/bin".to_string(),
            remote_qe_cpu_bin_dir: Some("~/qe/bin".to_string()),
            remote_qe_gpu_bin_dir: Some("~/qe-gpu/bin".to_string()),
            remote_epw_path: None,
            remote_wannier90_path: Some("wannier90.x".to_string()),
            remote_postw90_path: Some("postw90.x".to_string()),
            remote_pseudo_dir: "~/qe/pseudo".to_string(),
            remote_workspace_root: "~/qcortado/work".to_string(),
            remote_project_root: "~/qcortado/projects".to_string(),
            resource_mode: HpcResourceMode::Both,
            launcher: HpcLauncher::Srun,
            launcher_cpu_extra_args: None,
            launcher_gpu_extra_args: None,
            launcher_extra_args: None,
            default_cpu_resources: default_cpu_resources(),
            default_gpu_resources: default_gpu_resources(),
            credential_persisted: false,
            warning_acknowledged: false,
            created_at: String::new(),
            updated_at: String::new(),
        }
    }

    #[test]
    fn parses_gpu_counts_from_various_gres_formats() {
        assert_eq!(parse_gpu_count("gpu:4"), 4);
        assert_eq!(parse_gpu_count("gpu:a100:2(S:0-1)"), 2);
        assert_eq!(parse_gpu_count("gpu"), 1);
        assert_eq!(parse_gpu_count("gpu:2,mps:100"), 2);
        assert_eq!(parse_gpu_count("(null)"), 0);
    }

    #[test]
    fn parses_nodes_and_queue_lines_with_markers() {
        let output = "\
__QCORTADO_CLUSTER_NODES_BEGIN__
g003|short*,medium|alloc|32|257000|128000|gpu:a100:4|a100,avx2|intel_xeon_w12|(null)
g004|short|idle|32|257000|250000|(null)|avx2|intel_xeon_silver|(null)
__QCORTADO_CLUSTER_NODES_END__
__QCORTADO_CLUSTER_QUEUE_BEGIN__
12345|jake|PENDING|short|1|00:00|12:00:00|Priority|(null)|qe-scf
12344|jake|RUNNING|short|1|00:03:12|12:00:00|None|g003|qe-dos
__QCORTADO_CLUSTER_QUEUE_END__";

        let parsed = parse_cluster_snapshot_output(output, true).expect("parse should succeed");
        assert_eq!(parsed.warnings.len(), 0);
        assert_eq!(parsed.nodes.len(), 2);
        assert_eq!(parsed.queue.len(), 2);
        assert_eq!(parsed.nodes[0].node_name, "g003");
        assert_eq!(parsed.nodes[0].gpus, 4);
        assert_eq!(parsed.nodes[0].gpu_name.as_deref(), Some("A100"));
        assert_eq!(parsed.nodes[0].cpu_name.as_deref(), Some("Intel Xeon W12"));
        assert_eq!(parsed.nodes[0].state, "ALLOC");
        assert!(parsed.nodes[0]
            .partitions
            .iter()
            .any(|entry| entry == "short"));
        assert!(parsed.nodes[0]
            .partitions
            .iter()
            .any(|entry| entry == "medium"));
        assert_eq!(parsed.queue[0].job_id, "12345");
        assert_eq!(parsed.queue[0].state, "PENDING");
        assert_eq!(parsed.queue[0].reason, "Priority");
        assert_eq!(parsed.queue[0].nodelist, "unknown");
    }

    #[test]
    fn infers_cpu_name_from_andromeda_hardware_profile_when_features_missing() {
        let output = "\
__QCORTADO_CLUSTER_NODES_BEGIN__
c032|short|idle|64|262144|252000|gpu:a100:4|avx2|(null)|(null)
__QCORTADO_CLUSTER_NODES_END__
__QCORTADO_CLUSTER_QUEUE_BEGIN__
__QCORTADO_CLUSTER_QUEUE_END__";
        let parsed = parse_cluster_snapshot_output(output, true).expect("parse should succeed");
        assert_eq!(parsed.nodes.len(), 1);
        assert_eq!(
            parsed.nodes[0].cpu_name.as_deref(),
            Some("Intel Xeon Platinum 8358/8362")
        );
        assert_eq!(parsed.nodes[0].gpu_name.as_deref(), Some("A100"));
    }

    #[test]
    fn tolerates_missing_queue_fields() {
        let output = "\
__QCORTADO_CLUSTER_NODES_BEGIN__
g010|short|mixed|16|64000|32000|(null)|avx2|intel_xeon_gold|(null)
__QCORTADO_CLUSTER_NODES_END__
__QCORTADO_CLUSTER_QUEUE_BEGIN__
99999|alice|RUNNING|medium|2|00:10:00|2-00:00:00|||my-job
__QCORTADO_CLUSTER_QUEUE_END__";
        let parsed = parse_cluster_snapshot_output(output, true).expect("parse should succeed");
        assert_eq!(parsed.queue.len(), 1);
        assert_eq!(parsed.queue[0].reason, "unknown");
        assert_eq!(parsed.queue[0].nodelist, "unknown");
    }

    #[test]
    fn keeps_partial_snapshot_when_queue_query_fails() {
        let output = "\
__QCORTADO_CLUSTER_NODES_BEGIN__
c001|short|idle|16|64000|62000|(null)|avx2|intel_xeon_platinum|(null)
__QCORTADO_CLUSTER_NODES_END__
__QCORTADO_CLUSTER_QUEUE_BEGIN__
__QCORTADO_QUEUE_ERROR__:squeue failed with permission denied
__QCORTADO_CLUSTER_QUEUE_END__";
        let parsed = parse_cluster_snapshot_output(output, true).expect("parse should succeed");
        assert_eq!(parsed.nodes.len(), 1);
        assert!(parsed.queue.is_empty());
        assert_eq!(parsed.warnings.len(), 1);
    }

    #[test]
    fn normalizes_queue_scope_inputs() {
        assert_eq!(normalize_queue_scope(None), QueueScope::All);
        assert_eq!(
            normalize_queue_scope(Some("all".to_string())),
            QueueScope::All
        );
        assert_eq!(
            normalize_queue_scope(Some("mine".to_string())),
            QueueScope::Mine
        );
    }

    #[test]
    fn andromeda_guard_accepts_cluster_or_host_match() {
        assert!(is_andromeda_profile(&test_profile(
            "andromeda",
            "something"
        )));
        assert!(is_andromeda_profile(&test_profile(
            "custom",
            "andromeda.bc.edu"
        )));
        assert!(!is_andromeda_profile(&test_profile(
            "custom",
            "cluster.example.edu"
        )));
    }
}
