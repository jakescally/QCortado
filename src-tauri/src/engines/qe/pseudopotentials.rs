//! QE pseudopotential, UPF, PSLibrary repair, and SSSP helpers.
//!
//! Tauri command handlers stay in crate::lib so command names remain stable;
//! this module owns the engine-specific implementation.

use base64::{engine::general_purpose::STANDARD as BASE64_STANDARD, Engine as _};
use flate2::read::GzDecoder;
use regex::Regex;
use std::collections::{BTreeSet, HashMap};
use std::io::Read;
use std::path::{Path, PathBuf};

/// Metadata extracted from a pseudopotential header.
#[derive(serde::Serialize, Clone)]
pub struct PseudopotentialMetadata {
    pub filename: String,
    pub supports_soc: bool,
    pub pseudo_type: Option<String>,
    pub relativistic: Option<String>,
    pub cutoff_wfc: Option<f64>,
    pub cutoff_rho: Option<f64>,
    pub cutoff_wfc_source: Option<String>,
    pub cutoff_rho_source: Option<String>,
    #[serde(default)]
    pub available_angular_momenta: Vec<u8>,
    pub available_angular_momenta_source: Option<String>,
    pub max_angular_momentum: Option<u8>,
}

#[derive(serde::Serialize, Clone)]
pub struct PseudopotentialInventoryEntry {
    pub filename: String,
    pub size_bytes: u64,
    pub modified_at_epoch: u64,
}

#[derive(serde::Deserialize)]
struct DjrepoCutoffHint {
    ecut: Option<f64>,
}

#[derive(serde::Deserialize)]
struct DjrepoHints {
    low: Option<DjrepoCutoffHint>,
    normal: Option<DjrepoCutoffHint>,
    high: Option<DjrepoCutoffHint>,
}

#[derive(serde::Deserialize)]
struct DjrepoMetadata {
    hints: Option<DjrepoHints>,
}

fn parse_upf_attr_map(attrs: &str) -> HashMap<String, String> {
    let mut map = HashMap::new();
    let attr_re =
        Regex::new(r#"([A-Za-z_][A-Za-z0-9_-]*)\s*=\s*(?:"([^"]*)"|'([^']*)'|([^\s>]+))"#)
            .expect("valid UPF attribute regex");

    for cap in attr_re.captures_iter(attrs) {
        let key = cap.get(1).map(|m| m.as_str().to_string());
        let value = cap
            .get(2)
            .or_else(|| cap.get(3))
            .or_else(|| cap.get(4))
            .map(|m| m.as_str().trim().to_string());
        if let (Some(key), Some(value)) = (key, value) {
            map.insert(key, value);
        }
    }

    map
}

fn parse_upf_bool(value: Option<&String>) -> bool {
    matches!(
        value.map(|v| v.trim().to_lowercase()),
        Some(ref v) if v == "true" || v == ".true." || v == "t" || v == "1"
    )
}

fn parse_upf_f64(value: Option<&String>) -> Option<f64> {
    value
        .and_then(|raw| raw.trim().parse::<f64>().ok())
        .filter(|value| *value > 0.0)
}

fn capture_xml_value(text: &str, tag: &str) -> Option<String> {
    let pattern = format!(
        r"(?is)<{}\b[^>]*>\s*([^<]+?)\s*</{}>",
        regex::escape(tag),
        regex::escape(tag)
    );
    capture_group(text, &pattern).map(|value| value.trim().to_string())
}

fn capture_group(text: &str, pattern: &str) -> Option<String> {
    let re = Regex::new(pattern).ok()?;
    re.captures(text)
        .and_then(|caps| caps.get(1))
        .map(|m| m.as_str().to_string())
}

fn parse_upf_cutoff(text: &str, label: &str) -> Option<f64> {
    let pattern = format!(
        r"(?im){}\s*:\s*([-+]?\d*\.?\d+(?:[Ee][-+]?\d+)?)",
        regex::escape(label),
    );
    capture_group(text, &pattern)
        .and_then(|value| value.parse::<f64>().ok())
        .filter(|value| *value > 0.0)
}

fn insert_valid_angular_channel(channels: &mut BTreeSet<u8>, raw: &str) {
    if let Ok(value) = raw.trim().parse::<i32>() {
        if (0..=6).contains(&value) {
            channels.insert(value as u8);
        }
    }
}

#[derive(Clone, Copy)]
struct AngularMomentumParseDetails {
    source: Option<&'static str>,
    max_angular_momentum: Option<u8>,
}

fn parse_upf_max_angular_momentum(content: &str) -> Option<u8> {
    capture_group(
        content,
        r#"(?is)<PP_HEADER\b[^>]*\bl_max\s*=\s*["']?([0-9]+)"#,
    )
    .or_else(|| {
        capture_group(
            content,
            r"(?im)^\s*([0-9]+)\s+Max angular momentum component\s*$",
        )
    })
    .and_then(|value| value.trim().parse::<u8>().ok())
    .filter(|value| *value <= 6)
}

fn parse_upf_angular_momentum_details(
    content: &str,
    info_text: &str,
) -> (Vec<u8>, AngularMomentumParseDetails) {
    let mut channels = BTreeSet::new();
    let mut source: Option<&'static str> = None;
    let max_angular_momentum = parse_upf_max_angular_momentum(content);

    if let Ok(tag_re) = Regex::new(r#"(?is)<PP_(?:BETA|CHI|PSWFC|AEWFC)\b([^>]*)>"#) {
        for caps in tag_re.captures_iter(content) {
            let attrs = caps
                .get(1)
                .map(|m| parse_upf_attr_map(m.as_str()))
                .unwrap_or_default();
            if let Some(value) = attrs
                .get("angular_momentum")
                .or_else(|| attrs.get("l"))
                .or_else(|| attrs.get("lll"))
            {
                insert_valid_angular_channel(&mut channels, value);
            }
        }
        if !channels.is_empty() {
            source = Some("upf_tags");
        }
    }

    if channels.is_empty() {
        if let Ok(info_re) = Regex::new(r"(?im)\bl(?:\(\d+\))?\s*=\s*([0-9]+)") {
            for caps in info_re.captures_iter(info_text) {
                if let Some(value) = caps.get(1) {
                    insert_valid_angular_channel(&mut channels, value.as_str());
                }
            }
        }
        if !channels.is_empty() {
            source = Some("upf_info_l_lines");
        }
    }

    if channels.is_empty() {
        if let Ok(orbital_row_re) = Regex::new(
            r"(?i)^\s*(\d+[spdfghi])(?:\s+(\d+))?(?:\s+(\d+))?(?:\s+[-+0-9.eed]+){1,}.*$",
        ) {
            for line in info_text.lines().chain(content.lines().take(240)) {
                if let Some(caps) = orbital_row_re.captures(line) {
                    if let Some(third) = caps.get(3) {
                        insert_valid_angular_channel(&mut channels, third.as_str());
                        continue;
                    }
                    if let Some(second) = caps.get(2) {
                        insert_valid_angular_channel(&mut channels, second.as_str());
                    }
                }
            }
        }
        if !channels.is_empty() {
            source = Some("upf_info_orbital_rows");
        }
    }

    if channels.is_empty() {
        if let Some(max_l) = max_angular_momentum {
            for channel in 0..=max_l {
                channels.insert(channel);
            }
            if !channels.is_empty() {
                source = Some("upf_l_max_fallback");
            }
        }
    }

    (
        channels.into_iter().collect(),
        AngularMomentumParseDetails {
            source,
            max_angular_momentum,
        },
    )
}

#[cfg(test)]
fn parse_upf_available_angular_momenta(content: &str, info_text: &str) -> Vec<u8> {
    parse_upf_angular_momentum_details(content, info_text).0
}

fn parse_djrepo_wavefunction_cutoff_ry(text: &str) -> Option<f64> {
    let metadata: DjrepoMetadata = serde_json::from_str(text).ok()?;
    let hints = metadata.hints?;
    let ecut_ha = hints
        .normal
        .as_ref()
        .and_then(|hint| hint.ecut)
        .or_else(|| hints.high.as_ref().and_then(|hint| hint.ecut))
        .or_else(|| hints.low.as_ref().and_then(|hint| hint.ecut))
        .filter(|value| *value > 0.0)?;
    Some(ecut_ha * 2.0)
}

fn is_hamann_oncv_lone_rho_wavefunction_hint(
    generated: Option<&str>,
    pseudo_type: Option<&str>,
    has_wfc: bool,
    cutoff_wfc: Option<f64>,
    cutoff_rho: Option<f64>,
    info_text: &str,
) -> bool {
    let generated_lower = generated.unwrap_or_default().to_lowercase();
    let info_lower = info_text.to_lowercase();
    let pseudo_type = pseudo_type.unwrap_or_default();

    cutoff_wfc.is_none()
        && cutoff_rho.is_some()
        && !has_wfc
        && pseudo_type.eq_ignore_ascii_case("nc")
        && (generated_lower.contains("oncvpsp code by d. r. hamann")
            || info_lower.contains("oncvpsp")
            || info_lower.contains("d. r. hamann"))
}

pub(crate) fn parse_pseudopotential_metadata_from_content(
    filename: String,
    content: &str,
) -> PseudopotentialMetadata {
    let header_block =
        capture_group(content, r"(?is)<PP_HEADER\b[^>]*>(.*?)</PP_HEADER>").unwrap_or_default();
    let header_attrs = capture_group(content, r"(?is)<PP_HEADER\b([^>]*)>")
        .map(|attrs| parse_upf_attr_map(&attrs))
        .unwrap_or_default();
    let info_text = capture_group(content, r"(?is)<PP_INFO>(.*?)</PP_INFO>")
        .unwrap_or_else(|| content.lines().take(200).collect::<Vec<_>>().join("\n"));
    let info_lower = info_text.to_lowercase();
    let lower_content = content.to_lowercase();
    let header_has_so_xml = capture_xml_value(&header_block, "has_so");
    let header_has_wfc_xml = capture_xml_value(&header_block, "has_wfc");
    let generated = header_attrs
        .get("generated")
        .cloned()
        .or_else(|| capture_xml_value(&header_block, "generated"));

    let pseudo_type = header_attrs
        .get("pseudo_type")
        .cloned()
        .or_else(|| capture_xml_value(&header_block, "pseudo_type"))
        .or_else(|| capture_group(&info_text, r"(?im)Pseudopotential type:\s*([^\r\n<]+)"))
        .map(|value| value.trim().to_string())
        .filter(|value| !value.is_empty());

    let relativistic = header_attrs
        .get("relativistic")
        .cloned()
        .or_else(|| capture_xml_value(&header_block, "relativistic"))
        .or_else(|| {
            capture_group(
                &info_text,
                r"(?im)((?:non-|scalar-|fully-)?relativistic pseudopotential)",
            )
        });
    let (available_angular_momenta, angular_momentum_details) =
        parse_upf_angular_momentum_details(content, &info_text);

    let has_so = parse_upf_bool(header_attrs.get("has_so"))
        || parse_upf_bool(header_has_so_xml.as_ref())
        || lower_content.contains("<pp_spin_orb")
        || info_lower.contains("fully-relativistic")
        || info_lower.contains("spin-orbit");
    let has_wfc =
        parse_upf_bool(header_attrs.get("has_wfc")) || parse_upf_bool(header_has_wfc_xml.as_ref());

    let (raw_cutoff_wfc, raw_cutoff_wfc_source) =
        if let Some(value) = parse_upf_f64(header_attrs.get("wfc_cutoff")) {
            (Some(value), Some("upf"))
        } else if let Some(value) = capture_xml_value(&header_block, "wfc_cutoff")
            .and_then(|value| value.parse::<f64>().ok())
            .filter(|value| *value > 0.0)
        {
            (Some(value), Some("upf"))
        } else if let Some(value) = parse_upf_f64(header_attrs.get("ecutwfc")) {
            (Some(value), Some("upf"))
        } else if let Some(value) =
            parse_upf_cutoff(&info_text, "Suggested minimum cutoff for wavefunctions")
        {
            (Some(value), Some("upf_info"))
        } else {
            (None, None)
        };
    let (raw_cutoff_rho, raw_cutoff_rho_source) =
        if let Some(value) = parse_upf_f64(header_attrs.get("rho_cutoff")) {
            (Some(value), Some("upf"))
        } else if let Some(value) = capture_xml_value(&header_block, "rho_cutoff")
            .and_then(|value| value.parse::<f64>().ok())
            .filter(|value| *value > 0.0)
        {
            (Some(value), Some("upf"))
        } else if let Some(value) = parse_upf_f64(header_attrs.get("ecutrho")) {
            (Some(value), Some("upf"))
        } else if let Some(value) =
            parse_upf_cutoff(&info_text, "Suggested minimum cutoff for charge density")
        {
            (Some(value), Some("upf_info"))
        } else {
            (None, None)
        };
    let lone_rho_is_hartree_wfc = is_hamann_oncv_lone_rho_wavefunction_hint(
        generated.as_deref(),
        pseudo_type.as_deref(),
        has_wfc,
        raw_cutoff_wfc,
        raw_cutoff_rho,
        &info_text,
    );
    let cutoff_wfc = if lone_rho_is_hartree_wfc {
        raw_cutoff_rho.map(|value| value * 2.0)
    } else {
        raw_cutoff_wfc
    };
    let cutoff_wfc_source = if lone_rho_is_hartree_wfc {
        raw_cutoff_rho.map(|_| "upf_fallback".to_string())
    } else {
        raw_cutoff_wfc_source.map(str::to_string)
    };
    let cutoff_rho = if lone_rho_is_hartree_wfc {
        None
    } else {
        raw_cutoff_rho
    };
    let cutoff_rho_source = if lone_rho_is_hartree_wfc {
        None
    } else {
        raw_cutoff_rho_source.map(str::to_string)
    };

    PseudopotentialMetadata {
        filename,
        supports_soc: has_so
            || matches!(relativistic.as_deref(), Some(value) if value.eq_ignore_ascii_case("full")),
        pseudo_type,
        relativistic,
        cutoff_wfc,
        cutoff_rho,
        cutoff_wfc_source,
        cutoff_rho_source,
        available_angular_momenta,
        available_angular_momenta_source: angular_momentum_details.source.map(str::to_string),
        max_angular_momentum: angular_momentum_details.max_angular_momentum,
    }
}

pub(crate) fn parse_pseudopotential_metadata_from_sources(
    filename: String,
    content: &str,
    djrepo_text: Option<&str>,
) -> PseudopotentialMetadata {
    let mut metadata = parse_pseudopotential_metadata_from_content(filename, content);
    if let Some(djrepo_text) = djrepo_text {
        if let Some(cutoff_wfc) = parse_djrepo_wavefunction_cutoff_ry(djrepo_text) {
            metadata.cutoff_wfc = Some(cutoff_wfc);
            metadata.cutoff_rho = None;
            metadata.cutoff_wfc_source = Some("djrepo".to_string());
            metadata.cutoff_rho_source = None;
        }
    }
    metadata
}

#[cfg(test)]
mod upf_angular_momentum_tests {
    use super::parse_upf_available_angular_momenta;

    #[test]
    fn parses_old_style_upf_info_tables() {
        let content = r#"
<PP_INFO>
Generated using Vanderbilt code
nl pn  l   occ               Rcut            Rcut US             E pseu
3S  3  0  2.00      0.00000000000      1.44000000000     -9.60941084200
3P  3  1  6.00      0.00000000000      1.55000000000     -6.69550237500
3D  3  2  8.00      0.00000000000      1.50000000000     -2.08482978100
4S  4  0  0.00      0.00000000000      1.44000000000     -1.59107856400
4P  4  1  0.00      0.00000000000      1.55000000000     -1.10274995100
</PP_INFO>
<PP_HEADER>
    2                  Max angular momentum component
 Wavefunctions         nl  l   occ
                       3S  0  2.00
                       3P  1  6.00
                       3D  2  8.00
</PP_HEADER>
"#;

        let parsed = parse_upf_available_angular_momenta(content, content);
        assert_eq!(parsed, vec![0, 1, 2]);
    }

    #[test]
    fn parses_modern_upf_info_generation_tables() {
        let content = r#"
<UPF version="2.0.1">
  <PP_INFO>
    Valence configuration:
    nl pn  l   occ       Rcut    Rcut US       E pseu
    3S  1  0  2.00      1.600      1.800    -0.794728
    3P  2  1  2.00      1.600      1.800    -0.299965
    Generation configuration:
    3S  1  0  0.00      1.600      1.800     6.000000
    3P  2  1  0.00      1.600      1.800     6.000000
    3D  3  2  0.00      1.600      1.800     0.100000
  </PP_INFO>
  <PP_HEADER l_max="2" />
</UPF>
"#;

        let parsed = parse_upf_available_angular_momenta(content, content);
        assert_eq!(parsed, vec![0, 1, 2]);
    }

    #[test]
    fn falls_back_to_lmax_when_only_header_metadata_exists() {
        let content = r#"<PP_HEADER l_max="1" />"#;
        let parsed = parse_upf_available_angular_momenta(content, "");
        assert_eq!(parsed, vec![0, 1]);
    }

    #[test]
    fn metadata_tracks_confident_vs_fallback_angular_momentum_sources() {
        let si_like = r#"
<UPF version="2.0.1">
  <PP_INFO>
    3S  1  0  2.00      1.600      1.800    -0.794728
    3P  2  1  2.00      1.600      1.800    -0.299965
    3D  3  2  0.00      1.600      1.800     0.100000
  </PP_INFO>
  <PP_HEADER l_max="2" />
  <PP_BETA.1 angular_momentum="0" />
  <PP_BETA.2 angular_momentum="1" />
  <PP_BETA.3 angular_momentum="2" />
</UPF>
"#;
        let si_metadata =
            super::parse_pseudopotential_metadata_from_content("Si.UPF".to_string(), si_like);
        assert_eq!(
            si_metadata.available_angular_momenta_source.as_deref(),
            Some("upf_tags")
        );
        assert_eq!(si_metadata.max_angular_momentum, Some(2));

        let fallback_only = r#"<PP_HEADER l_max="2" />"#;
        let fallback_metadata = super::parse_pseudopotential_metadata_from_content(
            "Fallback.UPF".to_string(),
            fallback_only,
        );
        assert_eq!(
            fallback_metadata
                .available_angular_momenta_source
                .as_deref(),
            Some("upf_l_max_fallback")
        );
        assert_eq!(fallback_metadata.available_angular_momenta, vec![0, 1, 2]);
    }
}

fn decode_remote_metadata_payload(kind: &str, payload: &str) -> Result<String, String> {
    if !matches!(kind, "upf_b64" | "djrepo_b64" | "djrepo_gz_b64") {
        return Ok(payload.to_string());
    }

    let compact_payload = payload.lines().map(str::trim).collect::<String>();
    let decoded = BASE64_STANDARD
        .decode(compact_payload.as_bytes())
        .map_err(|e| format!("Failed to decode remote {} payload: {}", kind, e))?;
    Ok(String::from_utf8_lossy(&decoded).into_owned())
}

pub fn parse_remote_pseudopotential_metadata_output(
    output: &str,
) -> Result<Vec<PseudopotentialMetadata>, String> {
    let mut upf_contents: HashMap<String, String> =
        HashMap::new();
    let mut djrepo_contents: HashMap<String, String> =
        HashMap::new();
    let mut current_kind: Option<String> = None;
    let mut current_filename: Option<String> = None;
    let mut current_content: Vec<String> = Vec::new();

    let flush_current = |kind: &mut Option<String>,
                         filename: &mut Option<String>,
                         content: &mut Vec<String>,
                         upf_contents: &mut HashMap<String, String>,
                         djrepo_contents: &mut HashMap<String, String>|
     -> Result<(), String> {
        if let (Some(kind), Some(name)) = (kind.take(), filename.take()) {
            let joined = content.join("\n");
            let decoded = decode_remote_metadata_payload(&kind, &joined)?;
            match kind.as_str() {
                "upf" | "upf_b64" => {
                    upf_contents.insert(name, decoded);
                }
                "djrepo" | "djrepo_gz" | "djrepo_b64" | "djrepo_gz_b64" => {
                    djrepo_contents.insert(name, decoded);
                }
                _ => {}
            }
        }
        content.clear();
        Ok(())
    };

    for line in output.lines() {
        let trimmed = line.trim();
        if let Some(path) = trimmed.strip_prefix("__QCORTADO_PSEUDO_DIR_MISSING__:") {
            return Err(format!(
                "Remote pseudopotential directory not found: {}",
                path
            ));
        }
        if let Some(rest) = trimmed.strip_prefix("__QCORTADO_REMOTE_METADATA_FILE__:") {
            flush_current(
                &mut current_kind,
                &mut current_filename,
                &mut current_content,
                &mut upf_contents,
                &mut djrepo_contents,
            )?;
            if let Some((kind, name)) = rest.split_once(':') {
                current_kind = Some(kind.to_string());
                current_filename = Some(name.to_string());
            }
            continue;
        }
        if trimmed == "__QCORTADO_REMOTE_METADATA_FILE_END__" {
            flush_current(
                &mut current_kind,
                &mut current_filename,
                &mut current_content,
                &mut upf_contents,
                &mut djrepo_contents,
            )?;
            continue;
        }
        if current_kind.is_some() {
            current_content.push(line.to_string());
        }
    }

    flush_current(
        &mut current_kind,
        &mut current_filename,
        &mut current_content,
        &mut upf_contents,
        &mut djrepo_contents,
    )?;

    let mut pseudos = Vec::new();
    let mut upf_paths = upf_contents.keys().cloned().collect::<Vec<_>>();
    upf_paths.sort();
    for upf_path in upf_paths {
        let Some(content) = upf_contents.get(&upf_path) else {
            continue;
        };
        let stem = upf_path
            .rsplit_once('.')
            .map(|(prefix, _)| prefix)
            .unwrap_or(upf_path.as_str());
        let djrepo_path = format!("{}.djrepo", stem);
        let djrepo_gz_path = format!("{}.djrepo.gz", stem);
        let djrepo_text = djrepo_contents
            .get(&djrepo_path)
            .or_else(|| djrepo_contents.get(&djrepo_gz_path))
            .map(String::as_str);
        pseudos.push(parse_pseudopotential_metadata_from_sources(
            upf_path,
            content,
            djrepo_text,
        ));
    }
    pseudos.sort_by(|a, b| a.filename.cmp(&b.filename));
    Ok(pseudos)
}

fn read_optional_djrepo_companion(path: &Path) -> Result<Option<String>, String> {
    let Some(parent) = path.parent() else {
        return Ok(None);
    };
    let Some(stem) = path.file_stem().and_then(|value| value.to_str()) else {
        return Ok(None);
    };

    let plain_path = parent.join(format!("{}.djrepo", stem));
    if plain_path.exists() {
        return std::fs::read_to_string(&plain_path).map(Some).map_err(|e| {
            format!(
                "Failed to read companion metadata {}: {}",
                plain_path.display(),
                e
            )
        });
    }

    let gz_path = parent.join(format!("{}.djrepo.gz", stem));
    if gz_path.exists() {
        let file = std::fs::File::open(&gz_path).map_err(|e| {
            format!(
                "Failed to open companion metadata {}: {}",
                gz_path.display(),
                e
            )
        })?;
        let mut decoder = GzDecoder::new(file);
        let mut content = String::new();
        decoder.read_to_string(&mut content).map_err(|e| {
            format!(
                "Failed to read companion metadata {}: {}",
                gz_path.display(),
                e
            )
        })?;
        return Ok(Some(content));
    }

    Ok(None)
}

pub(crate) fn read_pseudopotential_metadata(path: &Path) -> Result<PseudopotentialMetadata, String> {
    let filename = path
        .file_name()
        .and_then(|value| value.to_str())
        .ok_or_else(|| format!("Invalid pseudopotential file name: {}", path.display()))?
        .to_string();
    let file = std::fs::File::open(path)
        .map_err(|e| format!("Failed to open pseudopotential {}: {}", filename, e))?;
    let mut buffer = Vec::new();
    file.take(131_072)
        .read_to_end(&mut buffer)
        .map_err(|e| format!("Failed to read pseudopotential {}: {}", filename, e))?;
    let content = String::from_utf8_lossy(&buffer).into_owned();
    let djrepo_content = read_optional_djrepo_companion(path)?;
    Ok(parse_pseudopotential_metadata_from_sources(
        filename,
        &content,
        djrepo_content.as_deref(),
    ))
}

pub fn is_upf_name(name: &str) -> bool {
    name.ends_with(".UPF") || name.ends_with(".upf")
}

#[derive(Debug, Clone, serde::Serialize)]
pub struct PslibraryPseudoRepairResult {
    pub pseudo_dir: String,
    pub scanned: usize,
    pub candidates: usize,
    pub patched: usize,
    pub already_clean: usize,
    pub patched_files: Vec<String>,
    pub clean_files: Vec<String>,
}

fn is_pslibrary_upf_name(name: &str) -> bool {
    is_upf_name(name) && name.to_ascii_lowercase().contains("_psl.")
}

fn repair_pslibrary_pp_chi10_blocks(content: &str, malformed_re: &Regex) -> (String, usize) {
    let replacement_count = malformed_re.find_iter(content).count();
    if replacement_count == 0 {
        return (content.to_string(), 0);
    }

    let repaired = malformed_re
        .replace_all(content, "<PP_CHI.10${1}index=\"10\"${2}>${3}</PP_CHI.10>")
        .to_string();
    (repaired, replacement_count)
}

#[cfg(test)]
mod pseudopotential_metadata_tests {
    use base64::Engine as _;

    use super::{
        parse_djrepo_wavefunction_cutoff_ry, parse_pseudopotential_metadata_from_content,
        parse_pseudopotential_metadata_from_sources, parse_remote_pseudopotential_metadata_output,
        repair_pslibrary_pp_chi10_blocks, Regex, BASE64_STANDARD,
    };

    #[test]
    fn parses_v201_header_attributes() {
        let content = r#"
<UPF version="2.0.1">
  <PP_INFO>
    Suggested minimum cutoff for wavefunctions: 60 Ry
    Suggested minimum cutoff for charge density: 480 Ry
  </PP_INFO>
  <PP_HEADER element="Bi" pseudo_type="US" relativistic="full" has_so="T" wfc_cutoff="60.0" rho_cutoff="480.0"></PP_HEADER>
</UPF>
"#;

        let metadata =
            parse_pseudopotential_metadata_from_content("Bi.rel-pbe.UPF".to_string(), content);
        assert!(metadata.supports_soc);
        assert_eq!(metadata.pseudo_type.as_deref(), Some("US"));
        assert_eq!(metadata.relativistic.as_deref(), Some("full"));
        assert_eq!(metadata.cutoff_wfc, Some(60.0));
        assert_eq!(metadata.cutoff_rho, Some(480.0));
        assert_eq!(metadata.cutoff_wfc_source.as_deref(), Some("upf"));
        assert_eq!(metadata.cutoff_rho_source.as_deref(), Some("upf"));
    }

    #[test]
    fn parses_pslibrary_cutoffs_from_modern_upf_headers() {
        let content = r#"
<UPF version="2.0.1">
  <PP_INFO>
    Generated using "atomic" code by A. Dal Corso  v.7.5
    Author: ADC
    Pseudopotential type: PAW
    Element:  B
    Suggested minimum cutoff for wavefunctions:  43. Ry
    Suggested minimum cutoff for charge density: 325. Ry
    The Pseudo was generated with a Fully-Relativistic Calculation
  </PP_INFO>
  <PP_HEADER generated="Generated using 'atomic' code by A. Dal Corso  v.7.5"
             element=" B"
             pseudo_type="PAW"
             relativistic="full"
             has_so="true"
             has_wfc="true"
             wfc_cutoff="42.557957626222603"
             rho_cutoff="324.58942206554065"/>
</UPF>
"#;

        let metadata = parse_pseudopotential_metadata_from_content(
            "B.rel-pbesol-n-kjpaw_psl.1.0.0.UPF".to_string(),
            content,
        );
        assert!(metadata.supports_soc);
        assert_eq!(metadata.pseudo_type.as_deref(), Some("PAW"));
        assert_eq!(metadata.relativistic.as_deref(), Some("full"));
        assert_eq!(metadata.cutoff_wfc, Some(42.557957626222603));
        assert_eq!(metadata.cutoff_rho, Some(324.58942206554065));
        assert_eq!(metadata.cutoff_wfc_source.as_deref(), Some("upf"));
        assert_eq!(metadata.cutoff_rho_source.as_deref(), Some("upf"));
    }

    #[test]
    fn parses_xml_subtags_for_newer_header_style() {
        let content = r#"
<upf version="3.0.0">
  <pp_header>
    <pseudo_type>PAW</pseudo_type>
    <relativistic>full</relativistic>
    <has_so>true</has_so>
    <wfc_cutoff>45</wfc_cutoff>
    <rho_cutoff>360</rho_cutoff>
  </pp_header>
  <pp_spin_orb />
</upf>
"#;

        let metadata =
            parse_pseudopotential_metadata_from_content("Pt.rel-pbe.UPF".to_string(), content);
        assert!(metadata.supports_soc);
        assert_eq!(metadata.pseudo_type.as_deref(), Some("PAW"));
        assert_eq!(metadata.relativistic.as_deref(), Some("full"));
        assert_eq!(metadata.cutoff_wfc, Some(45.0));
        assert_eq!(metadata.cutoff_rho, Some(360.0));
        assert_eq!(metadata.cutoff_wfc_source.as_deref(), Some("upf"));
        assert_eq!(metadata.cutoff_rho_source.as_deref(), Some("upf"));
    }

    #[test]
    fn interprets_lone_hamann_oncv_rho_cutoff_as_hartree_wavefunction_hint() {
        let content = r#"
<UPF version="2.0.1">
  <PP_INFO>
    This pseudopotential file has been produced using the code
    ONCVPSP  (Optimized Norm-Conservinng Vanderbilt PSeudopotential)
    fully-relativistic version 3.3.0 08/16/2017 by D. R. Hamann
  </PP_INFO>
  <PP_HEADER
    generated="Generated using ONCVPSP code by D. R. Hamann"
    element="Ag"
    pseudo_type="NC"
    relativistic="full"
    has_so="T"
    has_wfc="F"
    rho_cutoff="1.53700000000E+01">
  </PP_HEADER>
</UPF>
"#;

        let metadata = parse_pseudopotential_metadata_from_content("Ag.upf".to_string(), content);
        assert!(metadata.supports_soc);
        assert_eq!(metadata.pseudo_type.as_deref(), Some("NC"));
        assert_eq!(metadata.relativistic.as_deref(), Some("full"));
        let cutoff_wfc = metadata
            .cutoff_wfc
            .expect("expected converted wavefunction cutoff");
        assert!((cutoff_wfc - 30.74).abs() < 1e-9);
        assert_eq!(metadata.cutoff_rho, None);
        assert_eq!(metadata.cutoff_wfc_source.as_deref(), Some("upf_fallback"));
        assert_eq!(metadata.cutoff_rho_source, None);
    }

    #[test]
    fn parses_pseudodojo_djrepo_normal_hint_as_ry_wavefunction_cutoff() {
        let djrepo = r#"
{
  "hints": {
    "high": { "ecut": 47.0 },
    "low": { "ecut": 37.0 },
    "normal": { "ecut": 41.0 }
  }
}
"#;

        let cutoff_wfc = parse_djrepo_wavefunction_cutoff_ry(djrepo)
            .expect("expected djrepo wavefunction cutoff");
        assert!((cutoff_wfc - 82.0).abs() < 1e-9);
    }

    #[test]
    fn prefers_djrepo_cutoffs_over_embedded_upf_cutoff_hints() {
        let content = r#"
<UPF version="2.0.1">
  <PP_INFO>
    This pseudopotential file has been produced using the code
    ONCVPSP  (Optimized Norm-Conservinng Vanderbilt PSeudopotential)
    fully-relativistic version 3.3.0 08/16/2017 by D. R. Hamann
  </PP_INFO>
  <PP_HEADER
    generated="Generated using ONCVPSP code by D. R. Hamann"
    element="Ag"
    pseudo_type="NC"
    relativistic="full"
    has_so="T"
    has_wfc="F"
    rho_cutoff="1.53700000000E+01">
  </PP_HEADER>
</UPF>
"#;
        let djrepo = r#"
{
  "hints": {
    "normal": { "ecut": 41.0 }
  }
}
"#;

        let metadata = parse_pseudopotential_metadata_from_sources(
            "Ag.upf".to_string(),
            content,
            Some(djrepo),
        );
        let cutoff_wfc = metadata
            .cutoff_wfc
            .expect("expected djrepo wavefunction cutoff");
        assert!((cutoff_wfc - 82.0).abs() < 1e-9);
        assert_eq!(metadata.cutoff_rho, None);
        assert_eq!(metadata.cutoff_wfc_source.as_deref(), Some("djrepo"));
        assert_eq!(metadata.cutoff_rho_source, None);
    }

    #[test]
    fn pairs_remote_djrepo_companions_even_with_banner_noise() {
        let output = r#"
Welcome to the cluster.
Authorized use only.
__QCORTADO_REMOTE_METADATA_FILE__:djrepo:Ag.djrepo
{
  "hints": {
    "normal": { "ecut": 41.0 }
  }
}
__QCORTADO_REMOTE_METADATA_FILE_END__
__QCORTADO_REMOTE_METADATA_FILE__:upf:Ag.upf
<UPF version="2.0.1">
  <PP_INFO>
    This pseudopotential file has been produced using the code
    ONCVPSP  (Optimized Norm-Conservinng Vanderbilt PSeudopotential)
    fully-relativistic version 3.3.0 08/16/2017 by D. R. Hamann
  </PP_INFO>
  <PP_HEADER
    generated="Generated using ONCVPSP code by D. R. Hamann"
    element="Ag"
    pseudo_type="NC"
    relativistic="full"
    has_so="T"
    has_wfc="F"
    rho_cutoff="1.53700000000E+01">
  </PP_HEADER>
</UPF>
__QCORTADO_REMOTE_METADATA_FILE_END__
"#;

        let pseudos =
            parse_remote_pseudopotential_metadata_output(output).expect("expected parsed metadata");
        assert_eq!(pseudos.len(), 1);
        assert_eq!(pseudos[0].filename, "Ag.upf");
        assert!((pseudos[0].cutoff_wfc.expect("expected djrepo cutoff") - 82.0).abs() < 1e-9);
        assert_eq!(pseudos[0].cutoff_rho, None);
        assert_eq!(pseudos[0].cutoff_wfc_source.as_deref(), Some("djrepo"));
        assert_eq!(pseudos[0].cutoff_rho_source, None);
    }

    #[test]
    fn parses_base64_remote_metadata_payloads_for_djrepo_companions() {
        let djrepo = r#"{
  "hints": {
    "normal": { "ecut": 41.0 }
  }
}"#;
        let upf = r#"<UPF version="2.0.1">
  <PP_INFO>
    This pseudopotential file has been produced using the code
    ONCVPSP  (Optimized Norm-Conservinng Vanderbilt PSeudopotential)
    fully-relativistic version 3.3.0 08/16/2017 by D. R. Hamann
  </PP_INFO>
  <PP_HEADER
    generated="Generated using ONCVPSP code by D. R. Hamann"
    element="Ag"
    pseudo_type="NC"
    relativistic="full"
    has_so="T"
    has_wfc="F"
    rho_cutoff="1.53700000000E+01">
  </PP_HEADER>
</UPF>"#;
        let output = format!(
            "Welcome to the cluster.\n__QCORTADO_REMOTE_METADATA_FILE__:djrepo_b64:Ag.djrepo\n{}\n__QCORTADO_REMOTE_METADATA_FILE_END__\n__QCORTADO_REMOTE_METADATA_FILE__:upf_b64:Ag.upf\n{}\n__QCORTADO_REMOTE_METADATA_FILE_END__\n",
            BASE64_STANDARD.encode(djrepo.as_bytes()),
            BASE64_STANDARD.encode(upf.as_bytes()),
        );

        let pseudos = parse_remote_pseudopotential_metadata_output(&output)
            .expect("expected parsed metadata");
        assert_eq!(pseudos.len(), 1);
        assert_eq!(pseudos[0].filename, "Ag.upf");
        assert!((pseudos[0].cutoff_wfc.expect("expected djrepo cutoff") - 82.0).abs() < 1e-9);
        assert_eq!(pseudos[0].cutoff_rho, None);
        assert_eq!(pseudos[0].cutoff_wfc_source.as_deref(), Some("djrepo"));
        assert_eq!(pseudos[0].cutoff_rho_source, None);
    }

    #[test]
    fn repairs_only_malformed_tenth_pslibrary_atomic_wfc_tag() {
        let malformed_re =
            Regex::new(r#"(?s)<PP_CHI\.1\b([^>]*)\bindex="10"([^>]*)>(.*?)</PP_CHI\.1>"#)
                .expect("valid regex");
        let content = r#"
<PP_PSWFC>
<PP_CHI.1 type="real" size="2" columns="4" index="1" label="5S" l="0">
 0.0
</PP_CHI.1>
<PP_CHI.1 type="real" size="2" columns="4" index="10" label="4F" l="3">
 1.0
</PP_CHI.1>
</PP_PSWFC>
"#;

        let (repaired, replacements) = repair_pslibrary_pp_chi10_blocks(content, &malformed_re);

        assert_eq!(replacements, 1);
        assert!(repaired.contains(r#"<PP_CHI.1 type="real" size="2" columns="4" index="1""#));
        assert!(repaired.contains(r#"<PP_CHI.10 type="real" size="2" columns="4" index="10""#));
        assert!(repaired.contains("</PP_CHI.10>"));
    }
}

/// Lists available pseudopotentials in a directory.
pub fn list_pseudopotentials(pseudo_dir: String) -> Result<Vec<String>, String> {
    let path = PathBuf::from(&pseudo_dir);
    if !path.exists() {
        return Err(format!("Directory not found: {}", pseudo_dir));
    }

    let mut pseudos = Vec::new();
    for entry in std::fs::read_dir(&path).map_err(|e| e.to_string())? {
        let entry = entry.map_err(|e| e.to_string())?;
        if !entry.path().is_file() {
            continue;
        }

        let name = entry.file_name().to_string_lossy().to_string();
        if is_upf_name(&name) {
            pseudos.push(name);
        }
    }
    pseudos.sort();
    Ok(pseudos)
}

/// Lists local pseudopotentials with cheap change-detection metadata.
pub fn list_pseudopotential_inventory(
    pseudo_dir: String,
) -> Result<Vec<PseudopotentialInventoryEntry>, String> {
    let path = PathBuf::from(&pseudo_dir);
    if !path.exists() {
        return Err(format!("Directory not found: {}", pseudo_dir));
    }

    let mut entries = Vec::new();
    for entry in std::fs::read_dir(&path).map_err(|e| e.to_string())? {
        let entry = entry.map_err(|e| e.to_string())?;
        let file_path = entry.path();
        if !file_path.is_file() {
            continue;
        }

        let filename = entry.file_name().to_string_lossy().to_string();
        if !is_upf_name(&filename) {
            continue;
        }
        let metadata = entry
            .metadata()
            .map_err(|e| format!("Failed to inspect pseudopotential {}: {}", filename, e))?;
        let modified_at_epoch = metadata
            .modified()
            .ok()
            .and_then(|time| time.duration_since(std::time::UNIX_EPOCH).ok())
            .map(|duration| duration.as_secs())
            .unwrap_or(0);
        entries.push(PseudopotentialInventoryEntry {
            filename,
            size_bytes: metadata.len(),
            modified_at_epoch,
        });
    }
    entries.sort_by(|left, right| left.filename.cmp(&right.filename));
    Ok(entries)
}

pub fn repair_local_pslibrary_pseudopotentials_sync(
    pseudo_dir: String,
) -> Result<PslibraryPseudoRepairResult, String> {
    let path = PathBuf::from(&pseudo_dir);
    if !path.exists() {
        return Err(format!("Directory not found: {}", pseudo_dir));
    }

    let malformed_re =
        Regex::new(r#"(?s)<PP_CHI\.1\b([^>]*)\bindex="10"([^>]*)>(.*?)</PP_CHI\.1>"#)
            .map_err(|e| format!("Failed to build pseudopotential repair matcher: {}", e))?;
    let clean_re = Regex::new(r#"<PP_CHI\.10\b[^>]*\bindex="10""#)
        .map_err(|e| format!("Failed to build pseudopotential repair matcher: {}", e))?;

    let mut result = PslibraryPseudoRepairResult {
        pseudo_dir,
        scanned: 0,
        candidates: 0,
        patched: 0,
        already_clean: 0,
        patched_files: Vec::new(),
        clean_files: Vec::new(),
    };

    let mut entries = std::fs::read_dir(&path)
        .map_err(|e| format!("Failed to read pseudo directory {}: {}", path.display(), e))?
        .collect::<Result<Vec<_>, _>>()
        .map_err(|e| e.to_string())?;
    entries.sort_by_key(|entry| entry.file_name());

    for entry in entries {
        let file_path = entry.path();
        if !file_path.is_file() {
            continue;
        }

        let filename = entry.file_name().to_string_lossy().to_string();
        if !is_upf_name(&filename) {
            continue;
        }
        result.scanned += 1;
        if !is_pslibrary_upf_name(&filename) {
            continue;
        }

        let content = std::fs::read_to_string(&file_path)
            .map_err(|e| format!("Failed to read pseudopotential {}: {}", filename, e))?;
        let (repaired, replacements) = repair_pslibrary_pp_chi10_blocks(&content, &malformed_re);
        if replacements > 0 {
            std::fs::write(&file_path, repaired).map_err(|e| {
                format!(
                    "Failed to write repaired pseudopotential {}: {}",
                    filename, e
                )
            })?;
            result.candidates += 1;
            result.patched += 1;
            result.patched_files.push(filename);
        } else if clean_re.is_match(&content) {
            result.candidates += 1;
            result.already_clean += 1;
            result.clean_files.push(filename);
        }
    }

    Ok(result)
}

/// Lists pseudopotentials and parses SOC/cutoff metadata from their headers.
pub fn list_pseudopotential_metadata_sync(
    pseudo_dir: String,
) -> Result<Vec<PseudopotentialMetadata>, String> {
    let path = PathBuf::from(&pseudo_dir);
    if !path.exists() {
        return Err(format!("Directory not found: {}", pseudo_dir));
    }

    let mut pseudos = Vec::new();
    for entry in std::fs::read_dir(&path).map_err(|e| e.to_string())? {
        let entry = entry.map_err(|e| e.to_string())?;
        let file_path = entry.path();
        if !file_path.is_file() {
            continue;
        }

        let name = entry.file_name().to_string_lossy().to_string();
        if is_upf_name(&name) {
            pseudos.push(read_pseudopotential_metadata(&file_path)?);
        }
    }
    pseudos.sort_by(|a, b| a.filename.cmp(&b.filename));
    Ok(pseudos)
}

/// Parses SOC/cutoff metadata for one local pseudopotential.
pub fn get_pseudopotential_metadata(
    pseudo_dir: String,
    filename: String,
) -> Result<PseudopotentialMetadata, String> {
    if filename.contains('/') || filename.contains('\\') || !is_upf_name(&filename) {
        return Err(format!("Invalid pseudopotential file name: {}", filename));
    }
    let path = PathBuf::from(&pseudo_dir);
    if !path.exists() {
        return Err(format!("Directory not found: {}", pseudo_dir));
    }
    let file_path = path.join(&filename);
    if !file_path.is_file() {
        return Err(format!("Pseudopotential not found: {}", filename));
    }
    read_pseudopotential_metadata(&file_path)
}

/// SSSP element data from JSON
#[derive(serde::Serialize, serde::Deserialize, Clone)]
pub struct SSSPElementData {
    pub filename: String,
    pub md5: Option<String>,
    pub pseudopotential: Option<String>,
    pub cutoff_wfc: f64,
    pub cutoff_rho: f64,
}

/// Loads SSSP JSON data from the pseudo directory.
/// Looks for any file matching SSSP*.json pattern.
pub fn load_sssp_data_sync(
    pseudo_dir: String,
) -> Result<HashMap<String, SSSPElementData>, String> {
    let path = PathBuf::from(&pseudo_dir);
    if !path.exists() {
        return Err(format!("Directory not found: {}", pseudo_dir));
    }

    // Find SSSP JSON file
    let mut sssp_file: Option<PathBuf> = None;
    for entry in std::fs::read_dir(&path).map_err(|e| e.to_string())? {
        let entry = entry.map_err(|e| e.to_string())?;
        let name = entry.file_name().to_string_lossy().to_string();
        if name.starts_with("SSSP") && name.ends_with(".json") {
            sssp_file = Some(entry.path());
            break;
        }
    }

    let sssp_path = sssp_file.ok_or("No SSSP JSON file found in pseudo directory")?;
    let content = std::fs::read_to_string(&sssp_path)
        .map_err(|e| format!("Failed to read SSSP file: {}", e))?;

    let data: HashMap<String, SSSPElementData> =
        serde_json::from_str(&content).map_err(|e| format!("Failed to parse SSSP JSON: {}", e))?;

    Ok(data)
}
