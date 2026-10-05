//! QE point decomposition for generated batch commands. Use MPI within one
//! allocation: never run independent SCF jobs against the same charge density.
use super::profile::{ResourceValidation, SlurmResourceRequest};
use regex::Regex;

pub fn validate(ranks: u32, pools: u32, images: u32) -> Result<(), String> {
    let groups = pools
        .checked_mul(images)
        .filter(|n| *n > 0)
        .ok_or("QE pools and images must be positive integers without overflow.")?;
    if ranks == 0 || ranks % groups != 0 {
        return Err(format!("QE point parallelization requires MPI tasks ({ranks}) to be divisible by pools × images ({pools} × {images})."));
    }
    Ok(())
}

fn reduced_launcher(line: &str, ranks: u32, nodes: u32, executable: &Regex) -> Option<String> {
    if line.starts_with("srun ") {
        // Place rank overrides after profile launcher arguments so a profile's
        // task/node flags cannot override the collection step's smaller layout.
        Some(
            executable
                .replace(
                    line,
                    format!(
                        "--ntasks={ranks} --nodes=1-{} ${{exe}}${{space}}",
                        nodes.max(1).min(ranks)
                    ),
                )
                .into_owned(),
        )
    } else if line.starts_with("mpirun -np \"${SLURM_NTASKS:-1}\" ") {
        Some(line.replacen(
            "mpirun -np \"${SLURM_NTASKS:-1}\"",
            &format!("mpirun -np {ranks}"),
            1,
        ))
    } else {
        None
    }
}

pub fn commands(
    job_name: &str,
    lines: &[String],
    resources: &SlurmResourceRequest,
    validation: &mut ResourceValidation,
) -> Vec<String> {
    let pools = resources.qe_kpoint_pools.unwrap_or(1);
    let images = resources.qe_qpoint_images.unwrap_or(1);
    if pools == 1 && images == 1 {
        return lines.to_vec();
    }
    let phonon = job_name == "qcortado-phonon";
    let electronic = matches!(job_name, "qcortado-scf" | "qcortado-bands" | "qcortado-dos");
    if !phonon && !electronic || !phonon && images != 1 {
        validation.errors.push("QE point parallelization is supported for SCF, bands, DOS and phonon workflows; images require phonons.".into());
        return lines.to_vec();
    }
    let ranks = resources.ntasks.unwrap_or(1);
    if let Err(error) = validate(ranks, pools, images) {
        validation.errors.push(error);
        return lines.to_vec();
    }
    let executable = if phonon { "ph" } else { "pw" };
    // Match only the executable token, including path and module installations.
    let pattern = Regex::new(&format!(
        r#"(?P<exe>(?:"[^"\n]*/{executable}\.x"|(?:[^\s"']*/)?{executable}\.x))(?P<space>\s)"#
    ))
    .unwrap();
    let mut found = false;
    let mut result = Vec::new();
    for line in lines {
        if !pattern.is_match(line) {
            result.push(line.clone());
            continue;
        }
        found = true;
        if Regex::new(r"(?:^|\s)-(?:nk|npool|npools|ni|nimage|nimages)(?:\s|=)")
            .unwrap()
            .is_match(line)
        {
            validation.errors.push("Remove explicit QE pool/image flags before using the point parallelization controls.".into());
            result.push(line.clone());
            continue;
        }
        let flags = if phonon {
            format!(" -nk {pools} -nimage {images}")
        } else {
            format!(" -nk {pools}")
        };
        result.push(
            pattern
                .replace_all(line, format!("${{exe}}{flags}${{space}}"))
                .into_owned(),
        );
        if phonon && images > 1 {
            // Generated phonon workflows use ph.in/ph.out. QE needs a single
            // image recovery with the SAME ranks per image and pool layout.
            let per_image = ranks / images;
            let mut collect = line
                .replace("-in ph.in", "-in ph_collect.in")
                .replace(" > ph.out", " >> ph.out");
            if collect == *line || !line.contains("-in ph.in") {
                validation
                    .errors
                    .push("Phonon images require the generated ph.in/ph.out workflow.".into());
                continue;
            }
            if let Some(reduced) =
                reduced_launcher(&collect, per_image, resources.nodes.unwrap_or(1), &pattern)
            {
                collect = reduced;
            } else {
                validation.errors.push(
                    "Phonon image collection requires a supported srun or mpirun launcher.".into(),
                );
                continue;
            }
            result.push("# Collect image dynamical matrices before q2r; retain all image scratch on shared storage.".into());
            result.push("awk 'tolower($0) ~ /^[[:space:]]*recover[[:space:]]*=/ {next} {print} tolower($0) ~ /^[[:space:]]*&inputph/ {print \"  recover = .true.,\"}' ph.in > ph_collect.in".into());
            result.push(
                pattern
                    .replace_all(
                        &collect,
                        format!("${{exe}} -nk {pools} -nimage 1${{space}}"),
                    )
                    .into_owned(),
            );
        }
    }
    if !found {
        validation
            .errors
            .push("No supported QE executable was found for point parallelization.".into());
    }
    validation.warnings.push(format!(
        "QE decomposition: {images} image(s), {pools} pool(s) per image, {} MPI ranks per pool.",
        ranks / images / pools
    ));
    result
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn decomposition_rejects_invalid_and_overflowing_groups() {
        assert!(validate(64, 2, 8).is_ok());
        assert!(validate(63, 2, 8).is_err());
        assert!(validate(64, 0, 8).is_err());
        assert!(validate(64, u32::MAX, 2).is_err());
    }
    #[test]
    fn images_collect_before_postprocessing_for_both_launchers() {
        for launcher in ["srun", "mpirun -np \"${SLURM_NTASKS:-1}\""] {
            let resources = SlurmResourceRequest {
                ntasks: Some(64),
                nodes: Some(8),
                qe_kpoint_pools: Some(2),
                qe_qpoint_images: Some(8),
                ..Default::default()
            };
            let mut validation = ResourceValidation::default();
            let lines = vec![
                format!("{launcher} \"$QE_BIN/ph.x\" -in ph.in > ph.out 2>&1"),
                "srun q2r.x -in q2r.in > q2r.out".into(),
            ];
            let result = commands("qcortado-phonon", &lines, &resources, &mut validation);
            assert!(validation.errors.is_empty(), "{:?}", validation.errors);
            assert!(result[0].contains("-nk 2 -nimage 8"));
            assert!(result[3].contains("-nk 2 -nimage 1 -in ph_collect.in >> ph.out"));
            assert!(result[3].contains(if launcher == "srun" {
                "--ntasks=8"
            } else {
                "-np 8"
            }));
            assert_eq!(result[4], lines[1]);
        }
    }
    #[test]
    fn pools_target_pw_only_and_preserve_defaults() {
        let lines = vec![
            "srun pw.x -in nscf.in > nscf.out".into(),
            "srun dos.x -in dos.in > dos.out".into(),
        ];
        let mut validation = ResourceValidation::default();
        assert_eq!(
            commands("qcortado-dos", &lines, &Default::default(), &mut validation),
            lines
        );
        let resources = SlurmResourceRequest {
            ntasks: Some(8),
            qe_kpoint_pools: Some(4),
            ..Default::default()
        };
        let result = commands("qcortado-dos", &lines, &resources, &mut validation);
        assert!(result[0].contains("pw.x -nk 4"));
        assert_eq!(result[1], lines[1]);
    }

    #[test]
    fn generated_image_pipeline_collects_and_stops_on_failure() {
        let directory =
            std::env::temp_dir().join(format!("qcortado-point-test-{}", uuid::Uuid::new_v4()));
        std::fs::create_dir(&directory).unwrap();
        std::fs::write(
            directory.join("ph.in"),
            "phonons\n&INPUTPH\n  recover = .false.,\n  prefix = 'source',\n/\n",
        )
        .unwrap();
        for launcher in ["srun", "mpirun -np \"${SLURM_NTASKS:-1}\""] {
            let resources = SlurmResourceRequest {
                ntasks: Some(64),
                nodes: Some(8),
                qe_kpoint_pools: Some(2),
                qe_qpoint_images: Some(8),
                ..Default::default()
            };
            let mut validation = ResourceValidation::default();
            let lines = vec![
                format!("{launcher} ph.x -in ph.in > ph.out 2>&1"),
                format!("{launcher} q2r.x -in q2r.in > q2r.out"),
            ];
            let result = commands("qcortado-phonon", &lines, &resources, &mut validation);
            assert!(validation.errors.is_empty());
            let mock = r#"
set -euo pipefail
export SLURM_NTASKS=64
srun() { while [[ "$1" == --* ]]; do shift; done; "$@"; }
mpirun() { shift 2; "$@"; }
ph.x() {
  local images=1 input=''
  while (($#)); do
    case "$1" in -nimage) images="$2"; shift 2;; -in) input="$2"; shift 2;; *) shift;; esac
  done
  if [[ "$images" == 1 ]]; then
    [[ "${FAIL_COLLECTION:-0}" == 0 ]] || return 9
    [[ "$(awk '/recover = .true./ {n++} END {print n+0}' "$input")" == 1 ]]
    grep -q "prefix = 'source'" "$input"
    printf 'collected\n'
  else printf 'image work\n'; fi
}
q2r.x() { grep -q collected ph.out; printf 'postprocessed\n'; }
"#;
            let script = format!("{mock}\n{}", result.join("\n"));
            let success = std::process::Command::new("bash")
                .arg("-c")
                .arg(&script)
                .current_dir(&directory)
                .output()
                .unwrap();
            assert!(
                success.status.success(),
                "{}",
                String::from_utf8_lossy(&success.stderr)
            );
            assert_eq!(
                std::fs::read_to_string(directory.join("ph.out")).unwrap(),
                "image work\ncollected\n"
            );
            std::fs::remove_file(directory.join("q2r.out")).unwrap();
            let failure = std::process::Command::new("bash")
                .arg("-c")
                .arg(&script)
                .env("FAIL_COLLECTION", "1")
                .current_dir(&directory)
                .output()
                .unwrap();
            assert!(!failure.status.success());
            assert!(!directory.join("q2r.out").exists());
        }
        std::fs::remove_dir_all(directory).unwrap();
    }
}
