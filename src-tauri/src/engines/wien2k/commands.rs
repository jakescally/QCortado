//! Hidden WIEN2k command-plan skeleton.
//!
//! These helpers produce data-only plans for future remote execution. They do
//! not run WIEN2k and are not exposed through Tauri commands.

use super::case_state::{core_case_artifacts, initialized_case_artifacts};
use super::types::{
    Wien2kCaseReference, Wien2kCommandPlan, Wien2kCommandProgram, Wien2kInitializationSettings,
    Wien2kScfRunSettings, Wien2kSpinMode,
};
use super::Wien2kCasePhase;

pub fn build_init_lapw_plan(
    case_ref: &Wien2kCaseReference,
    settings: &Wien2kInitializationSettings,
) -> Wien2kCommandPlan {
    let mut argv = vec![
        "-b".to_string(),
        "-rkmax".to_string(),
        format_float(settings.rkmax),
        "-gmax".to_string(),
        format_float(settings.gmax),
        "-lmax".to_string(),
        settings.lmax.to_string(),
        "-numk".to_string(),
        "0".to_string(),
        settings.k_mesh[0].to_string(),
        settings.k_mesh[1].to_string(),
        settings.k_mesh[2].to_string(),
        "-vxc".to_string(),
        settings.exchange_correlation.to_string(),
        "-ecut".to_string(),
        format_float(settings.lstart_energy_cutoff_ry),
    ];

    if let Some(reduction) = settings.rmt_reduction_percent {
        argv.push("-red".to_string());
        argv.push(format_float(reduction));
    }
    if settings.spin_mode == Wien2kSpinMode::SpinPolarized {
        argv.push("-sp".to_string());
    }

    command_plan(
        Wien2kCommandProgram::InitLapw,
        argv,
        case_ref,
        Wien2kCasePhase::InitializationRunning,
        initialized_case_artifacts(&case_ref.case_name),
    )
}

pub fn build_scf_run_plan(
    case_ref: &Wien2kCaseReference,
    settings: &Wien2kScfRunSettings,
) -> Wien2kCommandPlan {
    let program = match settings.spin_mode {
        Wien2kSpinMode::NonSpinPolarized => Wien2kCommandProgram::RunLapw,
        Wien2kSpinMode::SpinPolarized => Wien2kCommandProgram::RunspLapw,
    };
    let mut argv = vec![
        "-ec".to_string(),
        format_float(settings.energy_convergence_ry),
        "-cc".to_string(),
        format_float(settings.charge_convergence),
        "-i".to_string(),
        settings.max_iterations.to_string(),
    ];
    if let Some(force_cutoff) = settings.force_convergence_mry_bohr {
        argv.push("-fc".to_string());
        argv.push(format_float(force_cutoff));
    }
    if settings.parallel {
        argv.push("-p".to_string());
    }

    command_plan(
        program,
        argv,
        case_ref,
        Wien2kCasePhase::ScfRunning,
        core_case_artifacts(&case_ref.case_name),
    )
}

fn command_plan(
    program: Wien2kCommandProgram,
    argv: Vec<String>,
    case_ref: &Wien2kCaseReference,
    phase: Wien2kCasePhase,
    expected_artifacts: Vec<super::types::Wien2kCaseArtifact>,
) -> Wien2kCommandPlan {
    let mut environment = Vec::new();
    if let Some(scratch_dir) = &case_ref.remote_scratch_dir {
        environment.push(("SCRATCH".to_string(), scratch_dir.clone()));
    }

    Wien2kCommandPlan {
        program,
        argv,
        working_directory: case_ref.remote_case_dir.clone(),
        environment,
        phase,
        expected_artifacts,
    }
}

fn format_float(value: f64) -> String {
    let formatted = format!("{value:.8}");
    formatted
        .trim_end_matches('0')
        .trim_end_matches('.')
        .to_string()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn case_ref() -> Wien2kCaseReference {
        Wien2kCaseReference {
            case_name: "Si".to_string(),
            remote_case_dir: "/scratch/qcortado/projects/project-1/Si".to_string(),
            remote_scratch_dir: Some("/scratch/qcortado/work/Si".to_string()),
            remote_archive_dir: None,
            local_shadow_dir: None,
            project_id: Some("project-1".to_string()),
            cif_id: Some("cif-1".to_string()),
        }
    }

    #[test]
    fn init_plan_is_case_directory_remote_plan() {
        let plan = build_init_lapw_plan(&case_ref(), &Wien2kInitializationSettings::default());

        assert_eq!(plan.program.script_name(), "init_lapw");
        assert_eq!(
            plan.working_directory,
            "/scratch/qcortado/projects/project-1/Si"
        );
        assert!(plan.argv.iter().any(|arg| arg == "-rkmax"));
        assert!(plan
            .expected_artifacts
            .iter()
            .any(|artifact| artifact.basename == "Si.struct"));
        assert!(plan
            .expected_artifacts
            .iter()
            .any(|artifact| artifact.basename == "Si.clmsum"));
    }

    #[test]
    fn spin_polarized_scf_uses_runsp_lapw() {
        let settings = Wien2kScfRunSettings {
            spin_mode: Wien2kSpinMode::SpinPolarized,
            parallel: true,
            ..Wien2kScfRunSettings::default()
        };
        let plan = build_scf_run_plan(&case_ref(), &settings);

        assert_eq!(plan.program.script_name(), "runsp_lapw");
        assert!(plan.argv.iter().any(|arg| arg == "-p"));
        assert_eq!(
            plan.environment,
            vec![(
                "SCRATCH".to_string(),
                "/scratch/qcortado/work/Si".to_string()
            )]
        );
    }
}
