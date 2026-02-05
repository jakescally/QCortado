//! Execution management for Quantum ESPRESSO programs.

use std::path::PathBuf;
use std::process::Stdio;
use tokio::io::{AsyncBufReadExt, BufReader};
use tokio::process::{Child, Command};
use tokio::sync::mpsc;

use super::output::parse_pw_output;
use super::types::QEResult;

/// Errors that can occur during QE execution.
#[derive(Debug, thiserror::Error)]
pub enum RunnerError {
    #[error("QE binary not found: {0}")]
    BinaryNotFound(String),

    #[error("Failed to start process: {0}")]
    ProcessStart(#[from] std::io::Error),

    #[error("Failed to write input: {0}")]
    InputWrite(String),

    #[error("Process failed with exit code: {0}")]
    ProcessFailed(i32),

    #[error("Process was killed")]
    ProcessKilled,
}

/// Configuration for the QE runner.
#[derive(Debug, Clone)]
pub struct QERunner {
    /// Path to the QE bin directory (containing pw.x, bands.x, etc.)
    pub qe_bin_dir: PathBuf,
    /// Number of MPI processes (0 = no MPI)
    pub nprocs: u32,
    /// MPI command (e.g., "mpirun", "mpiexec", "srun")
    pub mpi_command: Option<String>,
}

impl QERunner {
    /// Creates a new QE runner with the given bin directory.
    pub fn new(qe_bin_dir: impl Into<PathBuf>) -> Self {
        Self {
            qe_bin_dir: qe_bin_dir.into(),
            nprocs: 0,
            mpi_command: None,
        }
    }

    /// Sets the number of MPI processes.
    pub fn with_mpi(mut self, nprocs: u32, command: Option<String>) -> Self {
        self.nprocs = nprocs;
        self.mpi_command = command;
        self
    }

    /// Gets the path to a QE executable.
    pub fn get_executable(&self, name: &str) -> PathBuf {
        self.qe_bin_dir.join(name)
    }

    /// Checks if a QE executable exists.
    pub fn executable_exists(&self, name: &str) -> bool {
        self.get_executable(name).exists()
    }

    /// Runs pw.x with the given input and returns the parsed result.
    ///
    /// This is a blocking operation that waits for completion.
    pub async fn run_pw(
        &self,
        input: &str,
        working_dir: &PathBuf,
    ) -> Result<QEResult, RunnerError> {
        let (output, _) = self.run_executable("pw.x", input, working_dir, None).await?;
        Ok(parse_pw_output(&output))
    }

    /// Runs pw.x with streaming output.
    ///
    /// Returns a channel receiver that yields output lines as they're produced.
    pub async fn run_pw_streaming(
        &self,
        input: &str,
        working_dir: &PathBuf,
    ) -> Result<(mpsc::Receiver<String>, tokio::task::JoinHandle<Result<QEResult, RunnerError>>), RunnerError> {
        let (tx, rx) = mpsc::channel(100);

        let (output, _) = self.run_executable("pw.x", input, working_dir, Some(tx.clone())).await?;

        // Parse the final output
        let handle = tokio::spawn(async move {
            Ok(parse_pw_output(&output))
        });

        Ok((rx, handle))
    }

    /// Runs a QE executable with input from stdin.
    ///
    /// Returns (stdout, stderr).
    async fn run_executable(
        &self,
        exe_name: &str,
        input: &str,
        working_dir: &PathBuf,
        output_tx: Option<mpsc::Sender<String>>,
    ) -> Result<(String, String), RunnerError> {
        let exe_path = self.get_executable(exe_name);

        if !exe_path.exists() {
            return Err(RunnerError::BinaryNotFound(exe_name.to_string()));
        }

        // Build the command
        let mut cmd = if self.nprocs > 0 {
            let mpi = self.mpi_command.as_deref().unwrap_or("mpirun");
            let mut c = Command::new(mpi);
            c.args(["-np", &self.nprocs.to_string()]);
            c.arg(&exe_path);
            c
        } else {
            Command::new(&exe_path)
        };

        cmd.current_dir(working_dir)
            .stdin(Stdio::piped())
            .stdout(Stdio::piped())
            .stderr(Stdio::piped());

        let mut child = cmd.spawn()?;

        // Write input to stdin
        if let Some(mut stdin) = child.stdin.take() {
            use tokio::io::AsyncWriteExt;
            stdin.write_all(input.as_bytes()).await
                .map_err(|e| RunnerError::InputWrite(e.to_string()))?;
            // stdin is dropped here, closing the pipe
        }

        // Read output and wait for completion
        if let Some(tx) = output_tx {
            // Streaming mode: send lines as they come
            let stdout = child.stdout.take().expect("stdout not captured");
            let stderr = child.stderr.take().expect("stderr not captured");

            let tx_clone = tx.clone();
            let stdout_handle = tokio::spawn(async move {
                let mut reader = BufReader::new(stdout);
                let mut line = String::new();
                let mut full_output = String::new();

                loop {
                    line.clear();
                    match reader.read_line(&mut line).await {
                        Ok(0) => break, // EOF
                        Ok(_) => {
                            full_output.push_str(&line);
                            let _ = tx_clone.send(line.clone()).await;
                        }
                        Err(_) => break,
                    }
                }
                full_output
            });

            let stderr_handle = tokio::spawn(async move {
                let mut reader = BufReader::new(stderr);
                let mut output = String::new();
                let _ = reader.read_line(&mut output).await;
                output
            });

            let stdout_result = stdout_handle.await.unwrap_or_default();
            let stderr_result = stderr_handle.await.unwrap_or_default();

            // Wait for process to complete in streaming mode
            let status = child.wait().await?;
            if !status.success() {
                if let Some(code) = status.code() {
                    return Err(RunnerError::ProcessFailed(code));
                } else {
                    return Err(RunnerError::ProcessKilled);
                }
            }

            Ok((stdout_result, stderr_result))
        } else {
            // Non-streaming mode: wait for completion (consumes child)
            let output = child.wait_with_output().await?;

            if !output.status.success() {
                if let Some(code) = output.status.code() {
                    return Err(RunnerError::ProcessFailed(code));
                } else {
                    return Err(RunnerError::ProcessKilled);
                }
            }

            Ok((
                String::from_utf8_lossy(&output.stdout).to_string(),
                String::from_utf8_lossy(&output.stderr).to_string(),
            ))
        }
    }

    /// Runs bands.x for band structure post-processing.
    pub async fn run_bands(
        &self,
        input: &str,
        working_dir: &PathBuf,
    ) -> Result<String, RunnerError> {
        let (output, _) = self.run_executable("bands.x", input, working_dir, None).await?;
        Ok(output)
    }

    /// Runs dos.x for DOS calculation.
    pub async fn run_dos(
        &self,
        input: &str,
        working_dir: &PathBuf,
    ) -> Result<String, RunnerError> {
        let (output, _) = self.run_executable("dos.x", input, working_dir, None).await?;
        Ok(output)
    }

    /// Runs projwfc.x for projected DOS calculation.
    pub async fn run_projwfc(
        &self,
        input: &str,
        working_dir: &PathBuf,
    ) -> Result<String, RunnerError> {
        let (output, _) = self.run_executable("projwfc.x", input, working_dir, None).await?;
        Ok(output)
    }

    /// Runs ph.x for phonon calculation.
    pub async fn run_phonon(
        &self,
        input: &str,
        working_dir: &PathBuf,
    ) -> Result<String, RunnerError> {
        let (output, _) = self.run_executable("ph.x", input, working_dir, None).await?;
        Ok(output)
    }
}

/// Represents a running QE calculation that can be monitored/cancelled.
pub struct RunningCalculation {
    child: Child,
    output_rx: Option<mpsc::Receiver<String>>,
}

impl RunningCalculation {
    /// Kills the running calculation.
    pub async fn kill(&mut self) -> Result<(), std::io::Error> {
        self.child.kill().await
    }

    /// Takes the output receiver (can only be called once).
    pub fn take_output_rx(&mut self) -> Option<mpsc::Receiver<String>> {
        self.output_rx.take()
    }
}
