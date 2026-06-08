// QE-specific HPC adapter exports.
//
// The implementations still live in the legacy HPC config module while profile
// serialization remains unchanged. Importing through this file makes QE path,
// pseudopotential, and executable-command ownership explicit at call sites.

export {
  buildHpcQeRuntimeSetupLines,
  buildHpcQeInputCommandLine,
  getRemotePseudopotentialMetadata,
  listRemotePseudopotentialInventory,
  listRemotePseudopotentialMetadata,
  listRemotePseudopotentials,
  loadRemoteSsspData,
  qeEngineUsesModuleMode,
  resolveProfileRemoteQeAuxiliaryExecutable,
  resolveProfileRemotePseudoDir,
  resolveProfileRemoteQeBinDir,
} from "../../hpcConfig";
