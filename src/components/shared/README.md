# Shared Components

This folder marks platform and viewer ownership for engine-neutral UI. Existing components remain in `src/components` for now; `index.ts` provides an explicit shared boundary for future imports without moving files.

## Current Shared Owners

- Project and storage shell: `ProjectBrowser.tsx`, `ProjectDashboard.tsx`, `StorageManagerPage.tsx`, project dialogs.
- Task and execution shell: `TaskManagerDrawer.tsx`, `LiveOutputPanel.tsx`, `ProgressBar.tsx`, timers.
- HPC shell: `HpcActivityPanel.tsx`, `HpcNodeActivityPage.tsx`, `HpcRunSettings.tsx`, `RemoteUtilizationPanel.tsx`, `HpcSetupWizard.tsx`.
- Structure viewers: `UnitCellViewer.tsx`, `BrillouinZoneViewer.tsx`, `MagnetismViewer.tsx`.
- Result viewers: `BandPlot.tsx`, `ElectronicDOSPlot.tsx`, `PhononPlot.tsx`, `TransportPlot.tsx`, `BandsMultiview.tsx`.
- Small shared UI: `InfoTooltip.tsx`.

Some shared candidates still accept legacy QE-shaped payloads, such as `BandData`. Keep those payloads stable until adapter phases switch viewers to normalized datasets.

## Rules

- Shared components should not create or validate engine-native input configs.
- Viewer components should prefer normalized datasets once adapters exist.
- HPC shell components can remain here even when they still expose QE path fields; engine-specific profile details should be extracted in later PRs.
- Do not rename user-facing QE labels during organization-only work.
