import { useCallback, useEffect, useMemo, useState } from "react";
import {
  TransportAxisSettings,
  TransportResult,
  TransportViewerState,
  clampTransportIndex,
  createDefaultTransportAxisSettings,
  formatTransportNumber,
  getDefaultTransportMetric,
  getDefaultTransportMuIndex,
  getDefaultTransportTemperatureIndex,
  getNextTransportSliceSelection,
  getTransportComponentOptions,
} from "../../lib/transport";

function buildInitialState(result: TransportResult): TransportViewerState {
  const metric = getDefaultTransportMetric(result);
  const componentOptions = getTransportComponentOptions(result[metric]);
  const tdfOptions = getTransportComponentOptions(result.tdf ?? null);
  return {
    metric,
    component: componentOptions[0]?.value ?? "value",
    axisMode: "mu",
    selectedTemperatureIndex: getDefaultTransportTemperatureIndex(result.temperature_values_k),
    selectedMuIndex: getDefaultTransportMuIndex(result.mu_offsets_ev),
    heatmapOpen: true,
    settingsOpen: false,
    runContextOpen: false,
    tdfOpen: false,
    tdfComponent: tdfOptions[0]?.value ?? "value",
    axisSettings: createDefaultTransportAxisSettings(),
  };
}

function syncRangeInputs(settings: TransportAxisSettings): TransportAxisSettings {
  if (!settings.customYRange) {
    return {
      ...settings,
      manualYMinInput: "",
      manualYMaxInput: "",
    };
  }

  return {
    ...settings,
    manualYMinInput: formatTransportNumber(settings.customYRange[0], "manual", 6),
    manualYMaxInput: formatTransportNumber(settings.customYRange[1], "manual", 6),
  };
}

export function useTransportViewerState(result: TransportResult) {
  const [state, setState] = useState<TransportViewerState>(() => buildInitialState(result));
  const [rangeError, setRangeError] = useState<string | null>(null);

  const fallbackMetric = useMemo(() => getDefaultTransportMetric(result), [result]);
  const dataset = result[state.metric];
  const componentOptions = useMemo(
    () => getTransportComponentOptions(dataset),
    [dataset],
  );
  const tdfOptions = useMemo(
    () => getTransportComponentOptions(result.tdf ?? null),
    [result.tdf],
  );

  useEffect(() => {
    setState((prev) => {
      let next = prev;
      if (!result[prev.metric]) {
        next = { ...next, metric: fallbackMetric };
      }

      const safeTemperatureIndex = clampTransportIndex(
        next.selectedTemperatureIndex,
        result.temperature_values_k.length,
      );
      const safeMuIndex = clampTransportIndex(next.selectedMuIndex, result.mu_offsets_ev.length);

      if (
        safeTemperatureIndex !== next.selectedTemperatureIndex
        || safeMuIndex !== next.selectedMuIndex
      ) {
        next = {
          ...next,
          selectedTemperatureIndex: safeTemperatureIndex,
          selectedMuIndex: safeMuIndex,
        };
      }

      return next;
    });
  }, [fallbackMetric, result]);

  useEffect(() => {
    if (componentOptions.some((option) => option.value === state.component)) {
      return;
    }
    setState((prev) => ({
      ...prev,
      component: componentOptions[0]?.value ?? "value",
      axisSettings: syncRangeInputs({
        ...prev.axisSettings,
        customYRange: null,
      }),
    }));
  }, [componentOptions, state.component]);

  useEffect(() => {
    if (tdfOptions.some((option) => option.value === state.tdfComponent)) {
      return;
    }
    setState((prev) => ({
      ...prev,
      tdfComponent: tdfOptions[0]?.value ?? "value",
    }));
  }, [state.tdfComponent, tdfOptions]);

  const safeTemperatureIndex = clampTransportIndex(
    state.selectedTemperatureIndex,
    result.temperature_values_k.length,
  );
  const safeMuIndex = clampTransportIndex(state.selectedMuIndex, result.mu_offsets_ev.length);

  const setMetric = useCallback((metric: TransportViewerState["metric"]) => {
    setState((prev) => ({
      ...prev,
      metric,
      component: getTransportComponentOptions(result[metric])[0]?.value ?? "value",
      axisSettings: syncRangeInputs({
        ...prev.axisSettings,
        customYRange: null,
      }),
    }));
    setRangeError(null);
  }, [result]);

  const setComponent = useCallback((component: string) => {
    setState((prev) => ({
      ...prev,
      component,
      axisSettings: syncRangeInputs({
        ...prev.axisSettings,
        customYRange: null,
      }),
    }));
    setRangeError(null);
  }, []);

  const setAxisMode = useCallback((axisMode: TransportViewerState["axisMode"]) => {
    setState((prev) => ({
      ...prev,
      axisMode,
      axisSettings: syncRangeInputs({
        ...prev.axisSettings,
        customYRange: null,
      }),
    }));
    setRangeError(null);
  }, []);

  const setSelectedTemperatureIndex = useCallback((selectedTemperatureIndex: number) => {
    setState((prev) => ({
      ...prev,
      selectedTemperatureIndex: clampTransportIndex(
        selectedTemperatureIndex,
        result.temperature_values_k.length,
      ),
    }));
  }, [result.temperature_values_k.length]);

  const setSelectedMuIndex = useCallback((selectedMuIndex: number) => {
    setState((prev) => ({
      ...prev,
      selectedMuIndex: clampTransportIndex(selectedMuIndex, result.mu_offsets_ev.length),
    }));
  }, [result.mu_offsets_ev.length]);

  const applyHeatmapSelection = useCallback((temperatureIndex: number, muIndex: number) => {
    setState((prev) => ({
      ...prev,
      ...getNextTransportSliceSelection(prev.axisMode, temperatureIndex, muIndex, prev),
    }));
    setRangeError(null);
  }, []);

  const setCustomYRange = useCallback((customYRange: [number, number] | null) => {
    setState((prev) => ({
      ...prev,
      axisSettings: syncRangeInputs({
        ...prev.axisSettings,
        customYRange,
      }),
    }));
    setRangeError(null);
  }, []);

  const updateAxisSettings = useCallback((patch: Partial<TransportAxisSettings>) => {
    setState((prev) => ({
      ...prev,
      axisSettings: {
        ...prev.axisSettings,
        ...patch,
      },
    }));
  }, []);

  const applyManualYRange = useCallback(() => {
    const minValue = Number(state.axisSettings.manualYMinInput.trim());
    const maxValue = Number(state.axisSettings.manualYMaxInput.trim());
    if (!Number.isFinite(minValue) || !Number.isFinite(maxValue)) {
      setRangeError("Enter finite Y-range values.");
      return false;
    }
    if (maxValue <= minValue) {
      setRangeError("Y max must be greater than Y min.");
      return false;
    }
    setCustomYRange([minValue, maxValue]);
    return true;
  }, [setCustomYRange, state.axisSettings.manualYMaxInput, state.axisSettings.manualYMinInput]);

  const resetView = useCallback(() => {
    setState((prev) => ({
      ...prev,
      axisSettings: syncRangeInputs({
        ...prev.axisSettings,
        customYRange: null,
      }),
    }));
    setRangeError(null);
  }, []);

  const resetAxisSettings = useCallback(() => {
    setState((prev) => ({
      ...prev,
      axisSettings: createDefaultTransportAxisSettings(),
    }));
    setRangeError(null);
  }, []);

  const resetAll = useCallback(() => {
    setState(buildInitialState(result));
    setRangeError(null);
  }, [result]);

  const toggleSettingsOpen = useCallback(() => {
    setState((prev) => ({ ...prev, settingsOpen: !prev.settingsOpen }));
  }, []);

  const toggleHeatmapOpen = useCallback(() => {
    setState((prev) => ({ ...prev, heatmapOpen: !prev.heatmapOpen }));
  }, []);

  const toggleRunContextOpen = useCallback(() => {
    setState((prev) => ({ ...prev, runContextOpen: !prev.runContextOpen }));
  }, []);

  const toggleTdfOpen = useCallback(() => {
    setState((prev) => ({ ...prev, tdfOpen: !prev.tdfOpen }));
  }, []);

  const setTdfComponent = useCallback((tdfComponent: string) => {
    setState((prev) => ({ ...prev, tdfComponent }));
  }, []);

  return {
    state,
    dataset,
    componentOptions,
    tdfOptions,
    safeTemperatureIndex,
    safeMuIndex,
    rangeError,
    setMetric,
    setComponent,
    setAxisMode,
    setSelectedTemperatureIndex,
    setSelectedMuIndex,
    applyHeatmapSelection,
    setCustomYRange,
    updateAxisSettings,
    applyManualYRange,
    resetView,
    resetAxisSettings,
    resetAll,
    toggleHeatmapOpen,
    toggleSettingsOpen,
    toggleRunContextOpen,
    toggleTdfOpen,
    setTdfComponent,
  };
}
