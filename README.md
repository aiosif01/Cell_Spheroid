# Mitosis_O2 — Cell spheroid oxygen-driven growth model

This repository contains a CompuCell3D simulation modeling oxygen-driven growth, fate and mitosis of a single-cell-initiated spheroid. The model demonstrates a two-threshold oxygen-driven phenotype switching (Normoxic → Hypoxic → Necrotic), deterministic growth of normoxic cells, and optional radiotherapy using a linear-quadratic (LQ) kill model.

## Key components

- `mitosis_O2.cc3d` — CompuCell3D project file (entry point for the GUI runner).
- `Simulation/mitosis_O2.py` — Main Python runner that registers steppables and starts the simulation.
- `Simulation/mitosis_O2.xml` — CompuCell3D XML configuration containing the Potts model settings and `UserParameters` used by steppables.
- `Simulation/mitosis_O2Steppables.py` — Python steppables implementing initialization, oxygen-driven fate and growth, mitosis, radiotherapy, compaction, and light analysis/plotting.

## Model overview

- Geometry: 3D Potts model with domain size defined in `mitosis_O2.xml` (default 100x100x100 voxels).

- Cell types:
  - `Medium` (TypeId=0)
  - `Normoxic` (TypeId=1) — actively growing, can divide
  - `Hypoxic` (TypeId=2) — reduced activity, does not grow
  - `Necrotic` (TypeId=3) — dead, shrinks and may be removed after a lifetime

- Oxygen field:
  - DiffusionSolverFE provides a global `Oxygen` field with boundary planes held at concentration 1.0.
  - Cells consume oxygen via Michaelis–Menten uptake (per-type uptake parameters in the XML).
  - Two threshold rules drive phenotype transitions (parameters in `UserParameters`):
    - `O2_Thresh_NormoxicHypoxic` (default 0.15)
    - `O2_Thresh_HypoxicNecrotic` (default 0.05)

- Fate and growth:
  - `O2DrivenFateSteppable` samples oxygen at each cell's center-of-mass and switches types according to thresholds.
  - Normoxic cells grow their `targetVolume` deterministically by `GrowthRateNormoxic` per MCS and update `targetSurface` accordingly.
  - Hypoxic cells do not grow but keep targets unchanged.
  - Necrotic cells shrink at `NecroticShrinkageRate` and are removed probabilistically after `NecroticLifetime` MCS.

- Mitosis:
  - Implemented in `O2MitosisSteppable`: when a cell's `targetVolume` reaches `FinalTargetVolume` (default 63), it may divide with probability `DivProbNormoxic`.
  - Post-mitosis the parent and child targets are split equally and surfaces recomputed.

- Radiotherapy (optional):
  - `RadiotherapySteppable` applies external-beam fractions using a linear-quadratic survival model controlled by `RT_*` parameters in XML (`RT_Enable`, `RT_DoseGy`, `RT_Alpha`, `RT_Beta`, `RT_StartMCS`, `RT_PeriodMCS`, `RT_Fractions`).

- Center compaction:
  - `CenterCompactionSteppable` applies an inward force (`CenterPushStrength`) to compact the spheroid toward the domain center.

- Analysis:
  - `LightAnalysisSteppable` provides simple plots for total volume and cell counts (Normoxic/Hypoxic/Necrotic) and logs oxygen statistics periodically.

## Important parameters (found in `mitosis_O2.xml` `UserParameters`)

- `InitialCellRadius` (voxels) — starting sphere radius for the single seed cell. Default: 1.8
- `FinalTargetVolume` — target volume to trigger mitosis. Default: 63
- `GrowthRateNormoxic` — targetVolume increment per MCS for normoxic cells. Default: 1.95
- `LambdaVolume*`, `LambdaSurface*` — per-type constraint strengths
- `O2_Thresh_NormoxicHypoxic`, `O2_Thresh_HypoxicNecrotic` — oxygen thresholds for phenotype switching
- `NecroticShrinkageRate`, `NecroticLifetime` — necrotic cell dynamics
- `CenterPushStrength` — inward compaction force magnitude
- `RT_*` — radiotherapy controls (enable, dose, timing, alpha/beta)
- `OutputFrequency` — analysis/logging frequency


## How to run

Notes:
- Parameters can be tweaked directly in `mitosis_O2.xml` under `<UserParameters>`; the steppables read those values at start.
- Enable radiotherapy by setting `RT_Enable` to `1` and adjust `RT_*` parameters as needed.

## Suggested experiments

- Vary `InitialCellRadius` and `GrowthRateNormoxic` to observe different spheroid growth rates.
- Toggle `RT_Enable` and test different fractionation schemes using `RT_Fractions` and `RT_PeriodMCS`.
- Adjust `O2_Thresh_*` values to probe sensitivity of necrotic core formation.

## Notes and limitations

- This model uses spherical approximations for surface calculations and simplified per-cell mechanics. It is designed for demonstration and exploratory simulations, not for direct clinical predictions.
- The code assumes a 3D domain and uses COM-based sampling which may be noisy for very small cells; adjust `InitialCellRadius` accordingly.
