# Oxygen Field Parameters Reference

## Key Oxygen Parameters in mitosis_O2.xml

### Diffusion Parameters
| Parameter | Value | Description |
|-----------|-------|-------------|
| `DiffusionConstant` | 0.1 | Rate of oxygen diffusion through tissue |
| `DecayConstant` | 0.0 | Spontaneous oxygen decay (disabled - only cell consumption) |
| `InitialConcentration` | 1.0 | Starting oxygen level throughout domain |

### Cell Oxygen Consumption Rates
| Cell Type | MaxUptake | RelativeUptakeRate | Description |
|-----------|-----------|-------------------|-------------|
| Normoxic | 0.05 | 0.02 | Active metabolism, high O2 consumption |
| Hypoxic | 0.03 | 0.015 | Reduced metabolism, lower O2 consumption |
| Necrotic | 0.01 | 0.005 | Minimal metabolism, very low consumption |

### Boundary Conditions
- **All faces**: Constant oxygen supply = 1.0 (simulates vascular delivery)

### Chemotaxis (Oxygen-Seeking Movement)
| Cell Type | Lambda | Behavior |
|-----------|--------|----------|
| Normoxic | 10.0 | Moderate movement toward higher oxygen |
| Hypoxic | 15.0 | Strong movement toward higher oxygen |
| Necrotic | 0.0 | No movement (frozen) |

## Oxygen Thresholds (UserParameters)

### Cell Fate Transition Thresholds
| Parameter | Value | Transition |
|-----------|-------|------------|
| `O2_Thresh_Normoxic` | 0.90 | Hypoxic → Normoxic (recovery) |
| `O2_Thresh_Hypoxic` | 0.08 | Normoxic → Hypoxic (stress) |
| `O2_Thresh_Necrotic` | 0.03 | Any living cell → Necrotic (death) |

### Early Oxygen Reset
| Parameter | Value | Purpose |
|-----------|-------|---------|
| `O2_Reset_Threshold` | 0.5 | Reset O2 field if center drops below this early in simulation |
| `EarlyO2CheckWindowFactor` | 0.2 | Fraction of OutputFrequency defining early reset window |

### Division Probability Thresholds
| Parameter | Value | Division Probability | Cell Type |
|-----------|-------|---------------------|-----------|
| `DivO2High` | 0.30 | 99% (`DivProbHighO2`) | Normoxic |
| `DivO2Med` | 0.15 | 75% (`DivProbMedO2`) | Normoxic |
| `DivO2Lower` | 0.10 | 67.5% (75% × 0.9) | Normoxic |
| `DivO2LowEdge` | 0.08 | 60% (75% × 0.8) | Normoxic |
| `DivHypoxicMin` | 0.07 | 2% (5% × 0.4) | Hypoxic |

## Growth Rate Parameters

### Base Growth Rates (per MCS)
| Parameter | Value | Cell Type | Conditions |
|-----------|-------|-----------|------------|
| `GrowthRateNormoxic` | 15.0 | Normoxic | When O2 > 0.08 |
| `GrowthRateHypoxic` | 1.0 | Hypoxic | When O2 > 0.03 |

### Growth Multipliers
| Parameter | Value | Effect |
|-----------|-------|--------|
| `GrowthBoostNormoxic` | 2.0 | Normoxic growth: 15.0 × 2.0 = 30.0/MCS |
| `GrowthBoostHypoxic` | 2.0 | Hypoxic growth: 1.0 × 2.0 = 2.0/MCS |
| `ForcedGrowthIncrement` | 1.0 | Additional growth for all living cells |
| `PostDivisionGrowthFactor` | 0.3 | Post-division boost factor |

### Volume Constraints
| Parameter | Value | Purpose |
|-----------|-------|---------|
| `TV_Max` | 80.0 | Maximum target volume for any cell |
| `HypoxicMaxVolumeFrac` | 0.8 | Hypoxic cells limited to 80% of TV_Max (64.0) |

## Analysis Parameters

### Oxygen Sampling Radii
| Parameter | Value | Purpose |
|-----------|-------|---------|
| `AnalysisRadius1` | 10 | Inner radius for oxygen gradient monitoring |
| `AnalysisRadius2` | 20 | Outer radius for oxygen gradient monitoring |

## Example Oxygen Gradient Profile

In a typical simulation, you might observe:

| Location | Distance from Center | Expected O2 Level | Cell Type Expected |
|----------|---------------------|-------------------|-------------------|
| Boundary | 0 (at edge) | 1.0 | Normoxic |
| Periphery | 5-15 voxels | 0.8-0.9 | Normoxic |
| Mid-layer | 15-25 voxels | 0.3-0.6 | Mixed Normoxic/Hypoxic |
| Inner region | 25-35 voxels | 0.05-0.15 | Hypoxic |
| Core | >35 voxels | <0.05 | Necrotic |

## Modifying Oxygen Behavior

### To increase oxygen penetration:
- Increase `DiffusionConstant` (>0.1)
- Decrease cell uptake rates
- Increase boundary oxygen supply (>1.0)

### To decrease oxygen penetration:
- Decrease `DiffusionConstant` (<0.1)
- Increase cell uptake rates
- Add bulk decay (`DecayConstant` > 0.0)

### To change cell sensitivity:
- Adjust threshold values (`O2_Thresh_*`)
- Modify growth rate parameters
- Change division probability thresholds