# Oxygen Field Implementation in Cell Spheroid Simulation

## Overview

The oxygen field in this Cell Spheroid simulation is a key component that drives cell fate transitions, growth patterns, and spatial organization. This document explains how the oxygen field is implemented and how it controls cellular behavior.

## Implementation Architecture

### 1. XML Configuration (`mitosis_O2.xml`)

The oxygen field is configured through several XML components:

#### A. DiffusionSolverFE Steppable (Lines 69-105)

```xml
<Steppable Type="DiffusionSolverFE">
  <DiffusionField Name="Oxygen">
    <DiffusionData>
      <FieldName>Oxygen</FieldName>
      <DiffusionConstant>0.1</DiffusionConstant>        <!-- Realistic tissue O2 diffusion -->
      <DecayConstant>0.0</DecayConstant>                <!-- No bulk decay, only cell consumption -->
      <InitialConcentration>1.0</InitialConcentration>  <!-- Uniform O2 at start -->
    </DiffusionData>
```

**Key Parameters:**
- **DiffusionConstant (0.1)**: Controls how quickly oxygen spreads through tissue
- **DecayConstant (0.0)**: No spontaneous oxygen decay - only consumed by cells
- **InitialConcentration (1.0)**: Entire domain starts oxygen-rich

#### B. Cell Oxygen Consumption

```xml
<SecretionData>
  <Uptake Type="Normoxic" MaxUptake="0.05" RelativeUptakeRate="0.02"/>
  <Uptake Type="Hypoxic"  MaxUptake="0.03" RelativeUptakeRate="0.015"/>
  <Uptake Type="Necrotic" MaxUptake="0.01" RelativeUptakeRate="0.005"/>
</SecretionData>
```

**Cell-Type Dependent Consumption:**
- **Normoxic cells**: Highest consumption (active metabolism)
- **Hypoxic cells**: Reduced consumption (stressed metabolism)
- **Necrotic cells**: Minimal consumption (dying/dead cells)

#### C. Boundary Conditions

```xml
<BoundaryConditions>
  <!-- Oxygen supply from all boundaries (simulates vasculature) -->
  <Plane Axis="X">
    <ConstantValue PlanePosition="Min" Value="1.0"/>
    <ConstantValue PlanePosition="Max" Value="1.0"/>
  </Plane>
  <!-- Similar for Y and Z axes -->
</BoundaryConditions>
```

This creates **constant oxygen supply** from all 6 faces of the simulation domain, mimicking vascular oxygen delivery.

#### D. Chemotaxis Plugin

```xml
<Plugin Name="Chemotaxis">
  <ChemicalField Name="Oxygen" Source="DiffusionSolverFE">
    <ChemotaxisByType Type="Normoxic" Lambda="10.0"/>  <!-- Moderate O2 seeking -->
    <ChemotaxisByType Type="Hypoxic" Lambda="15.0"/>   <!-- Strong O2 seeking -->
    <!-- Necrotic cells: NO chemotaxis (frozen) -->
  </ChemicalField>
</Plugin>
```

Living cells actively migrate toward higher oxygen concentrations.

### 2. Python Implementation (`mitosis_O2Steppables.py`)

#### A. Oxygen Field Access

```python
# In O2DrivenFateSteppable.step()
oxygen_field = self.field.Oxygen

# Get oxygen level at cell's center of mass
x = max(0, min(int(round(cell.xCOM)), self.dim.x - 1))
y = max(0, min(int(round(cell.yCOM)), self.dim.y - 1))
z = max(0, min(int(round(cell.zCOM)), self.dim.z - 1))
o2_level = float(oxygen_field[x, y, z])
```

**Key Points:**
- Oxygen field accessed as `self.field.Oxygen`
- Oxygen concentration sampled at each cell's center of mass
- Coordinates are clamped to valid lattice bounds

#### B. Oxygen-Driven Cell Fate Transitions

```python
# Critical oxygen thresholds (from XML UserParameters)
O2_Thresh_Normoxic = 0.90   # Hypoxic → Normoxic transition
O2_Thresh_Hypoxic = 0.08    # Normoxic → Hypoxic transition  
O2_Thresh_Necrotic = 0.03   # Any cell → Necrotic transition

# State transition logic
if cell.type == self.NORMOXIC:
    if o2_level < O2_Thresh_Necrotic:
        self._to_necrotic(cell)      # Direct normoxic → necrotic
    elif o2_level < O2_Thresh_Hypoxic:
        self._to_hypoxic(cell)       # Normoxic → hypoxic
elif cell.type == self.HYPOXIC:
    if o2_level > O2_Thresh_Normoxic:
        self._to_normoxic(cell)      # Recovery: hypoxic → normoxic
    elif o2_level < O2_Thresh_Necrotic:
        self._to_necrotic(cell)      # Hypoxic → necrotic
```

#### C. Oxygen-Modulated Growth

```python
# Growth rates depend on cell type and oxygen availability
if cell.type == self.NORMOXIC and o2_level > O2_Thresh_Hypoxic:
    # Fast growth in high oxygen
    dv = GrowthRateNormoxic * GrowthBoostNormoxic  # 15.0 * 2.0 = 30.0
    cell.targetVolume = min(TV_Max, cell.targetVolume + dv)

elif cell.type == self.HYPOXIC and o2_level > O2_Thresh_Necrotic:
    # Slow growth in low oxygen  
    dv = GrowthRateHypoxic * GrowthBoostHypoxic    # 1.0 * 2.0 = 2.0
    max_vol = TV_Max * HypoxicMaxVolumeFrac        # Limited max size
    cell.targetVolume = min(max_vol, cell.targetVolume + dv)
```

#### D. Oxygen-Dependent Division (O2MitosisSteppable)

```python
# Division probability scales with oxygen availability
div_prob = 0.0
if cell.type == self.NORMOXIC:
    if o2_level > 0.30:           # High oxygen
        div_prob = 0.99           # 99% division chance
    elif o2_level > 0.15:         # Medium oxygen  
        div_prob = 0.75           # 75% division chance
    elif o2_level > 0.08:         # Low oxygen
        div_prob = 0.75 * 0.9     # Reduced chance
elif cell.type == self.HYPOXIC and o2_level > 0.07:
    div_prob = 0.05 * 0.4         # Very low hypoxic division
```

#### E. Oxygen Monitoring (LightAnalysisSteppable)

```python
# Sample oxygen at different distances from center
oxygen_field = self.field.Oxygen
center_x, center_y, center_z = self.dim.x // 2, self.dim.y // 2, self.dim.z // 2

o2_center = oxygen_field[center_x, center_y, center_z]
o2_r1 = oxygen_field[min(center_x + r1, self.dim.x - 1), center_y, center_z]  
o2_r2 = oxygen_field[min(center_x + r2, self.dim.x - 1), center_y, center_z]

print(f"O2: center={o2_center:.3f} r={r1}:{o2_r1:.3f} r={r2}:{o2_r2:.3f}")
```

### 3. Emergent Behavior

#### A. Spatial Oxygen Gradients

The implementation creates natural oxygen gradients:

1. **High O2 at boundaries** (constant supply = 1.0)
2. **Decreasing O2 toward center** (due to cell consumption)
3. **Lowest O2 in spheroid core** (maximum distance from oxygen source)

#### B. Cell Fate Zonation

This leads to realistic spheroid architecture:

- **Outer rim**: Normoxic cells (high O2, active growth/division)
- **Middle zone**: Hypoxic cells (low O2, reduced activity)  
- **Inner core**: Necrotic cells (very low O2, cell death)

#### C. Self-Regulating Growth

The oxygen field creates negative feedback:

1. More cells → more O2 consumption → lower O2 levels
2. Lower O2 → reduced growth/division → fewer new cells
3. Necrotic cells → less consumption → O2 recovery possible

### 4. Key Implementation Features

#### A. Early Oxygen Reset Mechanism

```python
# Prevent early oxygen depletion in small spheroids
if mcs < early_window:
    center_o2 = float(oxygen_field[cx, cy, cz])
    if center_o2 < O2_Reset_Threshold:    # 0.5
        oxygen_field[:, :, :] = 1.0       # Reset entire field
        print(f"Oxygen field reset (center={center_o2:.3f})")
```

#### B. Robust Coordinate Handling

```python
# Ensure coordinates stay within lattice bounds
x = max(0, min(int(round(cell.xCOM)), self.dim.x - 1))
y = max(0, min(int(round(cell.yCOM)), self.dim.y - 1))
z = max(0, min(int(round(cell.zCOM)), self.dim.z - 1))
```

#### C. Parameter-Driven Design

All oxygen thresholds and rates are loaded from XML UserParameters, making the model highly configurable without code changes.

## Summary

The oxygen field implementation in this Cell Spheroid simulation creates a realistic, dynamic microenvironment that:

- **Drives cell fate decisions** through oxygen-dependent state transitions
- **Controls growth patterns** via oxygen-modulated growth rates  
- **Regulates cell division** through oxygen-dependent proliferation
- **Enables spatial organization** via chemotaxis and gradient formation
- **Provides system stability** through negative feedback mechanisms

This creates emergent tumor-like behavior with characteristic zones of proliferation, hypoxia, and necrosis that match experimental observations of real spheroid cultures.