# How the Oxygen Field is Implemented in Cell_Spheroid

## Quick Answer

The oxygen field in this Cell Spheroid simulation is implemented using **CompuCell3D's DiffusionSolverFE** and drives all major cellular behaviors including growth, division, and fate transitions. Here's how it works:

## 🔬 **Core Implementation**

### 1. **Oxygen Field Setup** (in `mitosis_O2.xml`)
```xml
<Steppable Type="DiffusionSolverFE">
  <DiffusionField Name="Oxygen">
    <!-- Creates a 3D oxygen concentration field -->
    <!-- Oxygen diffuses through tissue (DiffusionConstant=0.1) -->
    <!-- Consumed by cells at different rates based on cell type -->
    <!-- Constant supply from all boundaries (Value=1.0) -->
  </DiffusionField>
</Steppable>
```

### 2. **Oxygen Access in Python** (in `mitosis_O2Steppables.py`)
```python
# Access oxygen field
oxygen_field = self.field.Oxygen

# Get oxygen level at any cell's location
x, y, z = int(cell.xCOM), int(cell.yCOM), int(cell.zCOM) 
o2_level = float(oxygen_field[x, y, z])
```

## 🎯 **What Oxygen Controls**

| Oxygen Level | Cell Fate | Growth Rate | Division Chance |
|-------------|-----------|-------------|----------------|
| **> 0.90** | Normoxic | Fast (30/MCS) | High (99%) |
| **0.08-0.90** | Normoxic→Hypoxic | Slowing | Medium (75%) |
| **0.03-0.08** | Hypoxic | Slow (2/MCS) | Low (2%) |
| **< 0.03** | Necrotic | None | None |

## 📍 **Spatial Organization**

The oxygen field creates realistic **tumor-like zonation**:

```
🔵 BOUNDARY (O2=1.0) → Normoxic cells (high growth/division)
    ↓ Oxygen gradients form due to:
    • Constant supply from boundaries  
    • Cell consumption (varies by type)
    • Diffusion through tissue

🟡 PERIPHERY (O2~0.5-0.8) → Mixed Normoxic/Hypoxic
    ↓ Lower oxygen as distance from boundary increases

🟠 MIDDLE (O2~0.1-0.3) → Hypoxic cells (stressed, slow growth)
    ↓ Oxygen further depleted by outer cell consumption

🔴 CORE (O2<0.05) → Necrotic cells (dying, no activity)
```

## ⚙️ **Key Files**

1. **`mitosis_O2.xml`** - Oxygen field configuration
2. **`mitosis_O2Steppables.py`** - Oxygen-driven cell behaviors  
3. **`OXYGEN_IMPLEMENTATION.md`** - Detailed technical documentation
4. **`OXYGEN_PARAMETERS.md`** - Parameter reference guide

## 🧬 **Dynamic Feedback Loop**

1. **Cells consume oxygen** → Creates gradients (high periphery, low center)
2. **Low oxygen** → Cells become hypoxic → Reduced growth/division  
3. **Very low oxygen** → Cells become necrotic → Minimal consumption
4. **Necrotic regions** → Less oxygen consumption → Possible recovery
5. **High oxygen regions** → Active proliferation → More consumption

This creates a **self-regulating system** that naturally forms spheroid architecture!

## 🔄 **Real-Time Monitoring**

The simulation tracks oxygen levels at different distances:
```
[MCS  1000] Cells: N=45  H=12  D=3  | O2: center=0.045 r=10:0.234 r=20:0.678
```
- **center**: Oxygen at spheroid center (often necrotic zone)
- **r=10**: Oxygen 10 voxels from center (hypoxic zone) 
- **r=20**: Oxygen 20 voxels from center (normoxic periphery)

## 📝 **Summary**

The oxygen field is the **central control mechanism** that:
- ✅ **Drives realistic cell fate transitions**
- ✅ **Creates spatial tumor organization** 
- ✅ **Controls growth and division patterns**
- ✅ **Provides automatic system regulation**
- ✅ **Mimics real tumor physiology**

For detailed implementation details, see `OXYGEN_IMPLEMENTATION.md`  
For parameter tuning, see `OXYGEN_PARAMETERS.md`