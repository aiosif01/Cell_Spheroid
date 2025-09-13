from cc3d.core.PySteppables import *
import random, os, math
from xml.dom import minidom

# Pre-computed constant for spherical surface calculations
SPHERE_SURF_COEFF = (36.0 * math.pi) ** (1.0 / 3.0)

def surface_from_volume(vol: float) -> float:
    """Return spherical surface estimate for a given volume."""
    return SPHERE_SURF_COEFF * (vol ** (2.0 / 3.0))

# ------------------------- PARAMETER HANDLING ---------------------------- #
def _load_user_parameters(xml_filename: str) -> dict:
    """Parses <UserParameters><Param Name= Value= /></UserParameters> into a dict."""
    xml_path = xml_filename
    if not os.path.exists(xml_path):
        xml_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), xml_filename)
    if not os.path.exists(xml_path):
        raise FileNotFoundError(f"Cannot find XML file '{xml_filename}' relative to simulation script or steppable directory")
    doc = minidom.parse(xml_path)
    user_nodes = doc.getElementsByTagName('UserParameters')
    if not user_nodes:
        raise RuntimeError('UserParameters block missing in XML')
    params = {}
    for node in user_nodes[0].getElementsByTagName('Param'):
        name = node.getAttribute('Name')
        val = node.getAttribute('Value')
        if not name:
            continue
        try:
            params[name] = float(val) if ('.' in val or 'e' in val.lower()) else int(val)
        except ValueError:
            raise RuntimeError(f"Parameter {name} has non-numeric value '{val}'")
    return params

PARAMS = _load_user_parameters('mitosis_O2.xml')
print(f"[PARAM] Loaded {len(PARAMS)} parameters from XML")

def P(name: str):
    if name not in PARAMS:
        raise KeyError(f"Missing required parameter '{name}' in XML <UserParameters>")
    return PARAMS[name]


# ------------------------- OXYGEN INITIALIZATION ---------------------------- #
class OxygenInitSteppable(SteppableBasePy):
    def start(self):
        # Ensure entire oxygen field starts at 1.0
        self.field.Oxygen[:, :, :] = 1.0
        print(f"[OXYGEN] Initialized entire domain to concentration 1.0")


# ------------------------- SINGLE CELL INITIALIZATION ---------------------------- #
class SingleCellInitSteppable(SteppableBasePy):
    def start(self):
        # Create exactly one cell at center
        cx, cy, cz = self.dim.x // 2, self.dim.y // 2, self.dim.z // 2
        initial_radius = P('InitialCellRadius')
        cell = self.new_cell(self.NORMOXIC)
        
        # Fill spherical region
        radius_int = int(initial_radius) + 1
        for x in range(cx - radius_int, cx + radius_int + 1):
            for y in range(cy - radius_int, cy + radius_int + 1):
                for z in range(cz - radius_int, cz + radius_int + 1):
                    if 0 <= x < self.dim.x and 0 <= y < self.dim.y and 0 <= z < self.dim.z:
                        if (x-cx)**2 + (y-cy)**2 + (z-cz)**2 <= initial_radius*initial_radius:
                            self.cell_field[x, y, z] = cell

        # Set per-type lambdas
        cell.lambdaVolume = P('LambdaVolumeNormoxic')
        cell.lambdaSurface = P('LambdaSurfaceNormoxic')

        # Calculate initial volume from radius directly: V = (4/3)πr³
        initial_volume = (4.0/3.0) * math.pi * (initial_radius ** 3)
        cell.targetVolume = initial_volume
        cell.targetSurface = surface_from_volume(initial_volume)
        print(f"[INIT] Created single cell {cell.id} vol={cell.volume:.1f} target={cell.targetVolume:.1f}")

        # Ensure only one cell
        for other_cell in list(self.cell_list):
            if other_cell.id != cell.id:
                try:
                    self.delete_cell(other_cell)
                    print(f"[INIT] Removed extra cell {other_cell.id}")
                except Exception:
                    pass
        print(f"[INIT] Final cell count: {len(self.cell_list)}")


# ------------------------- OXYGEN-DRIVEN FATE & GROWTH ---------------------------- #
class O2DrivenFateSteppable(SteppableBasePy):
    def __init__(self, frequency=1):
        super().__init__(frequency)
        # Cache frequently used parameters to avoid repeated lookups
        self.o2_thresh_hypoxic_necrotic = P('O2_Thresh_HypoxicNecrotic')
        self.o2_thresh_normoxic_hypoxic = P('O2_Thresh_NormoxicHypoxic')
        self.lambda_volume_normoxic = P('LambdaVolumeNormoxic')
        self.lambda_surface_normoxic = P('LambdaSurfaceNormoxic')
        self.lambda_volume_hypoxic = P('LambdaVolumeHypoxic')
        self.lambda_surface_hypoxic = P('LambdaSurfaceHypoxic')
        self.lambda_volume_necrotic = P('LambdaVolumeNecrotic')
        self.lambda_surface_necrotic = P('LambdaSurfaceNecrotic')
        self.growth_rate_normoxic = P('GrowthRateNormoxic')
        self.necrotic_shrinkage_rate = P('NecroticShrinkageRate')
        self.necrotic_lifetime = int(P('NecroticLifetime'))
        self.output_frequency = int(P('OutputFrequency'))
    def step(self, mcs):
        oxy = self.field.Oxygen
        for cell in self.cell_list:
            if cell.type == 0:  # Medium
                continue
            if cell.type == self.NECROTIC:
                self._process_necrotic(cell)
                continue

            # Sample oxygen at COM safely
            x = min(max(int(round(cell.xCOM)), 0), self.dim.x - 1)
            y = min(max(int(round(cell.yCOM)), 0), self.dim.y - 1)
            z = min(max(int(round(cell.zCOM)), 0), self.dim.z - 1)
            o2 = float(oxy[x, y, z])

            # FIELD-DRIVEN PHENOTYPE CHANGES - Two-Threshold System
            if cell.type == self.NORMOXIC:
                if o2 < self.o2_thresh_hypoxic_necrotic:
                    self._to_necrotic(cell); continue
                elif o2 < self.o2_thresh_normoxic_hypoxic:
                    self._to_hypoxic(cell)
                    if mcs % self.output_frequency == 0:
                        print(f"[FATE] Cell {cell.id} NORMOXIC->HYPOXIC at O2={o2:.3f}")
            elif cell.type == self.HYPOXIC:
                if o2 >= self.o2_thresh_normoxic_hypoxic:
                    self._to_normoxic(cell)
                    if mcs % self.output_frequency == 0:
                        print(f"[FATE] Cell {cell.id} HYPOXIC->NORMOXIC at O2={o2:.3f}")
                elif o2 < self.o2_thresh_hypoxic_necrotic:
                    self._to_necrotic(cell); continue

            self._grow(cell)

    def _to_normoxic(self, cell):
        # Preserve/restore per-cell targets on type switch to avoid XML target resets
        prev_tv = getattr(cell, 'targetVolume', cell.volume)
        cell.type = self.NORMOXIC
        # Assign per-type lambdas under Python control
        cell.lambdaVolume = self.lambda_volume_normoxic
        cell.lambdaSurface = self.lambda_surface_normoxic
        # Keep previous target volume (do not reset to current volume)
        cell.targetVolume = prev_tv
        # Keep surface coherent with target volume
        cell.targetSurface = surface_from_volume(cell.targetVolume)

    def _to_hypoxic(self, cell):
        prev_tv = getattr(cell, 'targetVolume', cell.volume)
        cell.type = self.HYPOXIC
        cell.lambdaVolume = self.lambda_volume_hypoxic
        cell.lambdaSurface = self.lambda_surface_hypoxic
        # Keep previous target volume (do not reset to current volume)
        cell.targetVolume = prev_tv
        # For hypoxic cells, don't use mathematical surface equations - keep current surface
        try:
            cell.targetSurface = cell.surface
        except Exception:
            pass  # Keep whatever surface value was there

    def _to_necrotic(self, cell):
        cell.type = self.NECROTIC
        # Assign per-type lambdas under Python control
        cell.lambdaVolume = self.lambda_volume_necrotic
        cell.lambdaSurface = self.lambda_surface_necrotic
        # Freeze volume when becoming necrotic - NO mathematical formulas
        cell.targetVolume = cell.volume
        # Keep current surface, no mathematical calculations for necrotic cells
        try:
            cell.targetSurface = cell.surface
        except Exception:
            pass  # Don't calculate surface for necrotic cells
        cell.dict['necrotic_mcs'] = self.mcs
        print(f"[FATE] Cell {cell.id} -> Necrotic at MCS {self.mcs} volume={cell.volume:.1f}")

    def _grow(self, cell):
        # Growth only for normoxic cells - use mathematical surface formula only for growing cells
        if cell.type == self.NORMOXIC:
            growth_rate = self.growth_rate_normoxic

            if growth_rate <= 0:
                # Keep previous targets unchanged to maintain a stable volume constraint
                # Update surface to remain consistent with the current target volume for normoxic only
                cell.targetSurface = surface_from_volume(cell.targetVolume)
                return

            # Accumulate target volume growth deterministically
            old_tv = cell.targetVolume
            old_ts = cell.targetSurface
            cell.targetVolume = old_tv + growth_rate
            # Linearized update for surface: dS ≈ (2/3)*(S/V)*dV
            cell.targetSurface = old_ts + (2.0 / 3.0) * (old_ts / old_tv) * growth_rate

            if self.mcs % self.output_frequency == 0:
                print(
                    f"[GROW] MCS {self.mcs} cell {cell.id} type={cell.type} "
                    f"vol={cell.volume:.1f} targetVol {old_tv:.1f}->{cell.targetVolume:.1f} "
                    f"targetSurf->{cell.targetSurface:.1f}"
                )
        
        elif cell.type == self.HYPOXIC:
            # Hypoxic cells don't grow but keep their current volume/surface targets
            # No mathematical surface calculations for hypoxic cells
            pass

    def _process_necrotic(self, cell):
        if 'necrotic_mcs' not in cell.dict:
            cell.dict['necrotic_mcs'] = self.mcs
        age = self.mcs - cell.dict['necrotic_mcs']
        lifetime_limit = self.necrotic_lifetime

        if age >= lifetime_limit:
            if random.random() < 0.75:
                self.delete_cell(cell)
                return
            else:
                cell.dict['necrotic_mcs'] = self.mcs
                return

        reduction_ratio = self.necrotic_shrinkage_rate
        if reduction_ratio > 0:
            cell.lambdaVolume = self.lambda_volume_necrotic
            cell.targetVolume = max(1, int(cell.volume * (1.0 - reduction_ratio)))


# ------------------------- MITOSIS ---------------------------- #
class O2MitosisSteppable(MitosisSteppableBase):
    def __init__(self, frequency=1):
        super().__init__(frequency)
        self.set_parent_child_position_flag(0)
        # Cache parameters used during division checks
        self.final_target_volume = P('FinalTargetVolume')
        self.growth_rate_normoxic = P('GrowthRateNormoxic')
        self.div_prob_normoxic = P('DivProbNormoxic')

    def step(self, mcs):
        # Division when target volume reaches final target volume
        cells_to_divide = []
        
        for cell in self.cell_list:
            if cell.type in (0, self.NECROTIC):
                continue
            # Only normoxic cells divide (hypoxic growth rate is 0)
            if cell.type == self.NORMOXIC and self.growth_rate_normoxic > 0:
                # Divide when target volume reaches final target volume
                if cell.targetVolume >= self.final_target_volume:
                    if random.random() < self.div_prob_normoxic:
                        cells_to_divide.append(cell)

        for cell in cells_to_divide:
            print(f"[MITOSIS-TRIGGER] MCS {mcs} cell {cell.id} vol={cell.volume:.1f} target={cell.targetVolume:.1f}")
            self.divide_cell_random_orientation(cell)

    def update_attributes(self):
        # Post-mitosis: split target volumes and keep surfaces coherent.
        self.parent_cell.targetVolume /= 2.0
        self.clone_parent_2_child()
        self.child_cell.type = self.parent_cell.type

        # Recompute target surfaces using spherical estimate (cheap)
        for c in (self.parent_cell, self.child_cell):
            c.targetSurface = surface_from_volume(c.targetVolume)

        print(
            f"[MITOSIS] Parent {self.parent_cell.id} child {self.child_cell.id} "
            f"parentTV={self.parent_cell.targetVolume:.1f} childTV={self.child_cell.targetVolume:.1f}"
        )

# ------------------------- CENTRAL COMPACTION ---------------------------- #
class CenterCompactionSteppable(SteppableBasePy):
    """Applies a polar radial force to ALL cell types for spherical compaction toward center."""

    def step(self, mcs):
        cx, cy, cz = self.dim.x / 2.0, self.dim.y / 2.0, self.dim.z / 2.0
        strength = P('CenterPushStrength')
        
        for cell in self.cell_list:
            # Skip Medium cells (type 0) but compact all living and necrotic cells
            if cell.type == 0:  # Medium
                continue
                
            # Calculate radial distance vector from cell to center
            dx = cx - cell.xCOM
            dy = cy - cell.yCOM
            dz = cz - cell.zCOM
            
            # Calculate radial distance
            r = math.sqrt(dx*dx + dy*dy + dz*dz)
            
            # Avoid division by zero for cells exactly at center
            if r < 1e-6:
                cell.lambdaVecX = 0.0
                cell.lambdaVecY = 0.0
                cell.lambdaVecZ = 0.0
                continue
            
            # Make force stronger and distance-dependent to overcome contact energies
            # Use quadratic scaling: farther cells get stronger inward force
            force_magnitude = strength * r * r * 0.1  # Scale factor to prevent excessive force
            
            # Normalize direction vector and apply radial force
            cell.lambdaVecX = force_magnitude * (dx / r)
            cell.lambdaVecY = force_magnitude * (dy / r)
            cell.lambdaVecZ = force_magnitude * (dz / r)
            
            # Debug output every 50 MCS for a few cells
            if mcs % 50 == 0 and cell.id <= 3:
                print(f"[COMPACT] MCS={mcs} Cell={cell.id} r={r:.1f} force_mag={force_magnitude:.2f} "
                      f"lambdaVec=({cell.lambdaVecX:.2f},{cell.lambdaVecY:.2f},{cell.lambdaVecZ:.2f})")




# ------------------------- LIGHT ANALYSIS / PLOTTING ---------------------------- #
class LightAnalysisSteppable(SteppableBasePy):
    def start(self):
        self.plot_total = self.add_new_plot_window(title='Total Volume', x_axis_title='MCS', y_axis_title='Volume',
                                                   x_scale_type='linear', y_scale_type='linear', grid=True)
        self.plot_total.add_plot("Volume", style='Lines', color='blue', size=2)
        self.plot_counts = self.add_new_plot_window(title='Cell Counts', x_axis_title='MCS', y_axis_title='Count',
                                                    x_scale_type='linear', y_scale_type='linear', grid=True)
        for name, color in (('Normoxic','green'), ('Hypoxic','orange'), ('Necrotic','red')):
            self.plot_counts.add_plot(name, style='Lines', color=color, size=2)

    def step(self, mcs):
        if mcs % int(P('OutputFrequency')) != 0:
            return
        counts = {'Normoxic':0,'Hypoxic':0,'Necrotic':0}
        total_vol = 0
        o2_samples = []  # Track oxygen levels
        
        for cell in self.cell_list:
            total_vol += cell.volume
            if cell.type == self.NORMOXIC:
                counts['Normoxic'] += 1
            elif cell.type == self.HYPOXIC:
                counts['Hypoxic'] += 1
            elif cell.type == self.NECROTIC:
                counts['Necrotic'] += 1

            # Sample oxygen at cell location for monitoring
            if cell.type != 0:  # Not medium
                x = min(max(int(round(cell.xCOM)), 0), self.dim.x - 1)
                y = min(max(int(round(cell.yCOM)), 0), self.dim.y - 1)
                z = min(max(int(round(cell.zCOM)), 0), self.dim.z - 1)
                o2 = float(self.field.Oxygen[x, y, z])
                o2_samples.append(o2)
                
        self.plot_total.add_data_point('Volume', mcs, total_vol)
        for k in counts:
            self.plot_counts.add_data_point(k, mcs, counts[k])
        
        # Oxygen monitoring
        if o2_samples:
            o2_min, o2_max, o2_avg = min(o2_samples), max(o2_samples), sum(o2_samples)/len(o2_samples)
            print(f"[O2-MONITOR] MCS {mcs} O2: min={o2_min:.3f} avg={o2_avg:.3f} max={o2_max:.3f} | "
                  f"Thresholds: N/H={P('O2_Thresh_NormoxicHypoxic'):.3f} H/Nec={P('O2_Thresh_HypoxicNecrotic'):.3f}")
        
        # Simple debug info
        if mcs % (int(P('OutputFrequency')) * 2) == 0:
            print(f"[STAT] MCS {mcs} Vol={total_vol:.1f} Cells={len(self.cell_list)} "
                  f"N={counts['Normoxic']} H={counts['Hypoxic']} Nec={counts['Necrotic']}")
