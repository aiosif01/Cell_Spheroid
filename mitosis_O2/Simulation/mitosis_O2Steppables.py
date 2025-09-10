from cc3d.core.PySteppables import *
import random
import math
from xml.dom import minidom

# Spherical surface approximation
SPHERE_SURFACE_FACTOR = (36 * math.pi) ** (1.0 / 3.0)

def spherical_surface_from_volume(volume: float) -> float:
    if volume <= 0:
        return 0.0
    return SPHERE_SURFACE_FACTOR * (volume ** (2.0 / 3.0))

def _load_user_parameters(xml_path: str):
    """Load all parameters from XML UserParameters/Param elements."""
    import os
    # If relative path fails, try in same directory as this script
    if not os.path.exists(xml_path):
        script_dir = os.path.dirname(os.path.abspath(__file__))
        xml_path = os.path.join(script_dir, xml_path)
    doc = minidom.parse(xml_path)
    up_nodes = doc.getElementsByTagName('UserParameters')
    if not up_nodes:
        raise RuntimeError('UserParameters block missing in XML')
    params = {}
    for p in up_nodes[0].getElementsByTagName('Param'):
        name = p.getAttribute('Name')
        val = p.getAttribute('Value')
        if not name:
            continue
        try:
            if '.' in val or 'e' in val.lower():
                params[name] = float(val)
            else:
                params[name] = int(val)
        except ValueError:
            raise RuntimeError(f'Invalid numeric value for parameter {name}: {val}')
    return params

# Load parameters safely with error handling
try:
    PARAMS = _load_user_parameters('mitosis_O2.xml')
    print(f"Successfully loaded {len(PARAMS)} parameters from XML")
except Exception as e:
    print(f"Warning: Failed to load XML parameters: {e}")
    # Use safe defaults to prevent crash
    PARAMS = {
        'InitRadius': 2, 'SeedTargetVolume': 35, 'SeedLambdaVolume': 10.0, 'SeedLambdaSurface': 15.0,
        'EarliestDivisionMCS': 50, 'TV_Min': 15.0, 'TV_Max': 80.0, 'GrowthRateNormoxic': 15.0,
        'GrowthRateHypoxic': 1.0, 'GrowthBoostNormoxic': 2.0, 'GrowthBoostHypoxic': 2.0,
        'ForcedGrowthIncrement': 1.0, 'PostDivisionGrowthFactor': 0.3, 'HypoxicVolumeReduceFactor': 0.95,
        'O2_Thresh_Normoxic': 0.90, 'O2_Thresh_Hypoxic': 0.08, 'O2_Thresh_Necrotic': 0.03,
        'O2_Reset_Threshold': 0.5, 'DivisionVolume': 25.0, 'DivProbHighO2': 0.99, 'DivProbMedO2': 0.75,
        'DivProbLowO2': 0.05, 'DivO2High': 0.30, 'DivO2Med': 0.15, 'DivO2Lower': 0.10,
        'DivO2LowEdge': 0.08, 'DivHypoxicMin': 0.07, 'NecroticFreezeLambda': 50.0,
        'OutputFrequency': 50, 'ActivityFrequency': 10, 'AnalysisRadius1': 10, 'AnalysisRadius2': 20,
        'DivHypoxicScale': 0.4, 'DivMedScale1': 0.9, 'DivMedScale2': 0.8, 'HypoxicMaxVolumeFrac': 0.8,
        'ActivityPerturbation': 0.1, 'ActivityStableDeltaThreshold': 2, 'EarlyO2CheckWindowFactor': 0.2
    }
    print("Using fallback default parameters")

def P(name):
    """Convenience accessor for parameters with error checking."""
    if name not in PARAMS:
        print(f"Warning: Missing parameter {name}, using default")
        return 1.0  # Safe default
    return PARAMS[name]


class InitializationSteppable(SteppableBasePy):
    def start(self):
        """Initialize exactly one spherical seed cell at lattice center."""
        cx, cy, cz = self.dim.x // 2, self.dim.y // 2, self.dim.z // 2
        seed_cell = self.new_cell(self.NORMOXIC)
        
        r = int(P('InitRadius'))
        for x in range(cx - r, cx + r + 1):
            for y in range(cy - r, cy + r + 1):
                for z in range(cz - r, cz + r + 1):
                    if 0 <= x < self.dim.x and 0 <= y < self.dim.y and 0 <= z < self.dim.z:
                        if (x - cx)**2 + (y - cy)**2 + (z - cz)**2 <= r**2:
                            self.cell_field[x, y, z] = seed_cell
        
        # Set initial cell properties from XML parameters
        seed_cell.targetVolume = P('SeedTargetVolume')
        seed_cell.lambdaVolume = P('SeedLambdaVolume')
        seed_cell.targetSurface = spherical_surface_from_volume(seed_cell.targetVolume)
        seed_cell.lambdaSurface = P('SeedLambdaSurface')
        
        print(f"Seed cell initialized: radius={r}, targetV={seed_cell.targetVolume}")


class O2DrivenFateSteppable(SteppableBasePy):
    def step(self, mcs):
        """Handle cell fate transitions and growth based on oxygen levels."""
        oxygen_field = self.field.Oxygen
        
        # Early safety re-initialization of oxygen if central value too low
        early_window = int(P('OutputFrequency') * P('EarlyO2CheckWindowFactor'))
        if mcs < early_window:
            cx, cy, cz = self.dim.x // 2, self.dim.y // 2, self.dim.z // 2
            center_o2 = float(oxygen_field[cx, cy, cz])
            if center_o2 < P('O2_Reset_Threshold'):
                oxygen_field[:, :, :] = 1.0
                print(f"[MCS {mcs}] Oxygen field reset (center={center_o2:.3f})")
        
        for cell in self.cell_list:
            if cell.type == 0:  # skip medium
                continue
            
            # Necrotic cells are frozen
            if cell.type == self.NECROTIC:
                continue
            
            x = max(0, min(int(round(cell.xCOM)), self.dim.x - 1))
            y = max(0, min(int(round(cell.yCOM)), self.dim.y - 1))
            z = max(0, min(int(round(cell.zCOM)), self.dim.z - 1))
            o2_level = float(oxygen_field[x, y, z])
            
            # State transitions
            if cell.type == self.NORMOXIC:
                if o2_level < P('O2_Thresh_Necrotic'):
                    self._to_necrotic(cell)
                    continue
                elif o2_level < P('O2_Thresh_Hypoxic'):
                    self._to_hypoxic(cell)
            elif cell.type == self.HYPOXIC:
                if o2_level > P('O2_Thresh_Normoxic'):
                    self._to_normoxic(cell)
                elif o2_level < P('O2_Thresh_Necrotic'):
                    self._to_necrotic(cell)
                    continue
            
            # Growth
            self._grow(cell, o2_level)
    
    def _to_normoxic(self, cell):
        cell.type = self.NORMOXIC
        cell.lambdaVolume = P('SeedLambdaVolume')
        cell.lambdaSurface = P('SeedLambdaSurface')
        cell.targetSurface = spherical_surface_from_volume(cell.targetVolume)
    
    def _to_hypoxic(self, cell):
        cell.type = self.HYPOXIC
        cell.lambdaVolume = P('SeedLambdaVolume')
        cell.lambdaSurface = P('SeedLambdaSurface')
        cell.targetVolume = max(P('TV_Min'), cell.targetVolume * P('HypoxicVolumeReduceFactor'))
        cell.targetSurface = spherical_surface_from_volume(cell.targetVolume)
    
    def _to_necrotic(self, cell):
        cell.type = self.NECROTIC
        cell.lambdaVolume = P('NecroticFreezeLambda')
        cell.lambdaSurface = P('NecroticFreezeLambda')
        cell.targetVolume = cell.volume
        cell.targetSurface = spherical_surface_from_volume(cell.targetVolume)
    
    def _grow(self, cell, o2_level: float):
        if cell.type == self.NECROTIC:
            return
        
        if cell.type == self.NORMOXIC and o2_level > P('O2_Thresh_Hypoxic'):
            dv = P('GrowthRateNormoxic') * P('GrowthBoostNormoxic')
            cell.targetVolume = min(P('TV_Max'), cell.targetVolume + dv)
        elif cell.type == self.HYPOXIC and o2_level > P('O2_Thresh_Necrotic'):
            dv = P('GrowthRateHypoxic') * P('GrowthBoostHypoxic')
            max_vol = P('TV_Max') * P('HypoxicMaxVolumeFrac')
            cell.targetVolume = min(max_vol, cell.targetVolume + dv)
        
        # Forced mild growth for all living cells
        if cell.type != self.NECROTIC:
            cell.targetVolume = min(P('TV_Max'), cell.targetVolume + P('ForcedGrowthIncrement'))
        
        cell.targetSurface = spherical_surface_from_volume(cell.targetVolume)


class O2MitosisSteppable(MitosisSteppableBase):
    def __init__(self, frequency=1):
        super().__init__(frequency)
        self.set_parent_child_position_flag(0)
    
    def step(self, mcs):
        """Division with initial delay so only one cell exists early on."""
        if mcs < P('EarliestDivisionMCS'):
            return
        
        try:
            oxygen_field = self.field.Oxygen
            cells_to_divide = []
            
            for cell in self.cell_list:
                if cell.type == self.NECROTIC or cell.type == 0:
                    continue
                
                if cell.volume > P('DivisionVolume'):
                    x = max(0, min(int(round(cell.xCOM)), self.dim.x - 1))
                    y = max(0, min(int(round(cell.yCOM)), self.dim.y - 1))
                    z = max(0, min(int(round(cell.zCOM)), self.dim.z - 1))
                    o2_level = float(oxygen_field[x, y, z])
                    
                    div_prob = 0.0
                    if cell.type == self.NORMOXIC:
                        if o2_level > P('DivO2High'):
                            div_prob = P('DivProbHighO2')
                        elif o2_level > P('DivO2Med'):
                            div_prob = P('DivProbMedO2')
                        elif o2_level > P('DivO2Lower'):
                            div_prob = P('DivProbMedO2') * P('DivMedScale1')
                        elif o2_level > P('DivO2LowEdge'):
                            div_prob = P('DivProbMedO2') * P('DivMedScale2')
                    elif cell.type == self.HYPOXIC and o2_level > P('DivHypoxicMin'):
                        div_prob = P('DivProbLowO2') * P('DivHypoxicScale')
                    
                    if random.random() < div_prob:
                        cells_to_divide.append(cell)
            
            for cell in cells_to_divide:
                self.divide_cell_random_orientation(cell)
            
            # Post-division forced growth
            for cell in self.cell_list:
                if cell.type != 0 and cell.type != self.NECROTIC:
                    inc = P('GrowthRateNormoxic') * P('PostDivisionGrowthFactor')
                    cell.targetVolume = min(cell.targetVolume + inc, P('TV_Max'))
                    cell.targetSurface = spherical_surface_from_volume(cell.targetVolume)
        
        except Exception as e:
            print(f"[ERROR MCS {mcs}] O2MitosisSteppable: {e}")
            pass
    
    def update_attributes(self):
        """Update parent and child cell attributes after division."""
        # Split volume equally
        self.parent_cell.targetVolume *= 0.5
        self.clone_parent_2_child()
        self.child_cell.targetVolume = self.parent_cell.targetVolume
        
        # Update spherical surfaces for both daughter cells
        self.parent_cell.targetSurface = spherical_surface_from_volume(self.parent_cell.targetVolume)
        self.child_cell.targetSurface = spherical_surface_from_volume(self.child_cell.targetVolume)
        self.parent_cell.lambdaSurface = P('SeedLambdaSurface')
        self.child_cell.lambdaSurface = P('SeedLambdaSurface')
        
        # Ensure both cells have same type
        self.child_cell.type = self.parent_cell.type


class LightAnalysisSteppable(SteppableBasePy):
    def step(self, mcs):
        """Analyze and report simulation progress."""
        if mcs % P('OutputFrequency') != 0:
            return
        
        try:
            # Count cells by type
            cell_counts = {'Normoxic': 0, 'Hypoxic': 0, 'Necrotic': 0}
            total_volume = 0
            center_x, center_y, center_z = self.dim.x // 2, self.dim.y // 2, self.dim.z // 2
            
            for cell in self.cell_list:
                total_volume += cell.volume
                if cell.type == self.NORMOXIC:
                    cell_counts['Normoxic'] += 1
                elif cell.type == self.HYPOXIC:
                    cell_counts['Hypoxic'] += 1
                elif cell.type == self.NECROTIC:
                    cell_counts['Necrotic'] += 1
            
            # Get oxygen levels at different distances from center
            oxygen_field = self.field.Oxygen
            o2_center = oxygen_field[center_x, center_y, center_z]
            
            # Sample oxygen at different radii
            r1 = int(min(P('AnalysisRadius1'), self.dim.x // 2))
            r2 = int(min(P('AnalysisRadius2'), self.dim.x // 2))
            
            o2_r1 = oxygen_field[min(center_x + r1, self.dim.x - 1), center_y, center_z]
            o2_r2 = oxygen_field[min(center_x + r2, self.dim.x - 1), center_y, center_z]
            
            print(f"[MCS {mcs:5d}] Cells: N={cell_counts['Normoxic']:3d} H={cell_counts['Hypoxic']:3d} "
                  f"D={cell_counts['Necrotic']:3d} | TotalVol={total_volume:6.0f} | "
                  f"O2: center={o2_center:.3f} r={r1}:{o2_r1:.3f} r={r2}:{o2_r2:.3f}")
            
            # Calculate spheroid radius (approximate)
            if len(self.cell_list) > 0:
                max_dist = 0
                for cell in self.cell_list:
                    dist = math.sqrt((cell.xCOM - center_x)**2 +
                                   (cell.yCOM - center_y)**2 +
                                   (cell.zCOM - center_z)**2)
                    max_dist = max(max_dist, dist)
                
                if mcs % (P('OutputFrequency') * 4) == 0:
                    print(f"     Spheroid radius ~ {max_dist:.1f} voxels")
        
        except Exception as e:
            print(f"[ERROR MCS {mcs}] LightAnalysisSteppable: {e}")
            pass


class VisualizationSteppable(SteppableBasePy):
    """Configure visualization settings for spherical cell display."""
    
    def start(self):
        """Set up visualization preferences."""
        print("Visualization: Using CompuCell3D default cell rendering")
        print("Note: Cells will appear with spherical constraints from Surface plugin")
    
    def step(self, mcs):
        # Do nothing - visualization is handled by CompuCell3D's default renderer
        pass


class ActivityForcerSteppable(SteppableBasePy):
    """Force continuous activity to prevent simulation from pausing."""
    
    def __init__(self):
        super().__init__(frequency=int(P('ActivityFrequency')))
        self.last_cell_count = 0
    
    def step(self, mcs):
        """Force activity by slightly perturbing cell properties."""
        try:
            current_cell_count = len(self.cell_list)
            
            # If cell count hasn't changed much, force some activity
            threshold = P('ActivityStableDeltaThreshold')
            if abs(current_cell_count - self.last_cell_count) < threshold:
                for cell in self.cell_list:
                    if cell.type != 0:  # Skip medium
                        # Slightly perturb target volume to force movement
                        perturbation = P('ActivityPerturbation') * (random.random() - 0.5)
                        cell.targetVolume += perturbation
                        # Keep within reasonable bounds
                        cell.targetVolume = max(P('TV_Min'), min(cell.targetVolume, P('TV_Max')))
            
            self.last_cell_count = current_cell_count
            
            # Print activity status
            if mcs % (P('OutputFrequency') * 2) == 0:
                print(f"[ACTIVITY] MCS {mcs}: activity enforcement check")
        
        except Exception as e:
            print(f"[ERROR MCS {mcs}] ActivityForcerSteppable: {e}")
            pass
