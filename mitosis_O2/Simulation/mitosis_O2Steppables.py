from cc3d.core.PySteppables import *
import random, os
from xml.dom import minidom

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
        radius = 2.2
        cell = self.new_cell(self.NORMOXIC)
        
        # Create spherical cell manually to ensure single cell
        # Convert radius to integer for range, but use float for distance calculation
        radius_int = int(radius) + 1  # Add 1 to ensure we cover the full radius
        for x in range(cx - radius_int, cx + radius_int + 1):
            for y in range(cy - radius_int, cy + radius_int + 1):
                for z in range(cz - radius_int, cz + radius_int + 1):
                    if 0 <= x < self.dim.x and 0 <= y < self.dim.y and 0 <= z < self.dim.z:
                        if (x-cx)**2 + (y-cy)**2 + (z-cz)**2 <= radius*radius:
                            self.cell_field[x, y, z] = cell

        # Initialize target volume and a persistent division threshold per cell
        # Each cell will divide when it reaches ~2x its "birth" targetVolume
        cell.targetVolume = cell.volume
        cell.dict['V_div'] = max(2.0 * cell.targetVolume, cell.volume + 30.0)
        print(f"[INIT] Created single cell {cell.id} at center ({cx},{cy},{cz}) vol={cell.volume:.1f} targetVol={cell.targetVolume:.1f} V_div={cell.dict['V_div']:.1f}")

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

            if cell.type == self.NORMOXIC:
                if o2 < P('O2_Thresh_Necrotic'):
                    self._to_necrotic(cell); continue
                elif o2 < P('O2_Thresh_Hypoxic'):
                    self._to_hypoxic(cell)
            elif cell.type == self.HYPOXIC:
                if o2 > P('O2_Thresh_Normoxic'):
                    self._to_normoxic(cell)
                elif o2 < P('O2_Thresh_Necrotic'):
                    self._to_necrotic(cell); continue

            self._grow(cell)

    def _to_normoxic(self, cell):
        cell.type = self.NORMOXIC

    def _to_hypoxic(self, cell):
        cell.type = self.HYPOXIC

    def _to_necrotic(self, cell):
        cell.type = self.NECROTIC
        cell.targetVolume = cell.volume  # freeze current size
        cell.dict['necrotic_mcs'] = self.mcs
        print(f"[FATE] Cell {cell.id} -> Necrotic at MCS {self.mcs} volume={cell.volume:.1f}")

    def _grow(self, cell):
        # Simple growth like in CC3D examples - just increment targetVolume directly
        if cell.type in (self.NORMOXIC, self.HYPOXIC):
            rate = P('GrowthRateNormoxic') if cell.type == self.NORMOXIC else P('GrowthRateHypoxic')
            if rate <= 0:
                return
            old = cell.targetVolume
            cell.targetVolume += rate
            # Keep a reasonable cap so target doesn't run away far beyond division threshold
            if 'V_div' in cell.dict:
                cell.targetVolume = min(cell.targetVolume, cell.dict['V_div'])
            if self.mcs % int(P('OutputFrequency')) == 0:
                print(f"[GROW] MCS {self.mcs} cell {cell.id} type={cell.type} vol={cell.volume:.1f} targetVol {old:.1f}->{cell.targetVolume:.1f} V_div={cell.dict.get('V_div','n/a')}")

    def _process_necrotic(self, cell):
        if 'necrotic_mcs' not in cell.dict:
            cell.dict['necrotic_mcs'] = self.mcs
        age = self.mcs - cell.dict['necrotic_mcs']
        if age >= P('NecroticLifetime'):
            self._safe_delete(cell)
            return
        shrink = P('NecroticShrinkageRate')
        if shrink > 0:
            cell.targetVolume = max(0, cell.targetVolume - shrink)

    def _safe_delete(self, cell):
        try:
            self.delete_cell(cell)
            print(f"[DELETE] Necrotic cell {cell.id} removed after lifetime")
        except Exception as e:
            print(f"[ERROR] Failed to delete cell {cell.id}: {e}")


# ------------------------- MITOSIS ---------------------------- #
class O2MitosisSteppable(MitosisSteppableBase):
    def __init__(self, frequency=1):
        super().__init__(frequency)
        self.set_parent_child_position_flag(0)

    def step(self, mcs):
        # Divide when actual volume reaches a persistent per-cell threshold
        for cell in self.cell_list:
            if cell.type in (0, self.NECROTIC):
                continue
            v_div = cell.dict.get('V_div', None)
            if v_div is None:
                # Initialize if missing (e.g., cells loaded from old state)
                cell.dict['V_div'] = max(2.0 * max(cell.targetVolume, 1.0), cell.volume + 30.0)
                v_div = cell.dict['V_div']
            if cell.volume >= v_div:
                prob = P('DivProbNormoxic') if cell.type == self.NORMOXIC else P('DivProbHypoxic')
                if random.random() < prob:
                    print(f"[MITOSIS-TRIGGER] MCS {mcs} cell {cell.id} vol={cell.volume:.1f} target={cell.targetVolume:.1f} V_div={v_div:.1f}")
                    self.divide_cell_random_orientation(cell)

    def update_attributes(self):
        # After division, split targets and set a new division threshold for each daughter
        # First, clone parent's attributes to child so both share phenotype and CC3D state
        self.clone_parent_2_child()
        # Halve target volumes so each daughter will grow again
        self.parent_cell.targetVolume = max(1.0, self.parent_cell.volume)
        self.child_cell.targetVolume = max(1.0, self.child_cell.volume)
        # Keep the same phenotype
        self.child_cell.type = self.parent_cell.type
        # Set new division thresholds for both daughters
        self.parent_cell.dict['V_div'] = max(2.0 * self.parent_cell.targetVolume, self.parent_cell.volume + 30.0)
        self.child_cell.dict['V_div'] = max(2.0 * self.child_cell.targetVolume, self.child_cell.volume + 30.0)
        print(f"[MITOSIS] Parent {self.parent_cell.id} child {self.child_cell.id} targets=({self.parent_cell.targetVolume:.1f},{self.child_cell.targetVolume:.1f}) next V_div=({self.parent_cell.dict['V_div']:.1f},{self.child_cell.dict['V_div']:.1f})")


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
        for cell in self.cell_list:
            total_vol += cell.volume
            if cell.type == self.NORMOXIC: counts['Normoxic'] += 1
            elif cell.type == self.HYPOXIC: counts['Hypoxic'] += 1
            elif cell.type == self.NECROTIC: counts['Necrotic'] += 1
        self.plot_total.add_data_point('Volume', mcs, total_vol)
        for k in counts:
            self.plot_counts.add_data_point(k, mcs, counts[k])
        
        # Debug info every 10*OutputFrequency
        if mcs % (int(P('OutputFrequency')) * 10) == 0:
            print(f"[STAT] MCS {mcs} Vol={total_vol:.1f} Cells={len(self.cell_list)} N={counts['Normoxic']} H={counts['Hypoxic']} Nec={counts['Necrotic']}")
            # Show individual cell details
            for cell in self.cell_list:
                print(f"  Cell {cell.id}: type={cell.type} vol={cell.volume:.1f} target={cell.targetVolume:.1f}")

