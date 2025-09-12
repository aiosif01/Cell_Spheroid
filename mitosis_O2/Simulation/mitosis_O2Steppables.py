from cc3d.core.PySteppables import *
import random, os, math
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
        # Create exactly one cell at center - simple approach like CC3D examples
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

        # Set per-type lambdas under Python control
        cell.lambdaVolume = P('LambdaVolumeNormoxic')
        cell.lambdaSurface = P('LambdaSurfaceNormoxic')

        # Simple initialization - set targets to current measured values to avoid XML defaults
        cell.targetVolume = cell.volume
        try:
            cell.targetSurface = cell.surface
        except Exception:
            # Fallback to spherical estimate if surface not available
            cell.targetSurface = (36.0 * math.pi) ** (1.0 / 3.0) * (cell.targetVolume ** (2.0 / 3.0))
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

            # FIELD-DRIVEN PHENOTYPE CHANGES (this is what you want to keep!)
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
        # Preserve/restore per-cell targets on type switch to avoid XML target resets
        prev_tv = getattr(cell, 'targetVolume', cell.volume)
        cell.type = self.NORMOXIC
        # Assign per-type lambdas under Python control
        cell.lambdaVolume = P('LambdaVolumeNormoxic')
        cell.lambdaSurface = P('LambdaSurfaceNormoxic')
        # Keep previous target volume (do not reset to current volume)
        cell.targetVolume = prev_tv
        # Keep surface coherent with target volume
        cell.targetSurface = (36.0 * math.pi) ** (1.0 / 3.0) * (cell.targetVolume ** (2.0 / 3.0))

    def _to_hypoxic(self, cell):
        prev_tv = getattr(cell, 'targetVolume', cell.volume)
        cell.type = self.HYPOXIC
        cell.lambdaVolume = P('LambdaVolumeHypoxic')
        cell.lambdaSurface = P('LambdaSurfaceHypoxic')
        # Keep previous target volume (do not reset to current volume)
        cell.targetVolume = prev_tv
        # Keep surface coherent with target volume
        cell.targetSurface = (36.0 * math.pi) ** (1.0 / 3.0) * (cell.targetVolume ** (2.0 / 3.0))

    def _to_necrotic(self, cell):
        cell.type = self.NECROTIC
        # Assign per-type lambdas under Python control
        cell.lambdaVolume = P('LambdaVolumeNecrotic')
        cell.lambdaSurface = P('LambdaSurfaceNecrotic')
        # Freeze volume when becoming necrotic
        cell.targetVolume = cell.volume
        try:
            cell.targetSurface = cell.surface
        except Exception:
            cell.targetSurface = (36.0 * math.pi) ** (1.0 / 3.0) * (cell.targetVolume ** (2.0 / 3.0))
        cell.dict['necrotic_mcs'] = self.mcs
        print(f"[FATE] Cell {cell.id} -> Necrotic at MCS {self.mcs} volume={cell.volume:.1f}")

    def _grow(self, cell):
        # Absolute control over target volume/surface
        if cell.type in (self.NORMOXIC, self.HYPOXIC):
            growth_rate = P('GrowthRateNormoxic') if cell.type == self.NORMOXIC else P('GrowthRateHypoxic')

            if growth_rate <= 0:
                # Keep previous targets unchanged to maintain a stable volume constraint
                # Update surface to remain consistent with the current target volume
                cell.targetSurface = (36.0 * math.pi) ** (1.0 / 3.0) * (cell.targetVolume ** (2.0 / 3.0))
                return

            # Accumulate target volume growth deterministically
            old_tv = cell.targetVolume
            cell.targetVolume = old_tv + growth_rate
            cell.targetSurface = (36.0 * math.pi) ** (1.0 / 3.0) * (cell.targetVolume ** (2.0 / 3.0))

            if self.mcs % int(P('OutputFrequency')) == 0:
                print(
                    f"[GROW] MCS {self.mcs} cell {cell.id} type={cell.type} "
                    f"vol={cell.volume:.1f} targetVol {old_tv:.1f}->{cell.targetVolume:.1f} "
                    f"targetSurf->{cell.targetSurface:.1f}"
                )

    def _process_necrotic(self, cell):
        if 'necrotic_mcs' not in cell.dict:
            cell.dict['necrotic_mcs'] = self.mcs
        
        age = self.mcs - cell.dict['necrotic_mcs']
        
        if age >= int(P('NecroticLifetime')):
            try:
                # Force cell to shrink to zero before deletion
                cell.targetVolume = 1
                cell.lambdaVolume = 100.0  # Very strong constraint
                # Try deletion
                self.delete_cell(cell)
                print(f"[DELETE-SUCCESS] Cell {cell.id} removed after {age} MCS")
                return
            except Exception as e:
                # If deletion fails, just make it invisible by shrinking
                cell.targetVolume = 1
                cell.lambdaVolume = 100.0
                print(f"[DELETE-FALLBACK] Cell {cell.id} shrunk to near-zero, age={age}")
                return
        
        # Normal shrinkage
        shrink_rate = P('NecroticShrinkageRate')
        if shrink_rate > 0:
            old_target = cell.targetVolume
            cell.targetVolume = max(0.1, cell.targetVolume - shrink_rate)
            
            # Debug shrinkage
            if age % 10 == 0:
                print(f"[SHRINK] Cell {cell.id} age={age} target: {old_target:.1f}->{cell.targetVolume:.1f} actual={cell.volume:.1f}")

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
        # Simple volume-based division like CC3D examples
        cells_to_divide = []
        for cell in self.cell_list:
            if cell.type in (0, self.NECROTIC):
                continue
            # Do not allow division if growth is disabled for this phenotype
            if (cell.type == self.NORMOXIC and P('GrowthRateNormoxic') <= 0) or \
               (cell.type == self.HYPOXIC and P('GrowthRateHypoxic') <= 0):
                continue
            # Simple rule: divide when volume > 2x initial volume
            if cell.volume > P('DivisionVolumeThreshold'):
                prob = P('DivProbNormoxic') if cell.type == self.NORMOXIC else P('DivProbHypoxic')
                if random.random() < prob:
                    cells_to_divide.append(cell)

        for cell in cells_to_divide:
            print(f"[MITOSIS-TRIGGER] MCS {mcs} cell {cell.id} vol={cell.volume:.1f} target={cell.targetVolume:.1f}")
            self.divide_cell_random_orientation(cell)

    def update_attributes(self):
        # Post-mitosis: split target volumes and keep surfaces coherent.
        self.parent_cell.targetVolume /= 2.0
        self.clone_parent_2_child()
        self.child_cell.type = self.parent_cell.type

        # Recompute target surfaces. Do not reset targets to current volume even if growth is zero.
        for c in (self.parent_cell, self.child_cell):
            c.targetSurface = (36.0 * math.pi) ** (1.0 / 3.0) * (c.targetVolume ** (2.0 / 3.0))

        print(
            f"[MITOSIS] Parent {self.parent_cell.id} child {self.child_cell.id} "
            f"parentTV={self.parent_cell.targetVolume:.1f} childTV={self.child_cell.targetVolume:.1f}"
        )


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
        
        # Simple debug info
        if mcs % (int(P('OutputFrequency')) * 2) == 0:
            print(f"[STAT] MCS {mcs} Vol={total_vol:.1f} Cells={len(self.cell_list)} "
                  f"N={counts['Normoxic']} H={counts['Hypoxic']} Nec={counts['Necrotic']}")
