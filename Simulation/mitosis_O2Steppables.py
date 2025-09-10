"""Clean, demo-style steppables for mitosis_O2 simulation.

Implements a moderate-growth spheroid with oxygen-driven fate:
 - Normoxic -> Hypoxic -> Necrotic transitions based on local oxygen.
 - Volume growth tied to oxygen.
 - Simple mitosis when above threshold (with global cell cap safeguard).
"""

import math
from cc3d.core.PySteppables import SteppableBasePy, MitosisSteppableBase

# Parameters (balanced)
INITIAL_RADIUS = 3
TARGET_VOL_START = 40.0
GROWTH_RATE_HIGH_O2 = 2.0      # volume units per MCS in good oxygen
GROWTH_RATE_LOW_O2 = 0.5        # reduced growth
O2_THRESH_HYPOXIC = 0.4
O2_THRESH_NECROTIC = 0.15
DIVISION_VOLUME = 55.0
MAX_CELLS = 5000

def _bounded(val, low, high):
    return low if val < low else high if val > high else val


class SpheroidInitializationSteppable(SteppableBasePy):
    def __init__(self, frequency=1):
        SteppableBasePy.__init__(self, frequency)

    def start(self):
        # Create single cell in center of domain
        cx, cy, cz = self.dim.x // 2, self.dim.y // 2, self.dim.z // 2
        cell = self.new_cell(self.cell_type.Normoxic)
        
        # Place cell at center pixel
        self.cell_field[cx, cy, cz] = cell
        
        cell.targetVolume = TARGET_VOL_START
        cell.lambdaVolume = 5.0
        cell.targetSurface = 2.0 * (4.0 * math.pi * (TARGET_VOL_START ** (2.0/3.0))) ** 0.5  
        cell.lambdaSurface = 2.0
        print(f"[Init] Single seed cell created at center ({cx}, {cy}, {cz})")


class OxygenRegulationSteppable(SteppableBasePy):
    def __init__(self, frequency=1):
        SteppableBasePy.__init__(self, frequency)

    def step(self, mcs):
        try:
            oxygen_field = self.field.Oxygen
        except Exception:
            return

        for cell in self.cell_list:
            if cell.type == self.cell_type.Medium:  # medium
                continue

            x = max(0, min(int(round(cell.xCOM)), self.dim.x - 1))
            y = max(0, min(int(round(cell.yCOM)), self.dim.y - 1))
            z = max(0, min(int(round(cell.zCOM)), self.dim.z - 1))
            o2 = float(oxygen_field[x, y, z])

            # Fate transitions
            if cell.type == self.cell_type.Normoxic and o2 < O2_THRESH_NECROTIC:
                cell.type = self.cell_type.Necrotic
            elif cell.type == self.cell_type.Normoxic and o2 < O2_THRESH_HYPOXIC:
                cell.type = self.cell_type.Hypoxic
            elif cell.type == self.cell_type.Hypoxic and o2 >= O2_THRESH_HYPOXIC:
                cell.type = self.cell_type.Normoxic
            elif cell.type == self.cell_type.Hypoxic and o2 < O2_THRESH_NECROTIC:
                cell.type = self.cell_type.Necrotic

            # Growth (necrotic cells stop changing volume)
            if cell.type != self.cell_type.Necrotic:
                if o2 >= O2_THRESH_HYPOXIC:
                    cell.targetVolume += GROWTH_RATE_HIGH_O2
                elif o2 >= O2_THRESH_NECROTIC:
                    cell.targetVolume += GROWTH_RATE_LOW_O2
                cell.targetVolume = _bounded(cell.targetVolume, 10.0, 120.0)


class SpheroidMitosisSteppable(MitosisSteppableBase):
    def __init__(self, frequency=1):
        MitosisSteppableBase.__init__(self, frequency)

    def start(self):
        self.set_parent_child_position_flag(0)  # random orientation

    def step(self, mcs):
        try:
            oxygen_field = self.field.Oxygen
        except Exception:
            oxygen_field = None

        if len(self.cell_list) >= MAX_CELLS:
            if len(self.cell_list) == MAX_CELLS:
                print(f"[Mitosis] Cap {MAX_CELLS} reached; divisions halted")
            return

        for cell in list(self.cell_list):
            if cell.type == self.cell_type.Medium or cell.type == self.cell_type.Necrotic:
                continue
            if cell.volume > DIVISION_VOLUME:
                # Optional oxygen gate: require minimal oxygen to divide
                if oxygen_field is not None:
                    x = max(0, min(int(round(cell.xCOM)), self.dim.x - 1))
                    y = max(0, min(int(round(cell.yCOM)), self.dim.y - 1))
                    z = max(0, min(int(round(cell.zCOM)), self.dim.z - 1))
                    if float(oxygen_field[x, y, z]) < O2_THRESH_NECROTIC:
                        continue
                self.divide_cell_random_orientation(cell)

    def update_attributes(self):
        self.parent_cell.targetVolume *= 0.5
        self.clone_parent_2_child()
        self.child_cell.targetVolume = self.parent_cell.targetVolume
        self.child_cell.type = self.parent_cell.type


class SpheroidAnalysisSteppable(SteppableBasePy):
    def __init__(self, frequency=50):
        SteppableBasePy.__init__(self, frequency)

    def step(self, mcs):
        if mcs % self.frequency:
            return
        counts = {"N": 0, "H": 0, "D": 0}
        for c in self.cell_list:
            if c.type == self.cell_type.Normoxic:
                counts["N"] += 1
            elif c.type == self.cell_type.Hypoxic:
                counts["H"] += 1
            elif c.type == self.cell_type.Necrotic:
                counts["D"] += 1
        print(f"[MCS {mcs:5d}] Cells Normoxic={counts['N']} Hypoxic={counts['H']} Necrotic={counts['D']}")