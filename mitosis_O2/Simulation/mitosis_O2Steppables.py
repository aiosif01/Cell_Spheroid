from cc3d.core.PySteppables import *
import random
import math

# ---- Improved Parameters for 3D Spheroid Formation ----

# Initial cell parameters
INIT_CELL_RADIUS = 3           # Start with slightly larger seed cell
INIT_TARGET_VOLUME = 45.0      # Target volume for initial growth

# Volume constraints - MUCH MORE PERMISSIVE
TV_MIN = 15.0                  # Minimum cell volume
TV_MAX = 80.0                  # MUCH HIGHER: Allow cells to grow big enough to divide
TV_NECROTIC_SHRINK = 20.0      # Necrotic cells shrink

# Growth rates - EVEN MORE EXPLOSIVE FOR NORMOXIC
GROW_RATE_NORMOXIC = 15.0      # MASSIVE: Even faster growth for normoxic
GROW_RATE_HYPOXIC = 1.0        # MUCH SLOWER: Minimal growth for hypoxic

# Oxygen thresholds for cell fate decisions - EXTREMELY LOW TRANSITION
O2_THRESH_NORMOXIC = 0.90      # EXTREMELY HIGH: Stay normoxic almost always
O2_THRESH_HYPOXIC = 0.08       # EXTREMELY LOW: Only become hypoxic in severely depleted oxygen
O2_THRESH_NECROTIC = 0.03      # EXTREMELY LOW: Almost never necrotic

# Division parameters - MUCH SMALLER THRESHOLD
DIVISION_VOLUME = 25.0         # MUCH SMALLER: Divide at very small size

# Division probabilities (per MCS) - MUCH MORE EXPLOSIVE FOR NORMOXIC
DIV_PROB_HIGH_O2 = 0.99       # ALMOST GUARANTEED: 99% division probability for normoxic
DIV_PROB_MED_O2 = 0.75        # VERY HIGH: 95% division in medium oxygen
DIV_PROB_LOW_O2 = 0.05        # MUCH LOWER: Only 5% division for hypoxic

class InitializationSteppable(SteppableBasePy):
    def start(self):
        """Initialize exactly one spherical seed cell at lattice center."""
        cx, cy, cz = self.dim.x // 2, self.dim.y // 2, self.dim.z // 2

        seed_cell = self.new_cell(self.NORMOXIC)

        # Sphere radius chosen so voxel count ~ targetVolume (r=2 => ~33 voxels)
        r = 2
        for x in range(cx - r, cx + r + 1):
            for y in range(cy - r, cy + r + 1):
                for z in range(cz - r, cz + r + 1):
                    if 0 <= x < self.dim.x and 0 <= y < self.dim.y and 0 <= z < self.dim.z:
                        if (x - cx)**2 + (y - cy)**2 + (z - cz)**2 <= r**2:
                            self.cell_field[x, y, z] = seed_cell

        # Align with XML plugin parameters for smooth relaxation
        seed_cell.targetVolume = 35.0
        seed_cell.lambdaVolume = 5.0
        seed_cell.targetSurface = 40.0
        seed_cell.lambdaSurface = 1.5

        print("Initialized single spherical seed cell (r=2)")

class O2DrivenFateSteppable(SteppableBasePy):
    def step(self, mcs):
        """Handle cell fate transitions and growth based on oxygen levels"""

        try:
            oxygen_field = self.field.Oxygen

            # SAFETY CHECK: Ensure oxygen field is properly initialized in early steps
            if mcs < 10:
                # Check oxygen at a few sample points instead of whole field
                cx, cy, cz = self.dim.x // 2, self.dim.y // 2, self.dim.z // 2
                sample_o2 = float(oxygen_field[cx, cy, cz])
                if sample_o2 < 0.5:  # If center oxygen is too low, reinitialize
                    oxygen_field[:, :, :] = 1.0
                    print(f"[MCS {mcs}] Reinitialized oxygen field - center was too low (o2={sample_o2:.3f})")

            for cell in self.cell_list:
                if cell.type == 0:  # Skip medium
                    continue

                # Get oxygen concentration at cell center with bounds checking
                x = max(0, min(int(round(cell.xCOM)), self.dim.x - 1))
                y = max(0, min(int(round(cell.yCOM)), self.dim.y - 1))
                z = max(0, min(int(round(cell.zCOM)), self.dim.z - 1))

                o2_level = float(oxygen_field[x, y, z])

                # --- Cell Fate Transitions ---
                if cell.type == self.NORMOXIC:
                    if o2_level < O2_THRESH_NECROTIC:
                        self._transition_to_necrotic(cell)
                    elif o2_level < O2_THRESH_HYPOXIC:
                        self._transition_to_hypoxic(cell)

                elif cell.type == self.HYPOXIC:
                    if o2_level > O2_THRESH_NORMOXIC:
                        self._transition_to_normoxic(cell)
                    elif o2_level < O2_THRESH_NECROTIC:
                        self._transition_to_necrotic(cell)

                elif cell.type == self.NECROTIC:
                    # NECROTIC CELLS ARE PERMANENTLY FROZEN
                    # No phenotype changes, no growth, no movement
                    # Skip all processing for necrotic cells
                    continue

                # --- Growth Based on Oxygen ---
                self._update_cell_growth(cell, o2_level)

        except Exception as e:
            print(f"[ERROR MCS {mcs}] O2DrivenFateSteppable: {e}")
            # Continue simulation even if this step fails
            pass

    def _transition_to_normoxic(self, cell):
        """Transition cell to normoxic state"""
        cell.type = self.NORMOXIC
        cell.lambdaVolume = 5.0
        cell.lambdaSurface = 1.5

    def _transition_to_hypoxic(self, cell):
        """Transition cell to hypoxic state"""
        cell.type = self.HYPOXIC
        cell.lambdaVolume = 5.0
        cell.lambdaSurface = 1.5
        # Slightly reduce target volume
        cell.targetVolume = max(TV_MIN, cell.targetVolume * 0.95)

    def _transition_to_necrotic(self, cell):
        """Transition cell to necrotic state - PERMANENTLY FROZEN"""
        cell.type = self.NECROTIC

        # FREEZE the cell completely - no more changes allowed
        cell.lambdaVolume = 50.0      # VERY HIGH: Prevent any volume changes
        cell.lambdaSurface = 50.0     # VERY HIGH: Prevent any surface changes

        # Set final volume - no more shrinking
        cell.targetVolume = cell.volume  # Keep current volume forever

    def _update_necrotic_cell(self, cell):
        """Update necrotic cell properties - COMPLETELY FROZEN"""
        # DO NOTHING - necrotic cells are frozen forever
        # No shrinkage, no growth, no changes at all
        pass

    def _update_cell_growth(self, cell, o2_level):
        """Update cell growth based on oxygen availability - MUCH MORE AGGRESSIVE"""

        # NECROTIC CELLS ARE COMPLETELY FROZEN - NO GROWTH AT ALL
        if cell.type == self.NECROTIC:
            return  # Exit immediately - no changes allowed

        if cell.type == self.NORMOXIC and o2_level > O2_THRESH_HYPOXIC:
            # EXPLOSIVE healthy growth
            growth = GROW_RATE_NORMOXIC * 2.0  # Double the growth rate
            cell.targetVolume = min(TV_MAX, cell.targetVolume + growth)

        elif cell.type == self.HYPOXIC and o2_level > O2_THRESH_NECROTIC:
            # EXPLOSIVE growth even in hypoxic conditions
            growth = GROW_RATE_HYPOXIC * 2.0  # Double the growth rate
            cell.targetVolume = min(TV_MAX * 0.8, cell.targetVolume + growth)

        # FORCE growth for LIVING cells only (exclude necrotic)
        if cell.type != self.NECROTIC and cell.type != 0:
            cell.targetVolume += 1.0  # Additional forced growth every step

class O2MitosisSteppable(MitosisSteppableBase):
    EARLIEST_DIVISION_MCS = 50  # allow initial single cell to stabilize

    def __init__(self, frequency=1):
        super().__init__(frequency)
        self.set_parent_child_position_flag(0)

    def step(self, mcs):
        """Division with initial delay so only one cell exists early on."""
        if mcs < self.EARLIEST_DIVISION_MCS:
            return

        try:
            oxygen_field = self.field.Oxygen
            cells_to_divide = []

            for cell in self.cell_list:
                if cell.type == self.NECROTIC or cell.type == 0:
                    continue

                if cell.volume > DIVISION_VOLUME:
                    x = max(0, min(int(round(cell.xCOM)), self.dim.x - 1))
                    y = max(0, min(int(round(cell.yCOM)), self.dim.y - 1))
                    z = max(0, min(int(round(cell.zCOM)), self.dim.z - 1))
                    o2_level = float(oxygen_field[x, y, z])

                    div_prob = 0.0
                    if cell.type == self.NORMOXIC:
                        if o2_level > 0.3:
                            div_prob = DIV_PROB_HIGH_O2
                        elif o2_level > 0.15:
                            div_prob = DIV_PROB_MED_O2
                        elif o2_level > 0.1:
                            div_prob = DIV_PROB_MED_O2 * 0.9
                        elif o2_level > 0.08:
                            div_prob = DIV_PROB_MED_O2 * 0.8
                    elif cell.type == self.HYPOXIC and o2_level > 0.07:
                        div_prob = DIV_PROB_LOW_O2 * 0.4

                    if random.random() < div_prob:
                        cells_to_divide.append(cell)

            for cell in cells_to_divide:
                self.divide_cell_random_orientation(cell)

            # Moderate forced growth to avoid instant size explosion
            for cell in self.cell_list:
                if cell.type != 0 and cell.type != self.NECROTIC:
                    cell.targetVolume = min(cell.targetVolume + GROW_RATE_NORMOXIC * 0.3, TV_MAX)

        except Exception as e:
            print(f"[ERROR MCS {mcs}] O2MitosisSteppable: {e}")
            pass

    def update_attributes(self):
        """Update parent and child cell attributes after division"""
        # Split volume equally
        self.parent_cell.targetVolume *= 0.5
        self.clone_parent_2_child()
        self.child_cell.targetVolume = self.parent_cell.targetVolume

        # Ensure both cells have same type
        self.child_cell.type = self.parent_cell.type

class LightAnalysisSteppable(SteppableBasePy):
    def step(self, mcs):
        """Analyze and report simulation progress"""
        if mcs % self.frequency != 0:
            return

        try:
            # Count cells by type
            cell_counts = {
                'Normoxic': 0,
                'Hypoxic': 0,
                'Necrotic': 0
            }

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
            r1 = min(10, self.dim.x // 4)
            r2 = min(20, self.dim.x // 3)

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

                if mcs % (self.frequency * 4) == 0:  # Less frequent detailed output
                    print(f"     Spheroid radius \u2248 {max_dist:.1f} voxels")

        except Exception as e:
            print(f"[ERROR MCS {mcs}] LightAnalysisSteppable: {e}")
            # Continue simulation even if analysis fails
            pass


class ActivityForcerSteppable(SteppableBasePy):
    """Force continuous activity to prevent simulation from pausing"""

    def __init__(self, frequency=10):
        super().__init__(frequency)
        self.last_cell_count = 0

    def step(self, mcs):
        """Force activity by slightly perturbing cell properties"""
        try:
            current_cell_count = len(self.cell_list)

            # If cell count hasn't changed much, force some activity
            if abs(current_cell_count - self.last_cell_count) < 2:
                for cell in self.cell_list:
                    if cell.type != 0:  # Skip medium
                        # Slightly perturb target volume to force movement
                        perturbation = 0.1 * (random.random() - 0.5)
                        cell.targetVolume += perturbation
                        # Keep within reasonable bounds
                        cell.targetVolume = max(20.0, min(cell.targetVolume, 80.0))

            self.last_cell_count = current_cell_count

            # Print activity status
            if mcs % 100 == 0:
                print(f"[ACTIVITY] MCS {mcs}: Forcing continuous simulation activity")

        except Exception as e:
            print(f"[ERROR MCS {mcs}] ActivityForcerSteppable: {e}")
            pass