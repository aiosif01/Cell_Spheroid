from cc3d.core.PySteppables import SteppableBasePy
import os
import cc3d

class VTKOutputSteppable(SteppableBasePy):
    """
    Export VTK files for ParaView visualization using CompuCell3D's built-in output system.
    Based on the CC3D repository code analysis.
    """
    
    def __init__(self, frequency=50):
        SteppableBasePy.__init__(self, frequency)

    def start(self):
        """Initialize VTK output using CompuCell3D's persistent globals system"""
        print(f"[VTK] Initializing VTK output system")
        
        # Get the persistent globals object (this is how CC3D manages output)
        pg = cc3d.CompuCellSetup.persistent_globals
        
        # The output directory should already be set by the run_script command line
        # But we can check if it exists
        if hasattr(pg, 'output_dir') and pg.output_dir:
            print(f"[VTK] Using output directory: {pg.output_dir}")
        else:
            print(f"[VTK] Warning: No output directory set in persistent globals")
        
        # Write a reference file for cell type mappings
        self._write_cell_type_mapping()
    
    def step(self, mcs):
        """The actual VTK output is handled automatically by CompuCell3D when output_frequency is set"""
        if mcs % self.frequency == 0:
            print(f"[VTK] MCS {mcs}: VTK files being written automatically by CC3D")
    
    def _write_cell_type_mapping(self):
        """Write cell type mapping file for ParaView reference in LOCAL output directory"""
        try:
            # Ensure output directory exists LOCALLY in Cell_Spheroid working directory
            output_dir = "output"  # Local output folder in Cell_Spheroid
            
            # Create output directory if it doesn't exist
            if not os.path.exists(output_dir):
                os.makedirs(output_dir)
                print(f"[VTK] Created local output directory: {output_dir}")
            
            # Write mapping file in LOCAL output directory
            mapping_file = os.path.join(output_dir, "CellTypes_Reference.txt")
            
            with open(mapping_file, 'w') as f:
                f.write("# Cell Type Mappings for ParaView Visualization\n")
                f.write("# VTK files are saved LOCALLY in Cell_Spheroid/output/ directory\n")
                f.write("# Use these values to interpret the CellType scalar in VTK files\n")
                f.write("# In ParaView: Color by 'CellType' and use these mappings:\n")
                f.write("#\n")
                f.write("TypeId=0, TypeName=Medium\n")
                f.write("TypeId=1, TypeName=Normoxic\n") 
                f.write("TypeId=2, TypeName=Hypoxic\n")
                f.write("TypeId=3, TypeName=Necrotic\n")
                f.write("#\n")
                f.write("# For Volume Rendering in ParaView:\n")
                f.write("# - Set Medium (0) to transparent (opacity=0)\n")
                f.write("# - Set Normoxic (1) to green (opacity=0.8)\n")
                f.write("# - Set Hypoxic (2) to yellow (opacity=0.8)\n")
                f.write("# - Set Necrotic (3) to red (opacity=0.8)\n")
            
            print(f"[VTK] Created LOCAL cell type mapping reference: {mapping_file}")
            
        except Exception as e:
            print(f"[VTK] Could not write cell type mapping: {e}")