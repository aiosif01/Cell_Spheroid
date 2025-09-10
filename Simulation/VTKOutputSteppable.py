from cc3d.core.PySteppables import SteppableBasePy
from cc3d.core.CMLFieldHandler import CMLFieldHandler
import os

class VTKOutputSteppable(SteppableBasePy):
    """
    Export VTK files for ParaView visualization of 3-phenotype spheroid.
    This creates Step_XXXXX.vtk files that can be directly opened in ParaView.
    """
    
    def __init__(self, frequency=50, output_dir="vtk_output", core_name="Step"):
        SteppableBasePy.__init__(self, frequency)
        self.output_dir = output_dir
        self.core_name = core_name
        self.handler = None
        
    def start(self):
        """Initialize VTK output handler"""
        print(f"[VTK] Initializing VTK output to {self.output_dir}/")
        
        # Create the CMLFieldHandler
        self.handler = CMLFieldHandler()
        
        # Configure output paths
        self.handler.output_dir_name = self.output_dir
        self.handler.output_file_core_name = self.core_name
        self.handler.output_freq = self.frequency
        
        # Create output directory if it doesn't exist
        if not os.path.exists(self.output_dir):
            os.makedirs(self.output_dir)
            print(f"[VTK] Created output directory: {self.output_dir}")
        
        # Initialize the handler's storage directory
        try:
            self.handler.create_storage_dir()
        except AttributeError:
            # If this method doesn't exist, we'll handle directory creation above
            pass
        
        # Write XML description file with cell type mappings
        # This creates a .dml file that maps TypeId to TypeName for ParaView
        try:
            self.handler.write_xml_description_file()
            print(f"[VTK] Created XML description file for cell type mappings")
        except AttributeError:
            # If this method doesn't exist, we'll create our own mapping file
            self._write_custom_cell_type_mapping()
    
    def step(self, mcs):
        """Export VTK files at specified frequency"""
        if mcs % self.frequency == 0:
            print(f"[VTK] Writing VTK files at MCS {mcs}")
            
            try:
                # Write the VTK files using CMLFieldHandler
                # This should create Step_XXXXXX.vtk files
                self.handler.write_fields(mcs=mcs)
                
                # Also write oxygen field if available
                self._write_oxygen_field(mcs)
                
                print(f"[VTK] Successfully wrote {self.core_name}_{mcs:06d}.vtk")
                
            except Exception as e:
                print(f"[VTK ERROR] Failed to write VTK at MCS {mcs}: {e}")
                # Try alternative method
                self._write_vtk_alternative(mcs)
    
    def _write_custom_cell_type_mapping(self):
        """Write custom cell type mapping file if CMLFieldHandler doesn't provide one"""
        mapping_file = os.path.join(self.output_dir, f"{self.core_name}_CellTypes.txt")
        
        with open(mapping_file, 'w') as f:
            f.write("# Cell Type Mappings for ParaView\n")
            f.write("# Use this to interpret the scalar values in VTK files\n")
            f.write("TypeId=0, TypeName=Medium\n")
            f.write("TypeId=1, TypeName=Normoxic\n") 
            f.write("TypeId=2, TypeName=Hypoxic\n")
            f.write("TypeId=3, TypeName=Necrotic\n")
        
        print(f"[VTK] Created custom cell type mapping: {mapping_file}")
    
    def _write_oxygen_field(self, mcs):
        """Export oxygen field as separate VTK file"""
        try:
            # Get the oxygen field
            oxygen_field = self.field.Oxygen
            
            # Create filename for oxygen field
            oxygen_filename = os.path.join(
                self.output_dir, 
                f"Oxygen_{mcs:06d}.vtk"
            )
            
            # Write oxygen field to VTK format
            # This is a simplified approach - you might need to use the field handler
            # for more complex scenarios
            self._write_scalar_field_vtk(oxygen_field, oxygen_filename, "Oxygen")
            
        except Exception as e:
            print(f"[VTK] Could not write oxygen field at MCS {mcs}: {e}")
    
    def _write_scalar_field_vtk(self, field, filename, field_name):
        """Write a scalar field to VTK format manually"""
        try:
            with open(filename, 'w') as f:
                # VTK header
                f.write("# vtk DataFile Version 3.0\n")
                f.write(f"{field_name} field from CompuCell3D\n")
                f.write("ASCII\n")
                f.write("DATASET STRUCTURED_POINTS\n")
                
                # Dimensions
                f.write(f"DIMENSIONS {self.dim.x} {self.dim.y} {self.dim.z}\n")
                f.write("ORIGIN 0 0 0\n")
                f.write("SPACING 1 1 1\n")
                
                # Point data
                total_points = self.dim.x * self.dim.y * self.dim.z
                f.write(f"POINT_DATA {total_points}\n")
                f.write(f"SCALARS {field_name} float 1\n")
                f.write("LOOKUP_TABLE default\n")
                
                # Write field values
                for z in range(self.dim.z):
                    for y in range(self.dim.y):
                        for x in range(self.dim.x):
                            value = float(field[x, y, z])
                            f.write(f"{value}\n")
            
            print(f"[VTK] Wrote {field_name} field to {filename}")
            
        except Exception as e:
            print(f"[VTK] Error writing {field_name} field: {e}")
    
    def _write_vtk_alternative(self, mcs):
        """Alternative VTK writing method if CMLFieldHandler fails"""
        try:
            # Create cell field VTK file manually
            cell_filename = os.path.join(
                self.output_dir,
                f"{self.core_name}_{mcs:06d}.vtk"
            )
            
            with open(cell_filename, 'w') as f:
                # VTK header
                f.write("# vtk DataFile Version 3.0\n")
                f.write("Cell field from CompuCell3D spheroid simulation\n")
                f.write("ASCII\n")
                f.write("DATASET STRUCTURED_POINTS\n")
                
                # Dimensions
                f.write(f"DIMENSIONS {self.dim.x} {self.dim.y} {self.dim.z}\n")
                f.write("ORIGIN 0 0 0\n")
                f.write("SPACING 1 1 1\n")
                
                # Point data
                total_points = self.dim.x * self.dim.y * self.dim.z
                f.write(f"POINT_DATA {total_points}\n")
                f.write("SCALARS CellType int 1\n")
                f.write("LOOKUP_TABLE default\n")
                
                # Write cell types (0=Medium, 1=Normoxic, 2=Hypoxic, 3=Necrotic)
                for z in range(self.dim.z):
                    for y in range(self.dim.y):
                        for x in range(self.dim.x):
                            cell = self.cell_field[x, y, z]
                            if cell:
                                cell_type = cell.type
                            else:
                                cell_type = 0  # Medium
                            f.write(f"{cell_type}\n")
            
            print(f"[VTK] Alternative method: wrote {cell_filename}")
            
        except Exception as e:
            print(f"[VTK] Alternative VTK writing also failed: {e}")
