from cc3d import CompuCellSetup

from mitosis_O2Steppables import (
    InitializationSteppable,
    O2DrivenFateSteppable,
    O2MitosisSteppable,
    LightAnalysisSteppable,
    ActivityForcerSteppable
)

# Import the VTK output steppable
from VTKOutputSteppable import VTKOutputSteppable

CompuCellSetup.register_steppable(InitializationSteppable(frequency=1))
CompuCellSetup.register_steppable(O2DrivenFateSteppable(frequency=1))
CompuCellSetup.register_steppable(O2MitosisSteppable(frequency=1))
CompuCellSetup.register_steppable(LightAnalysisSteppable(frequency=50))
CompuCellSetup.register_steppable(ActivityForcerSteppable(frequency=10))

# Register VTK output steppable - exports every 50 MCS for ParaView visualization
CompuCellSetup.register_steppable(VTKOutputSteppable(frequency=50, output_dir="vtk_output"))

CompuCellSetup.run()