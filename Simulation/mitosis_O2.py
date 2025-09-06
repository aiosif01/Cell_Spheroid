from cc3d import CompuCellSetup

from mitosis_O2Steppables import (
    InitializationSteppable,
    O2DrivenFateSteppable,
    O2MitosisSteppable,
    LightAnalysisSteppable,
    ActivityForcerSteppable
)

CompuCellSetup.register_steppable(InitializationSteppable(frequency=1))
CompuCellSetup.register_steppable(O2DrivenFateSteppable(frequency=1))
CompuCellSetup.register_steppable(O2MitosisSteppable(frequency=1))
CompuCellSetup.register_steppable(LightAnalysisSteppable(frequency=50))
CompuCellSetup.register_steppable(ActivityForcerSteppable(frequency=10))

CompuCellSetup.run()