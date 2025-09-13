from cc3d import CompuCellSetup

from mitosis_O2Steppables import (
	OxygenInitSteppable,
	SingleCellInitSteppable,
	O2DrivenFateSteppable,
	O2MitosisSteppable,
        CenterCompactionSteppable,
	LightAnalysisSteppable
)

# All steppables now get frequencies from XML parameters internally
# Initialization runs once at start
CompuCellSetup.register_steppable(OxygenInitSteppable(frequency=1))
CompuCellSetup.register_steppable(SingleCellInitSteppable(frequency=1))
CompuCellSetup.register_steppable(O2DrivenFateSteppable(frequency=5))
CompuCellSetup.register_steppable(O2MitosisSteppable(frequency=1))
CompuCellSetup.register_steppable(CenterCompactionSteppable(frequency=1))
CompuCellSetup.register_steppable(LightAnalysisSteppable())

CompuCellSetup.run()