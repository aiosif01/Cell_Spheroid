from cc3d import CompuCellSetup

from mitosis_O2Steppables import (
	OxygenInitSteppable,
	SingleCellInitSteppable,
	O2DrivenFateSteppable,
	O2MitosisSteppable,
	LightAnalysisSteppable
)

# All steppables now get frequencies from XML parameters internally
# Initialization runs once at start
CompuCellSetup.register_steppable(OxygenInitSteppable(frequency=1))
CompuCellSetup.register_steppable(SingleCellInitSteppable(frequency=1))
# Core simulation steppables run every step (frequency=1)
CompuCellSetup.register_steppable(O2DrivenFateSteppable(frequency=1))
CompuCellSetup.register_steppable(O2MitosisSteppable(frequency=1))
# Analysis steppable uses internal XML-controlled frequency
CompuCellSetup.register_steppable(LightAnalysisSteppable(frequency=1))

CompuCellSetup.run()