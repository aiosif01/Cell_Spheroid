from cc3d import CompuCellSetup

from mitosis_O2Steppables import (
	InitializationSteppable,
	O2DrivenFateSteppable,
	O2MitosisSteppable,
	LightAnalysisSteppable,
	VisualizationSteppable,
	ActivityForcerSteppable
)

# All steppables now get frequencies from XML parameters internally
# Core simulation steppables run every step (frequency=1)
CompuCellSetup.register_steppable(InitializationSteppable(frequency=1))
CompuCellSetup.register_steppable(O2DrivenFateSteppable(frequency=1))
CompuCellSetup.register_steppable(O2MitosisSteppable(frequency=1))
# Visualization setup runs once at start
CompuCellSetup.register_steppable(VisualizationSteppable(frequency=1))
# Analysis and Activity steppables use internal XML-controlled frequencies
CompuCellSetup.register_steppable(LightAnalysisSteppable(frequency=1))
CompuCellSetup.register_steppable(ActivityForcerSteppable())

CompuCellSetup.run()