"""Entry point for mitosis_O2 simulation (clean demo-style)."""

from cc3d import CompuCellSetup

from mitosis_O2Steppables import (
    SpheroidInitializationSteppable,
    OxygenRegulationSteppable,
    SpheroidMitosisSteppable,
    SpheroidAnalysisSteppable,
)

# Register steppables (order similar to demo conventions)
CompuCellSetup.register_steppable(steppable=SpheroidInitializationSteppable(frequency=1))
CompuCellSetup.register_steppable(steppable=OxygenRegulationSteppable(frequency=1))
CompuCellSetup.register_steppable(steppable=SpheroidMitosisSteppable(frequency=1))
CompuCellSetup.register_steppable(steppable=SpheroidAnalysisSteppable(frequency=50))

CompuCellSetup.run()