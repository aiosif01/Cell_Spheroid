
from cc3d import CompuCellSetup
        


from HeikonanoHUBDemoSteppables import ConstraintInitializerSteppable

CompuCellSetup.register_steppable(steppable=ConstraintInitializerSteppable(frequency=1))




from HeikonanoHUBDemoSteppables import GrowthSteppable

CompuCellSetup.register_steppable(steppable=GrowthSteppable(frequency=1))




from HeikonanoHUBDemoSteppables import MitosisSteppable

CompuCellSetup.register_steppable(steppable=MitosisSteppable(frequency=1))




from HeikonanoHUBDemoSteppables import DeathSteppable

CompuCellSetup.register_steppable(steppable=DeathSteppable(frequency=1))



        
from HeikonanoHUBDemoSteppables import LineageTrackerSteppable
CompuCellSetup.register_steppable(steppable=LineageTrackerSteppable(frequency=100))

CompuCellSetup.run()
