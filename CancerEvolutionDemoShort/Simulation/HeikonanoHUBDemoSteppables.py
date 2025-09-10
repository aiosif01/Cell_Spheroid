from cc3d.cpp.PlayerPython import *
from cc3d import CompuCellSetup

from cc3d.core.PySteppables import *
import numpy as np
import random

global simulation_type

# Probability that a stem cell has two stem offspring. Default is 0.5 (neutral)
P_init = 0.75

# Number of generations for a somatic cell. Default is 5.
# Bigger makes the system less evolvable since it is more likely to freeze
Rho_init = 5.0

# Default is 0.1, Probability of mutation in stemness coefficient per division
Prob_P_Mutation = 0.1

# Default is 0.1, Probability of mutation of senescence counter per division
Prob_Rho_Mutation = 0.1

# Default is 0.05. Typical amplitude of mutation of stemness coefficient on mutation
P_Mutation_Amplitude = 0.05

# Default is 1. Typical amplitude of mutation of generation lifetime coefficient on mutation
Rho_Mutation_Amplitude = 1.0

# If true then stem cells become somatic, if false, stem only
Somatic_Flag = True

# Controls compressibility--default is 4.0 bigger makes the system stiffer and more freezable
lambdaVolume = 8.0 #Was 4

# Default is 60. Controls the degree of overpressure allowed in the simulation.
# Bigger makes the system more collective and evolve faster.
maxTargetVolume = 60.0

# Set Pressure threshold for growing cells
Pressure_Threshold = 30.00 # Was 50 

# Information on cell death and simulated treatments

# if constant killing reduce kill rate by 1E4  a scale factor that we have to play with

kill_reduction=1.0E3

# Default is 1000. in MCS Simulate Chemo time of exposure
chemo_time = 1000

# Default is 1000. in MCS Simulate Radiation time of exposure
radiation_time = 1000

# Default is 2000. in MCS Simulate Chemo time of exposure
surgery_time = 2000

# Relative killing efficiency of radiation for stem cells
radiation_protection = 0.5

# Relative killing efficiency of radiation for quiescent protected cells
radiation_quiescence_protection = 0.5

# killing efficacy of chemo for stem cells Default is 0.8
chemo_efficacy_stem = 0.8

# killing efficacy of chemo for somatic cells
chemo_efficacy_somatic = 0.9

mutation_count = 0

# time over which to track back cell lineages
lineage_lookback = 200


class ConstraintInitializerSteppable(SteppableBasePy):
    def __init__(self, frequency=1):
        SteppableBasePy.__init__(self, frequency)

    def start(self):
        for cell in self.cell_list:
            cell.targetVolume = 25
            cell.lambdaVolume = lambdaVolume
            # Add stemness property and Number of generations of somatic cells
            # Probability that a stem cell has two stem-like offspring
            cell.dict['Pstem'] = P_init

            # Number of generations the somatic cells can divide
            cell.dict['Rho'] = Rho_init

            # to keep somatic counter separate from parent count
            cell.dict['Generations_left'] = cell.dict['Rho']

            # Tracks the ancestry of a cell
            cell.dict['Lineage'] = []
            cell.dict['Pressure'] = 2.0 *max(0., (cell.targetVolume - cell.volume) * cell.lambdaVolume)

            # Tracks mutations, for each mutation include mutation_count, time, Pstem, Rho
            # and any other evolving parameters
            cell.dict['Mutations'] = []

            # and any other evolving parameters
            #cell.dict['Mutations'].append([mutation_count, 0, cell.dict['Pstem'], cell.dict['Rho']])

            # Allows display of current heterogeneity of mutation history
            cell.dict['Current_Mutation'] = mutation_count


class GrowthSteppable(SteppableBasePy):
    def __init__(self, frequency=1):
        SteppableBasePy.__init__(self, frequency)

    def add_steering_panel(self):
        #Selects which type of simulation--Neutral evolution, normal evolution, more elaborate models
        self.add_steering_param(name='simulation_type', val="Select Simulation to Run", enum=["Neutral Variation","Simple Somatic Simulation"],
                                widget_name='combobox')
        # Selects whether to have any cell killing mechanisms, options, default=none, constant, surgical one time, surgical periodic, radiation one time, radiation periodic, chemo periodic                       
        self.add_steering_param(name='killing_type', val="Select Type of Treatment", enum=["None","Constant Killing All Cells","Constant Killing Somatic Only","Surgical","Radiation Simple",
            "Radiation Stem Protected","Radiation Quiescent Protected"],widget_name='combobox')
                                
        self.add_steering_param(name='Killing_efficacy', val=0.5, min_val=0.0, max_val=1.0, decimal_precision=2,
                                widget_name='slider')
                                
        self.add_steering_param(name='growth_type', val="Growth Law", enum=["Exponential","Pressure Limited","Contact Inhibited"],
                                widget_name='combobox')
        
        self.add_steering_param(name='P_init', val=0.5, min_val=0.0, max_val=1.0, decimal_precision=2,
                                widget_name='slider')

        self.add_steering_param(name='Rho_init', val=5, min_val=1, max_val=20,
                                widget_name='slider')

        self.add_steering_param(name='Prob_P_Mutation', val=0.1, min_val=0.0, max_val=1.0, decimal_precision=2,
                                widget_name='slider')

        self.add_steering_param(name='Prob_Rho_Mutation', val=0.1, min_val=0.0, max_val=1.0, decimal_precision=2,
                                widget_name='slider')

        self.add_steering_param(name='P_Mutation_Amplitude', val=0.05, min_val=0.0, max_val=1.0, decimal_precision=2,
                                widget_name='slider')

        self.add_steering_param(name='Rho_Mutation_Amplitude', val=1.0, min_val=0.0, max_val=3.0, decimal_precision=2,
                                widget_name='slider')

        self.add_steering_param(name='chemo_time', val=2000, min_val=100, max_val=20000,
                                widget_name='slider')

        self.add_steering_param(name='chemo_efficacy_somatic', val=0.9, min_val=0.0, max_val=1.0, decimal_precision=2,
                                widget_name='slider')

    def process_steering_panel_data(self):
        global simulation_type
        simulation_type = self.get_steering_param('simulation_type')
        
        global killing_type
        killing_type = self.get_steering_param('killing_type')
        
        global Killing_efficacy
        Killing_efficacy = self.get_steering_param('Killing_efficacy')
        
        global growth_type
        growth_type = self.get_steering_param('growth_type')
        
        global P_init
        P_init = self.get_steering_param('P_init')

        global Rho_init
        Rho_init = self.get_steering_param('Rho_init')

        global Prob_P_Mutation
        Prob_P_Mutation = self.get_steering_param('Prob_P_Mutation')

        global Prob_Rho_Mutation
        Prob_Rho_Mutation = self.get_steering_param('Prob_Rho_Mutation')

        global P_Mutation_Amplitude
        P_Mutation_Amplitude = self.get_steering_param('P_Mutation_Amplitude')

        global Rho_Mutation_Amplitude
        Rho_Mutation_Amplitude = self.get_steering_param('Rho_Mutation_Amplitude')

        global chemo_time
        chemo_time = self.get_steering_param('chemo_time')

        global chemo_efficacy_somatic
        chemo_efficacy_somatic = self.get_steering_param('chemo_efficacy_somatic')

    def start(self):

        # Add plot for mean of Stemness

        self.plot_win2 = self.add_new_plot_window(title='MeanStemness',
                                                  x_axis_title='MonteCarlo Step (MCS)',
                                                  y_axis_title='Stemness', x_scale_type='linear', y_scale_type='linear',
                                                  grid=False)
        self.plot_win2.add_plot("Stemness", style='Lines', color='red', size=5)

        self.plot_win4 = self.add_new_plot_window(title='MeanGenerations',
                                                  x_axis_title='MonteCarlo Step (MCS)',
                                                  y_axis_title='Generations', x_scale_type='linear',
                                                  y_scale_type='linear',
                                                  grid=False)
        self.plot_win4.add_plot("Generations", style='Lines', color='green', size=5)

    def step(self, mcs):
        #Hold until the simulation type is determined
        while self.get_steering_param('simulation_type')=="Select Simulation to Run":
            pass
        # make sure simulation type is allocated
        simulation_type=self.get_steering_param('simulation_type')
        #At the start of the simulation read the sliders for initialization
        if mcs == 0:
            
            # Default is 0.1, Probability of mutation in stemness coefficient per division
            Prob_P_Mutation = self.get_steering_param("Prob_P_Mutation")

            # Default is 0.1, Probability of mutation of senescence counter per division
            Prob_Rho_Mutation = self.get_steering_param("Prob_Rho_Mutation")

            # Default is 0.05. Typical amplitude of mutation of stemness coefficient on mutation
            P_Mutation_Amplitude = self.get_steering_param("P_Mutation_Amplitude")

            # Default is 1. Typical amplitude of mutation of generation lifetime coefficient on mutation
            Rho_Mutation_Amplitude = self.get_steering_param("Rho_Mutation_Amplitude")


            # self.add_steering_param(name='simulation_type', val="Select Simulation to Run", enum=["Neutral Variation","Simple Somatic Simulation","Competition Model"],
                                # widget_name='combobox')
            if simulation_type=="Neutral Variation":
                Somatic_Flag=False
                print("neutral variation simulation version")
                
            for cell in self.cell_list:
                cell.targetVolume = 25
                cell.lambdaVolume = lambdaVolume
                # Add stemness property and Number of generations of somatic cells
                # Probability that a stem cell has two stem-like offspring
                cell.dict['Pstem'] = self.get_steering_param('P_init')

                # Number of generations the somatic cells can divide
                cell.dict['Rho'] = self.get_steering_param('Rho_init')

                # to keep somatic counter separate from parent count
                cell.dict['Generations_left'] = cell.dict['Rho']

                # Tracks the ancestry of a cell
                cell.dict['Lineage'] = []
                cell.dict['Pressure'] = 2.0 *max(0.,(cell.targetVolume - cell.volume) * cell.lambdaVolume)

                # Tracks mutations, for each mutation include mutation_count, time, Pstem, Rho
                # and any other evolving parameters
                cell.dict['Mutations'] = []

                # and any other evolving parameters
                #cell.dict['Mutations'].append([mutation_count, 0, cell.dict['Pstem'], cell.dict['Rho']])

                # Allows display of current heterogeneity of mutation history
                cell.dict['Current_Mutation'] = mutation_count


 
            self.plot_win = self.add_new_plot_window(title='Histogram of Stemness', x_axis_title='P_stem',
                                                 y_axis_title='Probability')
            self.plot_win3 = self.add_new_plot_window(title='Histogram of Senescence Generations',
                                                  x_axis_title='Generations',
                                                  y_axis_title='Probability')
                                                  # self.plot_win.add_histogram_plot(plot_name='Hist 2', color='red', alpha=100)
        # self.plot_win.add_histogram_plot(plot_name='Hist 3', color='blue')
            self.plot_win.add_histogram_plot(plot_name='Stemness', color='green', alpha=100)
            self.plot_win3.add_histogram_plot(plot_name='Generations', color='blue', alpha=100)
        # _alpha is transparency 0 is transparent, 255 is opaque
        # Create empty list for the stemmess attribute
        pstem = []

        # Create empty list for number of generations of somatic cells
        generations = []
        
        # Read Growth law and put loop inside to speed things up
        growth_type=self.get_steering_param('growth_type')
        if growth_type=="Exponential":
            for cell in self.cellListByType(self.STEM, self.SOMATIC):
                cell.targetVolume += 0.2
                cell.targetVolume = min(cell.targetVolume, 60)
        elif growth_type=="Pressure Limited" or growth_type=="Growth Law" : # This is the default if none selected
            for cell in self.cellListByType(self.STEM, self.SOMATIC):
                if cell.dict['Pressure']< Pressure_Threshold: 
                    cell.targetVolume += 0.2
        elif growth_type=="Contact Inhibited":
            for cell in self.cellListByType(self.STEM, self.SOMATIC):
                cell.targetVolume += 2.0/(1.0+cell.dict['Pressure']/Pressure_Threshold)
       
        for cell in self.cellListByType(self.STEM, self.SOMATIC):
            # #print("pressure = ",cell.dict['Pressure'])
                       
            # Following is for lineage tracking probably slows down simulation check)
            if cell.type == self.STEM:
                pstem.append(cell.dict['Pstem'])
                generations.append(cell.dict['Rho'])

        if not mcs % 10:
            self.plot_win.add_histogram(plot_name='Stemness', value_array=pstem, number_of_bins=10)
            self.plot_win3.add_histogram(plot_name='Generations', value_array=generations, number_of_bins=10)
            pstem_mean = np.mean(pstem)
            generations_mean = np.mean(generations)
            # arguments are (name of the data series, x, y)
            self.plot_win2.add_data_point("Stemness", mcs, pstem_mean)
            self.plot_win4.add_data_point("Generations", mcs, generations_mean)


class MitosisSteppable(MitosisSteppableBase):
    def __init__(self, frequency=1):
        MitosisSteppableBase.__init__(self, frequency)

    def start(self):
        self.shared_steppable_vars['number_of_mutations_per_step'] = 0
        

    def step(self, mcs):
        # Default is 0.1, Probability of mutation in stemness coefficient per division
        Prob_P_Mutation = self.get_steering_param("Prob_P_Mutation")

        # Default is 0.1, Probability of mutation of senescence counter per division
        Prob_Rho_Mutation = self.get_steering_param("Prob_Rho_Mutation")

        # Default is 0.05. Typical amplitude of mutation of stemness coefficient on mutation
        P_Mutation_Amplitude = self.get_steering_param("P_Mutation_Amplitude")

        # Default is 1. Typical amplitude of mutation of generation lifetime coefficient on mutation
        Rho_Mutation_Amplitude = self.get_steering_param("Rho_Mutation_Amplitude")
        
        # set counter to keep track of number of mutations per time step
        
       

        if (mcs == 0):
                
                # Plot window 5 is for number of stem cells--needed in all simulations
                self.plot_win5 = self.add_new_plot_window(title='Number of Stem Cells',
                                                  x_axis_title='MonteCarlo Step (MCS)',
                                                  y_axis_title='Number of Stem Cells', x_scale_type='linear',
                                                  y_scale_type='linear',
                                                  grid=False)
                self.plot_win5.add_plot("Number of Stem Cells", style='Lines', color='green', size=5)
                
                # Plot window 10 is for number of mutations per time step--needed in all simulations
                self.plot_win10 = self.add_new_plot_window(title='Number of Mutations per Step',
                                                  x_axis_title='MonteCarlo Step (MCS)',
                                                  y_axis_title='Number of Mutations per Step', x_scale_type='linear',
                                                  y_scale_type='linear',
                                                  grid=False)
                self.plot_win10.add_plot("Number of Mutations per Step", style='Lines', color='green', size=5)
                
                #FOllowing plots not used in Neutral Variation case
                if self.get_steering_param('simulation_type') != "Neutral Variation":
                

                    # self.plot_win6 = self.add_new_plot_window(title='Number of Somatic Cells',
                                                      # x_axis_title='MonteCarlo Step (MCS)',
                                                      # y_axis_title='Number of Somatic Cells', x_scale_type='linear',
                                                      # y_scale_type='linear',
                                                      # grid=False)
                    # self.plot_win6.add_plot("Number of Somatic Cells", style='Lines', color='green', size=5)
                    # self.plot_win7 = self.add_new_plot_window(title='Number of Dividing Somatic Cells',
                                                              # x_axis_title='MonteCarlo Step (MCS)',
                                                              # y_axis_title='Number of Dividing Somatic Cells',
                                                              # x_scale_type='linear', y_scale_type='linear',
                                                              # grid=False)
                    # self.plot_win7.add_plot("Number of Dividing Somatic Cells", style='Lines', color='green', size=5)

                    self.plot_win8 = self.add_new_plot_window(title='Number of Dividing Stem Cells',
                                                              x_axis_title='MonteCarlo Step (MCS)',
                                                              y_axis_title='Number of Dividing Stem Cells', x_scale_type='linear',
                                                              y_scale_type='linear',
                                                              grid=False)
                    self.plot_win8.add_plot("Number of Dividing Stem Cells", style='Lines', color='green', size=5)

                    # self.plot_win9 = self.add_new_plot_window(title='Fraction of Stem Cells',
                                                              # x_axis_title='MonteCarlo Step (MCS)',
                                                              # y_axis_title='Fraction of Stem Cells', x_scale_type='linear',
                                                              # y_scale_type='linear',
                                                              # grid=False)
                    # self.plot_win9.add_plot("Fraction of Stem Cells", style='Lines', color='green', size=5)
               





        cells_to_divide = []
        stem_cells_to_divide = []
        stem_cell_list = []
        somatic_cell_list = []
        somatic_cells_to_divide = []
        for cell in self.cell_list:
            cell.dict['Pressure'] = 2.0 * (cell.targetVolume - cell.volume) * cell.lambdaVolume
            if cell.type == self.STEM:  # only stem cells mutate
                stem_cell_list.append(cell)
            elif cell.type == self.SOMATIC:
                somatic_cell_list.append(cell)
            if cell.volume > 25:
                cells_to_divide.append(cell)
                if cell.type == self.STEM:  # only stem cells mutate
                    stem_cells_to_divide.append(cell)
                if cell.type == self.SOMATIC:  # Somatic Cells lose a generation on division
                    cell.dict['Generations_left'] -= 1.0  # count down the generations
                    somatic_cells_to_divide.append(cell)

        self.plot_win5.add_data_point("Number of Stem Cells", mcs, len(stem_cell_list))
        if self.get_steering_param('simulation_type') != "Neutral Variation":
            #self.plot_win6.add_data_point("Number of Somatic Cells", mcs, len(somatic_cell_list))
            #self.plot_win7.add_data_point("Number of Dividing Somatic Cells", mcs, len(somatic_cells_to_divide))
            self.plot_win8.add_data_point("Number of Dividing Stem Cells", mcs, len(stem_cells_to_divide))
            # if float(len(somatic_cell_list) + len(stem_cell_list)) >0.0 : # Trap divide by zero error
                # self.plot_win9.add_data_point("Fraction of Stem Cells", mcs,
                                      # float(len(stem_cell_list)) / float(len(somatic_cell_list) + len(stem_cell_list)))

        # Better way to do mutation to ensure it is unbiassed
        # Pick delta of change
        # Pick pair of cells from division list
        # Check if you can add delta to one and subtract from the other without going out of range
        # If so, do so
        # If not, have to decide whether to do nothing or to repick the values from the distribution
        # There is a special case if there is only one cell to divide
        # Remember P_Mutation_Amplitude Typical amplitude of mutation of stemness coefficient on mutation
        # Rho_Mutation_Amplitude is typical mutation amlitude of lifetime of somatic cells
        print("number of stem cells ", len(stem_cell_list), "number of stem cells to divide", len(stem_cells_to_divide))
        if (len(stem_cells_to_divide) == 1) and (len(stem_cell_list) == 1):
            print("Failed to find pair to change")
        else:
            for cell in stem_cells_to_divide:
                if random.random() < Prob_P_Mutation:  # check if this one is going to mutate
                    delta = abs(random.gauss(0.0, P_Mutation_Amplitude))  # Pick mutation amplitude
                    #print("Mutating stemness for cell.id", cell.id, " by ", delta)
                    # If there are multiple stem cells to divide pick a second one from the list
                    if len(stem_cells_to_divide) > 1:
                        # pick second cell from list
                        cell_two = random.choice(stem_cells_to_divide)
                        # if the two cells are the same pick again
                        while cell_two.id == cell.id:
                            # pick second cell from list
                            cell_two = random.choice(cells_to_divide)
                            # now have two different cells
                    else:
                        # if only one stem cell to divide need to mutate a nondividing cell
                        # Check that there are enough stem cells to choose from
                        if len(stem_cell_list) > 1:
                            # pick second cell from list
                            cell_two = random.choice(stem_cell_list)

                            # if the two cells are the same pick again
                            while cell_two.id == cell.id:
                                # pick second cell from list
                                cell_two = random.choice(stem_cell_list)

                                # now have two different cells

                            print("Cell 1:", cell.id, "Pstem:", cell.dict['Pstem'], "Cell 2:", cell_two.id, " Pstem:",
                            cell_two.dict['Pstem'])
                    # Check that the change is legal
                    if cell.dict['Pstem'] > delta and cell_two.dict['Pstem'] < 1 - delta:
                        # if change is illegal do nothing (alternatively could repick)
                        # mutatio count is total number of nutations in simulation
                        # number_of_mutations_per_step is number in this step only
                        global mutation_count
                        mutation_count += 1
                        self.shared_steppable_vars['number_of_mutations_per_step'] +=1
                        cell.dict['Pstem'] -= delta
                        
                        # and any other evolving parameters
                        # I am commenting out the tracking data for the simplifed simulation
                        #cell.dict['Mutations'].append([mutation_count, mcs, cell.dict['Pstem'], cell.dict['Rho']])
                        cell_two.dict['Pstem'] += delta
                        #cell_two.dict['Mutations'].append([mutation_count, mcs, cell_two.dict['Pstem'],
                         #                                  cell_two.dict['Rho']])  # and any other evolving parameters
                        #print("Made Mutation", "Cell 1:", cell.id, "Pstem:", cell.dict['Pstem'], "Cell 2:", cell_two.id,
                              #" Pstem:", cell_two.dict['Pstem'])
                    else:
                        #print("Vetoed Mutation")
                        pass

                if random.random() < Prob_Rho_Mutation:
                    delta2 = abs(random.gauss(0.0, Rho_Mutation_Amplitude))
                    #print("Mutating number of generations for cell.id", cell.id, " by ", delta2)

                    if len(stem_cells_to_divide) > 1:
                        # If there are multiple stem cells to divide pick a second one from the list
                        # pick second cell from list
                        cell_two = random.choice(stem_cells_to_divide)

                        # if the two cells are the same pick again
                        while cell_two.id == cell.id:
                            # pick second cell from list
                            cell_two = random.choice(cells_to_divide)

                    else:
                        # if only one stem cell to divide need to mutate a nondividing cell

                        # Check that there are enough stem cells to choose from
                        if len(stem_cell_list) > 1:
                            # pick second cell from list
                            cell_two = random.choice(stem_cell_list)

                            # if the two cells are the same pick again
                            while cell_two.id == cell.id:
                                # pick second cell from list
                                cell_two = random.choice(stem_cell_list)

                    # now have two different cells if found_cell=True
                    #print("Cell 1:", cell.id, "Rho:", cell.dict['Rho'], "Cell 2:", cell_two.id, "Rho:",
                          #cell_two.dict['Rho'], )

                    # now have two different cells
                    # Check that the change is legal
                    if cell.dict['Rho'] > delta2:
                        # if change is illegal do nothing (alternatively could repick)
                        cell.dict['Rho'] -= delta2
                        # Keep track of number of mutations mutation_count is number for whole simulation
                        # number_of_mutations_per_step is number per step
                        mutation_count += 1
                        self.shared_steppable_vars['number_of_mutations_per_step'] +=1
                        # Update for the somatic cells
                        cell.dict['Generations_left'] = cell.dict['Rho']
                        cell_two.dict['Rho'] += delta2

                        # Update for the somatic cells
                        cell_two.dict['Generations_left'] = cell_two.dict['Rho']
                        #print("Made Mutation", "Cell 1:", cell.id, "Rho:", cell.dict['Rho'], "Cell 2:", cell_two.id,
                              #"Rho:", cell_two.dict['Rho'], )
                    else:
                        #print("Vetoed Mutation")
                        pass
        # add number of mutations per step to plot
        if not mcs%10: # Once per 10 MCS to save time                
            self.plot_win10.add_data_point("Number of Mutations per Step", mcs, self.shared_steppable_vars['number_of_mutations_per_step'])
            self.shared_steppable_vars['number_of_mutations_per_step']=0            
        # Actually do the division of cells                
        for cell in cells_to_divide:
            #Keeping track of lineage below may be expensive, check
            #cell.dict['Lineage'].append(cell.id)
            self.divide_cell_random_orientation(cell)

    def update_attributes(self):
        self.parentCell.targetVolume /= 2.0  # reducing parent target volume

        if self.get_steering_param('simulation_type')!="Neutral Variation":
            # in this case need to handle stem->somatic and somatic->necrotic cells

            if self.parentCell.type == self.SOMATIC:

                if self.parentCell.dict['Generations_left'] < 0:
                    # If somatic cell has run out of generations turn it to necrotic cells type
                    self.parentCell.type = self.NECROTIC
                    self.parentCell.targetVolume = 0.0
                    self.parentCell.lambdaVolume = 100.0

            self.clone_parent_2_child()

            if self.parentCell.type == self.STEM:
                # Need to differentiate offspring according to stemness
                self.childCell.type = self.SOMATIC
                # normally daughter is somatic. Maybe should make stemness 'Pstem' =0
                # If the parent is stem cell see what the offspring are
                # if stemness  'P' > 0.5 have chance both are stem, if < 0.5 both could be somatic

                if self.parentCell.dict['Pstem'] > 0.5:
                    # chance of daughter being stem
                    if random.random() < self.parentCell.dict['Pstem'] - 0.5:
                        self.childCell.type = self.STEM
                        print("Added Stem Cell")
                if self.parentCell.dict['Pstem'] < 0.5:
                    # chance of parent being somatic
                    if random.random() < 0.5 - self.parentCell.dict['Pstem']:
                        # Maybe should make stemness 'Pstem' =0
                        self.parentCell.type = self.SOMATIC
                        print("Eliminated Stem Cell")
        else:
            # If not doing somatic cells just clone parent
            self.clone_parent_2_child()


class DeathSteppable(SteppableBasePy):
    def __init__(self, frequency=1):
        SteppableBasePy.__init__(self, frequency)
    
    def start(self):   
        self.shared_steppable_vars['number_of_cells_killed_per_step'] = 0
        
        # Note--should have separate plot of the number dying due to senescence
        self.plot_win11 = self.add_new_plot_window(title='Number of Cells Killed Per Step',
                                                              x_axis_title='MonteCarlo Step (MCS)',
                                                              y_axis_title='Number of Cells Killed Per Step', x_scale_type='linear',
                                                              y_scale_type='linear',
                                                              grid=False)
        self.plot_win11.add_plot("Number of Cells Killed Per Step", style='Lines', color='green', size=5)

    def step(self, mcs):
        self.shared_steppable_vars['number_of_cells_killed_per_step'] = 0
        
        # ["None","Constant Killing All Cells","Constant Killing Somatic Only","Surgical","Radiation Simple",
        #    "Radiation Stem Protected","Radiation Quiescent Protected","Chemotherapy Simple","Chemotherapy Somatic Protected"]
                          
                                
        global killing_type
        killing_type = self.get_steering_param('killing_type')
       
        
        global Killing_efficacy
        Killing_efficacy = self.get_steering_param('Killing_efficacy')
        
        if killing_type=="Constant Killing Somatic Only" or killing_type=="Constant Killing All Cells":
            # if constant killing reduce kill rate by 1E4  a scale factor that we have to play with
            
            Killing_efficacy/=kill_reduction
            
        print("Killing type ",killing_type," Killing Efficacy ",Killing_efficacy)
            
        
             
        # Options--None
        if killing_type=="Select Type of Treatment" or killing_type=="None":
            pass # do nothing if no killing selected
        elif killing_type=="Constant Killing All Cells":
            # if constant killing loop over all cells
            for cell in self.cell_list_by_type(self.STEM, self.SOMATIC):
                if random.random()< Killing_efficacy: #Kill Cell
                    cell.type = self.NECROTIC
                    cell.targetVolume = 0.0
                    cell.lambdaVolume = 100.0
                    self.shared_steppable_vars['number_of_cells_killed_per_step']+=1
        elif killing_type=="Constant Killing Somatic Only":
            # if constant killing of somatic only loop over all somatic cells
            for cell in self.cell_list_by_type(self.SOMATIC):
                if random.random()< Killing_efficacy: #Kill Cell
                    cell.type = self.NECROTIC
                    cell.targetVolume = 0.0
                    cell.lambdaVolume = 100.0
                    self.shared_steppable_vars['number_of_cells_killed_per_step']+=1   
        elif killing_type=="Surgical":
            # Only do every 2000 MCS
            if not mcs%surgery_time: 
                # remove cells in left x percentage 
                for cell in self.cell_list:
                    if cell.xCOM < Killing_efficacy * 128: # Replace with x dim 
                        cell.type = self.NECROTIC
                        cell.targetVolume = 0.0
                        cell.lambdaVolume = 100.0
                        self.shared_steppable_vars['number_of_cells_killed_per_step']+=1
        elif killing_type=="Radiation Simple":
            # Only do every 2000 MCS
            if not mcs%radiation_time: 
                # remove x percentage of cells 
                for cell in self.cell_list_by_type(self.STEM, self.SOMATIC):
                    if random.random()< Killing_efficacy: #Kill Cell
                        cell.type = self.NECROTIC
                        cell.targetVolume = 0.0
                        cell.lambdaVolume = 100.0
                        self.shared_steppable_vars['number_of_cells_killed_per_step']+=1
        elif killing_type=="Radiation Stem Protected": # Could also add version where non-dividing cells are protected
            # Only do every 2000 MCS
            if not mcs%radiation_time: 
                # remove x percentage of cells favoring stem cells
                for cell in self.cell_list_by_type(self.SOMATIC):
                    if random.random()< Killing_efficacy: #Kill Somatic cells at higher rate
                        cell.type = self.NECROTIC
                        cell.targetVolume = 0.0
                        cell.lambdaVolume = 100.0 
                        self.shared_steppable_vars['number_of_cells_killed_per_step']+=1
                for cell in self.cell_list_by_type(self.SOMATIC):
                    if random.random()< Killing_efficacy*radiation_protection: #Kill Stem cells at lower rate
                        cell.type = self.NECROTIC
                        cell.targetVolume = 0.0
                        cell.lambdaVolume = 100.0 
                        self.shared_steppable_vars['number_of_cells_killed_per_step']+=1 
        elif killing_type=="Radiation Quiescent Protected": # Could also add version where non-dividing cells are protected
            # Only do every 2000 MCS
            if not mcs%radiation_time: 
                # remove x percentage of cells favoring stem cells
                for cell in self.cell_list_by_type(self.SOMATIC,self.STEM):
                    if cell.dict['Pressure']< Pressure_Threshold: # Kill cells which are dividing at higher rate
                        if random.random()< Killing_efficacy: #Kill Dividing cells at higher rate
                            cell.type = self.NECROTIC
                            cell.targetVolume = 0.0
                            cell.lambdaVolume = 100.0
                            self.shared_steppable_vars['number_of_cells_killed_per_step']+=1
                    else: # Kill cells which are quiescent at lower rate
                        if random.random()< Killing_efficacy*radiation_quiescence_protection: #Kill Quiescent cells at lower rate
                            cell.type = self.NECROTIC
                            cell.targetVolume = 0.0
                            cell.lambdaVolume = 100.0
                            self.shared_steppable_vars['number_of_cells_killed_per_step']+=1 
        else:
            print("Chemo Options not yet implemented")
            #Put the sources and diffusion fields in here
             
        self.plot_win11.add_data_point("Number of Cells Killed Per Step", mcs, self.shared_steppable_vars['number_of_cells_killed_per_step'])

class LineageTrackerSteppable(SteppableBasePy):
    def __init__(self, frequency=100):
        SteppableBasePy.__init__(self, frequency)

        #self.create_scalar_field_cell_level_py("LineageTracker")
        self.track_cell_level_scalar_attribute(field_name='Pressure', attribute_name='Pressure')
        self.track_cell_level_scalar_attribute(field_name='Generations_left', attribute_name='Generations_left')
        self.track_cell_level_scalar_attribute(field_name='Stemness', attribute_name='Pstem')
        self.track_cell_level_scalar_attribute(field_name='Generations_to_senescence', attribute_name='Rho')
        self.track_cell_level_scalar_attribute(field_name='Generations_left', attribute_name='Generations_left')
        #self.track_cell_level_scalar_attribute(field_name='Mutation_id', attribute_name='Current_Mutation')

        # and any other evolving parameters
        # cell.dict['Mutations'].append([mutation_count, mcs, cell.dict['Pstem'], cell.dict['Rho']])

    def start(self):
        print("LineageTrackerSteppable: This function is called once before simulation")

