from DFTGroundStateCalculation import DFTGroundStateCalculation
from EhrenfestStoppingSimulation import EhrenfestStoppingSimulation
from Trajectory import Trajectory, HyperchannellingTrajectory
import numpy as np


        #################################
        ## HYPERCHANNELLING TRAJECTORY ##
        #################################


# trajectory = HyperchannellingTrajectory()



        ############################
        ## PRESAMPLED1 TRAJECTORY ##
        ############################

projectile_starting_position = [3.83534014, 2.7002751, 0.44542105]
# TRANSLATE FROM PRESAMPLING CODE TO UNCENTERED SYSTEM
projectile_starting_position = [i - 1.0125 for i in projectile_starting_position]

initial_direction = [-0.39661179, -0.60274438, 0.69238594]

trajectory = Trajectory(projectile_starting_position, initial_direction)


        #####################
        ## DFT CALCULATION ##
        #####################

# supercell_size = (6, 2, 2)
supercell_size = (3, 3, 3)



calculation = DFTGroundStateCalculation(supercell_size, trajectory)
calculation.run_calculation()

        ###########################
        ## EHRENFEST CALCULATION ##
        ###########################

name = "Al_stopping"
# kinetic_energies = np.array([1, 20, 40, 60, 80, 100, 120]) * 1e3
kinetic_energies = np.array([400]) * 1e3

simulation = EhrenfestStoppingSimulation(name, trajectory, kinetic_energies)
simulation.run_simulation()
