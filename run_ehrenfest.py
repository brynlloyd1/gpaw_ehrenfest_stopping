from EhrenfestStoppingSimulation import EhrenfestStoppingSimulation
from Trajectory import Trajectory, HyperchannellingTrajectory
import numpy as np

name = "Al_stopping"
# kinetic_energies = np.array([1, 20, 40, 60, 80, 100, 120]) * 1e3
kinetic_energies = np.array([400]) * 1e3




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

        ############################
        ## PRESAMPLED2 TRAJECTORY ##
        ############################

#presampled2 trajectory
# projectile_starting_position = [1.8873831, 3.78791897, 0.14010488]
# initial_direction = [0.47792543, 0.84844362, 0.22744388]


        ################################
        ##  VOLUME CAPTURE TRAJECTORY ##
        ################################

# trajectory such that projectile makes an angle wrt channelling of 1 degree
# and it crosses the hyperchannelling trajectory after travelling a half-length of the supercell
# angle = 1 * np.pi/180
# offset = 3 * np.tan(angle)
# projectile_starting_position = [4.05/2, 4.05 + offset, 4.05]
# initial_direction = [1, np.tan(angle), 0]

trajectory = Trajectory(projectile_starting_position, initial_direction)

simulation = EhrenfestStoppingSimulation(name, trajectory, kinetic_energies)
simulation.run_simulation()
