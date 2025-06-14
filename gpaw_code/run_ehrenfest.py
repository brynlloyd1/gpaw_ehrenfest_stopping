from EhrenfestStoppingSimulation import EhrenfestStoppingSimulation
from Trajectory import Trajectory, HyperchannellingTrajectory
import numpy as np

name = "Al_stopping"
kinetic_energies = np.array([20, 40, 60, 80])
trajectory = HyperchannellingTrajectory()

simulation = EhrenfestStoppingSimulation(name, trajectory, kinetic_energies)
simulation.run_simulation()
