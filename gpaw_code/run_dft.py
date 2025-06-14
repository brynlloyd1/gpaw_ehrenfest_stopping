from DFTGroundStateCalculation import DFTGroundStateCalculation
from Trajectory import Trajectory, HyperchannellingTrajectory

trajectory = HyperchannellingTrajectory()
supercell_size = (6, 2, 2)
calculation = DFTGroundStateCalculation(supercell_size, trajectory)
calculation.run()
