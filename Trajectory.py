import numpy as np

class Trajectory:
    def __init__(self, starting_position, direction):
        self.starting_position = starting_position
        # normalise direction so that it is a unit vector
        direction = list(np.array(direction) / np.linalg.norm(np.array(direction)))
        self.direction = direction


class HyperchannellingTrajectory(Trajectory):
    def __init__(self):
        Al_lattice_constant = 4.05   # Angstroms
        hyperchannelling_starting_position = [Al_lattice_constant,
                                              Al_lattice_constant,
                                              Al_lattice_constant]
        hyperchannelling_direction = [1, 0, 0]
        super().__init__(hyperchannelling_starting_position,
                         hyperchannelling_direction)

