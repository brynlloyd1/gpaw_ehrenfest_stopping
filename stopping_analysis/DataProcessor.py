import numpy as np


class DataProcessor:
    def __init__(self):
        pass


    def calculate_stopping_powers(self, atoms_dict, calc_dict):
        """
        Parameters:
        atoms_dict (Dict[str, List[Atoms]])
        calc_dict (Dict[str, List[GPAW]])
        """




    # def get_projectile_kinetic_energies(self):
    #     kinetic_energies = {energy: [atoms.get_kinetic_energy()
    #                                  for atoms in atoms_list]
    #                         for energy, atoms_list in self.atoms_dict.items()}
    #     return kinetic_energies
    #
    # def get_projectile_positions(self):
    #     projectile_positions = {energy: [atoms.get_positions()[-1][0]
    #                                      for atoms in atoms_list]
    #                             for energy, atoms_list in self.atoms_dict.items()}
    #     return projectile_positions


