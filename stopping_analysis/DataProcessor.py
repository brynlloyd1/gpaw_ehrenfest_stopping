import numpy as np


class DataProcessor:
    def __init__(self):
        pass


    def calculate_stopping_powers(self, atoms_dict, something_about_cropping_data=None):
        """
        Parameters:
        atoms_dict (Dict[str, List[Atoms]])
        calc_dict (Dict[str, List[GPAW]])
        """

        # actually duplicate code from the plotting
        # can this be a function somewhere?
        projectile_positions = {energy: [atoms.get_positions()
                                         for atoms in atoms_list]
                                for energy, atoms_list in atoms_dict.items()}

        projectile_kinetic_energies = {energy: [atoms.get_kinetic_energy()
                                                for atoms in atoms_list]
                                       for energy, atoms_list in atoms_dict.items()}

        fits = {}
        covs = {}

        for i in range(len(projectile_kinetic_energies.keys())):
            # get raw data
            energy, kinetic_energies = list(projectile_kinetic_energies.items())[i]
            kinetic_energies = np.array(kinetic_energies) * 1e-3  # convert from eV to keV
            _, positions = list(projectile_positions.items())[i]
            positions = np.array(positions)[:,-1,0]   # for every timestep, want projectile x_pos

            fit, cov = np.polyfit(positions, kinetic_energies, 1, cov=True)
            fits[energy] = fit
            covs[energy] = cov

        return fits, covs
