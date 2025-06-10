from ase.geometry import distance
import numpy as np


class FitInformation:
    def __init__(self):
        self.fit = None
        self.cov = None
        self.crop = [None, None]



class DataProcessor:
    def __init__(self):
        pass


    def calculate_stopping_powers(self, atoms_dict, crop=[None, None]):
        """
        Parameters:
        atoms_dict (Dict[str, List[Atoms]])
        calc_dict (Dict[str, List[GPAW]])
        """

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

            initial_position = positions[0]
            distance_travelled = np.array([np.linalg.norm(position - initial_position) for position in positions])

            distance_per_timestep = np.diff(distance_travelled)[0]
            if crop[1] is not None:
                if i == 0:  # print debug info only first time
                    print(f"cropping {crop[1]} timesteps off the end = {crop[1] * distance_per_timestep} Angstroms")
                distance_traveleld = distance_travelled[:crop[1]]
                kinetic_energies = kinetic_energies[:crop[1]]
            if crop[0] is not None:
                if i == 0:
                    print(f" cropping {crop[0]} timesteps off the start = {crop[0]*distance_per_timestep} Angstroms")
                distance_travelled = distance_travelled[crop[0]:]
                kinetic_energies = kinetic_energies[crop[0]:]
            if i == 0:
                print(f"number of remaining timesteps: {len(distance_travelled)} = {len(distance_travelled) * distance_per_timestep} Angstroms")

            fit, cov = np.polyfit(distance_travelled, kinetic_energies, 1, cov=True)
            fits[energy] = fit
            covs[energy] = cov
        return fits, covs

