from DataHandler import DataHandler

import numpy as np

class EhrenfestAnalysis:
    def __init__(self):
        self.data_handlers = {}

    def initialise_analysis(self, directory):
        data_handler = DataHandler(directory)
        self.data_handlers[data_handler.trajectory_name] = data_handler

        return data_handler.trajectory_name

    def set_name(self, old_name, new_name):
        if new_name in self.data_handlers.keys():
            raise KeyError(f"{new_name} already exists")

        self.data_handlers[old_name].trajectory_name = new_name
        self.data_handlers[new_name] = self.data_handlers[old_name].pop()

    def set_energies(self, trajectory_name, which_energies):
        self.data_handlers[trajectory_name].data_loader.set_which_energies(which_energies)

    def set_timesteps(self, trajectory_name, which_timesteps):
        self.data_handlers[trajectory_name].data_loader.set_which_timesteps(which_timesteps)

    def load_gpw_data(self, trajectory_name):
        data_handler = self.data_handlers[trajectory_name]
        data_handler.atoms_dict, data_handler.calc_dict = data_handler.data_loader.load()

    def save_electron_density_to_npy(self, trajectory_name, energy):
        data_handler = self.data_handlers[trajectory_name]
        electron_density_array = np.array([calc.get_all_electron_density() for calc in data_handler.calc_dict[energy]])
        data_handler.data_loader.save_electron_density_to_npy(trajectory_name, energy, electron_density_array)

    def read_electron_density_from_npy(self, trajectory_name, energy):
        data_handler = self.data_handlers[trajectory_name]
        print("started loading")
        electron_density_array = data_handler.data_loader.read_electron_density_from_npy(trajectory_name, energy)
        print("finished loading")
        self.visualise_electron_density_change(trajectory_name, energy)

    def calculate_stopping_curve(self, trajectory_name):
        handler = self.data_handlers[trajectory_name]
        processor = handler.data_processor
        fits,covs = processor.calculate_stopping_powers(handler.atoms_dict)
        handler.fits[trajectory_name] = fits
        handler.covs[trajectory_name] = covs


    def view_fits(self, trajectory_name):
        handler = self.data_handlers[trajectory_name]
        visualiser = handler.data_visualiser

        # if DataProcessor.calculate_stopping_powers() has not been run for this trajectory
        # then the function can still be used to view the kinetic energy data
        if handler.fits[trajectory_name]:
            params = [handler.atoms_dict,
                     handler.fits[trajectory_name],
                     handler.covs[trajectory_name]]
        else:
            params = [handler.atoms_dict, None, None]

        visualiser.plot_all_fits(*params)

    def visualise_electron_density(self, trajectory_name, energy):
        handler = self.data_handlers[trajectory_name]
        visualiser = handler.data_visualiser
        electron_density_list = [calc.get_all_electron_density() for calc in handler.calc_dict[energy]]
        visualiser.visualise_electron_density(electron_density_list)

    def visualise_electron_density_change(self, trajectory_name, energy, electron_density_list=None):
        handler = self.data_handlers[trajectory_name]
        visualiser = handler.data_visualiser

        if not electron_density_list:
            electron_density_list = [calc.get_all_electron_density() for calc in handler.calc_dict[energy]]
        visualiser.visualise_electron_density_change(electron_density_list)





if __name__ == "__main__":
    directory = "/Users/brynlloyd/Developer/Coding/Python/dft/gpaw/my_own_stopping/data/322_supercell/"
    analysis = EhrenfestAnalysis()
    trajectory_name = analysis.initialise_analysis(directory)
    analysis.set_timesteps(trajectory_name, "::10")
    analysis.load_gpw_data(trajectory_name)

    # analysis.save_electron_density_to_npy(trajectory_name, "40 keV")

    analysis.read_electron_density_from_npy(trajectory_name, "40 keV")
    # analysis.visualise_electron_density_change(trajectory_name, "40 keV")














