from DataHandler import DataHandler

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

    def calculate_stopping_curve(self, trajectory_name, crop=[30, None]):
        handler = self.data_handlers[trajectory_name]
        processor = handler.data_processor
        fits_information = processor.calculate_stopping_powers(handler.atoms_dict, crop=crop)
        handler.fits = fits_information


    def view_fits(self, trajectory_name):
        handler = self.data_handlers[trajectory_name]
        visualiser = handler.data_visualiser
        visualiser.plot_all_fits(handler.atoms_dict, handler.fits)


    def compare_to_geant4(self, trajectory_names):
        stopping_power_data = {}
        for trajectory_name in trajectory_names:
            handler = self.data_handlers[trajectory_name]
            energies = []
            stopping_powers = []

            for energy, fit_info in handler.fits.items():
                energies.append(int(energy.rstrip(" keV")))
                stopping_powers.append(-fit_info.fit[0]*1e3)   # fits[0] is in keV/Angstrom so convert to eV/Ang

            stopping_power_data[trajectory_name] = {"energies": energies,
                                                    "stopping powers": stopping_powers}
        self.data_handlers[trajectory_names[0]].data_visualiser.geant4_comparison(stopping_power_data)


    def visualise_electron_density(self, trajectory_name, energy):
        handler = self.data_handlers[trajectory_name]
        loader = handler.data_loader
        visualiser = handler.data_visualiser

        if loader.check_for_npy(handler.directory_name, energy):
            electron_density_list = loader.load_from_npy(handler.directory_name, energy)
        else:
            electron_density_list = [calc.get_all_electron_density() for calc in handler.calc_dict[energy]]
            loader.save_to_npy(handler.directory_name, energy, electron_density_list)

        visualiser.visualise_electron_density(electron_density_list)




if __name__ == "__main__":
    directory = "/Users/brynlloyd/Developer/Coding/Python/dft/gpaw/my_own_stopping/data/322_presampled1/"
    analysis = EhrenfestAnalysis()
    trajectory_name = analysis.initialise_analysis(directory)
    analysis.set_timesteps(trajectory_name, "all")
    analysis.load_gpw_data(trajectory_name)
    analysis.visualise_electron_density(trajectory_name, "40 keV")











