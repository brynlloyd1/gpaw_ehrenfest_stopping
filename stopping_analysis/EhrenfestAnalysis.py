from DataLoader import DataLoader
from DataProcessor import DataProcessor
from DataVisualiser import DataVisualiser

import os

class DataHandler:
    def __init__(self, directory):
        self.trajectory_name = os.path.basename(directory.rstrip("/"))
        self.data_loader = DataLoader(directory)
        self.data_processor = DataProcessor()
        self.data_visualiser = DataVisualiser()

        self.atoms_dict = {}
        self.calc_dict = {}
        self.fits = {}

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

    def load_data(self, trajectory_name):
        data_handler = self.data_handlers[trajectory_name]
        data_handler.atoms_dict, data_handler.calc_dict = data_handler.data_loader.load()


    def calculate_stopping_curve(self, trajectory_name):
        handler = self.data_handlers[trajectory_name]
        processor = handler.data_processor
        handler.fits = processor.calculate_stopping_curve()

    def view_fits(self, trajectory_name):
        handler = self.data_handlers[trajectory_name]
        visualiser = handler.data_visualiser
        visualiser.plot_all_fits(handler.atoms_dict, handler.fits)

    def visualise_electron_density(self, trajectory_name, energy):
        handler = self.data_handlers[trajectory_name]
        visualiser = handler.data_visualiser
        electron_density = [calc.get_all_electron_density() for calc in handler.calc_dict[energy]]
        visualiser.visualise_electron_density(electron_density)









    # def load_trajectory(self, directory, energies=["all"], timesteps="all"):
    #     """
    #     uses gpaw.restart() to load Atoms and GPAW objects from .gpw files
    #     these objects are then stored in nested dictionaries/lists as follows
    #
    #     self.data {
    #         "trajectory 1" : [atoms_dict {
    #                             "20 keV" : [Atoms],
    #                             "40 keV" : [Atoms],
    #                          },
    #                          calc_dict {
    #                             "20 keV" : [GPAW],
    #                             "40 keV" : [GPAW],
    #                          }]
    #     }
    #
    #
    #     Parameters:
    #     directory (str): path to directory containing .gpw files for the ehrenfest stopping trajectory
    #     energies (List[str]): list of energies to be loaded, default is ["all"] to load in all energy files
    #     timesteps (str): which timesteps to load. Options are "all" (default), "::10" (for every 10'th)
    #     """
    #
    #     trajectory_name = os.path.basename(directory)
    #     atoms_dict, calc_dict = self.data_loader.load(directory, energies, timesteps)
    #
    #     self.data[trajectory_name] = [atoms_dict, calc_dict]
    #
    #
    # def calculate_stopping_powers(self, trajectory_name, which_energies=["all"]):
    #     """
    #     calculates stopping power curve for a given trajectory
    #     by default, calculates entire stopping curve,
    #     but can optionally only calculate stopping power at specified projectile energies
    #
    #     Parameters:
    #     trajectory_name (str): must be one of the keys of self.data
    #     which_energies (list[str]): list of energies at which to calculate stopping power. Default is all energies
    #
    #     Returns:
    #     """
    #     ##################
    #     # PROCESS INPUTS #
    #     ##################
    #
    #     # verify that trajectory_name is one of the keys of self.data
    #     if not trajectory_name in self.data.keys():
    #         raise KeyError(f"{trajectory_name} not found in self.data.keys()")
    #
    #     if which_energies == ["all"]:
    #         which_energies = self.data[trajectory_name].keys()
    #
    #     # check that all specified energies are in all_files.keys()
    #     for energy in which_energies:
    #         if energy in self.data[trajectory_name].keys():
    #             continue
    #         else:
    #             raise KeyError(f"{energy} not found in this self.data[{trajectory_name}]")
    #
    #     ############################
    #     # CALL SELF.DATA_PROCESSOR #
    #     ############################
    #
    #     # only want to pass specified energies into the DataProcessor
    #     atoms_list = {key:value for key, value in self.data[trajectory_name][0] if key in which_energies}
    #     calc_list = {key:value for key, value in self.data[trajectory_name][1] if key in which_energies}
    #
    #     initial_projectile_energies, stopping_powers = self.data_processor.calculate_stopping_powers(atoms_list, calc_list)
    #
    #     return initial_projectile_energies, stopping_powers
    #
    #


if __name__ == "__main__":
    directory = "/Users/brynlloyd/Developer/Coding/Python/dft/gpaw/my_own_stopping/data/322_supercell"
    analysis = EhrenfestAnalysis()
    trajectory_name = analysis.initialise_analysis(directory)
    analysis.load_data(trajectory_name)

    # analysis.view_fits(trajectory_name)
    # analysis.calculate_stopping_curve(trajectory_name)
    # analysis.visualise_electron_density(trajectory_name)














