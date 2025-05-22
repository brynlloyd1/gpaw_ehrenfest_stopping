from gpaw import restart

from ase.io import write

import numpy as np
import matplotlib.pyplot as plt

import utils

import os
import re
import json


class StoppingPowerAnalysis:
    def __init__(self, data_directory, energy=40):

        # get dictionary of filenames
        self.data_directory = data_directory
        self.filenames = self.extract_gpaw_files()
        print(json.dumps(self.filenames, indent=4))
        self.trajectory_file = "output.xyz"

        self.atoms_dict = {}
        for energy in self.filenames.keys():
            self.load_data(energy)

        write = True
        if write:
            energy_to_write = "40 keV"
            self.write_data(energy_to_write)

        self.kinetic_energies = {energy: [atoms.get_kinetic_energy()
                                          for atoms in atoms_list]
                                 for energy, atoms_list in self.atoms_dict.items()}

    def extract_gpaw_files(self):
        """
        gets .gpw files from a given directory

        Returns:
        dict: key is energy (in keV), value is a list of filenames, sorted by ascending timestep
        """

        # get all files, and throw any that arent .gpw
        all_files = os.listdir(self.data_directory)
        all_gpw_files = [f for f in all_files if f.endswith(".gpw")]

        filenames_temp = {}
        timesteps_temp = {}

        # regex to extract energy and timestep
        energy_timestep_pattern = re.compile(r"(\d+)k_step(\d+)")
        for filename in all_gpw_files:
            match = energy_timestep_pattern.search(filename)
            if not match:
                continue
            energy, timestep = match.group(1), match.group(2)

            filenames_temp = utils.append_to_dict(
                filenames_temp, energy, filename)
            timesteps_temp = utils.append_to_dict(
                timesteps_temp, energy, timestep)

        # sort files in ascending order according to the timestep
        # yes ik its sorting strings numerically but it works so idc
        for energy in filenames_temp.keys():
            paired = sorted(
                zip(timesteps_temp[energy], filenames_temp[energy]))

            _, filenames_sorted = zip(*paired)
            filenames_temp[energy] = list(filenames_sorted)

        # rename keys
        filenames = {f"{key} keV": value for key,
                     value in filenames_temp.items()}

        return filenames

    def load_data(self, energy):
        """
        Parameters:
        energy (string): "__ keV"
        """

        for filename in self.filenames[energy]:
            atoms, calc = restart(self.data_directory + filename)
            self.atoms_dict = utils.append_to_dict(
                self.atoms_dict, energy, atoms)

    def write_data(self, energy):
        """
        writes data to a .xyz file so that it can be visualised in ovito or something

        Parameters:
        energy (string): "__ keV"
        """

        if energy not in self.atoms_dict:
            raise ValueError("nothing to write. self.atoms is an empty list")

        for i, atoms in enumerate(self.atoms_dict[energy]):
            write(self.trajectory_file, atoms, append=i > 0)

    def plot_kinetic_energies(self):
        """
        extract kinetic energy from the atoms instances stored inside self.atoms_dict[energy]
        then plot, with a subplots for each energy
        """

        n_subplots = len(self.kinetic_energies.keys())
        fig, axs = plt.subplots(n_subplots, figsize=(10, 5*n_subplots))

        for i, (energy, kinetic_energies) in enumerate(self.kinetic_energies.items()):
            axs[i].plot(kinetic_energies)

        plt.show()


if __name__ == "__main__":

    data_directory = "/Users/brynlloyd/Developer/Coding/Python/dft/gpaw/data/small_unitcell/"
    # data_directory = "/Users/brynlloyd/Developer/Coding/Python/dft/gpaw/data/larger_unitcell/"

    analysis = StoppingPowerAnalysis(data_directory)
    analysis.plot_kinetic_energies()
