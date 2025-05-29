from gpaw import restart

from ase.io import write

import numpy as np
import matplotlib.pyplot as plt
# from mpl_toolkits.mplot3d import Axes3D
from matplotlib.widgets import Slider

import utils

import os
import re
import json


class StoppingPowerAnalysis:
    def __init__(self, data_directory):

        # get dictionary of filenames
        if not data_directory.endswith("/"):
            data_directory += "/"
        self.data_directory = data_directory
        all_gpw_files, self.filenames = self.extract_gpaw_files()

        print(self.data_directory)
        print(json.dumps(self.filenames, indent=4))

        # TODO: if theres multiple energies, and you write for each of them
        # TODO: atm it just overwrites the same file multiple times
        self.trajectory_file = "output.xyz"

        self.atoms_dict = {}
        self.calc_dict = {}
        for energy in self.filenames.keys():
            self.load_data(energy)

        write = True
        if write:
            energy_to_write = "40 keV"
            self.write_data(energy_to_write)

        self.kinetic_energies = self.get_projectile_kinetic_energies()
        self.projectile_positions = self.get_projectile_positions()

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
        for energy in filenames_temp.keys():
            paired = sorted(
                zip(map(int, timesteps_temp[energy]), filenames_temp[energy]))

            _, filenames_sorted = zip(*paired)
            filenames_temp[energy] = list(filenames_sorted)

        # rename keys
        filenames = {f"{key} keV": value for key, value in filenames_temp.items()}

        return all_gpw_files, filenames

    def load_data(self, energy):
        """
        Parameters:
        energy (string): "__ keV"
        """

        for filename in self.filenames[energy]:
            atoms, calc = restart(self.data_directory + filename)
            self.atoms_dict = utils.append_to_dict(self.atoms_dict, energy, atoms)

            self.calc_dict = utils.append_to_dict(self.calc_dict, energy, calc)

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

    def get_projectile_kinetic_energies(self):
        kinetic_energies = {energy: [atoms.get_kinetic_energy()
                                     for atoms in atoms_list]
                            for energy, atoms_list in self.atoms_dict.items()}
        return kinetic_energies

    def get_projectile_positions(self):
        projectile_positions = {energy: [atoms.get_positions()[-1][0]
                                         for atoms in atoms_list]
                                for energy, atoms_list in self.atoms_dict.items()}
        return projectile_positions

    def plot_kinetic_energies(self):
        """
        plots projectile kinetic energy over the projectile trajectory
        performs a sliding window straight line fit (r^2 minimisation) to fit
        a straight line to the region where the projectile passes through the supercell
        extracts stopping power from this fit
        """

        n_subplots = len(self.kinetic_energies.keys())
        fig, axs = plt.subplots(n_subplots, figsize=(10, 5*n_subplots), sharex=True)
        # so that is works for the =1 case
        if n_subplots == 1:
            axs = [axs]

        fig.suptitle("Fitting KE to extract stopping powers")

        axs[-1].set_xlabel(r"projectile position [$\AA$]")
        _ = [ax.set_ylabel("kinetic energy [keV]") for ax in axs]

        for i in range(n_subplots):
            energy, kinetic_energies = list(self.kinetic_energies.items())[i]
            kinetic_energies = np.array(
                kinetic_energies) * 1e-3  # convert from eV to keV
            _, projectile_positions = list(
                self.projectile_positions.items())[i]

            ##################################
            # PERFORM SLIDING WINDOW FITTING #
            #     and get stopping power
            ##################################

            # TODO: should be able to calculat the minimum window size from the size of the unit cell and from the distance travelled per timestep
            fit, cov, x_window, y_window = utils.sliding_fit(projectile_positions, kinetic_energies, 3, 5)
            stopping_power = -fit[0]
            stopping_power_uncertainty = np.sqrt(cov[0][0])

            axs[i].set_title(f"Projectile Initial energy: {energy}")
            axs[i].plot(projectile_positions, kinetic_energies, "x")
            axs[i].plot(x_window, y_window, "x", color="red")
            axs[i].plot(x_window, np.poly1d(fit)(x_window), color="red",
                        label=rf"$S_e$ = {1e3*stopping_power:.1f} $\pm$ {1e3*stopping_power_uncertainty:.1f} [eV/$\AA$]")

            ############################################
            # SHOWING ON THE PLOT WHERE THE LATTICE IS #
            ############################################

            # can get atom positions from any of the timesteps becuase lattice positions dont change
            # [0] -> first timestep, [:-1] -> exclude the projectile
            lattice_positions = self.atoms_dict[energy][0].get_positions()[:-1]
            lattice_start = min(lattice_positions[:, 0])
            lattice_end = max(lattice_positions[:, 0])
            axs[i].axvspan(lattice_start, lattice_end,
                           color="yellow", alpha=0.25, label="nuclei positions")

        _ = [ax.legend() for ax in axs]
        plt.show()

    def visualise_electron_density(self, energy, timestep):
        """
        produces an interactive 3D plot of electron density
        """

        timestep_index = int(timestep / 10 - 1)
        atoms = self.atoms_dict[energy][timestep_index]

        print(f"""

        DEBUG: number of atoms: {len(atoms)}

        """)

        calc = self.calc_dict[energy][timestep_index]
        electron_density = calc.get_all_electron_density()
        atom_positions = atoms.get_positions()
        lattice_positions = atom_positions[:-1]
        projectile_positions = atom_positions[-1]

        # setup density slice things
        nx, ny, nz = np.shape(electron_density)
        xx, yy = np.meshgrid(np.linspace(0, 1, nz), np.linspace(0, 1, ny))
        X = xx
        Y = yy

        # setup plotting atom position things
        # want to normalise atom positions to unit cube using atoms.get_cell()
        cell_size = atoms.get_cell()
        fcc_xs = lattice_positions[:, 0] / (cell_size[0, 0])
        fcc_ys = lattice_positions[:, 1] / (cell_size[1, 1])
        fcc_zs = lattice_positions[:, 2] / (cell_size[2, 2])

        proton_xs = projectile_positions[0] / (cell_size[0, 0])
        proton_ys = projectile_positions[1] / (cell_size[1, 1])
        proton_zs = projectile_positions[2] / (cell_size[2, 2])

        fig = plt.figure(figsize=(20, 10))
        ax = fig.add_subplot(111, projection='3d')
        ax.set_zlim((0, 1))

        # Create the slider
        ax_slider = plt.axes([0.15, 0.02, 0.7, 0.03],
                             facecolor='lightgoldenrodyellow')
        slice_location = 0
        slider = Slider(ax_slider, 'Frame', 0, nx - 1,
                        valinit=slice_location, valstep=1)

        vmax = np.log10(np.max(electron_density))
        vmin = vmax - 8

        def animate(f):
            ax.clear()
            slice = electron_density[f, :, :]
            ax.contourf(X, Y, np.log10(slice), zdir='z',offset=f/nx, cmap="jet", vmax=vmax)
            ax.scatter(fcc_ys, fcc_zs, fcc_xs, s=500, alpha=1)
            ax.scatter(proton_ys, proton_zs, proton_xs, s=300, alpha=1)
            ax.set_zlim((0, 1))

        # dont need to call FuncAnimation when you have the slider going
        slider.on_changed(animate)
        plt.show()


if __name__ == "__main__":
    data_directory = "/Users/brynlloyd/Developer/Coding/Python/dft/gpaw/my_own_stopping/data/larger_unitcell"
    analysis = StoppingPowerAnalysis(data_directory)
    # analysis.visualise_electron_density("40 keV", 40)
