from gpaw import restart

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, CheckButtons

import utils
import os
import re
import json


class StoppingPowerAnalysis:
    def __init__(self, data_directory):

        # get dictionary of filenames
        if not data_directory.endswith("/"): data_directory += "/"
        self.data_directory = data_directory
        self.filenames = self.extract_gpaw_files()
        print(self.data_directory)
        print(json.dumps(self.filenames, indent=4))

        self.atoms_dict = {}
        self.calc_dict = {}
        for energy in self.filenames.keys():
            self.load_data(energy)

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

            filenames_temp = utils.append_to_dict(filenames_temp, energy, filename)
            timesteps_temp = utils.append_to_dict(timesteps_temp, energy, timestep)

        # sort files in ascending order according to the timestep
        for energy in filenames_temp.keys():
            paired = sorted(zip(map(int, timesteps_temp[energy]), filenames_temp[energy]))
            _, filenames_sorted = zip(*paired)
            filenames_temp[energy] = list(filenames_sorted)

        # rename keys
        filenames = {f"{key} keV": value for key, value in filenames_temp.items()}
        return filenames

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


    def plot_stopping_data(self):
        """
        for every projectile energy, creates 2 plots
        1. projectile kinetic energy against distance travelled
        2. gradient of kinetic energy against distance travelled

        TODO: want to extract stopping powers
        """
        #################
        # CREATE FIGURE #
        #################

        n_subplots = len(self.kinetic_energies.keys())
        fig,axs = plt.subplots(n_subplots, 2, figsize=(15, 5*n_subplots), sharex="col")
        if n_subplots == 1:
            print(1)
            axs = np.array([axs])

        fig.suptitle("Stopping Power Things")
        _ = [ax.set_ylabel("Kinetic Energy [keV]") for ax in axs[:,0]]
        _ = [ax.set_ylabel(r"$\frac{dKE}{dx}$ [keV/$\AA$]") for ax in axs[:,1]]
        _ = [ax.set_xlabel(r"Projectile Position [$\AA$]") for ax in axs[n_subplots-1]]

        #######################################
        # ITERATE OVER EACH PROJECTILE ENERGY #
        #######################################

        for i in range(n_subplots):

            # get raw data
            energy, kinetic_energies = list(self.kinetic_energies.items())[i]
            kinetic_energies = np.array(kinetic_energies) * 1e-3  # convert from eV to keV
            _, projectile_positions = list(self.projectile_positions.items())[i]

            # get lattice positions, and highlight region where there are nuclei
            lattice_positions = self.atoms_dict[energy][0].get_positions()[:-1]
            lattice_start = min(lattice_positions[:, 0])
            lattice_end = max(lattice_positions[:, 0])
            _ = [ax.axvspan(lattice_start, lattice_end, color="yellow", alpha=0.25, label="nuclei positions") for ax in axs.flatten()]

            # get the data specifically for when the proton is inside the supercell
            timesteps_in_supercell = np.where((projectile_positions > lattice_start) & (projectile_positions < lattice_end))[0]
            projectile_positions_supercell = projectile_positions[timesteps_in_supercell[0] : timesteps_in_supercell[-1]]
            kinetic_energies_supercell = kinetic_energies[timesteps_in_supercell[0] : timesteps_in_supercell[-1]]

            # plot kinetic energy on left plot
            axs[i, 0].plot(projectile_positions, kinetic_energies, "x")
            axs[i, 0].plot(projectile_positions_supercell, kinetic_energies_supercell, "x", color="red")

            # plot gradient of kinetic energy on right plot
            axs[i, 1].plot(projectile_positions[:-1], np.diff(kinetic_energies))
            axs[i, 1].plot(projectile_positions_supercell[:-1], np.diff(kinetic_energies_supercell), "x", color="red")


    def visualise_ehrenfest(self, energy):
        """
        creates animation using matplotlib sliders to display 2d slice of electron density
        has sliders to control slice position, timestep, vmin&vmax
        can also toggle logarithmic plotting
        """
        atoms_list = self.atoms_dict[energy]
        # nuclei_positions_list = [atoms.get_positions() for atoms in atoms_list]
        calc_list = self.calc_dict[energy]
        electron_density_list = [calc.get_all_electron_density() for calc in calc_list]

        timesteps = len(atoms_list)
        shape = np.shape(electron_density_list)

        # Initial values
        init_timestep = 0
        init_slice = shape[1] // 2
        use_log = False

        # Compute global vmin/vmax for raw and log data
        all_data = np.stack(electron_density_list)
        vmin_raw = all_data.min()
        vmax_raw = all_data.max()
        log_data = np.log10(np.clip(all_data, 1e-10, None))
        vmin_log = log_data.min()
        vmax_log = log_data.max()

        # Setup plot
        fig, ax = plt.subplots(figsize=(8, 6))
        plt.subplots_adjust(left=0.25, bottom=0.45)  # Extra room for sliders

        def get_plot_data(t, slice, log=False):
            data = electron_density_list[t][:, slice, :]
            data = np.rot90(data)
            if log:
                data = np.log10(np.clip(data, 1e-10, None))
            return data

        # Initial image
        plot_data = get_plot_data(init_timestep, init_slice, use_log)
        im = ax.imshow(plot_data, aspect='auto', vmin=vmin_raw, vmax=vmax_raw)
        cb = fig.colorbar(im, ax=ax)
        title = ax.set_title(f"Timestep: {init_timestep}, Slice: {init_slice}, Log: {use_log}")

        # Axes for sliders and checkbox
        ax_timestep = plt.axes((0.25, 0.35, 0.65, 0.03))
        ax_slice = plt.axes((0.25, 0.3, 0.65, 0.03))
        ax_vmin = plt.axes((0.25, 0.2, 0.65, 0.03))
        ax_vmax = plt.axes((0.25, 0.15, 0.65, 0.03))
        ax_check = plt.axes((0.05, 0.5, 0.15, 0.1))

        # Create sliders and checkbox
        slider_timestep = Slider(ax_timestep, 'Timestep', 0, len(electron_density_list) - 1, valinit=init_timestep, valstep=1)
        slider_slice = Slider(ax_slice, 'Slice', 0, shape[1] - 1, valinit=init_slice, valstep=1)
        slider_vmin = Slider(ax_vmin, 'vmin', vmin_raw, vmax_raw, valinit=vmin_raw)
        slider_vmax = Slider(ax_vmax, 'vmax', vmin_raw, vmax_raw, valinit=vmax_raw)
        check = CheckButtons(ax_check, ['Log scale'], [use_log])

        # update called when a slider is moved
        def update(val):
            t = int(slider_timestep.val)
            s = int(slider_slice.val)
            log_flag = check.get_status()[0]

            # Determine color range
            if log_flag:
                current_vmin = slider_vmin.val
                current_vmax = slider_vmax.val
            else:
                current_vmin = slider_vmin.val
                current_vmax = slider_vmax.val

            # Update data
            new_data = get_plot_data(t, s, log=log_flag)
            im.set_data(new_data)
            im.set_clim(vmin=current_vmin, vmax=current_vmax)
            title.set_text(f"Timestep: {t}, Slice: {s}, Log: {log_flag}")
            fig.canvas.draw_idle()

        # Log scale toggle — resets vmin/vmax sliders
        def toggle_log(label):
            log_flag = check.get_status()[0]
            if log_flag:
                slider_vmin.valmin = vmin_log
                slider_vmin.valmax = vmax_log
                slider_vmax.valmin = vmin_log
                slider_vmax.valmax = vmax_log
                slider_vmin.set_val(vmin_log)
                slider_vmax.set_val(vmax_log)
            else:
                slider_vmin.valmin = vmin_raw
                slider_vmin.valmax = vmax_raw
                slider_vmax.valmin = vmin_raw
                slider_vmax.valmax = vmax_raw
                slider_vmin.set_val(vmin_raw)
                slider_vmax.set_val(vmax_raw)
            slider_vmin.ax.set_xlim(slider_vmin.valmin, slider_vmin.valmax)
            slider_vmax.ax.set_xlim(slider_vmax.valmin, slider_vmax.valmax)
            update(None)

        # runs update functions on slider movement
        slider_timestep.on_changed(update)
        slider_slice.on_changed(update)
        slider_vmin.on_changed(update)
        slider_vmax.on_changed(update)
        check.on_clicked(toggle_log)

        plt.show()

if __name__ == "__main__":
    data_directory = "/Users/brynlloyd/Developer/Coding/Python/dft/gpaw/my_own_stopping/data/322_supercell"
    analysis = StoppingPowerAnalysis(data_directory)
    analysis.visualise_ehrenfest("40 keV")
