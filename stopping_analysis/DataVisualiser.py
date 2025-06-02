import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, CheckButtons


class DataVisualiser:
    def __init__(self):
        pass

    ########################################
    # FOR PLOTTING STOPPING POWER ANALYSIS #
    ########################################

    def plot_single_fit(self, projectile_positions, projectile_kinetic_energies):
        
        #################
        # CREATE FIGURE #
        #################

        fig,axs = plt.subplots(1,2, figsize=(15,5))
        _ = [ax.set_xlabel(r"Projectile Position [$\AA$]") for ax in axs]
        axs[0].set_ylabel("Kinetic Energy [keV]")
        axs[1].set_ylabel(r"Projectile Position [$\AA$]")

        #############
        # PLOT DATA #
        #############

        # plot kinetic energy on left plot
        axs[0].plot(projectile_positions, projectile_kinetic_energies, "x")
        # plot gradient of kinetic energy on right plot
        axs[1].plot(projectile_positions[:-1], np.diff(projectile_kinetic_energies))

        #############
        # PLOT FITS #
        #############




    def plot_all_fits(self, atoms_dict, fits=None):
        """
        for every projectile energy, creates 2 subplots
        1. projectile kinetic energy against distance travelled
        2. gradient of kinetic energy against distance travelled
        """

        ######################
        # SORT OUT VARIABLES #
        ######################


        projectile_positions = {energy: [atoms.get_positions()
                                         for atoms in atoms_list]
                                for energy, atoms_list in atoms_dict.items()}

        projectile_kinetic_energies = {energy: [atoms.get_kinetic_energy()
                                                for atoms in atoms_list]
                                       for energy, atoms_list in atoms_dict.items()}


        # also need to sort out fits



        #################
        # CREATE FIGURE #
        #################

        n_subplots = len(projectile_kinetic_energies.keys())
        fig,axs = plt.subplots(n_subplots, 2, figsize=(15, 5*n_subplots), sharex="col")
        if n_subplots == 1:
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
            energy, kinetic_energies = list(projectile_kinetic_energies.items())[i]
            kinetic_energies = np.array(kinetic_energies) * 1e-3  # convert from eV to keV
            _, positions = list(projectile_positions.items())[i]
            positions = np.array(positions)[:,-1,0]   # for every timestep, want projectile x_pos

            # plot kinetic energy on left plot
            axs[i, 0].plot(positions, kinetic_energies, "x")

            # plot gradient of kinetic energy on right plot
            axs[i, 1].plot(positions[:-1], np.diff(kinetic_energies))


    ######################################
    # FOR VISUALISING SIMULATION RESULTS #
    ######################################

    def visualise_electron_density(self):
        pass
