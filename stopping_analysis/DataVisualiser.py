import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, CheckButtons


class DataVisualiser:
    def __init__(self):
        pass

    ########################################
    # FOR PLOTTING STOPPING POWER ANALYSIS #
    ########################################

    def plot_all_fits(self, atoms_dict, fits, covs):
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


            if fits:
                fit = fits[energy]
                pfit = np.poly1d(fit)
                cov = covs[energy]
                stopping_power = -fit[0]
                stopping_power_uncertainty = np.sqrt(cov[0][0])
                label = rf"$S_e$ = {(1e3*stopping_power):.1f} $\pm$ {(1e3*stopping_power_uncertainty)} [eV/$\AA$]"
                axs[i, 0].plot(positions, pfit(positions), color="red", label=label)
                axs[i, 1].plot(positions[:-1], np.diff(pfit(positions)), color="red")

                axs[i,0].legend()



    ######################################
    # FOR VISUALISING SIMULATION RESULTS #
    ######################################

    def visualise_electron_density(self, electron_density_list):
        """
        creates animation using matplotlib sliders to display 2d slice of electron density
        has sliders to control slice position, timestep, vmin&vmax
        can also toggle logarithmic plotting
        """

        shape = np.shape(electron_density_list)
        timesteps = shape[0]

        # Initial values
        init_timestep = 0
        init_slice = shape[1] // 2
        use_log = False

        #########################
        # CALCULATE VMIN/VMAX'S #
        #########################

        all_data = np.stack(electron_density_list)
        vmin_raw = all_data.min()
        vmax_raw = all_data.max()
        log_data = np.log10(np.clip(all_data, 1e-10, None))
        vmin_log = log_data.min()
        vmax_log = log_data.max()

        ##############
        # SETUP PLOT #
        ##############

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

    def visualise_electron_density_change(self, electron_density_list):
        """
        creates animation using matplotlib sliders to display 2d slice of electron density
        has sliders to control slice position, timestep, vmin&vmax
        can also toggle logarithmic plotting
        """

        electron_density_array = np.array(electron_density_list)
        electron_density_change = electron_density_array = electron_density_array[0]
        print("here...")
        shape = np.shape(electron_density_list)
        print(shape)
        timesteps = shape[0]

        # Initial values
        init_timestep = 0
        init_slice = shape[1] // 2
        use_log = False

        #########################
        # CALCULATE VMIN/VMAX'S #
        #########################

        vmin_raw = np.min(electron_density_change)
        vmax_raw = np.max(electron_density_change)
        log_data = np.log10(np.clip(electron_density_change, 1e-10, None))
        vmin_log = np.log10(vmin_raw)
        vmax_log = np.log10(vmax_raw)

        ##############
        # SETUP PLOT #
        ##############

        fig, ax = plt.subplots(figsize=(8, 6))
        plt.subplots_adjust(left=0.25, bottom=0.45)  # Extra room for sliders

        def get_plot_data(t, slice, log=False):
            print(np.shape(electron_density_change))
            data = electron_density_change[t, :, slice]
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


