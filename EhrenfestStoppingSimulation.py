from gpaw.tddft import TDDFT
from gpaw.tddft.ehrenfest import EhrenfestVelocityVerlet
from gpaw import setup_paths
from ase.units import Bohr, AUT, _amu, _me, Hartree
from ase.parallel import parprint
import numpy as np

import time
from datetime import timedelta


class EhrenfestStoppingSimulation:
    def __init__(self, name, trajectory, kinetic_energies):
        """
        Parameters:
        name (str)
        trajectory (Trajectory)
        kinetic_energies (np.ndarray[int])
        """

        # so can use "H.H+.LDA"
        setup_paths.insert(0, "./custom_setups")

        self.name = name
        self.ground_state_name = name + "_ground_state.gpw"
        self.trajectory = trajectory
        self.kinetic_energies = kinetic_energies / Hartree

        self.td_calc = None

    def initialise_calculation(self, kinetic_energy):
        """
        creates instance of gpaw.tddft.TDDFT

        Parameters:
        kinetic_energy (int): initial kinetic energy of projectile in eV
        """

        # grid spacing is only defined in the ground state calculation
        self.td_calc = TDDFT(filename = self.ground_state_name,                        # file containing the ground state
                            propagator = "EFSICN",                                     # time propagator for the Kohn-Sham equations - need this one for ehrenfest
                            solver = "BiCGStab",                                       # solver for the propagator
                            txt = self.name + '_log.txt')

        projectile_index = len(self.td_calc.atoms.positions) - 1
        initial_velocities = np.zeros_like(self.td_calc.atoms.get_velocities())
        projectile_mass = self.td_calc.atoms.get_masses()[projectile_index] * (_amu / _me)
        speed = np.sqrt((2 * kinetic_energy) / projectile_mass) * Bohr / AUT
        initial_velocities[projectile_index] = (1.0 * speed * np.array(self.trajectory.direction)).tolist()

        self.td_calc.atoms.set_velocities(initial_velocities)


    def run_ehrenfest(self, kinetic_energy):
        """
        creates instance of gpaw.tddft.EhrenfestVelocityVerlet
        and runs .propagate for 100 timesteps, saving .gpw file for every timestep

        Parameters:
        kinetic_energy (int): initial kinetic energy of projectile in eV
        """

        if not self.td_calc:
            raise ValueError("EhrenfestStoppingSimulation.td_calc not initialised.")

        # timestemp was previously * 16.0
        timestep = np.sqrt(1e3 / kinetic_energy) * 4.0
        parprint(timestep)
        niters = 400
        save_every = 1

        evv = EhrenfestVelocityVerlet(self.td_calc)
        for n in range(niters):

            if n != 0 and n % save_every == 0:
                self.td_calc.write(f"{self.name}_{str(round(kinetic_energy * Hartree * 1e-3))}k_step{str(n)}.gpw")
            start_time = time.time()
            evv.propagate(timestep)

            end_time = time.time()
            formatted_time = str(timedelta(seconds = int(end_time - start_time)))

            parprint(f"    completed timestep: {n+1}/{niters} in {formatted_time}")

    def run_simulation(self):
        for i, kinetic_energy in enumerate(self.kinetic_energies):
            parprint(f"Running Simulation {i+1}/{len(self.kinetic_energies)}: {kinetic_energy * Hartree * 1e-3}keV")
            self.initialise_calculation(kinetic_energy)
            self.run_ehrenfest(kinetic_energy)

