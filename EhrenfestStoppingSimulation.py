from gpaw.tddft import TDDFT
from gpaw.tddft.ehrenfest import EhrenfestVelocityVerlet
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

        self.td_calc = TDDFT(filename = self.ground_state_name,                        # file containing the ground state
                            propagator = "EFSICN",                                     # time propagator for the Kohn-Sham equations - need this one for ehrenfest
                            solver = "BiCGStab",                                       # solver for the propagator
                            txt = self.name + '_log.txt')

        projectile_index = len(self.td_calc.atoms.positions) - 1
        initial_velocities = np.zeros_like(self.td_calc.atoms.get_velocities())
        projectile_mass = self.td_calc.atoms.get_masses()[projectile_index] * (_amu / _me)
        speed = np.sqrt((2 * kinetic_energy) / projectile_mass) * Bohr / AUT
        initial_velocities[projectile_index] = (speed * np.array(self.trajectory.direction)).tolist()

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

        timestep = np.sqrt(1e3 / kinetic_energy) * 16.0 
        niters = 100
        save_every = 1

        evv = EhrenfestVelocityVerlet(self.td_calc)
        for i in range(niters):

            if i != 0 and i % save_every == 0:
                self.td_calc.write(f"{self.name}_{str(round(kinetic_energy*1e-3))}k_step{str(i)}.gpw")
            start_time = time.time()
            evv.propagate(timestep)
            end_time = time.time()
            formatted_time = str(timedelta(seconds = int(end_time - start_time)))

            parprint(f"    completed timestep: {i+1}/{niters} in {formatted_time}")

    def run_simulation(self):
        for i, kinetic_energy in enumerate(self.kinetic_energies):
            parprint(f"Running Simulation {i+1}/{len(self.kinetic_energies)}: {kinetic_energy}keV")
            self.initialise_calculation(kinetic_energy)
            self.run_ehrenfest(kinetic_energy)

