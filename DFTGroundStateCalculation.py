from gpaw import GPAW, FermiDirac
from gpaw.external import ConstantPotential

from ase.lattice.cubic import FaceCenteredCubic
from ase import Atoms
from ase.units import Bohr

import numpy as np


class DFTGroundStateCalculation:
    def __init__(self, supercell_size, trajectory):
        self.name = "Al_stopping_ground_state"
        self.supercell_size = supercell_size
        self.trajectory = trajectory

        self.system = None

    def create_supercell_projectile_system(self):

        AL_LATTICE_CONSTANT = 4.05
        al_supercell = FaceCenteredCubic("Al",
                                         size = self.supercell_size,
                                         latticeconstant = AL_LATTICE_CONSTANT)
        al_supercell.center()
        projectile = Atoms("H",
                           cell = al_supercell.cell,
                           positions = [self.trajectory.position])

        self.system = al_supercell  + projectile


    def set_up_calculation(self):
        """
        sets up the calculation of the electronic ground state of the system.
        Want to set an external potential so that the projectile can be added in as a charged particle.
        To do this, need to run an initial calculation (require only very poor convergence to make it as quick as possible),
        giving an external potential that can be modified.
        Actual ground state calculation is then set up.
        """

        if not self.system:
            raise ValueError("system not yet initialised")

        constant_potential = ConstantPotential(0.0)
        conv_fast = {"energy" : 1.0,
                    "density" : 1.0,
                    "eigenstates" : 1.0, }

        calc_initial = GPAW(mode = "fd",
                    occupations = FermiDirac(0.02),
                    symmetry = {"point_group": False},
                    external = constant_potential,
                    convergence = conv_fast,
                    txt = "output.txt")

        self.system.calc = calc_initial
        self.system.get_potential_energy()


        vext = self.system.calc.hamiltonian.vext
        gd = self.system.calc.hamiltonian.finegd
        n_c = gd.n_c
        h_c = gd.get_grid_spacings()
        b_c = gd.beg_c
        vext.vext_g.flags.writeable = True

        A = 2.0
        X0 = self.system.positions[len(self.system.positions) - 1] / Bohr   # projectile position in units of Bohr radii
        cutoff_radius = 3.0 / Bohr
        for i in range(n_c[0]):
            for j in range(n_c[1]):
                for k in range(n_c[2]):
                    x = h_c[0] * (b_c[0] + i)
                    y = h_c[1] * (b_c[1] + j)
                    z = h_c[2] * (b_c[2] + k)

                    X = np.array([x, y, z])
                    dist_to_proton = np.linalg.norm(X - X0)

                    if dist_to_proton < cutoff_radius:
                        vext.vext_g[i, j, k] += A * np.exp(- dist_to_proton**2)


        # set-up of actual dft calculation
        calc = GPAW(mode = "fd",
                    occupations = FermiDirac(0.02),
                    symmetry = {"point_group": False},
                    external = vext,
                    charge = 1,
                    txt = "output.txt")

        self.system.calc = calc

    def run_calculation(self):
        """
        runs the actual calculation of the electronic ground state of the system
        """

        if not self.system:
            raise ValueError("system not yet initialised")

        self.system.get_potential_energy()
        self.system.calc.write(self.name + '.gpw', mode='all')

    def run(self):
        self.create_supercell_projectile_system()
        self.set_up_calculation()
        self.run_calculation()
