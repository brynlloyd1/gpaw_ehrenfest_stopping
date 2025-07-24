from gpaw import GPAW, FermiDirac, setup_paths
from gpaw.utilities import h2gpts

from ase.lattice.cubic import FaceCenteredCubic
from ase import Atoms

class DFTGroundStateCalculation:
    def __init__(self, supercell_size, trajectory):
        # to allow gpaw to find the custom potentials
        setup_paths.insert(0, "./custom_setups")

        self.name = "Al_stopping_ground_state"
        self.supercell_size = supercell_size
        self.trajectory = trajectory

        self.system = None

    def create_supercell_projectile_system(self):

        AL_LATTICE_CONSTANT = 4.05
        al_supercell = FaceCenteredCubic("Al",
                                         size = self.supercell_size,
                                         latticeconstant = AL_LATTICE_CONSTANT)
        # al_supercell.center()

        projectile = Atoms("H",
                           cell = al_supercell.cell,
                           positions = [self.trajectory.starting_position])

        self.system = al_supercell  + projectile
        self.system.pbc = True
        self.system.wrap()


    def run_calculation(self):
        if not self.system:
            raise ValueError("System not initialised")

        calc = GPAW(mode = "fd",
                    gpts = h2gpts(0.16, self.system.get_cell(), idiv=8),
                    setups = {"H": "H+",
                              "Al": "Al_2frozen"},
                    occupations = FermiDirac(0.03),
                    symmetry = {"point_group": False},
                    txt = "output.txt")

        self.system.calc = calc
        self.system.get_potential_energy()
