from gpaw import GPAW, FermiDirac

from ase.lattice.cubic import FaceCenteredCubic
from ase import Atoms



name = "Al_stopping_ground_state"

##########################
# create al + projectile #
##########################

AL_LATTICE_CONSTANT = 4.05
al_supercell = FaceCenteredCubic("Al",
                                 size=(3, 2, 2),
                                 latticeconstant=AL_LATTICE_CONSTANT)
# al_supercell.center(vacuum=15.0, axis=0)
al_supercell.center()

# HYPERCHANNELLING STARTING POSITION
# projectile_starting_position = [0.0, AL_LATTICE_CONSTANT, AL_LATTICE_CONSTANT]

# PRESAMPLED TRAJECTORY STARTING POSITION
projectile_starting_position = [3.83534014, 2.7002751, 0.44542105]
projectile = Atoms("H",
                   cell=al_supercell.cell,
                   positions=[projectile_starting_position])

system = al_supercell  + projectile

#####################
# find ground state #
#####################

calc = GPAW(mode = "fd",
            occupations = FermiDirac(0.02),
            symmetry = {"point_group": False},
            txt = "output.txt")


system.calc = calc
system.get_potential_energy()
system.calc.write(name + '.gpw', mode='all')
