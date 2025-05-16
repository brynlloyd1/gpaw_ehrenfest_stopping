from gpaw import GPAW, RMMDIIS
from gpaw.external import ConstantPotential
from gpaw.utilities import h2gpts

from ase.lattice.cubic import FaceCenteredCubic
from ase import Atoms
from ase.units import Bohr

import numpy as np


name = "Al_stopping"

##########################
# create al + projectile #
##########################

AL_LATTICE_CONSTANT = 4.05
al_supercell = FaceCenteredCubic("Al",
                                 size=(1, 2, 2),
                                 latticeconstant=AL_LATTICE_CONSTANT)
al_supercell.center(vacuum=15.0, axis=0)
al_supercell.center()

projectile_starting_position = [5.0, AL_LATTICE_CONSTANT, AL_LATTICE_CONSTANT]
projectile = Atoms("H",
                   cell=al_supercell.cell,
                   positions=[projectile_starting_position])

system = al_supercell + projectile
# PREVIOUSLY HAD PBC'S ON (IN ALL THREE DIRECTIONS?? WHICH DOESNT MAKE SENSE BC I ADD VACUUM TO THE END...)
# system.set_pbc(True)


#####################
# find ground state #
#####################

# supposedly two-stage convergence is required for stopping of a charged particle
conv_fast = {"energy": 1.0, "density": 1.0, "eigenstates": 1.0}
conv_par = {"energy": 1e-3, "density": 1e-3, "eigenstates": 1e-7}


constant_potential = ConstantPotential(1.0)


# first ground_state calculation with fast convergence
calc = GPAW(mode='fd',
            gpts=h2gpts(0.2, al_supercell.cell, idiv=8),
            nbands=110,
            xc='LDA',
            charge=1,
            txt=name + 'ground_state.txt',
            convergence=conv_fast,
            external=constant_potential,
            symmetry={'point_group': False})


system.calc = calc
system.get_potential_energy()


# second ground state calculation with par?? convergence

# External potential used to prevent charge tranfer from graphene to ion.
A = 1.0
X0 = system.positions[-1] / Bohr
rcut = 3.0 / Bohr
vext = calc.hamiltonian.vext
gd = calc.hamiltonian.finegd
n_c = gd.n_c
h_c = gd.get_grid_spacings()
b_c = gd.beg_c
vext.vext_g.flags.writeable = True
vext.vext_g[:] = 0.0
for i in range(n_c[0]):
    for j in range(n_c[1]):
        for k in range(n_c[2]):
            x = h_c[0] * (b_c[0] + i)
            y = h_c[1] * (b_c[1] + j)
            z = h_c[2] * (b_c[2] + k)
            X = np.array([x, y, z])
            dist = np.linalg.norm(X - X0)
            if dist < rcut:
                vext.vext_g[i, j, k] += A * np.exp(-dist**2)

system.calc = calc.new(convergence=conv_par, eigensolver=RMMDIIS(5), external=vext, txt=name + '_vext_gs.txt')

system.get_potential_energy()
system.calc.write(name + '.gpw', mode='all')
