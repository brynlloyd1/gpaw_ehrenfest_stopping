from gpaw.tddft import TDDFT
from gpaw.tddft.ehrenfest import EhrenfestVelocityVerlet

from ase.units import Bohr, AUT, _amu, _me

import numpy as np

####################
# create filenames #
####################

name = "Al_stopping"
kinetic_energy = 40e3
filename_body = f"{name}_{str(kinetic_energy*1e-3):0f}k"


######################
# create calculation #
######################


# will need to add parallel argument to get parallelisation to work
td_calc = TDDFT(filename=name + "_ground_state.gpw",   # file containing the ground state from which to propagate
                propagator="EFSICN",                   # time propagator for the Kohn-Sham equations - need this one for ehrenfest
                solver="BiCGStab",                     # solver for the propagator
                text="output.txt")


#########################
# initialise projectile #
#########################

projectile_index = len(td_calc.atoms.positions) - 1
initial_velocities = np.zeros_like(td_calc.atoms.velocities)
projectile_mass = td_calc.atoms.get_masses()[projectile_index] * (_amu / _me)
initial_velocities[projectile_index, 0] = np.sqrt((2 * kinetic_energy) / projectile_mass) * Bohr / AUT

td_calc.atoms.set_velocities(initial_velocities)

##########################
# run ehrenfest dynamics #
##########################

timestep = np.sqrt(1e3 / kinetic_energy) * 16          # this is what is used by the example


niters = 100
save_every = 10

evv = EhrenfestVelocityVerlet(td_calc)
for i in range(niters):
    if i != 0 and i % save_every == 0:
        td_calc.write(filename_body + "_step" + i + ".gpw")

    evv.propagate(timestep)
