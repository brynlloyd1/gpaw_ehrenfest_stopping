from gpaw.tddft import TDDFT
from gpaw.tddft.ehrenfest import EhrenfestVelocityVerlet
from ase.units import Bohr, AUT, _amu, _me, Hartree
import numpy as np

####################
# create filenames #
####################

name = "Al_stopping"
kinetic_energy = 10e3
filename_body = f"{name}_{str(round(kinetic_energy*1e-3))}k"


######################
# create calculation #
######################


# will need to add parallel argument to get parallelisation to work
td_calc = TDDFT(filename=name + "_ground_state" + ".gpw",                # file containing the ground state
                propagator="EFSICN",                                     # time propagator for the Kohn-Sham equations - need this one for ehrenfest
                solver="BiCGStab")                                       # solver for the propagator


#########################
# initialise projectile #
#########################

projectile_index = len(td_calc.atoms.positions) - 1
initial_velocities = np.zeros_like(td_calc.atoms.get_velocities())
projectile_mass = td_calc.atoms.get_masses()[projectile_index] * (_amu / _me)
kinetic_energy *= 1 / Hartree

# UNCOMMENT FOR HYPERCHANNELLING TRAJECTORY
#initial_velocities[projectile_index, 0] = np.sqrt((2 * kinetic_energy) / projectile_mass) * Bohr / AUT

# UNCOMMENT FOR PRESAMPLED TRAJECTORY
speed = np.sqrt((2 * kinetic_energy) / projectile_mass) * Bohr / AUT
initial_direction = np.array([-0.39661179, -0.60274438, 0.69238594])
initial_velocities[projectile_index] = (speed * np.array(initial_direction)).tolist()

# check that the projectile kinetic energy is correct
projectile_KE_check = (np.linalg.norm(initial_velocities[projectile_index]) * AUT / Bohr)**2 * projectile_mass/2
print(f"""calculated projectile KE: {projectile_KE_check}
       expected: {kinetic_energy} """)

td_calc.atoms.set_velocities(initial_velocities)


##########################
# run ehrenfest dynamics #
##########################

timestep = np.sqrt(1e3 / kinetic_energy) * 16.0         # this is what is used by the example

niters = 100
save_every = 1

evv = EhrenfestVelocityVerlet(td_calc)
for i in range(niters):

    print(f"timestep: {i+1} / {niters}")
    if i != 0 and i % save_every == 0:
        td_calc.write(filename_body + "_step" + str(i) + ".gpw")

    evv.propagate(timestep)
