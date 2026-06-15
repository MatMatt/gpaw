# Script for reading previous structure and converging
from gpaw import GPAW, PW
from ase.constraints import FixAtoms
from ase.optimize import QuasiNewton
import numpy as np
from ase.io import read, write

for name in ['N2Ru-top.traj', '2Nads.traj']:
    slab = read(name)

    z = slab.positions[:, 2]
    constraint = FixAtoms(mask=(z < z.min() + 1.0))
    slab.set_constraint(constraint)

    calc = GPAW(xc='PBE',
                mode=PW(500),
                txt='tight-' + name[:-4] + 'txt',
                kpts={'size': (4, 4, 1), 'gamma': True},
                convergence={'eigenstates': 1e-8})
    slab.calc = calc

    dyn = QuasiNewton(slab,
                      trajectory='tight-' + name,
                      maxstep=0.02,
                      logfile='tight-' + name[:-4] + 'log')
    dyn.run(fmax=0.01)

