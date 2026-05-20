# basic imports
import numpy as np
from ase import Atoms # Import atoms
from ase.build import bulk

Si = bulk('Si', 'diamond', a=5.4) # student:
Ge = bulk('Ge','diamond', a=5.7) # student:
C = bulk('C', 'diamond', a=3.6) # student:
CdTe = bulk('CdTe', 'zincblende', 6.5) # student:
GaAs = bulk('GaAs','zincblende',a = 5.6) # student:
BN = Atoms('BN', pbc=[True, True, False], positions=[[0., 0., 7.],[0., 2.5/np.sqrt(3), 7.]], cell=[[2.50, 0., 0.],[-2.5/2., 2.5*np.sqrt(3.)/2., 0.],[0., 0., 14.]]) # student:

atoms = Si # student: atoms = ???
label = 'Si' # student: label = '???'

#view(atoms) # check your initial structure # student: view(atoms)
from gpaw import GPAW, PW, FermiDirac

# We create the calculator object
calc = GPAW(xc = 'PBE',
            txt = label + '_relax.txt',
            occupations = FermiDirac(0.01), # smearing the occupation # student: occupations = ???
            mode = PW(600), # plane wave basis # student: mode = ???
            kpts = {'size': (6,6,6)} # 661 for 2D # student: kpts = ???
)

atoms.calc = calc
from ase.filters import UnitCellFilter
from ase.optimize import BFGS

filt = UnitCellFilter(atoms) # student: filt = ???
op = BFGS(filt) # student: op = ???

# Run the optimization. This will take some time, do not get nervous.
# Only if it takes longer than 4-5 minutes or if it does not print anything
# contact us :)
op.run(fmax=0.05)

# save the results in a file
calc.write(label + '_gs.gpw')


