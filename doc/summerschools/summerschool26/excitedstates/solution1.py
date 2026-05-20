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

