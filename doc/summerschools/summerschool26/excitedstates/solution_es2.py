from ase.build import bulk
from gpaw import GPAW, FermiDirac
from gpaw import PW

atoms = bulk('Si', 'diamond', a=5.4) # student:

calc = GPAW(mode=PW(300),  # student: # energy cutoff for plane wave basis (in eV)
            kpts={'size': (3, 3, 3), 'gamma': True}, # student: kpts={'size': (?, ?, ?), 'gamma': True},
            xc='LDA',
            occupations=FermiDirac(0.001),
            parallel={'domain': 1},
            txt='Si_groundstate.txt') # student: txt='???_groundstate.txt'

atoms.calc = calc
atoms.get_potential_energy()

calc.diagonalize_full_hamiltonian()  # determine all bands
calc.write('Si_groundstate.gpw', 'all')  # student: '???_groundstate.gpw', 'all' # write out wavefunctions

from gpaw.response.g0w0 import G0W0

gw = G0W0(calc='Si_groundstate.gpw', # student: calc='???_groundstate.gpw',
          nbands=30,  # student: nbands=???,  # number of bands for calculation of self-energy
          bands=(3, 5), # student: bands=(?, ?), # VB and CB
          ecut=20.0,  # student: ecut=???, # plane-wave cutoff for self-energy (20-200)
          integrate_gamma='WS',  # Use supercell Wigner-Seitz truncation for W.
          filename='Si-g0w0') # student: filename='???-g0w0'
gw.calculate()

import pickle

results = pickle.load(open('Si-g0w0_results_GW.pckl', 'rb')) # student: open('???-g0w0_results_GW.pckl', 'rb')
direct_gap = results['qp'][0, 0, -1] - results['qp'][0, 0, -2]

print('Direct bandgap of Si:', direct_gap) # student: 'Direct bandgap of ???:', direct_gap
