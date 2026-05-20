from ase.parallel import paropen
from ase.io import read, write
# from ase.io.trajectory import Trajectory
from ase.dft.bee import BEEFEnsemble
from gpaw import GPAW, FermiDirac, Mixer, PW
fepo4 = read('fepo4.traj')

params_GPAW = {}
params_GPAW['mode']        = PW(500)                     #The used plane wave energy cutoff
params_GPAW['nbands']      = -40                           #The number on empty bands had the system been spin-paired
params_GPAW['kpts']        = {'size':  (2,4,5),            #The k-point mesh
                              'gamma': True}
params_GPAW['spinpol']     = True                          #Performing spin polarized calculations
params_GPAW['xc']          = 'BEEF-vdW'                    #The used exchange-correlation functional
params_GPAW['occupations'] = FermiDirac(width = 0.1,      #The smearing
                                        fixmagmom = True)  #Total magnetic moment fixed to the initial value
params_GPAW['convergence'] = {'eigenstates': 1.0e-4,       #eV^2 / electron
                              'energy':      2.0e-4,       #eV / electron
                              'density':     1.0e-3,}
params_GPAW['mixer']       = Mixer(0.1, 5, weight=100.0)   #The mixer used during SCF optimization
params_GPAW['setups']      = {'Fe': ':d,4.3'}              #U=4.3 applied to d orbitals

calc = GPAW(**params_GPAW)
fepo4.calc = calc
epot_fepo4_cell=fepo4.get_potential_energy()
print('E_Pot=', epot_fepo4_cell)

write('fepo4_out.traj', fepo4)

ens = BEEFEnsemble(calc)
dE = ens.get_ensemble_energies(2000)

with paropen('ensemble_fepo4.dat', 'a') as result:
    for e in dE:
        print(e, file=result)
