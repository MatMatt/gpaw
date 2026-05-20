calc_fm = GPAW('CrI3_fm.gpw', txt=None)      # Ferromagnetic calculation
calc_afm = GPAW('CrI3_afm.gpw', txt=None)    # Anti-ferromagnetic calculation
    
N = 3
E_fm = calc_fm.get_potential_energy() / 2    # Energy per magnetic atom
E_afm = calc_afm.get_potential_energy() / 2  # Energy per magnetic atom
dE = E_afm - E_fm

J = dE / S**2 / 3                            # Student version: J = ???
print(f'J = {J * 1000:1.3f} meV')

from ase.units import kB                     # Boltzmann's constant in eV/K
T_c = N * J * S * (S + 1) / 3 / kB           # Student version: T_c = ???
print(f'T_c(MF) = {T_c:1.1f} K')
