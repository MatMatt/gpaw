# Script for analysing energy contributions from vibrations
from ase.io import read
from ase.thermochemistry import HarmonicThermo
from matplotlib import pyplot as plt

energies = vib.get_energies()  # initial vib. energy levels
energiesfin = vib2.get_energies()  # final vib. energy levels

# The actual analysis
ads = HarmonicThermo(energies)
adsfin = HarmonicThermo(energiesfin)
print('DFT reaction energy: ', Ufin - Uini)
Sini = []
Sfin = []
U = []
U2 = []
temp = [10, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500]
for T in temp:
    Sini.append(-T * ads.get_entropy(T, verbose=False))
    U.append(ads.get_internal_energy(T, verbose=False))
    Sfin.append(-T * adsfin.get_entropy(T, verbose=False))
    U2.append(adsfin.get_internal_energy(T, verbose=False))
    Efree = ads.get_helmholtz_energy(T, verbose=False)
    Efreefin = adsfin.get_helmholtz_energy(T, verbose=False)
    print(f'Reaction free energy at {T} K: ',
          Ufin + Efreefin - (Uini + Efree))
plt.plot(temp, Sini, 'r', label='-T*S ini')
plt.plot(temp, U, 'b', label='U ini')
plt.plot(temp, Sfin, 'r--', label='-T*S fin')
plt.plot(temp, U2, 'b--', label='U fin')
plt.plot(temp, np.array(U) + np.array(Sini), 'm', label='Efree ini')
plt.plot(temp, np.array(U2) + np.array(Sfin), 'm--', label='Efree fin')
plt.legend()
plt.savefig('vib.png')
