ecut = 400.0
vacuum = 4.0
from ase import Atoms
from gpaw import GPAW, PW, Davidson
nitrogen = Atoms('N', magmoms=[3])
nitrogen.center(vacuum=4.0)
nitrogen.calc = GPAW(txt='N.txt',
                     mode=PW(ecut),
                     eigensolver=Davidson(niter=2))
en = nitrogen.get_potential_energy()

# %%
print(en, 'eV')

nlayers = 2
a = 2.72
c = 1.58 * a
from ase.build import hcp0001
slab = hcp0001('Ru', a=a, c=c, size=(1, 1, nlayers), vacuum=vacuum)

# %%
from ase.visualize import view
view(slab, repeat=(3, 3, 2))

# %%
nkpts = 7
slab.calc = GPAW(txt='Ru.txt',
                 mode=PW(ecut),
                 eigensolver=Davidson(niter=2),
                 kpts={'size': (nkpts, nkpts, 1), 'gamma': True},
                 xc='PBE')
eru = slab.get_potential_energy()
