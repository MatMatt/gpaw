# web-page: nru2.png
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


# %%
import numpy as np
height = 1.1
nslab = hcp0001('Ru', a=a, c=c, size=(1, 1, nlayers))
# Calculate the coordianates of the N-atoms:
z = slab.positions[:, 2].max() + height
x, y = np.dot([2 / 3, 2 / 3], slab.cell[:2, :2])
nslab.append('N')
nslab.positions[-1] = [x, y, z]
nslab.center(vacuum=vacuum, axis=2)  # 2: z-axis


# %%
height = 1.1
nslab = hcp0001('Ru', a=a, c=c, size=(1, 1, nlayers))
from ase.build import add_adsorbate
add_adsorbate(nslab, 'N', position='hcp', height=height)
nslab.center(vacuum=vacuum, axis=2)

# %%
view(nslab, repeat=(3, 3, 2))

# %%
from ase.io import write
write('nru2.png', nslab.repeat((3, 3, 1)))


# %%
nslab.calc = GPAW(txt='NRu.txt',
                  mode=PW(ecut),
                  eigensolver=Davidson(niter=2),
                  poissonsolver={'dipolelayer': 'xy'},
                  kpts={'size': (nkpts, nkpts, 1), 'gamma': True},
                  xc='PBE')
enru0 = nslab.get_potential_energy()

# %%
print('Unrelaxed adsoption energy:', enru0 - eru - en, 'eV')

f = nslab.get_forces()
print(f[-1])

# %%
# This cell will take a few minutes to finish ...
from ase.constraints import FixAtoms
from ase.optimize import BFGSLineSearch
nslab.constraints = FixAtoms(indices=list(range(nlayers)))
optimizer = BFGSLineSearch(nslab, trajectory='NRu.traj')
optimizer.run(fmax=0.01)
height = nslab.positions[-1, 2] - nslab.positions[:-1, 2].max()
print('Height:', height, 'Ang')

# %%
enru = nslab.get_potential_energy()
print('Relaxed adsorption energy:', enru - eru - en, 'eV')
