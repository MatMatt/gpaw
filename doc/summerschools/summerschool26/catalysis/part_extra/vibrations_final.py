from ase.io import read
from ase.vibrations import Vibrations
from gpaw import GPAW, PW

slab = read('tight-2Nads.traj')
calc = GPAW(xc='PBE',
            mode=PW(500),
            kpts={'size': my_kpts, 'gamma': True},  # student: kpts=?,
            symmetry={'point_group': False},
            txt='vib2.txt')
slab.calc = calc
Ufin = slab.get_potential_energy()

vib2 = Vibrations(slab,
                  name='vib_adsorbate',
                  indices=my_indices,
                  nfree=4)
vib2.run()
vib2.summary(log='vib2_summary.log')
for i in range(6):
    vib2.write_mode(i)

