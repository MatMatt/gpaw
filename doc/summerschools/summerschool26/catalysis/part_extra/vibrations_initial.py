# Intended values (hidden from student)
my_kpts = (4, 4, 1)
my_structure = 'tight-N2Ru-top.traj'
my_indices = [8, 9]

# snippet-vibrations-start
from ase.io import read
from ase.vibrations import Vibrations
from gpaw import GPAW, PW

slab = read(my_structure)  # student: slab = read('your_structure')
calc = GPAW(xc='PBE',
            mode=PW(500),
            kpts={'size': my_kpts, 'gamma': True},  # student: kpts=?,
            symmetry={'point_group': False},
            txt='vib.txt')
slab.calc = calc
Uini = slab.get_potential_energy()

vib = Vibrations(slab,
                 name=my_name,
                 indices=my_indices,  # student: indices=[?, ?],
                 nfree=4)
vib.run()
vib.summary(log='vib_summary.log')
for i in range(6):
    vib.write_mode(i)
# snippet-vibrations-end
