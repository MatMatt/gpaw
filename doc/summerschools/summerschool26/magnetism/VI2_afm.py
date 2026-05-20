from ase.io import read
from gpaw import GPAW

layer_afm = layer.repeat((2, 1, 1))
layer_afm.set_initial_magnetic_moments([m, 0, 0, -m, 0, 0])

calc = GPAW(mode=PW(400),
            xc='LDA',
            kpts=(2, 4, 1))  # student: kpts=???))
layer_afm.calc = calc
layer_afm.get_potential_energy()
calc.write('V2I4_afm.gpw')

layer_fm = layer.repeat((2, 1, 1))
layer_fm.set_initial_magnetic_moments([m, 0, 0, m, 0, 0])
calc = GPAW(mode=PW(400), xc='LDA', kpts=(2, 4, 1))
layer_fm.calc = calc
layer_fm.get_potential_energy()
calc.write('V2I4_fm.gpw')

