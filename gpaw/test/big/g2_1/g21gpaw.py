from ase.build import molecule
from ase.optimize.bfgs import BFGS
from gpaw import GPAW, PW


def relax(name):
    atoms = molecule(name)
    atoms.cell = [12, 13, 14]
    atoms.center()
    atoms.calc = GPAW(mode=PW(800),
                      xc='PBE',
                      txt=name + '.txt')
    atoms.get_forces()
    opt = BFGS(atoms, logfile=name + '.log')
    opt.run(0.01)
