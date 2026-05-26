from ase import Atoms
from ase.mep import NEB
from ase.optimize import BFGS
from ase.parallel import paropen
from gpaw import GPAW, FermiDirac, PW
from ase.constraints import FixAtoms

ccdist = 1.39
a = ccdist * 3**0.5
c = 3.7

initial = Atoms(
    'C6Li',
    positions=[(0, 0, 0),
               (0, ccdist, 0),
               (a, 0, 0),
               (-a, 0, 0),
               (-a / 2, -ccdist / 2, 0),
               (a / 2, -ccdist / 2, 0),
               (0, -ccdist, c / 2)],
    cell=([1.5 * a, -1.5 * ccdist, 0],
          [1.5 * a, 1.5 * ccdist, 0],
          [0, 0, c]),
    pbc=(1, 1, 1))
# snippet-final
final = initial.copy()
a, b = final.cell[:2]
final.positions[6] += a / 3 + b / 3
# snippet-neb1
# These will become the minimum energy path images.
images = [initial]
images += [initial.copy() for i in range(5)]
images += [final]
# snippet-neb2
neb = NEB(images, method='improvedtangent')
neb.interpolate()

