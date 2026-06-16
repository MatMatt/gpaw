from myqueue.workflow import run
from ase.data.g2_1 import atom_names, molecule_names
from gpaw.test.big.g2_1.g21gpaw import relax


def workflow():
    deps = []
    for name in molecule_names + atom_names:
        d = run(function=relax)
        deps.apend(d)
    run(script='analyse.py', deps=deps)
