import pytest
from ase import Atoms
from gpaw.dft import GPAW, PW
from gpaw.mpi import world


def par(size):
    kb = []
    for k in range(1, 3):
        if size % k != 0:
            continue
        for b in range(1, size // k + 1):
            if (size // k) % b == 0:
                kb.append((k, b))
    return kb


@pytest.mark.parametrize(
    'kb',
    [pytest.param((k, b), id=f'k{k}b{b}') for k, b in par(world.size)])
def test_all(kb):
    k, b = kb
    L = 2.6
    a = Atoms('Li2',
              [[0, 0, 0], [0.9, 0.9, 0]],
              cell=[L, L, 1.4],
              pbc=1)
    a.center()

    kwargs = dict(
        mode=PW(400),
        convergence={'density': 1e-7,
                     'forces': 1e-6},
        mixer={'beta': 0.25},
        xc='HSE06')
    a.calc = GPAW(
        kpts={'size': (1, 1, 4), 'gamma': True},
        parallel={'kpt': k, 'band': b},
        **kwargs)
    a.get_potential_energy()
    from gpaw.new.pw.bs import fixed
    import numpy as np
    print(a.calc.eigenvalues())
    kpts = np.zeros((21, 3))
    kpts[:, 2] = np.linspace(-0.25, 0.5, 21)
    i = fixed(a.calc.dft, kpts)
    import matplotlib.pyplot as plt
    plt.plot(kpts[:, 2], [wfs.eig_n[0] for wfs in i])
    plt.show()


if __name__ == '__main__':
    import sys
    test_all([int(x) for x in sys.argv[1:]])
    # print(par(2))
