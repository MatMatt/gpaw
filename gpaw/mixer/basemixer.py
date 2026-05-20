import numpy as np

from gpaw.new.density import Density
from gpaw.core.arrays import DomainType, XArray
from gpaw.core.atoms_arrays import AtomsArrayLayout, AtomsArray
from gpaw.typing import NDArray

class BaseMixer:
    name = 'no-mixing'

    def __init__(self,
                 desc: DomainType):
        self.desc = desc
        self.xp = self.desc.xp
        self.density_history = DensityHistory(
            nmaxold=1,
            nspins=desc.nspins,
            desc=desc,
        )

    def mix(self, density: Density) -> float:
        ntc_sR = add_compensation_charge(density)
        if self.density_history.nold == 0:
            self.density_history.add(ntc_sR)
            return np.inf

        Rc_sR = self.calculate_residual(ntc_sR)
        self.density_history.add(ntc_sR)
        return self.calculate_charge_sloshing(Rc_sR)

    def calculate_charge_sloshing(self, res_sX: XArray) -> float:
        slosh_sX = self.desc.empty()
        slosh_sX.data[:] = self.xp.abs(res_sX.data)
        return slosh_sX.integrate()

    def calculate_residual(self, ntc_sX: XArray, res_sX: XArray) -> XArray:
        assert ntc_sX.data is not out_sX.data
        # We do density history first, in case res_sX is the same as
        # self.density_history[-1].
        res_sX.data[:] = -self.density_history[-1].data
        res_sX.data += ntc_sX.data
        return res_sX

    def add_compensation_charge(self, density: Density) -> XArray:
        return density.nt_sR

    def __str__(self):
        return "No-mixing mixer"

class DensityHistory:
    def __init__(self,
                 nmaxold: int,
                 nspins: int,
                 desc: DomainType,
                 atoms_layout: AtomsArrayLayout | None = None):
        self.nmaxold = nmaxold
        self.nspins = nspins
        self.desc = desc
        self.atoms_layout = atoms_layout
        self._n_hsX = desc.zeroes((nmaxold, nspins))
        self._D_hasii = atoms_layout.zeroes((nmaxold, nspins)) \
            if atoms_layout is not None else None
        self.current_indicies = []

    def add_density(self, density: Density) -> None:
        self.add(density.n_sR, density.D_asii)

    def add(self, n_sX: XArray,
            D_asii: AtomsArray | None = None) -> None:
        if self.nold == self.nmaxold:
            self.delete_oldest()

        self._n_hsX.data[self.next_index] = n_sX.data
        if self._D_hasii is not None:
            self._D_hasii.data[self.next_index] = D_asii.data

        self.current_indicies.append(self.next_index)

    def delete_oldest(self) -> None:
        self.current_indicies.pop(0)

    @property
    def nold(self) -> int:
        return len(self.current_indicies)

    @property
    def next_index(self) -> int:
        for i in range(self.nmaxold):
            if i not in self.current_indicies:
                return i
        raise ValueError('No available index')

    def __getitem__(self, index: int) -> tuple[XArray, AtomsArray] | XArray:
        index = self.current_indicies[index]

        if self._D_hasii is not None:
            return self._n_hsX[index], self._D_hasii[index]
        else:
            return self._n_hsX[index]
