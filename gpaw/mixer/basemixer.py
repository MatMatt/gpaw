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
        self.density_history = DensityHistory(
            nmaxold=1,
            nspins=desc.nspins,
            desc=desc,
            atoms_layout=None,
        )

    def mix(self, density: Density) -> float:

        ntc_sR = add_compensation_charge(density)
        Rc_sR = self.calculate_residual(ntc_sR)
        self.prev_density = density
        return self.calculate_charge_sloshing(Rc_sR)

    def calculate_charge_sloshing(self, res_sX: XArray) -> float:
        return

    def calculate_residual(self, density: Density) -> Density:
        if self.prev_density is None:
            return density
        return density - self.prev_density

    def add_compensation_charge(self, density: Density) -> XArray:
        return density

    def __str__(self):
        return "No-mixing mixer"

class DensityHistory:
    def __init__(self,
                 nmaxold: int,
                 nspins: int,
                 desc: DomainType,
                 atoms_layout: AtomsArrayLayout | None):
        self.nmaxold = nmaxold
        self.nspins = nspins
        self.desc = desc
        self.atoms_layout = atoms_layout
        self._n_hsX = desc.zeroes((nmaxold, nspins))
        self._D_hasii = atoms_layout.zeroes((nmaxold, nspins)) \
            if atoms_layout is not None else None
        self.current_indicies = []

    def add_density(self, density: Density, copy=True) -> None:
        if copy:
            self.add(density.n_sR.copy(), density.D_asii.copy())
        else:
            self.add(density.n_sR, density.D_asii)

    def add(self, n_sX: XArray,
            D_asii: AtomsArray | None = None) -> None:
        if self.nold == self.nmaxold:
            self.delete_oldest()

        self._n_hsX.data[self.next_index] = n_sX
        if self._D_hasii is not None:
            self._D_hasii.data[self.next_index] = D_asii

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
