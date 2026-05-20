import numpy as np

from gpaw.new.density import Density
from gpaw.core.arrays import DomainType, XArray
from gpaw.core.atom_arrays import AtomArraysLayout, AtomArrays
from gpaw.new.density_history import DensityHistory


class PulayMixer(BaseMixer):
    def _initialize_history_(self):
        self.density_history = DensityHistory(
            nmaxold=16,
            nspins=self.nspins,
            desc=self.desc,
            atoms_layout=self.atoms_layout,
            xp=self.xp
        )

    def mix(self, density: Density) -> float:
        if self.density_history.nold == 0:
            self.density_history.add_density(density)
            return 0.0
