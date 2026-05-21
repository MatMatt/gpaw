import numpy as np

from gpaw.core.arrays import XArray
from gpaw.typing import ArrayND
from gpaw.core.plane_waves import PWDesc
from gpaw.core.uniform_grid import UGDesc


class BaseMetric:
    def __init__(self, g_ss: ArrayND | None = None):
        self.g_ss = g_ss

    def initialize(self,
                   ncomponents: int,
                   grid: UGDesc,
                   xp=np) -> None:
        self.ncomponents = ncomponents
        self.xp = xp
        if self.g_ss is None or len(self.g_ss) != ncomponents:
            self.g_ss = xp.eye(ncomponents)
        self.grid = grid
        self.pw = PWDesc(ecut=0.99 * grid.ekin_max(),
                         cell=grid.cell,
                         comm=grid.comm)

    def __call__(self, n_sX: XArray,
                 N_sX: XArray | None = None) -> XArray:
        if N_sX is None:
            N_sX = n_sX.new()
        N_sX.matrix.data[:] = self.g_ss @ n_sX.matrix.data
        return N_sX


class FFTMetric(BaseMetric):
    def __init__(self,
                 weight: float = 100,
                 sigma: float = 1.0,
                 g_ss: ArrayND | None = None):
        self.weight = weight
        self.sigma = sigma
        super().__init__(g_ss)

    def __call__(self, n_sX: XArray,
                 N_sX: XArray | None = None) -> XArray:
        if N_sX is None:
            N_sX = n_sX.new()
        # We will now assume that the density is represented
        # in real space. This may change in the future.

        ekin_G = self.pw.ekin_G
        w_G = self.weight * (self.sigma + ekin_G) \
            / (self.sigma + self.weight * ekin_G)

        # TODO: Fix non-pbc directions
        n_sG = n_sX.fft(pw=self.pw)
        n_sG.data *= w_G
        n_sG.ifft(grid=self.grid, out=N_sX)
        N_sX.matrix.data[:] = self.g_ss @ N_sX.matrix.data

        return N_sX
