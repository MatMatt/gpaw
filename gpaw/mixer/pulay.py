import numpy as np

from gpaw.new.density import Density
from gpaw.mixer.basemixer import DensityHistory, BaseMixer


class PulayMixer(BaseMixer):
    name = 'pulay'
    nmaxold = 16
    beta = 0.08

    def _initialize_history_(self):
        self.nt_hsX = DensityHistory(
            nmaxold=self.nmaxold,
            ncomponents=self.ncomponents,
            desc=self.desc,
            atom_layout=self.atom_layout,
            xp=self.xp
        )
        self.R_hsX = DensityHistory(
            nmaxold=self.nmaxold,
            ncomponents=self.ncomponents,
            desc=self.desc,
            atom_layout=self.atom_layout,
            xp=self.xp
        )
        self.Rc_hsX = DensityHistory(
            nmaxold=self.nmaxold,
            ncomponents=self.ncomponents,
            desc=self.desc,
            xp=self.xp
        )
        self.histories = [self.nt_hsX, self.R_hsX, self.Rc_hsX]

        # DIIS matrix:
        self.H_hh = self.xp.empty(
            (self.nmaxold, self.nmaxold), dtype=float)

    def mix(self, density: Density) -> float:
        # Step 1: Initialize
        nold = self.nt_hsX.nold
        if nold == 0:
            self.nt_hsX.add_density(density)
            return np.inf

        nt_sX = density.nt_sR
        D_asii = density.D_asii

        R_sX, R_asii = self.R_hsX.next()
        self.calculate_residual(nt_sX, self.nt_hsX[-1][0], R_sX)
        self.calculate_paw_residual(D_asii.data, self.nt_hsX[-1][1], R_asii)

        Rc_sX = self.Rc_hsX.next()
        self.add_compensation_charge(R_sX, R_asii, out=Rc_sX)

        # Step 2: Calculate the DIIS matrix
        # TODO: Apply Metric
        Rc_1sX = Rc_sX.new(data=Rc_sX.data[None], dims=(1,) + Rc_sX.dims)
        H_h = self.Rc_hsX.dotprod(Rc_1sX)
        self.H_hh[:nold, nold - 1] = H_h
        self.H_hh[nold - 1, :nold] = H_h

        # TODO: Regularize
        H_hh = self.xp.linalg.inv(self.H_hh[:nold, :nold])
        alpha_i = H_hh.sum(1)
        alpha_i /= alpha_i.sum()

        # Step 3: Mix the densities
        nt_sX.data[:] = 0
        D_asii.data[:] = 0
        for i, alpha in enumerate(alpha_i):
            nt_sX.data[:] += alpha * self.nt_hsX[i][0].data
            D_asii.data[:] += alpha * self.nt_hsX[i][1]
            nt_sX.data[:] += self.beta * alpha * self.R_hsX[i][0].data
            D_asii.data[:] += self.beta * alpha * self.R_hsX[i][1]

        # Step 4: Update the density history
        self.nt_hsX.add_density(density)

        # Step 5: Return the mixing error
        return self.calculate_charge_sloshing(Rc_sX)
