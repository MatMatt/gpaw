import numpy as np

from gpaw.new.density import Density
from gpaw.mixer.base import DensityHistory, BaseMixer


class MSR1Mixer(BaseMixer):
    def __init__(self,
                 nmaxold: int = 16,
                 beta: float = 0.08,
                 reg: float = 1e-4,
                 gb_scale: float = 1,
                 max_A: float = 0.45):
        self.nmaxold = nmaxold
        self.beta = beta
        self.reg = reg
        self.gb_scale = gb_scale
        self.A_lims = [0.02, max_A]
        self.B_lims = [0.6, 1.1]
        self.rate_ratio = [0.7, 1.3]
        self.A = 0
        sefl.B = 1
        self.B_boost = 0.1

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
        self.ntc_hsX = DensityHistory(
            nmaxold=self.nmaxold,
            ncomponents=self.ncomponents,
            desc=self.desc,
            xp=self.xp
        )
        self.Rc_hsX = DensityHistory(
            nmaxold=self.nmaxold,
            ncomponents=self.ncomponents,
            desc=self.desc,
            xp=self.xp
        )
        self.MRc_hsX = DensityHistory(
            nmaxold=self.nmaxold,
            ncomponents=self.ncomponents,
            desc=self.desc,
            xp=self.xp
        )
        self.histories = [self.nt_hsX, self.ntc_hsX, self.R_hsX,
                          self.Rc_hsX, self.MRc_hsX]

    def mix(self, density: Density) -> float:
        # Step 1: Initialize
        nold = self.nt_hsX.nold
        nt_sX = density.nt_sR
        D_asii = density.D_asii
        if nold == 0:
            self.nt_hsX.add_density(density)
            self.add_compensation_charge(nt_sX, D_asii, density,
                                         out=self.ntc_hsX.next())
            return np.inf

        R_sX, R_asii = self.R_hsX.next()
        self.calculate_residual(nt_sX, self.nt_hsX[-1][0], R_sX)
        self.calculate_paw_residual(D_asii.data, self.nt_hsX[-1][1], R_asii)

        Rc_sX = self.Rc_hsX.next()
        self.add_compensation_charge(R_sX, R_asii, density, out=Rc_sX)

        MRc_sX = self.MRc_hsX.next()
        self.metric(Rc_sX, out=MRc_sX)

        # Step 2: Do a pratt-step if nold <= 1
        if nold <= 1:
            return self.pratt_step(density)

        # Step 3: Calculate the multisecants
        s_hsX = self.ntc_hsX.to_multisecant()
        y_hsX = self.Rc_hsX.to_multisecant()
        my_hsX = self.MRc_hsX.to_multisecant()

        Ay_hh = y_hsX.matrix_elements(my_hsX).data
        As_hh = s_hsX.matrix_elements(my_hsX).data

        # Step 4: Decide on good broydenness
        good_broydenness = self.decide_good_broydenness(
            Ay_hh, As_hh)

        # Step 5: Rescale the multisecants and define t_hsX
        y_norm = self.xp.diag(Ay_hh)
        s_norm = self.xp.diag(As_hh)
        t_norm = y_norm + good_broydenness * s_norm
        if (t_norm < 0).any():
            # Somethings gone terribly wrong
            good_broydenness = 0
            t_norm = y_norm
        t_norm = 1 / self.xp.sqrt(t_norm)
        y_hsX.matrix.data *= t_norm[:, None]
        s_hsX.matrix.data *= t_norm[:, None]
        my_hsX.matrix.data *= t_norm[:, None]
        t_hsX = y_hsX.copy()
        t_hsX.data += good_broydenness * s_hsX.data

        # Step 6: Calculate the mixing weights
        A_hh = t_hsX.matrix_elements(my_hsX)
        B_hh = t_hsX.matrix_elements(s_hsX)

        A_diag = self.xp.diag(A_hh)
        B_diag = self.xp.diag(B_hh)

        S, V, D = self.xp.linalg.svd(A_hh)
        V /= V**2 + A_diag**2 * self.reg**2
        A_hh = D.T @ self.xp.diag(V) @ S.T

        S, V, D = self.xp.linalg.svd(B_hh)
        V /= V**2 + B_diag**2 * self.reg**2
        B_hh = D.T @ self.xp.diag(V) @ S.T

        MRc_1sX = MRc_sX.new(data=MRc_sX.data[None], dims=(1,) + MRc_sX.dims)
        alpha_h = t_isG.matrix_elements(MRc_1sX)[:, 0]
        alpha_h = A_hh @ alpha_h

        # TODO: Ensure ranks agree

        # Step 7: Predict mixing coefficients
        tmp_1sX = self.uk_1sX.copy()
        tmp_1sX.data -= self.Rc_sX.data
        tmp_1sX.data *= -1
        A1 = self.uk_1sX.norm2().sum()
        B1 = tmp_1sX.norm2().sum()

        A2_i = t_hsX.matrix_elements(self.uk_1sX)
        A2_j = self.uk_1sX.matrix_elements(y_hsX)
        A2 = A2_i @ B_hh @ A2_j

        B2_i = t_hsX.matrix_elements(self.pk_1sX)
        B2_j = tmp_1sX.matrix_elements(y_hsX)
        B2 = B2_i @ B_hh @ B2_j

        trig_fact = self.A_lims[-1] * 2 / np.pi
        A_target = np.clip(
            np.arctan(np.abs(A1 / (A2 * trig_fact))) * trig_fact,
            *self.A_lims
        )
        B_target = np.abs(B1 / B2) + self.B_boost
        if self.nold > 2:
            A_ratio = np.sqrt(A_target * self.A) / self.A
            self.A *= np.clip(A_ratio, *self.rate_ratio)
            self.A = np.clip(self.A, *self.A_lims)
            B_ratio = (self.B + B_target) / self.B
            self.B *= np.clip(B_ratio, *self.rate_ratio)
            self.B = np.clip(self.B, *self.B_lims)
        else:
            self.A = A_target

        A = self.A
        B = self.B

        # Step 8: Trust region control
        tmp_1sX.data[:] = A * self.uk_1sX.data
        tmp_1sX.data[:] += B * self.pk_1sX.data
        

        # Step 9: Mix the densities

        # Step 10: Update density history

        # Step 11: Return the mixing error

    def decide_good_broydenness(self,
                                Ay_hh: ArrayND,
                                As_hh: ArrayND) -> float:
        Ay_norm = self.xp.linalg.norm(Ay_hh, ord='fro')
        As_norm = self.xp.linalg.norm(As_hh, ord='fro')
        max_gb = max(Ay_norm / As_norm, 1)
        y_norm = self.xp.diag(Ay_hh)
        s_norm = self.xp.diag(As_hh)
        fracs_i = y_norm / s_norm
        if (fracs_i <= 0).any():
            # Handle negative values in fracs_i
            # Set max_gb to a safe value
            max_gb = min(
                max_gb,
                self.xp.min(-fracs_i[fracs_i <= 0])
            )

        good_broydenness = 0.5 * max_gb
        for iter in range(2, 12):
            norm_vec = 1 / self.xp.sqrt(
                y_norm + s_norm * good_broydenness
            )
            A_hh = Ay_hh + As_hh * good_broydenness
            A_hh *= norm_vec[None, :]
            A_hh *= norm_vec[:, None]
            try:
                eigs = self.xp.linalg.eigvals(A_ii)
                min_real = self.xp.min(eigs.real)
                max_imag = self.xp.max(self.xp.abs(eigs.imag))
            except np.linalg.LinAlgError:
                good_broydenness -= 2**(-iter) * max_gb
                continue
            if min_real > max_imag:
                good_broydenness += 2**(-iter) * max_gb
            else:
                good_broydenness -= 2**(-iter) * max_gb
        good_broydenness -= 2**(-iter) * max_gb
        good_broydenness *= self.gb_scale

    def pratt_step(self, density: Density) -> float:
        nt_sX = density.nt_sR
        D_asii = density.D_asii
        nt_sX.data[:] = self.nt_hsX[-1][0].data
        D_asii.data[:] = self.nt_hsX[-1][1]
        nt_sX.data += self.beta * self.R_hsX[-1][0].data
        D_asii.data += self.beta * self.R_hsX[-1][1]
        self.nt_hsX.add_density(density)
        self.add_compensation_charge(nt_sX, D_asii, density,
                                     out=self.ntc_hsX.next())
        return self.calculate_charge_sloshing(self.Rc_hsX[-1])
