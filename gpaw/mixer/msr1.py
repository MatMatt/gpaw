import numpy as np

from scipy.optimize import root_scalar

from gpaw.new.density import Density
from gpaw.mixer.base import DensityHistory, BaseMixer
from gpaw.typing import ArrayND


class MSR1Mixer(BaseMixer):
    def __init__(self,
                 nmaxold: int = 10,
                 beta: float = 0.05,
                 reg: float = 4e-3,
                 gb_scale: float = 1.0,
                 max_A: float = 0.5,
                 trust_scale: float = 3.0,
                 soft_lim: float = 2.0,
                 hard_lim: float = 3.0):
        self.nmaxold = nmaxold
        self.beta = beta
        self.reg = reg
        self.gb_scale = gb_scale
        self.trust_scale = trust_scale
        self.soft_lim = soft_lim
        self.hard_lim = hard_lim
        self.A_lims: list[float] = [0.02, max_A]
        self.B_lims: list[float] = [0.4, 1.0]
        self.rate_ratio: list[float] = [0.7, 1.3]
        self.A: float = 0.0
        self.B: float = 1.0
        self.B_boost = 0.1
        self.trust_radius: None | float = None

    def _initialize_history_(self):
        self.nt_hsX = DensityHistory(
            nmaxold=self.nmaxold,
            ncomponents=self.ncomponents,
            desc=self.desc,
            atom_layout=self.atom_layout,
            xp=self.xp)
        self.R_hsX = DensityHistory(
            nmaxold=self.nmaxold,
            ncomponents=self.ncomponents,
            desc=self.desc,
            atom_layout=self.atom_layout,
            xp=self.xp)
        self.ntc_hsX = DensityHistory(
            nmaxold=self.nmaxold,
            ncomponents=self.ncomponents,
            desc=self.desc,
            xp=self.xp)
        self.Rc_hsX = DensityHistory(
            nmaxold=self.nmaxold,
            ncomponents=self.ncomponents,
            desc=self.desc,
            xp=self.xp)
        self.MRc_hsX = DensityHistory(
            nmaxold=self.nmaxold,
            ncomponents=self.ncomponents,
            desc=self.desc,
            xp=self.xp)
        self.histories = [self.nt_hsX, self.ntc_hsX, self.R_hsX,
                          self.Rc_hsX, self.MRc_hsX]
        self.uk_1sX = self.desc.empty((1, self.ncomponents), xp=self.xp)
        self.pk_1sX = self.desc.empty((1, self.ncomponents), xp=self.xp)
        self.last_dNt = np.inf

    def mix(self, density: Density) -> float:
        # Step 1: Initialize
        nold = self.nt_hsX.nold
        nt_sX = density.nt_sR
        D_asii = density.D_asii
        if nold == 0:
            self.nt_hsX.add_density(density)
            self.add_compensation_charge(nt_sX, D_asii, density,
                                         out=self.ntc_hsX.next(),
                                         add_delta0=True)
            return np.inf

        _R_sX, _R_asii = self.R_hsX.next()
        self.calculate_residual(nt_sX, self.nt_hsX[-1][0], _R_sX)
        self.calculate_paw_residual(D_asii.data, self.nt_hsX[-1][1], _R_asii)

        self.add_compensation_charge(_R_sX, _R_asii, density,
                                     out=self.Rc_hsX.next())

        self.metric(self.Rc_hsX[-1], out=self.MRc_hsX.next())

        # Step 2: Do a pratt-step if nold <= 1, else discard bad steps
        if nold <= 1:
            return self.pratt_step(density)

        dNt = self.calculate_charge_sloshing(self.Rc_hsX[-1])
        increased_error = dNt / self.last_dNt
        if increased_error > self.soft_lim:
            dNt = self.last_dNt
            insert_pos = 0 if increased_error > self.hard_lim else -2
            for hist in self.histories:
                last_ind = hist.current_indicies.pop(-1)
                hist.current_indicies.insert(insert_pos, last_ind)

        # Step 3: Calculate the multisecants
        s_hsX = self.ntc_hsX.to_multisecant()
        y_hsX = self.Rc_hsX.to_multisecant()
        y_hsX.data *= -1
        my_hsX = self.MRc_hsX.to_multisecant()
        my_hsX.data *= -1

        Ay_hh = y_hsX.matrix_elements(my_hsX).data
        As_hh = s_hsX.matrix_elements(my_hsX).data

        # Step 4: Decide on good broydenness
        good_broydenness = self.decide_good_broydenness(
            Ay_hh, As_hh)
        good_broydenness *= (nold / self.nmaxold)**4

        # Step 5: Rescale the multisecants and define t_hsX
        y_norm = self.xp.diag(Ay_hh)
        s_norm = self.xp.diag(As_hh)
        t_norm = y_norm + good_broydenness * s_norm
        if (t_norm < 0).any():
            # Somethings gone terribly wrong
            good_broydenness = 0
            t_norm = y_norm
        t_norm = 1 / t_norm**0.5
        y_hsX.matrix.data *= t_norm[:, None]
        s_hsX.matrix.data *= t_norm[:, None]
        my_hsX.matrix.data *= t_norm[:, None]
        t_hsX = y_hsX.copy()
        t_hsX.data += good_broydenness * s_hsX.data

        # Step 6: Calculate the mixing weights
        A_hh = t_hsX.matrix_elements(my_hsX).data
        B_hh = t_hsX.matrix_elements(s_hsX).data

        A_diag = self.xp.mean(self.xp.diag(A_hh))
        B_diag = self.xp.mean(self.xp.diag(B_hh))

        S, V, D = self.xp.linalg.svd(A_hh)
        V /= V**2 + A_diag**2 * self.reg**2
        A_hh = D.T @ self.xp.diag(V) @ S.T

        S, V, D = self.xp.linalg.svd(B_hh)
        V /= V**2 + B_diag**2 * self.reg**2
        B_hh = D.T @ self.xp.diag(V) @ S.T

        MRc_sX = self.MRc_hsX[-1]
        MRc_1sX = MRc_sX.new(data=MRc_sX.data[None], dims=(1,) + MRc_sX.dims)
        BR_h = t_hsX.matrix_elements(MRc_1sX).data[:, 0]
        alpha_h = A_hh @ BR_h
        self.world.broadcast(alpha_h, 0)

        # Step 7: Predict mixing coefficients
        tmp_1sX = self.uk_1sX.copy()
        tmp_1sX.data -= self.Rc_hsX[-1].data
        tmp_1sX.data *= -1
        A1 = self.uk_1sX.norm2().sum()
        B1 = tmp_1sX.norm2().sum()
        A1 = A1 if self.xp is np else A1.get()
        B1 = B1 if self.xp is np else B1.get()

        A2_i = t_hsX.matrix_elements(self.uk_1sX).data[:, 0]
        A2_j = self.uk_1sX.matrix_elements(y_hsX).data[0, :]
        A2 = A2_i @ B_hh @ A2_j
        A2 = A2 if self.xp is np else A2.get()

        B2_i = t_hsX.matrix_elements(self.pk_1sX).data[:, 0]
        B2_j = tmp_1sX.matrix_elements(y_hsX).data[0, :]
        B2 = B2_i @ B_hh @ B2_j
        B2 = B2 if self.xp is np else B2.get()
        trig_fact = self.A_lims[-1] * 2 / np.pi
        A_target = np.clip(
            np.arctan(np.abs(A1 / (A2 * trig_fact))) * trig_fact,
            self.A_lims[0], self.A_lims[1]
        )
        if nold > 2:
            B_target = np.abs(B1 / B2) + self.B_boost
            A_ratio = np.sqrt(A_target * self.A) / self.A
            self.A *= np.clip(A_ratio, self.rate_ratio[0], self.rate_ratio[1])
            self.A = np.clip(self.A, self.A_lims[0], self.A_lims[1])
            B_ratio = (self.B + B_target) / (2 * self.B)
            self.B *= np.clip(B_ratio, self.rate_ratio[0], self.rate_ratio[1])
            self.B = np.clip(self.B, self.B_lims[0], self.B_lims[1])
        else:
            self.A = A_target

        A = self.A
        B = self.B

        # Step 8: Trust region control
        tmp_1sX.data[:] = self.uk_1sX.data * A
        tmp_1sX.data[:] += self.pk_1sX.data * B
        trust_radius = self.trust_scale * tmp_1sX.norm2().sum()**0.5
        if self.trust_radius is None:
            self.trust_radius = trust_radius
        else:
            self.trust_radius = (self.trust_radius + trust_radius) * 0.5

        self.pk_1sX.data[:] = 0
        self.uk_1sX.data[:] = self.Rc_hsX[-1].data

        for h, alpha in enumerate(alpha_h):
            self.pk_1sX.data[:] += alpha * s_hsX[h].data
            self.uk_1sX.data[:] -= alpha * y_hsX[h].data

        predicted_radius = B * self.pk_1sX.norm2().sum()**0.5
        if predicted_radius > self.trust_radius * 1.02:
            s_hh = s_hsX.matrix_elements(s_hsX).data
            s_hh *= B**2
            A_hh = self.xp.linalg.inv(A_hh)

            def err_fct(lamb):
                beta_h = self.xp.linalg.solve(
                    A_hh + lamb * self.xp.eye(nold - 1), BR_h
                )
                rtnval = (beta_h @ s_hh @ beta_h) - self.trust_radius**2
                return rtnval if self.xp is np else rtnval.get()

            try:
                upplimscale = A_diag if self.xp is np else A_diag.get()
                lamb = root_scalar(err_fct, bracket=[0, 1000 * upplimscale])
                root = lamb.root
            except ValueError:
                root = 1000 * A_diag
            beta_h = self.xp.linalg.solve(
                A_hh + root * self.xp.eye(nold - 1), BR_h
            )

            scale_factor = self.trust_radius / predicted_radius
            A *= np.clip(scale_factor, 0, 1)
            A = max(A, self.A_lims[0])
        else:
            beta_h = alpha_h
        self.world.broadcast(beta_h, 0)

        # Step 9: Finally mix the densities
        nt_sX.data[:] = self.nt_hsX[-1][0].data
        nt_sX.data[:] += self.R_hsX[-1][0].data * A
        D_asii.data[:] = self.nt_hsX[-1][1]
        D_asii.data[:] += self.R_hsX[-1][1] * A

        beta_h *= t_norm

        for h, beta in enumerate(beta_h):
            nt_sX.data += B * beta * (
                self.nt_hsX[h][0].data - self.nt_hsX[-1][0].data)
            D_asii.data += B * beta * (self.nt_hsX[h][1] - self.nt_hsX[-1][1])
            nt_sX.data += A * beta * (
                self.R_hsX[h][0].data - self.R_hsX[-1][0].data)
            D_asii.data += A * beta * (self.R_hsX[h][1] - self.R_hsX[-1][1])

        # Step 10: Update density history
        if increased_error > self.hard_lim:
            for hist in self.histories:
                hist.delete_oldest()
        self.nt_hsX.add_density(density)
        self.add_compensation_charge(nt_sX, D_asii, density,
                                     out=self.ntc_hsX.next(),
                                     add_delta0=True)

        # Step 11: Return the mixing error
        return dNt

    def decide_good_broydenness(self, Ay_hh: ArrayND,
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
            max_gb = min(max_gb,
                         self.xp.min(-fracs_i[fracs_i <= 0]))

        good_broydenness = 0.5 * max_gb
        for iter in range(2, 12):
            norm_vec = 1 / (y_norm + s_norm *
                            good_broydenness)**0.5
            A_hh = Ay_hh + As_hh * good_broydenness
            A_hh *= norm_vec[None, :]
            A_hh *= norm_vec[:, None]
            try:
                eigs = self.xp.linalg.eigvals(A_hh)
                min_real = self.xp.min(eigs.real)
                max_imag = self.xp.max(self.xp.abs(eigs.imag))
            except Exception:
                good_broydenness -= 2**(-iter) * max_gb
                continue
            if min_real > max_imag:
                good_broydenness += 2**(-iter) * max_gb
            else:
                good_broydenness -= 2**(-iter) * max_gb
        good_broydenness -= 2**(-iter) * max_gb
        good_broydenness *= self.gb_scale
        return good_broydenness

    def pratt_step(self, density: Density) -> float:
        nt_sX = density.nt_sR
        D_asii = density.D_asii
        nt_sX.data[:] = self.nt_hsX[-1][0].data
        D_asii.data[:] = self.nt_hsX[-1][1]
        nt_sX.data += self.beta * self.R_hsX[-1][0].data
        D_asii.data += self.beta * self.R_hsX[-1][1]

        self.nt_hsX.add_density(density)
        self.add_compensation_charge(nt_sX, D_asii, density,
                                     out=self.ntc_hsX.next(),
                                     add_delta0=True)
        self.uk_1sX.data[:] = self.Rc_hsX[-1].data
        self.pk_1sX.data[:] = 0
        self.last_dNt = self.calculate_charge_sloshing(self.Rc_hsX[-1])
        return self.last_dNt
