import warnings
from functools import cached_property
from math import pi

import numpy as np
from ase.units import Bohr, Ha
from scipy.sparse.linalg import LinearOperator, cg
from scipy.special import erf

from gpaw import get_scipy_version
from gpaw.core import PWArray, PWDesc, UGDesc
from gpaw.new.poisson import PoissonSolver

if get_scipy_version() >= [1, 14]:
    RTOL = 'rtol'
else:
    RTOL = 'tol'


def make_poisson_solver(pw: PWDesc,
                        grid: UGDesc,
                        charge: float,
                        strength: float = 1.0,
                        dipolelayer: bool = False,
                        **kwargs) -> PoissonSolver:
    if charge != 0.0 and not grid.pbc_c.any():
        return ChargedPWPoissonSolver(pw, grid, charge, strength, **kwargs)

    ps = PWPoissonSolver(pw, charge, strength)

    if dipolelayer:
        return DipoleLayerPWPoissonSolver(ps, grid, **kwargs)

    assert not kwargs

    return ps


class PWPoissonSolver(PoissonSolver):
    def __init__(self,
                 pw: PWDesc,
                 charge: float = 0.0,
                 strength: float = 1.0):
        self.pw = pw
        self.charge = charge
        self.strength = strength

        self.ekin_g = pw.ekin_G.copy()
        if pw.comm.rank == 0:
            # Avoid division by zero:
            self.ekin_g[0] = 1.0

    def __str__(self) -> str:
        txt = ('poisson solver:\n'
               f'  ecut: {self.pw.ecut * Ha}  # eV\n')
        if self.strength != 1.0:
            txt += f'  strength: {self.strength}\n'
        if self.charge != 0.0:
            txt += f'  uniform background charge: {self.charge}  # electrons\n'
        return txt

    def solve(self,
              vHt_g: PWArray,
              rhot_g: PWArray) -> float:
        """Solve Poisson equation.

        Places result in vHt_g ndarray.
        """
        epot = self._solve(vHt_g, rhot_g)
        return epot

    def _solve(self,
               vHt_g,
               rhot_g) -> float:
        vHt_g.data[:] = 2 * pi * self.strength * rhot_g.data
        if self.pw.comm.rank == 0:
            # Use uniform background charge in case we have a charged system:
            vHt_g.data[0] = 0.0
        if not isinstance(self.ekin_g, vHt_g.xp.ndarray):
            self.ekin_g = vHt_g.xp.array(self.ekin_g)
        vHt_g.data /= self.ekin_g
        epot = 0.5 * vHt_g.integrate(rhot_g)
        return epot


class ChargedPWPoissonSolver(PWPoissonSolver):
    def __init__(self,
                 pw: PWDesc,
                 grid: UGDesc,
                 charge: float,
                 strength: float = 1.0,
                 alpha: float = None,
                 eps: float = 1e-5):
        """Reciprocal-space Poisson solver for charged molecules.

        * Add a compensating Gaussian-shaped charge to the density
          in order to make the total charge neutral (placed in the
          middle of the unit cell

        * Solve Poisson equation.

        * Correct potential so that it has the correct 1/r
          asymptotic behavior

        * Correct energy to remove the artificial interaction with
          the compensation charge

        Parameters
        ----------
        pw: PWDesc
        grid: UGDesc
        charge: float
        strength: float
        alpha: float
        eps: float

        Attributes
        ----------
        alpha : float
        charge_g : np.ndarray
            Gaussian-shaped charge in reciprocal space
        potential_g : PWArray
             Potential in reciprocal space created by charge_g
        """
        super().__init__(pw, charge, strength)

        if alpha is None:
            # Shortest distance from center to edge of cell:
            rcut = 0.5 / (pw.icell**2).sum(axis=1).max()**0.5

            # Make sure e^(-alpha*rcut^2)=eps:
            alpha = -rcut**-2 * np.log(eps)

        self.alpha = alpha

        center_v = pw.cell_cv.sum(axis=0) / 2
        G2_g = 2 * pw.ekin_G
        G_gv = pw.G_plus_k_Gv
        self.charge_g = np.exp(-1 / (4 * alpha) * G2_g +
                               1j * (G_gv @ center_v))
        self.charge_g *= charge / pw.dv

        R_Rv = grid.xyz()
        d_R = ((R_Rv - center_v)**2).sum(axis=3)**0.5
        potential_R = grid.empty()

        # avoid division by 0
        zero_indx = d_R == 0
        d_R[zero_indx] = 1
        potential_R.data[:] = charge * erf(alpha**0.5 * d_R) / d_R
        # at zero we should have:
        # erf(alpha**0.5 * d_R) / d_R = alpha**0.5 * 2 / sqrt(pi)
        potential_R.data[zero_indx] = charge * alpha**0.5 * 2 / np.sqrt(pi)
        self.potential_g = potential_R.fft(pw=pw)

    def __str__(self) -> str:
        txt, x, _ = super().__str__().rsplit('\n', 2)
        assert x.startswith('  uniform background charge:')
        txt += (
            '\n  # using Gaussian-shaped compensation charge: e^(-alpha r^2)\n'
            f'  alpha: {self.alpha}   # bohr^-2')
        return txt

    def _solve(self,
               vHt_g,
               rhot_g) -> float:
        neutral_g = rhot_g.copy()
        neutral_g.data += self.charge_g

        if neutral_g.desc.comm.rank == 0:
            error = neutral_g.data[0]  # * self.pd.gd.dv
            assert error.imag == 0.0, error
            assert abs(error.real) < 0.00001, error
            neutral_g.data[0] = 0.0

        vHt_g.data[:] = 2 * pi * neutral_g.data
        vHt_g.data /= self.ekin_g
        epot = 0.5 * vHt_g.integrate(neutral_g)
        epot -= self.potential_g.integrate(rhot_g)
        epot -= self.charge**2 * (self.alpha / 2 / pi)**0.5
        vHt_g.data -= self.potential_g.data
        return epot


class DipoleLayerPWPoissonSolver(PoissonSolver):
    def __init__(self,
                 ps: PWPoissonSolver,
                 grid: UGDesc,
                 width: float = 1.0,  # Ångström
                 zero_vacuum=False):
        self.ps = ps
        self.grid = grid
        self.width = width / Bohr
        self.zero_vacuum = zero_vacuum
        assert grid.pbc_c.sum() == 2
        self.axis = np.where(~grid.pbc_c)[0][0]
        self.correction = np.nan
        self.pw = ps.pw

    def solve(self,
              vHt_g: PWArray,
              rhot_g: PWArray) -> float:
        epot = self.ps.solve(vHt_g, rhot_g)
        dip_v = -rhot_g.moment()
        c = self.axis
        L = self.grid.cell_cv[c, c]
        self.correction = 2 * np.pi * dip_v[c] * L / self.grid.volume
        vHt_g.data -= 2 * self.correction * self.sawtooth_g.data
        if self.zero_vacuum:
            v0 = vHt_g.boundary_value(self.axis)
            if vHt_g.desc.comm.rank == 0:
                vHt_g.data[0] += self.correction - v0
        return epot + 2 * np.pi * dip_v[c]**2 / self.grid.volume

    def dipole_layer_correction(self) -> float:
        return self.correction

    @cached_property
    def sawtooth_g(self) -> PWArray:
        grid = self.grid
        if grid.comm.rank == 0:
            c = self.axis
            L = grid.cell_cv[c, c]
            w = self.width / 2
            assert w < L / 2, (w, L, c)
            gc = int(w / L * grid.size_c[c])
            x = np.linspace(0, L, grid.size_c[c], endpoint=False)
            sawtooth = x / L - 0.5
            a = 1 / L - 0.75 / w
            b = 0.25 / w**3
            sawtooth[:gc] = x[:gc] * (a + b * x[:gc]**2)
            sawtooth[-gc:] = -sawtooth[gc:0:-1]
            sawtooth_r = grid.new(comm=None).empty()
            shape = [1, 1, 1]
            shape[c] = -1
            sawtooth_r.data[:] = sawtooth.reshape(shape)
            sawtooth_g = sawtooth_r.fft(pw=self.ps.pw.new(comm=None)).data
        else:
            sawtooth_g = None

        result_g = self.ps.pw.empty()
        result_g.scatter_from(sawtooth_g)
        return result_g


class FDPWsolver(PWPoissonSolver):
    def __init__(self,
                 pw: PWDesc,
                 grid,
                 dielectric,
                 charge: float = 0.0,
                 strength: float = 1.0,
                 eps=1e-9,
                 maxiter=1000,
                 real_space_solver=None,
                 dipolelayer: bool = True,
                 zero_vacuum=False):
        """Initialize the finite-difference Poisson solver.
        Parameters:
        """

        super().__init__(pw, charge, strength)
        if real_space_solver is None:
            from gpaw.solvation.poisson import WeightedFDPoissonSolver
            # I need to make eps 1e-9 to get convergence, equally I need
            # to set maxcharge to 1e-5 when solving. I think they are
            # connected
            real_space_solver = WeightedFDPoissonSolver(eps=eps,
                                                        maxiter=maxiter)

        self.dielectric = dielectric
        self.grid = grid
        self.pw0 = pw.new(comm=None)
        self.pwg0 = self.pw0
        self.grid0 = grid.new(comm=None)
        self.dipolelayer = dipolelayer
        self.correction = np.nan

        if pw.comm.rank == 0:
            self.ekin_g = self.pw0.ekin_G.copy()
            self.ekin_g[0] = 1.0

        if self.dipolelayer:
            self.corrterm = 1
            self.elcorr = None
            self.last_corrterm = None
            self.dirichlet = False #dirichlet

        self.eps = eps
        self.maxiter = maxiter

        self.eps0_R = None
        self.drho_g = None

        self.zero_vacuum = zero_vacuum
        self.real_space_solver = real_space_solver
        self.real_space_solver.set_dielectric(self.dielectric)
        self.real_space_solver.set_grid_descriptor(grid._gd)
        if zero_vacuum:
            self.drho_g = dipole_layer(grid).fft(pw=pw)

    def get_description(self):
        return 'FD Poisson Solver for PW mode'

    def __str__(self) -> str:
        txt = ('fd poissonsolver for pw mode:\n'
               f'  ecut: {self.pw.ecut * Ha}  # eV\n'
               f'  eps: {self.eps}\n'
               f'  real space solver: {self.real_space_solver}\n'
               f'  maxiter: {self.maxiter}\n')

        if self.strength != 1.0:
            txt += f'  strength: {self.strength}\n'
        if self.charge != 0.0:
            txt += f'  uniform background charge: {self.charge}  # electrons\n'
        return txt

    def dipole_layer_correction(self) -> float:
        return self.correction

    def _solve(self,
               vHt_g,
               rhot_g) -> float:
        rhot_r = self.grid.new(comm=self.grid.comm).empty()
        vHt_r = self.grid.new(comm=self.grid.comm).empty()
        rhot0_r, vHt0_r = None, None

        vHt0_g = vHt_g.gather()
        rhot0_g = rhot_g.gather()

        if rhot0_g is not None:
            rhot0_r = rhot0_g.ifft(grid=self.grid.new(comm=None))
            vHt0_r = vHt0_g.ifft(grid=self.grid.new(comm=None))

        rhot_r.scatter_from(rhot0_r)
        vHt_r.scatter_from(vHt0_r)

        if self.dipolelayer:
            gd = self.grid
            print(gd.__dict__.items())
            print(gd.cell_cv)
            print(gd.dv**(1/3))
            slope_lim = 1e-13
            slope = slope_lim * 10

            def calculate_dipole_moment(rho_r, center=False, origin_c=None):
                    """Calculate dipole moment of density."""
                    r_cz = [np.arange(gd.start_c[c], gd.end_c[c]) for c in range(3)]
                    if center:
                        assert origin_c is None
                        r_cz = [r_cz[c] - 0.5 * self.N_c[c] for c in range(3)]
                    elif origin_c is not None:
                        r_cz = [r_cz[c] - origin_c[c] for c in range(3)]

                    rho_01 = rho_r.data.sum(axis=2)
                    rho_02 = rho_r.data.sum(axis=1)
                    rho_cz = [rho_01.sum(axis=1), rho_01.sum(axis=0), rho_02.sum(axis=0)]
                    rhog_c = [np.dot(r_cz[c], rho_cz[c]) for c in range(3)]
                    #FIGURE OUT how to get h_cv or dipole moment directly
                    gd.h_cv = np.array([gd.cell_cv[c] / gd.size_c[c] for c in range(3)])
                    d_c = -np.dot(rhog_c, gd.h_cv) * gd.dv
                    gd.comm.sum(d_c)
                    return d_c

            dipmom = calculate_dipole_moment(rhot_r)[2]
            print(f'Initial dipole moment: {dipmom} e·Bohr')

            if self.elcorr is not None:
                vHt_r.data[:, :] -= self.elcorr

            self.real_space_solver.solve(vHt_r.data, rhot_r.data,
                                     maxcharge=1e-5)
#            iters2 = self.solve(vHt_g, rhot_g, **kwargs)

            from gpaw.new.sjm import modified_saw_tooth
            eps_r = self.grid.from_data(self.dielectric.eps_gradeps[0])
            eps0_r = eps_r.gather(broadcast=True)
            sawtooth_z = modified_saw_tooth(eps0_r) #sjm_sawtooth(dirichlet=self.dirichlet)
            L = gd.cell_cv[2, 2]
            count = 0
            while abs(slope) > slope_lim:
                count += 1
                vHt_r2 = vHt_r.copy()
                self.correction = 2 * np.pi * dipmom * L / \
                    gd.volume * self.corrterm

                elcorr = -2 * self.correction * sawtooth_z
                elcorr2 = elcorr[gd.start_c[2]:gd.end_c[2]]
                vHt_r2.data[:, :] += elcorr2

                vHt0_r = vHt_r2.gather(broadcast=True)
                if vHt0_r is not None:
                    vHt0_z = vHt0_r.data.mean(0).mean(0)

                    slope = (vHt0_z[3] - vHt0_z[8]) / (gd.h_cv[2][2] * Bohr)

                    print(f'FD Poisson dipole correction iteration {count}: '
                          f'slope = {slope}, corrterm = {self.corrterm}')
                    if abs(slope) > slope_lim:
                        if self.last_corrterm is None:
                            self.last_corrterm = self.corrterm
                            self.corrterm -= slope * 10.
                        else:
                            ds = (slope - self.last_slope) / \
                                (self.corrterm - self.last_corrterm)
                            con = slope - (ds * self.corrterm)
                            self.last_corrterm = self.corrterm
                            self.corrterm = -con / ds
                        self.last_slope = slope
                    else:
                        vHt_r.data[:, :] += elcorr2
                        self.elcorr = elcorr2
                self.corrterm
        else:
            self.real_space_solver.solve(vHt_r.data, rhot_r.data,
                                     maxcharge=1e-5)

        vHt0_r = vHt_r.gather()
        rhot0_r = rhot_r.gather()

        # The following still fails in parallel
        if vHt0_r is not None:
            vHt0_g = vHt0_r.fft(pw=self.pw.new(comm=None))
            rhot0_g = rhot0_r.fft(pw=self.pw.new(comm=None))

        rhot_g.scatter_from(rhot0_g)
        vHt_g.scatter_from(vHt0_g)

        epot = 0.5 * vHt_g.integrate(rhot_g)
        return epot

    def asolve(self,
              vHt_g: PWArray,
              rhot_g: PWArray) -> float:
        rhot_r = self.grid.new(comm=self.grid.comm).empty()
        vHt_r = self.grid.new(comm=self.grid.comm).empty()
        rhot0_r, vHt0_r = None, None

        vHt0_g = vHt_g.gather()
        rhot0_g = rhot_g.gather()

        if rhot0_g is not None:
            rhot0_r = rhot0_g.ifft(grid=self.grid.new(comm=None))
            vHt0_r = vHt0_g.ifft(grid=self.grid.new(comm=None))

        rhot_r.scatter_from(rhot0_r)
        vHt_r.scatter_from(vHt0_r)

        self.real_space_solver.solve(vHt_r.data, rhot_r.data,
                                     maxcharge=1e-5)

        vHt0_r = vHt_r.gather()
        rhot0_r = rhot_r.gather()

        # The following still fails in parallel
        if vHt0_r is not None:
            vHt0_g = vHt0_r.fft(pw=self.pw.new(comm=None))
            rhot0_g = rhot0_r.fft(pw=self.pw.new(comm=None))

        rhot_g.scatter_from(rhot0_g)
        vHt_g.scatter_from(vHt0_g)

        epot = 0.5 * vHt_g.integrate(rhot_g)
        return epot


class ConjugateGradientPoissonSolver(PWPoissonSolver):
    """Poisson solver using conjugate gradient method in reciprocal space.
    """

    def __init__(self,
                 pw: PWDesc,
                 grid,
                 dielectric,
                 charge: float = 0.0,
                 strength: float = 1.0,
                 eps=1e-4,
                 maxiter=15,
                 zero_vacuum=False):
        """Initialize the conjugate gradient Poisson solver.

        Parameters:
        -----------
        pw : PWDesc
            Plane wave descriptor
        charge : float, optional
            Total charge of the system
        strength : float, optional
            Scaling factor for the potential
        eps : float, optional
            Convergence threshold for conjugate gradient algorithm
        maxiter : int, optional
            Maximum number of iterations for the conjugate gradient algorithm
        """
        super().__init__(pw, charge, strength)
        self.dielectric = dielectric
        self.grid = grid
        self.pw0 = pw.new(comm=None)
        self.pwg0 = self.pw0
        self.grid0 = grid.new(comm=None)
        if pw.comm.rank == 0:
            self.ekin_g = self.pw0.ekin_G.copy()
            self.ekin_g[0] = 1.0

        self.eps = eps
        self.maxiter = maxiter

        self.eps0_R = None

        self.drho_g = None
        self.zero_vacuum = zero_vacuum
        if zero_vacuum:
            self.drho_g = dipole_layer(grid).fft(pw=pw)

    def __str__(self) -> str:
        txt = ('conjugate gradient poisson solver:\n'
               f'  ecut: {self.pw.ecut * Ha}  # eV\n'
               f'  eps: {self.eps}\n'
               f'  maxiter: {self.maxiter}\n')
        if self.strength != 1.0:
            txt += f'  strength: {self.strength}\n'
        if self.charge != 0.0:
            txt += f'  uniform background charge: {self.charge}  # electrons\n'
        return txt

    def get_description(self):
        return 'Conjugate Gradient Poisson Solver'

    def operator(self, phi_G):
        """Apply the generalized Poisson operator in reciprocal space.

        Parameters:
        -----------
        phi_q : ndarray
            Input potential in reciprocal space

        Returns:
        --------
        ndarray
            Result of operator application
        """
        pw = self.pw0
        grid = self.grid0

        G_vG = pw.G_plus_k_Gv.T

        ophi_G = np.zeros_like(phi_G)
        for G_G in G_vG:
            # Gradient in G-space is pw coefficient multiplied by potential?
            grad_G = pw.from_data(G_G * phi_G)
            # Transform to real space and multiply with dielectric function
            grad_R = grad_G.ifft(grid=grid)
            grad_R.data *= self.eps0_R.data
            # Transform back to G-space and multiply with G-vector again
            # Is this overal a ket-bra operation?
            ophi_G += grad_R.fft(pw=pw).data * G_G

        return ophi_G

    def _solve(self,
               vHt_g,
               rhot_g) -> float:
        vHt_g.data[:] = 4 * np.pi * self.strength * rhot_g.data

        eps_R = self.grid.from_data(self.dielectric.eps_gradeps[0])
        self.eps0_R = eps_R.gather()

        vHt0_g = vHt_g.gather()

        if self.pw.comm.rank == 0:
            vHt0_g.data[0] = 0.0

            N = len(vHt0_g.data)
            op = LinearOperator((N, N),
                                matvec=self.operator,
                                dtype=complex)

            # M is the preconditioner to help convergence
            M = LinearOperator((N, N),
                               matvec=lambda x: 0.5 * x / self.ekin_g,
                               dtype=complex)

            vHt0_g.data[:], info = cg(
                op, vHt0_g.data, maxiter=self.maxiter, M=M, **{RTOL: self.eps})
            if info != 0:
                warnings.warn(
                    f'Conjugate gradient did not converge (info={info})')

        vHt_g.scatter_from(vHt0_g)

        if self.zero_vacuum:
            self.zero_vacuum = False
            dphi_g = self.pw.zeros()
            self._solve(dphi_g, self.drho_g)
            v0s, v1s = xy_average_at_boundary(dphi_g)
            v0, v1 = xy_average_at_boundary(vHt_g)
            vHt_g.data -= dphi_g.data * (v1 / v1s)
            vHt_g.data[0] -= v0 - (v1 / v1s) * v0s
            self.zero_vacuum = True

        epot = 0.5 * vHt_g.integrate(rhot_g)
        return epot

    def correct_slope(self, vHt_g: PWArray):
        from gpaw.new.sjm import modified_saw_tooth
        eps_r = self.grid.from_data(self.dielectric.eps_gradeps[0])
        eps0_r = eps_r.gather()
        vHt0_g = vHt_g.gather()
        if eps0_r is not None:
            saw_tooth_z = modified_saw_tooth(eps0_r)
            vHt0_r = vHt0_g.ifft(grid=self.grid.new(comm=None))
            s1, s2 = saw_tooth_z[[2, 10]]
            v1, v2 = vHt0_r.data[:, :, [2, 10]].mean(axis=(0, 1))
            vHt0_r.data -= (v2 - v1) / (s2 - s1) * saw_tooth_z[np.newaxis,
                                                               np.newaxis]
            vHt0_r.data -= vHt0_r.data[:, :, -1].mean()
            vHt0_r.fft(out=vHt0_g)
        vHt_g.scatter_from(vHt0_g)


def dipole_layer(grid: UGDesc, z0: float = 2.0):
    a_r = grid.empty()
    z_r = grid.xyz()[:, :, :, 2]
    n = z_r.shape[2]
    h = grid.cell_cv[2, 2]
    d = h / n
    z0 = round(z0 / d) * d
    z_r += h / 2 - z0
    z_r %= h
    z_r -= h / 2
    alpha = 6.0 / z0**2
    a_r.data[:] = 4 * alpha**1.5 / np.pi**0.5 * z_r * np.exp(-alpha * z_r**2)
    return a_r


def xy_average_at_boundary(f_G: PWArray) -> np.ndarray:
    """Calculate average value and derivative at boundary of box."""
    pw = f_G.desc
    m0_G, m1_G = pw.indices_cG[:2, pw.ng1:pw.ng2] == 0
    mask_G = m0_G & m1_G
    f_q = f_G.data[mask_G]
    value = f_q.real.sum() * 2
    derivative = pw.G_plus_k_Gv[mask_G, 2] @ f_q.imag
    if pw.comm.rank == 0:
        value -= f_q[0].real
    result = np.array([value, derivative])
    pw.comm.sum(result)
    return result
