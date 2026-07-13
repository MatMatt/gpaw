"""Diagnose residual spatial distribution."""
import numpy as np
from gpaw.old.grid_descriptor import GridDescriptor
from gpaw.solvation.poisson import WeightedFDPoissonSolver
from gpaw.fd_operators import Laplace, Gradient

gd = GridDescriptor((24, 32, 128), (2.885, 4.997, 20.711), pbc_c=(True, True, False))
zv = np.linspace(0, gd.cell_cv[2,2], gd.n_c[2])
eps = gd.zeros(); dze = gd.zeros()
prof = 1.0 + 77.0 * 0.5 * (1.0 - np.tanh((zv-10.0)/0.67))
for k in range(gd.n_c[2]): eps[:,:,k] = prof[k]
dz = gd.cell_cv[2,2]/gd.n_c[2]
for k in range(1, gd.n_c[2]-1):
    dze[:,:,k] = (prof[k+1]-prof[k-1])/(2*dz)
class MD: eps_gradeps = [eps, gd.zeros(), gd.zeros(), dze]

rho = gd.zeros()
for k in range(gd.n_c[2]):
    rho[:,:,k] = np.exp(-(zv[k]-8.0)**2/0.5) - np.exp(-(zv[k]-12.0)**2/0.5)

solver = WeightedFDPoissonSolver(nn=4, relax='J', eps=1e-12, maxiter=500)
solver.set_grid_descriptor(gd); solver.set_dielectric(MD()); phi = gd.zeros()
try:
    niter = solver.solve(phi, rho)
except Exception as e:
    print(f"FAILED: {e}")
    import sys; sys.exit(1)

print(f"Solver: {niter} iters (maxiter=200)")

laplace = Laplace(gd, 1.0, 4)
grad_z = Gradient(gd, 2, 1.0, 4)
lap_phi = gd.zeros(); laplace.apply(phi, lap_phi)
gz_phi = gd.zeros(); grad_z.apply(phi, gz_phi)
eps_g, _, _, dz_e = MD().eps_gradeps
residual = eps_g * lap_phi + dz_e * gz_phi - rho

res_z = abs(residual).mean(axis=(0,1))
eps_z = eps_g.mean(axis=(0,1))

print(f"\n{'z':>7s}  {'|R|_avg':>10s}  {'eps':>7s}")
for k in range(gd.n_c[2]):
    if 8 < zv[k] < 13:
        print(f"{zv[k]:7.2f}  {res_z[k]:10.2e}  {eps_z[k]:7.1f}")

bm = (zv > 7) & (zv < 13)
total = np.sqrt((residual**2).sum())
boundary = np.sqrt((residual[bm]**2).sum())
print(f"\nTotal L2: {total:.2e}")
print(f"Boundary (±3Å): {boundary:.2e} ({100*boundary/total:.0f}%)")
print(f"Rest: {np.sqrt(total**2-boundary**2):.2e}")

eps_lap = eps_g * lap_phi
grad_term = dz_e * gz_phi
print(f"|eps*Lap|_L2 = {np.sqrt((eps_lap**2).sum()):.2e}")
print(f"|dz_eps*dz_phi|_L2 = {np.sqrt((grad_term**2).sum()):.2e}")