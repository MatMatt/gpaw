"""Quick residual diagnostic — reduced grid."""
import numpy as np
from gpaw.old.grid_descriptor import GridDescriptor
from gpaw.solvation.poisson import WeightedFDPoissonSolver
from gpaw.fd_operators import Laplace, Gradient

gd = GridDescriptor((16, 16, 64), (2.885, 4.997, 20.711), pbc_c=(True, True, False))
zv = np.linspace(0, 20.711, 64)
eps = gd.zeros(); dze = gd.zeros()
prof = 1.0 + 77.0 * 0.5 * (1.0 - np.tanh((zv-10.0)/0.67))
for k in range(64): eps[:,:,k] = prof[k]
dz = 20.711/64
for k in range(1, 63):
    dze[:,:,k] = (prof[k+1]-prof[k-1])/(2*dz)
class MD: eps_gradeps = [eps, gd.zeros(), gd.zeros(), dze]

rho = gd.zeros()
for k in range(64):
    rho[:,:,k] = np.exp(-(zv[k]-8.0)**2/0.5) - np.exp(-(zv[k]-12.0)**2/0.5)

# Run to see residual WITHOUT converging fully
solver = WeightedFDPoissonSolver(nn=4, relax='J', eps=1e-20, maxiter=50)
solver.set_grid_descriptor(gd); solver.set_dielectric(MD()); phi = gd.zeros()
try:
    niter = solver.solve(phi, rho)
except:
    niter = "maxiter"

# Compute residual
laplace = Laplace(gd, 1.0, 4)
grad_z = Gradient(gd, 2, 1.0, 4)
lap_phi = gd.zeros(); laplace.apply(phi, lap_phi)
gz_phi = gd.zeros(); grad_z.apply(phi, gz_phi)
eps_g, _, _, dz_e = MD().eps_gradeps
residual = eps_g * lap_phi + dz_e * gz_phi - rho

res_z = abs(residual).mean(axis=(0,1))
eps_z = eps_g.mean(axis=(0,1))

print(f"After {niter} iters:")
print(f"{'z':>7s}  {'|R|_avg':>10s}  {'eps':>7s}")
for k in range(64):
    if 7 < zv[k] < 14:
        bar = '*' * int(abs(residual[:,:,k]).max() * 10) if abs(residual[:,:,k]).max() > 0.1 else ''
        print(f"{zv[k]:7.2f}  {res_z[k]:10.2e}  {eps_z[k]:7.1f} {bar}")

bm = (zv > 7) & (zv < 13)
total = np.sqrt((residual**2).sum())
boundary = np.sqrt((residual[bm]**2).sum())
print(f"\nTotal L2: {total:.2e}")
print(f"Boundary (±3Å): {boundary:.2e} ({100*boundary/total:.0f}% of total)")
print(f"Rest: {np.sqrt(total**2-boundary**2):.2e}")

# Show what dominates
eps_lap = abs(eps_g * lap_phi).mean(axis=(0,1))
grad_term = abs(dz_e * gz_phi).mean(axis=(0,1))
print(f"\n{'z':>7s}  {'|eps*Lap|':>10s}  {'|dz_eps*dz_phi|':>12s}")
for k in range(64):
    if 7 < zv[k] < 14:
        print(f"{zv[k]:7.2f}  {eps_lap[k]:10.2e}  {grad_term[k]:12.2e}")