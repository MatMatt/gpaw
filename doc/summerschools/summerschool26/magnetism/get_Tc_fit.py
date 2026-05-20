from numpy import tanh, log

A = 0
J = 0

T0 = 1.52 # Ising limit
T_c = T0 * S**2 * J / kB * (tanh(6 / N * log(1 - 0.033 * A / J)))**(0.25)  # Student version: T_c = ???
print(f'T_c = {T_c:.1f} K')
