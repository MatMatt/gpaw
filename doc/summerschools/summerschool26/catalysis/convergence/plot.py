# web-page: layers.svg, kpts.svg, ecut.svg
from matplotlib import pyplot as plt
from pathlib import Path
import json

results = json.loads(Path('results.json').read_text())

N = range(1, 10)
h = [results[(n, 7, 400)][0] for n in N]
ea = [results[(n, 7, 400)][1] for n in N]
fig, axs = plt.subplots(2, 1, sharex=True)
fig.subplots_adjust(hspace=0)
axs[0].plot(N, h)
axs[1].plot(N, ea)
axs[0].set_ylabel('height [Å]')
axs[1].set_ylabel('ads. energy [eV]')
axs[1].set_xlabel('number of layers')
# plt.show()
plt.savefig('layers.svg')

# %%
K = range(4, 18)
h = [results[(2, k, 400)][0] for k in K]
ea = [results[(2, k, 400)][1] for k in K]
fig, axs = plt.subplots(2, 1, sharex=True)
fig.subplots_adjust(hspace=0)
axs[0].plot(K, h)
axs[1].plot(K, ea)
axs[0].set_ylabel('height [Å]')
axs[1].set_ylabel('ads. energy [eV]')
axs[1].set_xlabel('number of k-points')
plt.savefig('kpts.svg')

# %%
E = range(350, 801, 50)
h = [results[(2, 7, e)][0] for e in E]
ea = [results[(2, 7, e)][1] for e in E]
fig, axs = plt.subplots(2, 1, sharex=True)
fig.subplots_adjust(hspace=0)
axs[0].plot(E, h)
axs[1].plot(E, ea)
axs[0].set_ylabel('height [Å]')
axs[1].set_ylabel('ads. energy [eV]')
axs[1].set_xlabel('plane-wave cutoff energy [eV]')
plt.savefig('ecut.svg')
