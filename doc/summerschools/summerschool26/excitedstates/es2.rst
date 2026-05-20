.. _es2:

============================================
Quasiparticle bandgap using GW approximation
============================================


In this exercise we will calculate the quasiparticle band gap of the compound using the GW approximation. For a brief introduction to the GW theory and the details of its implementation in GPAW, see :ref:`GW theory <gw_theory>`

First, we need to do a regular groundstate calculation. We do this in plane wave mode and choose the LDA exchange-correlation functional. In order to keep the computational efforts small, you should start with reasonable k-points and plane wave basis.

.. literalinclude:: solution_es2.py
	:start-at: from ase.build
	:end-at: write out wavefunctions

Next, we set up the G0W0 calculator and calculate the quasi-particle spectrum for all the k-points present in the irreducible Brillouin zone from the ground state calculation and the specified bands. For example Carbon has 4 valence electrons and the bands are double occupied. Setting bands=(3,5) means including band index 3 and 4 which is the highest occupied band and the lowest unoccupied band.

.. literalinclude:: solution_es2.py
        :start-at: from gpaw.response.g0w0
        :end-at: gw.calculate()


The dictionary is stored in ``???-g0w0_results.pckl``. From the dict it is for example possible to extract the direct bandgap at the Gamma point.

.. literalinclude:: solution_es2.py
        :start-at: import pickle
        :end-at: 'Direct bandgap of ???:', direct_gap



Can we trust the calculated value of the direct bandgap? A check for convergence with respect to the plane wave cutoff energy and number of k points is necessary. This is done by changing the respective values in the groundstate calculation and restarting.

For more details on the convergence and other parameters you can look at this tutorial: :ref:`GW tutorial <gw tutorial>`

Typical convergence would require ``ecut=300 eV``, ``8x8x8 k-points`` and ``ecut_extrapolation=.True.`` for Carbon.

References:

   1. F. Hüser, T. Olsen, and K. S. Thygesen, “Quasiparticle GW calculations for solids, molecules, and two-dimensional materials”, Phys. Rev. B 87, 235132 (2013).
