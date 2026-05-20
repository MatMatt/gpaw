.. _es2:

============================================
Quasiparticle bandgap using GW approximation
============================================


In this exercise we will calculate the quasiparticle band gap of the compound using the GW approximation. For a brief introduction to the GW theory and the details of its implementation in GPAW, see :ref:`GW theory <gw_theory>`

First, we need to do a regular groundstate calculation. To start off, we do this in plane wave mode and choose the LDA exchange-correlation functional. In order to keep the computational efforts small, you should start with reasonable k-points and plane wave basis.

.. literalinclude:: solution_es2.py
	:start-after: from ase.build
	:end-before: write out wavefunctions

Next, we set up the G0W0 calculator and calculate the quasi-particle spectrum for all the k-points present in the irreducible Brillouin zone from the ground state calculation and the specified bands. For example Carbon has 4 valence electrons and the bands are double occupied. Setting bands=(3,5) means including band index 3 and 4 which is the highest occupied band and the lowest unoccupied band.

.. literalinclude:: solution_es2.py
        :start-afer: from gpaw.response.g0w0
        :end-before: gw.calculate()


The dictionary is stored in ``???-g0w0_results.pckl``. From the dict it is for example possible to extract the direct bandgap at the Gamma point.

.. literalinclude:: solution_es2.py
        :start-after: import pickle
        :end-before: 'Direct bandgap of ???:', direct_gap

Can we trust the calculated value of the direct bandgap? A check for convergence with respect to the plane wave cutoff energy and number of k points is necessary. This is done by changing the respective values in the groundstate calculation and restarting.

For more details on the convergence and other parameters you can look at this tutorial: :ref:`GW tutorial <gw tutorial>`

Typical convergence would require ``ecut=300 eV``, ``8x8x8 k-points`` and ``ecut_extrapolation=.True.`` for Silicon.

When doing convergence checks, you will notice the G0W0 computations becoming increasingly expensive with higher k point grid and plane wave cutoff energy. In order to decrease the computational cost of the G0W0 computations, we can start with a groundstate calculation in LCAO mode using the ``dzp`` basis set. Here, we the number of bands should reflect the number of atomic orbitals present in the computation.

.. literalinclude:: solution_es2.py
        :start-after: calc = GPAW(mode="lcao",  # student mode lcao
        :end-before: calc.write('Si_lcao_groundstate.gpw'

Before starting the G0W0 computation the LCAO groundstate needs to be converted to plane wave mode.

.. literalinclude:: solution_es2.py
        :start-after: from gpaw.new.calculation import DFTCalculation 
        :end-before: dft.write_gpw_file('Si_pw_from_lcao_groundstate.gpw

The converted groundstate is stored in ``???_pw_from_lcao_groundstate.gpw``. From this groundstate the G0W0 computation can be started.

.. literalinclude:: solution_es2.py
        :start-after: gw = G0W0(calc='Si_pw_from_lcao_groundstate.gpw 
        :end-before: gw.calculate

Next, the G0W0 bandgap can be computed.

.. literalinclude:: solution_es2.py
        :start-after: results = pickle.load(open('Si-from-lcao-g0w0_results_GW.pckl'
	:end-before: print('Direct bandgap of Si starting from LCAO groundstate:


Again we can check for convergence with respect to the number of k points. Typical convergence would require ``8x8x8 k-points`` for Silicon. Why does the calculated value of the direct bandgap change compared to the plane wave mode ?


References:

   1. F. Hüser, T. Olsen, and K. S. Thygesen, “Quasiparticle GW calculations for solids, molecules, and two-dimensional materials”, Phys. Rev. B 87, 235132 (2013).
