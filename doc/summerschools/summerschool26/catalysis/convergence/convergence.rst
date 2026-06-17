.. _catalysis_convergence:

Convergence checks
==================

In this notebook we look at the adsorption energy and height of a nitrogen
atom on a Ru(0001) surface in the hcp site.  We check for convergence with
respect to:

* number of layers
* number of k-points in the BZ
* plane-wave cutoff energy


Nitrogen atom
-------------

First step is an isolated nitrogen atom which has a magnetic moment of 3.
More information:
[Atoms][1] and [GPAW parameters][2].

[1]: https://ase-lib.org/ase/atoms.html#ase.Atoms
[2]: https://gpaw.readthedocs.io/documentation/basic.html#parameters

.. literalinclude:: eads.py


Clean slab
----------

We use the [ase.build.hcp0001()][3] function to build the Ru(0001) surface.

[3]: https://ase-lib.org/ase/build/surface.html#ase.build.hcp0001

.. literalinclude:: eads.py


N/Ru(0001)
----------

Now, let's add a nitrogen atom in the "HCP" site:

.. literalinclude:: eads.py

Alternatively, you can just use the [add_adsorbate()][4] function:

[4]: https://ase-lib.org/ase/build/surface.html#ase.build.add_adsorbate
"""

.. literalinclude:: eads.py

.. image:: nru2.png

Now, calculate the total energy and the unrelaxed adsorption energy:

.. literalinclude:: eads.py

If you also calculate for forces (``f = nslab.get_forces()``),
you will see that the force on the N-atom is quite big (``print(f[-1])``).
Let's freeze the surface and relax the
adsorbate.  We use
[ase.optimize.BFGSLineSearch][5] and
[ase.constraints.FixAtoms][6]
for this task.

[5]: https://ase-lib.org/ase/optimize.html#module-ase.optimize
[6]: https://ase-lib.org/ase/constraints.html#ase.constraints.FixAtoms

.. literalinclude:: eads.py


Convergence
-----------

In order to make it easy to check for convergence of the adsorption energy
and height we write functions (called ``adsorp()`` and ``atom_energy()``)
that does all of the stuff above taking
``nlayers``, ``nkpts`` and ``ecut`` as input parameters.

Take a look at the :download:`convergence.py` script and try to
understand what it does:

.. literalinclude:: convergence.py

If you run the script, the ``main()`` function will run and all the
convergence calculations will be run.  You can run the script yourself or
you can just download the results here: :download:`results.json`. We will
now analyse the results.

Try this:

.. literalinclude:: plot.py
   :begin-after: web-page
   :end-after: layers

This should produce this:

.. image:: layers.svg

Here are the plots for **k**-point and `E_{text{cut}}` convergence:

.. image:: kpts.svg
.. image:: ecut.svg


Conclusion
----------

For accurate calculations you would need:

* a plane-wave cutoff of 600 eV
* 5 layers of Ru
* 9x9 Monkhorst-Pack grid for BZ sampling (for a 1x1 unit cell)

For our quick'n'dirty calculations we will use 350 eV, 2 layers and a
4x4 `\Gamma`-centered Monkhorst-Pack grid (for a 2x2 unit cell).
