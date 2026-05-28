======================================
Transport barriers and Voltage profile
======================================

You will calculate the energy barriers for transport of Li intercalated in the
graphite anode. You will examine how sensitive this barrier is to the
interlayer distance in graphite. You will also examine the energy of
intermediate states during the charge/discharge process. This will allow some
basic discussion of the voltage profile of the battery.

* Create initial and final structures for a NEB calculation, that will
  determine the transition state

  - If time permits you can study the influence of changing the interlayer
    graphite distance on the energy barrier.

* Create structures for a Li vacancy in :mol:`LiFePO_4` and a single Li
  in :mol:`FePO_4`

* Calculate the Li vacancy/insertion energies and compare them to the
  equilibrium potential

  - What can they tell you about the charge/discharge potential curves?

You will in general be provided less code than yesterday, especially
towards the end of this notebook.  You will have to use what you have
already seen and learned so far.

There will be some natural pauses while you wait for calculations to
finish.  If you do not finish this entire notebook today, do not despair.


Transport barrier of Li in graphite
===================================

You will now calculate the energy barrier for Li diffusion in the graphite
anode.  You will do this using the [Nudged Elastic Band (NEB)
method](https://ase-lib.org/ase/neb.html#module-ase.neb)

You can use your work from Day 2, but for simplicity you are advised to
load in the initial atomic configuration from file.

.. literalinclude:: li_barrier.py
   :end-before: snippet-final

You will now make a final structure, where the Li atom has been moved to a
neighbouring equivalent site.  The
[`get_positions`](https://ase-lib.org/ase/atoms.html?highlight=get_positions#ase.Atoms.get_positions),
[`set_positions`](https://ase-lib.org/ase/atoms.html?highlight=get_positions#ase.Atoms.set_positions)
and
[`get_cell`](https://ase-lib.org/ase/atoms.html?highlight=get_positions#ase.Atoms.get_cell)
functions are highly useful for such a task.  HINT: Displace the Li atom
$\frac{1}{n} (\vec{a}+\vec{b})$

.. code::

    final = initial.copy()
    ...
    ...

Visualize that you have made the final strcuture correctly.

.. code::

    from ase.visualize import view
    view(final)

Make a band consisting of 7 images including the initial and final.

.. literalinclude:: li_barrier.py
   :start-after: snippet-neb1
   :end-before: snippet-neb2

It this point `images` consist of 6 copies of `initial` and one entry of
`final`. Use the `NEB` method to create an initial guess for the minimum
energy path (MEP). In the cell below a simple interpolation between the
`initial` and `final` image is used as initial guess.

.. literalinclude:: li_barrier.py
   :start-after: snippet-neb2
   :end-before: snippet-constraint-gpaw

Add a ``view(images)`` line to your script and run it so that you can
visualize the NEB images.

It turns out, that while running the NEB calculation, the largest amount
of resources will be spend translating the carbon layer without any
noticeable buckling.  You will thus
[constrain](https://ase-lib.org/ase/constraints.html#constraints) the
positions of the carbon atoms to save computational time.

Each image in the NEB requires a unique calculator.

This very simple case is highly symmetric.  To better illustrate how the
NEB method works, the symmetry is broken using the
[rattle](https://ase-lib.org/ase/atoms.html#ase.Atoms.rattle) function.

.. literalinclude:: li_barrier.py
   :start-after: snippet-constraint-gpaw
   :end-before: snippet-initial-final

Start by calculating the energy and forces of the first (`initial`) and
last (`final`) images as this is not done during the actual NEB
calculation.

Note, that this can take a while if you opt to do it inside the notebook.

.. literalinclude:: li_barrier.py
   :start-after: snippet-initial-final
   :end-before: snippet-optimize

You can run the NEB calculation by running an optimization on the NEB
object the same way you would on an atoms object.  Note the `fmax` is
larger for this tutorial example than you would normally use.

.. literalinclude:: li_barrier.py
   :start-after: snippet-optimize

Submit the calculation to the HPC cluster.  Do this by first building a
complete script in the cell below using the cells above (minus the
`view()` commands). Make sure the cell runs and then interrupt the kernel.

.. code:: bash

    $ mq submit NEB.py -R 8:1h  # submits the calculation to 8 cores, 1 hour
    $ ...
    $ mq ls
    $ ...
    $ tail neb.log

You can move on while you wait for the calculation to finish.

Once the maximum force (`fmax`) in the log is below 0.1, the calculation
is finished.  Load in the full trajectory.

You will use the `ase gui` to inspect the result.  The below line reads in
the last 7 images in the file.  In this case the MEP images.

.. code::

    $ ase gui neb.traj@-7:

In the GUI use :menuselection:`Tools --> NEB`.

Now inspect how the TS image has developed.

.. code::

   $ ase gui neb.traj@3::7

For more complicated MEP's, use the [climbing image
method](https://ase-lib.org/ase/neb.html?highlight=neb#climbing-image) to
determine the transition state.  Why is it not required here?


Bonus
=====

You will now study the influence of changing the interlayer graphite
distance on the energy barrier.  Due to the high degree of symmetry, this
can be done easily in this case.  Load in the initial state (IS) and
transition state (TS) images from the converged MEP.

.. literalinclude:: li_barrier_2.py
   :end-before: snippet-barrier

Now calculate the energy of the initial state (IS) image and the
transition state (TS) image using
[`get_potential_energy()`](https://ase-lib.org/ase/atoms.html?highlight=get_potential_energy#ase.Atoms.get_potential_energy)

.. literalinclude:: li_barrier_2.py
   :start-after: snippet-barrier
   :end-before: snippet-strain

Why does this not fully align with what you found before?

New change the graphite layer distance by changing the the size of the
unit cell in the *z* direction by ±3 %. and use the same calculator
object as you did above and calculate the potential energy of the
compressed initial and final state.

.. literalinclude:: li_barrier_2.py
   :start-after: snippet-strain
   :end-before: snippet-gpaw

How does the energy barrier change?



FePO$_4$ with one Li
====================

You will now calculate the energy gain of adding a single Li atom into
the FePO$_4$ cell you made on Day 3. This corresponds to a charge of 25
%. You can compare this energy to the equilibrium potential.

Start preparing a new Python script (say, ``fepo4_1li.py``) and load in the FePO$_4$ structure you wrote to file on in a previous
exercise and add Li.  Assume that the cell dimension remain unchanged.

.. literalinclude:: fepo4_1li.py
   :end-before: snippet-li

Add a Li atoms at `(x,y,z)=(0,0,0)`:

.. literalinclude:: fepo4_1li.py
   :start-after: snippet-li
   :end-before: snippet-gpaw

Now finish the script:

* calculate the energy with the BEEF functional
* write result to a ``.traj`` file:
  ``from ase.io import write`` and
  ``write('fepo4_1li_out.traj', fepo4_1li)``
* calculate a BEEF-ensemble

and submit the job to the queue.

Once the calculation is finished, you are ready to calculate the energy
gained by intercalating a single Li ion into the cathode.  Start by
loading in the relevant reference structures and obtain the potential
energies.  This should not require any new DFT calculations.

.. literalinclude:: fepo4_1li.py
   :start-after: snippet-results
   :end-before: snippet-results-end

Calculate the energy of intercalting a single Li in the FePO$_4$ cell.
How does this energy compare with the equilibirum potential?  What can it
tell you about the charge/discharge potential curves?


Bonus: LiFePO$_4$ with one vacancy
==================================


If time permits, you will now do a similar calculation but this time with LiFePO$_4$ contraining one vacancy. Once again you should assume that the cell dimension remain unchanged compaired to LiFePO$_4$.

There are numerous ways to obtain this structure. You can get inspiration from the way LiFePO$_4$ was made on Exercise day 3, use the [`del` or `pop()` methods](https://ase-lib.org/ase/atoms.html?highlight=pop#list-methods), or even use the GUI to delete an atom and save the structure afterwards.


.. code::

    # In this cell you create the vacancy in LiFePO4

    # lifepo4_vac = ...

    # ...

    # teacher
    lifepo4_wo_li=read('lifepo4_wo_li.traj')
    from numpy import identity
    cell=lifepo4_wo_li.get_cell()
    xyzcell = identity(3)
    lifepo4_wo_li.set_cell(xyzcell, scale_atoms=True)  # Set the unit cell and rescale
    #lifepo4_wo_li.append(Atom('Li', (0, 0, 0)))
    lifepo4_wo_li.append(Atom('Li', (0, 0.5, 0)))
    lifepo4_wo_li.append(Atom('Li', (0.5, 0.5, 0.5)))
    lifepo4_wo_li.append(Atom('Li', (0.5, 0, 0.5)))
    lifepo4_wo_li.set_cell(cell, scale_atoms=True)
    lifepo4_vac=lifepo4_wo_li.copy()


Visualize the structure


.. code::

    view(lifepo4_vac)


Now ensure that the total magnetic moment is equal to 17.


.. code::

    for atom in fepo4_1li:
        if atom.symbol == 'Fe':
            atom.magmom = 4.25

    print(sum(fepo4_1li.get_initial_magnetic_moments()))


Write your atoms object to file giving it the name `lifepo4_vac.traj`.


.. code::

    # ...

    # teacher
    write('lifepo4_vac.traj', lifepo4_vac)


Make a full script in the cell below similar to that you made above. Make sure the cell runs before interupting the notebook kernel.


.. code::

    # %%writefile 'lifepo4_vac.py'
    # from ase.parallel import paropen
    # from ase.io import read, write
    # from ase.dft.bee import BEEFEnsemble
    # from gpaw import GPAW, FermiDirac, Mixer, PW

    # Read in the structure you made and wrote to file above
    # lifepo4_vac = read('lifepo4_vac.traj')


    # ...

    # write('lifepo4_vac_out.traj', lifepo4_vac)

    # ens = BEEFEnsemble(calc)
    # dE = ens.get_ensemble_energies(2000)
    # with paropen('ensemble_lifepo4_vac.dat','a') as results:
    #     for e in dE:
    #         print(e, file=result)

    # teacher
    from ase.parallel import paropen
    from ase.io import read, write
    from ase.dft.bee import BEEFEnsemble
    from gpaw import GPAW, FermiDirac, Mixer, PW

    # Read in the structure you made and wrote to file above
    lifepo4_vac = read('lifepo4_vac.traj')

    params_GPAW = {}
    params_GPAW['mode']        = PW(500)                     #The used plane wave energy cutoff
    params_GPAW['nbands']      = -40                           #The number on empty bands had the system been spin-paired
    params_GPAW['kpts']        = {'size':  (2,4,5),            #The k-point mesh
                                  'gamma': True}
    params_GPAW['spinpol']     = True                          #Performing spin polarized calculations
    params_GPAW['xc']          = 'BEEF-vdW'                    #The used exchange-correlation functional
    params_GPAW['occupations'] = FermiDirac(width = 0.1,      #The smearing
                                            fixmagmom = True)  #Total magnetic moment fixed to the initial value
    params_GPAW['convergence'] = {'eigenstates': 1.0e-4,       #eV^2 / electron
                                  'energy':      2.0e-4,       #eV / electron
                                  'density':     1.0e-3,}
    params_GPAW['mixer']       = Mixer(0.1, 5, weight=100.0)   #The mixer used during SCF optimization
    params_GPAW['setups']      = {'Fe': ':d,4.3'}              #U=4.3 applied to d orbitals

    calc = GPAW(**params_GPAW)
    lifepo4_vac.calc = calc
    epot_lifepo4_vac_cell=lifepo4_vac.get_potential_energy()
    print('E_Pot=', epot_lifepo4_vac_cell)

    write('lifepo4_vac_out.traj', lifepo4_vac)

    ens = BEEFEnsemble(calc)
    dE = ens.get_ensemble_energies(2000)
    result = paropen('ensemble_lifepo4_vac.dat','a')
    for i in range(0,len(dE)):
        print(dE[i], file=result)
    result.close()


Once you have made sure the cell runs, submit it to the HPC cluster.


.. code::

    # magic: !mq submit lifepo4_vac.py -R 8:1h  # submits the calculation to 8 cores, 1 hour


Once the calculation has finished, load in the trajectory.


.. code::

    try:
        lifepo4_vac=read('lifepo4_vac_out.traj')
        print('Calculation finished')
    except FileNotFoundError:
        print('Calculation has not yet finished')


Once the calculation has finished you are ready to calculate the energy cost of creating a li vacancy in the fully lithiated LiFePO$_4$. Start by loading in the relevant reference structures and obtain the potential energies. This should not require any calculations.


.. code::

    # Loading in files from exercise day 3.
    li_metal = read('li_metal.traj')   # you should have already read this in above
    lifepo4 = read('lifepo4_out.traj')

    epot_li_metal = li_metal.get_potential_energy() / len(li_metal)


.. code::

    # epot_lifepo4 = ...
    # ...

    # teacher
    epot_lifepo4=lifepo4.get_potential_energy()
    epot_lifepo4_vac=lifepo4_vac.get_potential_energy()


.. code::

    # vac_cost = ...
    # print(vac_cost)

    # teacher
    vac_cost=epot_lifepo4_vac-epot_lifepo4+epot_li_metal
    print(vac_cost)


How does this energy compare with the equilibirum potential? What can it tell you about the charge/discharge potential curves?



Bonus
=====
Calculate the error estimates of the energy for the added Li atom and vacancy formation using the ensembles.


.. code::

    # Cell for bonus question


.. code::

    # Cell for bonus question


.. code::

    # Cell for bonus question
