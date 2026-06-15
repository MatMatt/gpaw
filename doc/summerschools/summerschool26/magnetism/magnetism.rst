.. _magnetism:

===============
Magnetism in 2D
===============

This exercise investigates magnetic order in two dimensions. While a collinear
spin-density functional theory calculation might reveal that it is
energetically favorable for a given 2D material to order magnetically, two
dimensional magnetic order at finite temperatures also requires presence of
magnetic anisotropy, usually arising from the spin-orbit coupling.

This exercise will teach you how to extract magnetic exchange and anisotropy
parameters for a localized spin model based on first principles calculations.
It will also touch upon the Mermin-Wagner theorem and show why anisotropy is
crucial in order to sustain magnetic order in two dimensions.

In the first part of the project, you will calculate the Curie temperature of
a :mol:`CrI_3` monolayer. Afterwards, you will investigate the magnetic order in
:mol:`VI_2`, which has antiferromagnetic coupling and noncollinear order. Finally,
you will finish the project by performing a search for new magnetic 2D materials
with high critical temperatures based on a database of hypothetical monolayers.

.. toctree::
   
   magnetism1
   magnetism2

=========================================================================
Part 3: Find new ferromagnetic monolayers with high critical temperatures
=========================================================================

Now that you know how to calculate exchange coupling constants and single-ion anisotropies, you are in a position to search for new ferromagnetic monolayers yourself.

First, we will try to find some suitable candidates based on the data, which is already in the  Computational 2D Materials Database (C2DB).

1.   Go to the C2DB website at https://cmrdb.fysik.dtu.dk/c2db and search for magnetic monolayers with finite band gaps (the theory developed so far is not applicable to metals, since these typically have long range interactions). (Hint: Toggle `Magnetic` to `yes` and set a minimum value for the `Band gap range` in the search menu, then press the search icon.)
2.   Refine the query to sort out the magnetic monolayers with a ferromagnetic ground state. (Hint: Key values such as the nearest neighbor exchange coupling can be queried by entering `J>0` into the search field.)
3.   Show the spin state, exchange coupling and magnetic anisotropy in the property overview. (Hint: You can add key values to the overview through the `Add Column` botton.)
4.   Decide on a material you find promising to have a high Curie temperature, based on what you learned in previous notebooks. You may also want to choose a material with a limited number of atoms in the unit cell in order for the calculations to run fast.

You are welcome to repeat the process for as many monolayers as you like.
