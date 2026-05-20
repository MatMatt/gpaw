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


Part 1: Critical temperature of :mol:`CrI_3`
============================================

In the turorial :ref:`magnetism1` you will set up and relax a monolayer of 
:mol:`CrI_3` after which you will calculate its exchange contants and critical 
temperature using different models

The procedure will be as follows:

1) Set up the atomic structure and optimize the geometry of :mol:`CrI_3`
2) Calculate the nearest neighbor Heisenberg exchange coupling based on a total
   energy mapping analysis
3) Show that the magnetic ground state is thermodynamically unstable when
   anisotropy is neglected (The Mermin-Wagner theorem)
4) Calculate the single-ion magnetic anisotropy and estimate the critical
   temperature


Part 2: Noncollinear magnetism in :mol:`CrI_3`
==============================================

In materials where the dominant magnetic exchange coupling is antiferromagnetic,
or in cases where different exchange couplings compete, the ground state may
have a complicated noncollinear magnetic order. In the tutorial :ref:`magnetism2`
you will examine a prototypical monolayer with a noncollinear ground state,
namely :mol:`VI_2`. Starting from the structure file :download:``VI2.xyz``,
you will:

1) Relax the atomic structure using LDA
2) Compare a collinear antiferromagnetic structure with the ferromagnetic state
3) Obtain the noncollinear ground state
4) Calculate the magnetic anisotropy and discuss whether or not the material
   will exhibit magnetic order at low temperatures


Part 3: Find new ferromagnetic monolayers with high critical temperatures
=========================================================================

In this last part of the project, you will try to find new ferromagnetic
monolayers that can preserve their magnetic ordering at elevated temperatures.
The procedure is outlined in the tutorial :ref:`magnetism3`. In particular you
will

1) Search through a database of monolayers to pick a material you might expect
   to have a high critical temperature
2) Carry out a total energy mapping analysis to obtain exchange coupling and
   anisotropy parameters
3) Calculate a first principles estimate of the critical temperature

You are welcome to repeat this procedure for as many monolayers as you like.
