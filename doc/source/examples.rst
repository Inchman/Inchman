Examples
========

.. _example-ab-reaction:

Simple Reaction Network
-----------------------

A simple network consisting of two species and four reactions. Details can be found in [Vigelius2010]_. The source code is available in :download:`ABReaction.xml <../examples/ABReaction.xml>`. A standard analysis can be done by::

   $python -m gpgmp.test.abreaction <output-file>

.. _example-ab-reaction-individual:

Simple Reaction Network with Individuals
----------------------------------------

The same reaction network as in :ref:`example-ab-reaction` is modelled with individual species. The source code is available in :download:`ABReactionIndividual.xml <../examples/ABReaction.xml>`. A standard analysis can be done by::

   $python -m gpgmp.test.abreaction <output-file> individual

.. _example-ab-reaction-individual-new:

Simple Reaction Network with Individuals (using New Individuals Method)
-----------------------------------------------------------------------

The same reaction network as in :ref:`example-ab-reaction` is modelled with individual species, where all newly created individuals are initialized with a random diffusivity and drift. The source code is available in :download:`ABReactionIndividualNew.xml <../examples/ABReaction.xml>`. A standard analysis can be done by::

   $python -m gpgmp.test.abreaction <output-file> individual

In addition, you can use any HDF5 viewer to check, if the individual properties were set correctly, i.e. uniformly-random distributed.

.. _example-fisher:

Fisher Problem
--------------

The famous Fisher problem, a non-linear reaction-diffusion problem  [Vigelius2010]_. The source code is in :download:`FisherProblem.xml <../examples/FisherProblem.xml>` and standard analysis can be done by::

   $python -m gpgmp.test.fisher_problem <output-file>

.. _example-fisher-individual:

Fisher Problem with Individuals
-------------------------------

The same as :ref:`example-fisher` except that the species are individuals. The source code is in :download:`FisherProblemIndividual.xml <../examples/FisherProblemIndividual.xml>` and standard analysis can be done by::

   $python -m gpgmp.test.fisher_problem <output-file> individual


.. _example-diffusion:

Homogeneous Diffusion
---------------------

A simple example to test the homogeneous solver involving one diffusing species (no drift) and no reactions [Vigelius2012a]_. The source code is in :download:`HomogeneousDiffusion.xml <../examples/HomogeneousDiffusion.xml>` and standard analysis can be done by::

   $python -m gpgmp.test.homogeneous_diffusion <output-file>

.. _example-drift:

Homogeneous Diffusivity and Drift
---------------------------------

The same example as :ref:`example-diffusion` with a drift field [Vigelius2012a]_. The source code is in :download:`HomogeneousDrift.xml <../examples/HomogeneousDrift.xml>` and standard analysis can be done by::

   $python -m gpgmp.test.homogeneous_drift_diffusion <output-file>

.. _example-drift-field:

Homogeneous Diffusivity and Drift with Field Parameter
------------------------------------------------------

Exactly the same as :ref:`example-drift` only that the diffusivity and drift field is now set using a field parameter. The source code is in :download:`HomogeneousDriftField.xml <../examples/HomogeneousDriftField.xml>` and standard analysis can be done by::

   $python -m gpgmp.test.homogeneous_drift_diffusion <output-file>

.. _example-drift-individual:

Homogeneous Diffusivity and Drift (Individuals)
------------------------------------------------------

Exactly the same as :ref:`example-drift` only that the diffusing species are individuals. The source code is in :download:`HomogeneousDriftIndividual.xml <../examples/HomogeneousDriftIndividual.xml>` and standard analysis can be done by::

   $python -m gpgmp.test.homogeneous_drift_diffusion <output-file>

.. _example-multiplicative-noise:

Multiplicative Noise (Geometric Brownian Motion)
------------------------------------------------

An example to test the inhomogeneous solver. It contains one species and position-dependent diffusivity and drift [Vigelius2012a]_. The source code is in :download:`MultiplicativeNoise.xml <../examples/MultiplicativeNoise.xml>` and standard analysis can be done by::

   $python -m gpgmp.test.multiplicative_noise <output-file>

.. _example-nonlinear:

Non-linear Drift Field
----------------------

A test problem involving a non-linear drift field [Vigelius2012a]_. The source code is in :download:`Nonlinear.xml <../examples/Nonlinear.xml>` and standard analysis can be done by::

   $python -m gpgmp.test.nonlinear <output-file>


.. _example-ornstein-uhlenbeck:

Ornstein-Uhlenbeck Process
--------------------------

An implementation of an Ornstein-Uhlenbeck process [Vigelius2012a]_. The source code is in :download:`OrnsteinUhlenbeck.xml <../examples/OrnsteinUhlenbeck.xml>` and standard analysis can be done by::

   $python -m gpgmp.test.ornstein_uhlenbeck <output-file>

.. _example-ab:

A+B Annihilation
----------------

This is a simple system with two species :math:`A` and :math:`B` and a corresponding reaction :math:`A+B\xrightarrow{k_1}\emptyset`. Details can be found in [Vigelius2010]_. You can analyze it using::

   $python -m gpgmp.test.annihilation_2d <output-file>

.. _example-ab-drift:

A+B Annihilation with Drift
---------------------------

This example is the same as example :ref:`example-ab` except that a constant drift field is present. The source code is found in :download:`ABAnnihilationDrift.xml <../examples/ABAnnihilationDrift.xml>`. Analysis can be done using::

   $python -m gpgmp.test.annihilation_2d_drift <output-file>

.. _example-random-drift:

Random Drift Field
------------------

In this example, a number of individuals is initialized with random diffusivity and drift, where each individual has a different diffusivity and drift. The source code is found in :download:`RandomDrift.xml <../examples/RandomDrift.xml>`. Analysis can be done using::

   $python -m gpgmp.test.random_drift <output-file>


.. _example-calcium:

Intracellular Calcium Distribution
----------------------------------

This example implements an actual biological model used to describe the intracellular distribution of Calcium ions (cf. [Vigelius2012b]_ for details). The model is found in :download:`Calcium.xml <../examples/Calcium.xml>` and can be analysed using::

   $python -m python -m gpgmp.models.calcium <output-file>

.. _example-oregonator:

Oregonator Model of the Belousov-Zhabotinsky Reaction
-----------------------------------------------------

This is the Oregonator model of the Belousov-Zhabotinsky reaction which is covered in the :ref:`tutorial <tutorial_bz_reaction>`. The model is found in :download:`Oregonator.xml <../examples/Oregonator.xml>`. Analysis is covered in the tutorial.

.. _example-slit:

Chemotaxis of Neurons in the Brain
----------------------------------

The model of migrating neurons in the brain which is also covered in the :ref:`tutorial <tutorial_slit>`. The model is found in :download:`Slit.xml <../examples/Slit.xml>`. Analysis is covered in the tutorial.

.. [Vigelius2010] Vigelius M, Lane A, Meyer B (2010):  `Accelerating reaction–diffusion simulations with general-purpose graphics processing units <http://bioinformatics.oxfordjournals.org/content/27/2/288>`_  Bioinformatics (2011) 27 (2): 288-290. 

.. [Vigelius2012a] Vigelius M, Meyer B (2012a) `Multi-Dimensional, Mesoscopic Monte Carlo Simulations of Inhomogeneous Reaction-Drift-Diffusion Systems on Graphics-Processing Units. <http://www.plosone.org/article/info%3Adoi%2F10.1371%2Fjournal.pone.0033384>`_ PLoS ONE 7(4): e33384. doi:10.1371/journal.pone.0033384

.. [Vigelius2012b]  Vigelius M, Meyer B (2012b): `Stochastic Simulations of Pattern Formation in Excitable Media. <http://www.plosone.org/article/info%3Adoi%2F10.1371%2Fjournal.pone.0042508>`_ PLoS ONE 7(8): e42508. doi:10.1371/journal.pone.0042508
