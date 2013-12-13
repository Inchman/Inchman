.. _section-hdf5:

Structure of the HDF5 File
==========================

This section describes the structure of the HDF5 output file. HDF5 files are structured hierarchically in groups, datasets and attributes.

Model Group
-----------

This group has the following attributes:

* ``/Model/numRuns`` The number of runs.
* ``/Model/numSpecies`` The number of species in the model.

The datasets are:

* ``/Model/InchmanFile (H5T_STRING)`` The model :ref:`SBML file <section-sbml>`.
* ``/Model/Species (H5T_STRING[..])`` List of species in the model.

Compartments Group
""""""""""""""""""

This group consists of one or more dataspaces which list the compartments in the model:

* ``/Model/Compartments/World (H5T_STD_I32LE[4])`` Array listing the position and size of the obligatory ``World`` compartment.


Parameters Group
----------------

This group lists all parameters if any, as attributes. For example:

* ``/Model/Parameters/k2 (H5T_IEEE_F32LE[1])`` Contains the value of the model parameter ``k2``.

.. note ::

  If you set up :ref:`parameter sweeping <parameters>` experiments, each parameter set will correspond to exactly one HDF5 output file. That means that the parameter values are always scalars and do not contain any information about the loop range.


Runs Group
----------

This group lists all runs, starting from run "0". 

.. note ::

  The runs are sorted lexically, i.e. run number "10" will be listed before run number "2".

Dumps Group
"""""""""""

Each run contains a number of output dumps, which are again sorted lexically.

* ``/Runs/<run-no>/Dump_<dump-no>/time (H5T_IEEE_F32LE[1])`` This attribute gives the simulation time corresponding to this dump.

* ``/Runs/<run-no>/Dump_<dump-no>/<species-name> (H5T_STD_I32LE[<width>,<height>])`` Eeach species is represented by a dataset with dimensions ``GridWidth x GridWidth``. The integers in this dataset are the number of particles in the corresponding grid cell.

Individuals Group
^^^^^^^^^^^^^^^^^

If the model contains any individuals, they are listed in this group. The individuals are stored as 64 bit unsigned integers representing the individual ID and the cell positon (see :ref:`analysis-individuals`).


