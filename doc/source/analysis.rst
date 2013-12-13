Analyzing the simulation output
===============================

.. _section-output-files:

Inchman output files
--------------------

All output files are stored in the current working directory (see :ref:`section-inchman`). If Inchman is invoked from the :ref:`section-editor` or using the :ref:`wrapper script <inchman-wrapper>`, the following files will be created in the working directory:

* ``stdout.txt.<job-id>`` The standard log file for the given job. The first place to look for errors.
* ``stderr.txt.<job-id>`` The standard error file for the given job.
* ``job_<job-id>.h5``     The :ref:`HDF5 File <section-hdf5>`.
* ``final_code_Gillespie.cl`` The Open-CL code for the kernels used by the Gillespie solver (only present if the model contains reactions).
* ``final_code_StochasticInhomogeneousDriftDiffusion.cl`` The Open-CL code for the inhomogeneous diffusion kernels (only present if the inhomogeneous solver is selected)
* ``final_code_StochasticHomogeneousDiffusion.cl`` The Open-CL code for the homogeneous diffusion kernels (only present if the homogeneous solver is selected)
* ``final_code_Deterministic.cl`` The Open-CL code for the kernels used by the deterministic solver (only present if the deterministic solver is selected)

The HDF5 File
-------------

The simulation output is stored in the corresponding HDF5 file (see :ref:`section-HDF5`). The files can be read out by the standard HDF5 tools. Inchman also comes with a variety of Python helpers that can be used to extract all information from the Inchman output file and analyse it using Numpy. 

Using Python for Analysis
^^^^^^^^^^^^^^^^^^^^^^^^^

Python in combination with Numpy is the recommended way to analyse Inchman output. All input-output functionality is provided by the model ``gpgmp.io`` which can be loaded by ::

  import gpgmp.io

Reading in the state arrays
"""""""""""""""""""""""""""

The state array can be read in using the function ``gpgmp.io.read_gmp_hdf5(file, deterministic=False)``::

  >>> n, times, species, nruns = gpgmp.io.read_gmp_hdf5('job_0', deterministic=False)

If the deterministic solver was used, you have to set the ``deterministic`` keyword to ``True``. In this case, the state array will be of type ``FLOAT`` instead of ``INT``. The variable ``nruns`` now contains the number of runs. The array ``species`` gives the species that were defined in the model:

  >>> species
  array(['A'], 
        dtype='<U1')

The array ``times`` gives the time stamps of the different output dumps:

  >>> times
  array([ 0.,  1.,  2.,  3.,  4.,  5.])

.. warning::

   Note that the output dumps are sorted lexically, i.e. the output at time ``10`` comes before the dump at time ``2``. That means that you might have to order the time stamps numerically first, for example through

     >>> tind = numpy.argsort(times)

Finally, the five-dimensional state array ``n`` gives the number count for each run (first index), each output number (second index), each species (third index) at each subvolume (fourth and fifth indices):

   >>> import numpy
   >>> numpy.shape(n)
   (10, 6, 1, 512, 512)

For example, to display species 'A' for run number '3' at time '4' you can use:

   >>> import matplotlib.pyplot
   >>> matplotlib.pyplot.imshow(n[3,4,0,:,:])
   <matplotlib.image.AxesImage object at 0x7fde0e424590>
   >>> matplotlib.pyplot.show()

Reading in the field arrays
"""""""""""""""""""""""""""

If the model contains :ref:`field parameters <field_parameter>`, the fields can be read using ``read_field(filename, fieldname)``:

  >>> tracefield = gpgmp.io.read_field('/home/matthias/xibalba/temp/results/job_0.h5', 'traceField')
  >>> numpy.shape(tracefield)
  (10, 6, 512, 512)

The first index denotes the run, the second index denotes the time stamp and the last two indices give the subvolume.

.. _analysis-individuals:

Simulations involving Individuals
"""""""""""""""""""""""""""""""""

If the simulations involves :ref:`individuals`, the position and ID of every individual can be read out using the function ``gpgmp.io.read_gmp_individuals_all(filename, species)``:

   >>> n = gpgmp.io.read_gmp_individuals_all('/home/matthias/xibalba/temp/results/job_0', 'A')

The result is a dictionary, where the keys denote the output dump number:

   >>> n.keys()
   dict_keys(['Dump_1', 'Dump_0', 'Dump_3', 'Dump_2', 'Dump_5', 'Dump_4'])

The elements are 64 bit unsigned integers which contain the cell position in the lower 32 bits and the individual id in the higher 32 bits:

   >>> n['Dump_1']
   array([108512348822306, 245749438833965, 413115724425990, ...,
          346741299916537, 179899000336154, 231030586001694], dtype=uint64)

The module ``gpgmp.individuals`` provides a function to separate the position/id pairs:

   >>> import gpgmp.individuals
   >>> pos, ids = gpgmp.individuals.split_keys(n['Dump_1'])
   >>> pos
   array([ 88866,  91437,  92934, ..., 175865, 175898, 182558], dtype=uint64)

.. note::

   If the user is only interested in the number of individuals per subvolume, she can read in the state array just normally as described above.

