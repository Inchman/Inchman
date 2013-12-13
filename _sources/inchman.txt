.. _section-inchman:

Running Inchman
===============

Invoke Inchman from the Editor
------------------------------

Inchman can be invoked on the local computer directly from the :ref:`section-editor` by pressing the |simulate| button. The editor then creates the :ref:`annotated SBML file <section-sbml>` corresponding to the model and sends a link back to the browser using the ``inchman://`` protocol.

.. note ::

   When the user submits the experiment by pressing the simulate button, the editor internally creates an SBML file and sends a link to this file back to the browser. The link is of the form ::

     inchman://simulate/<global id>

   and uses the specific protocol ``inchman://`` for communication between the browser and the editor. If you are using Inchman for the first time, or if you are using a different browser, you will get a message which asks you to specify the program for this protocol. You will then have to point the browser to the location of the :ref:`Inchman Python wrapper script <inchman-wrapper>` (for example ``/usr/local/bin/inchman``). This choice can be changed at any time by going to Preferences->Applications in your browser and look for the Inchman protocol. 

After invokation, a window pops up asking the user to specify the working directory. This is the directory where all :ref:`section-output-files` will be stored. If the user previously chose to perform parameter sweeps (see :ref:`Section  parameters <parameters>`), Inchman will perform the required number of simulations to cover all combinations of parameter values. These different simulations will be saved in different output files. In contrast, if the user runs a number of experiments for each parameter set (see :ref:`Section options <options>`), the results of all these experiments will be stored in the same file.

.. warning ::

   Always check the output log files for errors. Not all errors are caught by the minimal GUI. In particular, any errors which involve the custom :ref:`Diffusivity-Drift-method <drift_diffusivity_method>` will only reveal themselves at compile time.

Invoke Inchman from the Command Line
------------------------------------

Inchman can also be invoked from the command line. This can be useful, for example, if Inchman is run on a GPU-enabled cluster or if the user opts to generate the :ref:`annotated SBML files <section-sbml>` from a script. For this purpose, two executables are available to the user.

Inchman Python Wrapper
""""""""""""""""""""""
.. _inchman-wrapper:

This Python wrapper script is known by the name of ``inchman`` and, in a unix environment, by default resides in ``/usr/local/bin/``. This script is called by the editor if the simulation button is pressed by the user. On the command line, it can be invoked by::

  $ inchman [options] <SBML-FILE>|<GLOBAL_ID>

The wrapper looks for any parameter sweep options in the <SBML-file> and computes the required parameter combinations. For each combination, it then invokes the :ref:`Inchman binary <inchman-binary>`. The wrapper script can also be called with a reference to to an editor model of the form ``inchman://simulate/<global id>``. In this case, the wrapper attempts to fetch the model with the given id from the server first.

.. program:: inchman

.. option:: --no-display

   Specifically asks Inchman not to show the progress window. If this option is selected, a working directory must be given explicitly.

.. option:: --datadir

   This will be the working directory.

.. option:: --pbs

   This option is useful, if Inchman is run on cluster that supports the `PBS job scheduling system <http://en.wikipedia.org/wiki/Portable_Batch_System>`_. If this option is checked, Inchman will create a PBS script file for each experiment and submit it to the queue. Since the syntax of these files is very platform dependent, the user will most likely have to adapt the template file. The template can be found in the Python package ``gpgmp.common.pbs``.

.. option:: --pbs-walltime=<walltime>

   Can be used in conjunction with the ``--pbs`` option and allows the user to specify a maximum wall time for the PBS run.


Inchman Binary
""""""""""""""

.. _inchman-binary:

The binary file ``inchman-exec`` is the actual simulator. It can be invoked by::

  $ inchman-exec [options] <model-file>

and takes a <model-file> as an argument. Note that the binary cannot perform any parameter sweeps.

.. option:: --output=<output file>

   Tells Inchman to save the :ref:`simulation output <section-output-files>` to the file <output file>

.. option:: --p <name>=<value>

   Sets the parameter <name> to value <value>. This is used by the wrapper script to do parameter sweeps.

.. option:: --cl-source=<cl-directory>

   This option can be used to tell Inchman where to find the Open-CL template files.

.. |simulate| image:: images/simulate_button.png
