.. _section-installation:

Installation
============

Linux
-----

Prerequisites
"""""""""""""

Inchman for linux is provided as a source-code distribution which requires the user to build and install the package. The following prerequisites need to be met:

* `CMake <http://www.cmake.org/>`_ (Version 2.8.0+). CMake is required to create the platform-specific build script and is available as a package for all major distributions. On Ubuntu, type ::
    
    sudo apt-get install cmake

* `Open-CL <http://www.khronos.org/opencl/>`_ developer drivers and headers for your hardware. Unfortunately, the exact procedure to install these drivers still strongly depends on your hardware and your linux distribution. For NVIDIA cards on Ubuntu, installing the package "nvidia-opencl-dev" ::

    sudo apt-get install nvidia-opencl-dev

  should do the trick. More general information about setting up Open-CL for Intel cards on Ubuntu 13.04. can be found `here <http://mhr3.blogspot.com.au/2013/06/opencl-on-ubuntu-1304.html>`_. 

* `Python <http://www.python.org/>`_ (Version 2.7+ or Version 3.3+). We need the interpreter and the developer headers. For Ubuntu users ::

    sudo apt-get install python-dev

* `Numpy <http://www.numpy.org/>`_ (Version 1.7.1+) ::

    sudo apt-get install python-numpy

* `HDF5 <http://www.hdfgroup.org/HDF5/>`_ (Version 1.8.10+). Again, we need the library and the headers. ::
    
    sudo apt-get install libhdf5-dev

* `H5PY <http://www.h5py.org/>`_ (Version 2.1.3+). You only need H5PY if you would like to use Python for analyzing the Inchman output (recommended). Make sure that you H5PY version matches your HDF5 version - the latest H5PY version (2.2 at the time of writing) requires HDF5 1.8.4+ ::

    sudo apt-get install python-h5py

* `libSBML5 <http://sbml.org/Software/libSBML>`_ library and headers (Version 5.8.0+) ::

    sudo apt-get install libsbml5-dev

* `Boost <http://www.boost.org/>`_ (Version 1.54+) library and headers. We need the libraries "log", "chrono", "date_time", "filesystem", "log_setup", "system", "thread", "program_options", "python". ::

    sudo apt-get install libboost-all-dev

  Note that the default boost library in Ubuntu Saucy Salamander is still 1.53, which does not contain the program_options library. Until this has been updated, you will need to install it by hand.

Build and install
"""""""""""""""""

Once all pre-requisites are met, you can build and install Inchman.

1. Download the latest version from http://code.google.com/p/gpgmp/downloads/list

2. Open a shell and unpack it (assuming that the package is in the current directory)::

     $tar xzf inchman-2.0.0.tar.gz

3. Change into the new directory, create a build directory and run CMake in it::

     $cd inchman-2.0.0/inchman
     $mkdir build
     $cd build
     $cmake ..

4. Build the code and install it:::

     $make
     $sudo make install

5. You finally need to install the Python libraries. **The :ref:`Inchman wrapper script <inchman-wrapper>` will not work without access to the Python libraries!**. If you have root access to your machine::

     $cd ../python
     $sudo python setup.py install

   If you do not have root access to your machine, you can copy the ``python`` directory to any location and set the ``PYTHONPATH`` environment variable (see the `python documentation <http://docs.python.org/3/using/cmdline.html#envvar-PYTHONPATH>`_).

Mac OS X
--------

1. Download the disk image from http://code.google.com/p/gpgmp/downloads/detail?name=Inchman.app.dmg.
2. Open the image in Finder. You can then drag Inchman to your Application folder or start it by double-clicking.
3. The Inchman application will then listen to incoming simulation requests from the :ref:`editor <section-editor>`.
4. You might want to install the Python libraries. These are not required for the Mac version but might be helpful for analysis. The Python modules are available at http://code.google.com/p/gpgmp/downloads/detail?name=python.tar.gz. Download them, open a terminal window and unpack them. If you have root access to your machine::

     $cd python
     $sudo python setup.py install

   If you do not have root access to your machine, you can copy the ``python`` directory to any location and set the ``PYTHONPATH`` environment variable (see the `python documentation <http://docs.python.org/3/using/cmdline.html#envvar-PYTHONPATH>`_).
