# This python code snippet serves to initialize the species array
import numpy
import sys
import traceback

try:
    print 'Init script for Fisher problem called!'
    print 'Runtime information is: ', runtimeInformation
    print 'Optional parameters are: ', parameters
    print "Species list is:", species
    sys.stdout.flush()

    # extract some parameters we need
    nx = runtimeInformation['nx']
    ny = runtimeInformation['ny']
    length = runtimeInformation['length']
    nSpecies = runtimeInformation['nspecies']

    particles = parameters['particles']

    print "Parameter extraction successful!"

    # make a numpy array out of it
    if (runtimeInformation['solver']=='stochastic_inhomogeneous') or (runtimeInformation['solver']=='stochastic_homogeneous'):
        # we have a stochastic simulation going
        print "Detected stochastic simulation!"
        numpy_state = numpy.array(state, copy=False)
    else:
        # this is for the deterministic simulation
        print "Detected deterministic simulation!"
        numpy_state = numpy.array(deterministicState, copy=False)

    # reshape the array
    np_state_new = numpy_state.reshape(nSpecies, nx, ny)

    # and set up the initial fisher solution
    dx = length/float(nx)
    xc = (numpy.mgrid[0:nx]+0.5)*dx-3.*length/4.
    t = (1.+numpy.exp(-xc/numpy.sqrt(6.)))

    # and duplicate rows
    tcc = numpy.array([t]*nx)
    np_state_new[0,:,:] = numpy.round(1./(tcc*tcc)*particles)

except:
    (etype, value, traceback) = sys.exc_info()
    print "Unexpected error:", etype
    print traceback.print_tb(traceback)
    raise
