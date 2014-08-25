# This python code snippet serves to initialize the species array
import numpy
import sys
import math

# print out information about the model
print 'Init script for Gray-Scott (after Wang et al., 2007) model called!'
print 'Runtime information is: ', runtimeInformation
print 'Optional parameters are: ', parameters
print "Species list is:", species
sys.stdout.flush()

# extract some parameters we need
nx = runtimeInformation['nx']
ny = runtimeInformation['ny']
nSpecies = runtimeInformation['nspecies']

scale = parameters['omega']
k1 = parameters['k1']
k2 = parameters['k2']
k3 = parameters['k3']
k4 = parameters['k4']
k5 = parameters['k5']
k6 = parameters['k6']

# make a numpy array out of it
if (runtimeInformation['solver']=='stochastic_homogeneous'):
    # we have a stochastic simulation going
    print "Detected stochastic simulation!"
    numpy_state = numpy.array(state, copy=False)
else:
    # this is for the deterministic simulation
    print "Detected deterministic simulation!"
    numpy_state = numpy.array(deterministicState, copy=False)

# reshape the array
np_state_new = numpy_state.reshape(nSpecies_, dx_, dy_)

# zero the array
np_state_new[:,:,:] = 0

# We follow the notebook dimless_to_dimensional.nb
# (iii) Three-state dimensional with f = 3. (note that f is now different!!)
if (True):
    # parameters
    scalet = 3e-4
    k1t = 843.2
    k2t = 1.04472e-9/scale**3
    k3t = 4.216
    k4t = 5.2236e-15/scale**3
    k5t = 0.048
    k6t = 4.353e-6*1e7/scale**3
    
    print "Script finds (theo, real): k1=({0:g}, {1:g}), k2=({2:g}, {3:g}), k3=({4:g}, {5:g}), k4=({6:g}, {7:g}), k5=({8:g}, {9:g}), k6=({10:g}, {11:g})".format(k1t, k1, k2t, k2, k3t, k3, k4t, k4, k5t, k5, k6t, k6)

    # create coordinate system
    xc=range(-dx_/2, dy_/2)
    yc=range(-dx_/2, dy_/2)
    gxc=numpy.outer(xc, numpy.ones_like(xc))
    gyc=numpy.outer(numpy.ones_like(yc), yc)
    theta = numpy.arctan2(gyc, gxc) + numpy.pi

    # create float arrays (for stochastic sim)
    nsf = numpy.zeros(numpy.shape(np_state_new), dtype=float)

    # u = {0.8 (0<theta<0.5), uss (else)}
    q = 0.002;
    f = 1.5;
    uss = 0.5*(1.-(f+q)+numpy.sqrt((f+q-1.)**2+4.*q*(1.+f)))

    tindex = numpy.where(theta<0.5)
    nsf[species['X'], :, : ] = uss
    nsf[species['X'], tindex[0], tindex[1] ] = 0.8

    # v = vss + theta/(8 pi f)
    nsf[species['Z'], :, : ] = uss + theta/(8.*numpy.pi*f)

    # and for Y we take the dynamic equilibrium value
    nsf[species['Y'], :, :] = f * nsf[species['Z'], :, : ]/(q + nsf[species['X'], :, : ])

    # and convert all to dimensional values
    np_state_new[species['X'], :, : ] = nsf[species['X'], :, : ] / (2.*k4/k3)
    np_state_new[species['Y'], :, : ] = nsf[species['Y'], :, : ] * (k3/k2)
    np_state_new[species['Z'], :, : ] = nsf[species['Z'], :, : ] / (k5*k4/k3**2)


# (i) Stationary state/no diffusion for three-state dimensionless oregonator
if (False):
    np_state_new[species['X'], :, :] = 0.00977001
    np_state_new[species['Y'], :, :] = 1.24511
    np_state_new[species['Z'], :, :] = 0.00977001

# (ii) three-state dimensionless with/without diffusion and f = 1.5
if (False):
    # create coordinate system
    xc=range(-dx_/2, dy_/2)
    yc=range(-dx_/2, dy_/2)
    gxc=numpy.outer(xc, numpy.ones_like(xc))
    gyc=numpy.outer(numpy.ones_like(yc), yc)
    theta = numpy.arctan2(gyc, gxc) + numpy.pi

    # u = {0.8 (0<theta<0.5), uss (else)}
    q = 0.002;
    f = 1.5;
    uss = 0.5*(1.-(f+q)+numpy.sqrt((f+q-1.)**2+4.*q*(1.+f)))

    tindex = numpy.where(theta<0.5)
    np_state_new[species['X'], :, : ] = uss
    np_state_new[species['X'], tindex[0], tindex[1] ] = 0.8

    # v = vss + theta/(8 pi f)
    np_state_new[species['Z'], :, : ] = uss + theta/(8.*numpy.pi*f)

    # and for Y we take the dynamic equilibrium value
    np_state_new[species['Y'], :, :] = f * np_state_new[species['Z'], :, : ]/(q + np_state_new[species['X'], :, : ])

# dummy test problem to validate alpha-QSS
if (False):
    np_state_new[species['U'], :, :] = 0.
    np_state_new[species['V'], :, :] = 0.
