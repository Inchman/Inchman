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
omega = parameters['omega']
k = parameters['k']
F = parameters['F']
u0 = parameters['u0']
k1 = parameters['k1']
print "Parameter extraction successful!"
sys.stdout.flush()
k2 = F + k
kf = F

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

if (True):
    # set stable homogeneous background (assume A=2)
    np_state_new[species['U'],:,:] = 2.*k2*omega/math.sqrt(k1*kf)
    np_state_new[species['V'],:,:] = 0.

    # now set initial region to start wave
    np_state_new[species['U'], dx_/2-5:dx_/2, 0:dx_/4] = 0.
    np_state_new[species['V'], dx_/2-5:dx_/2, 0:dx_/4] = 4.*k2*omega/math.sqrt(k1*kf)

    # and refractory region
    np_state_new[species['U'], dx_/2:, 0:dx_/4] = k2*omega/math.sqrt(k1*kf)/2.
    np_state_new[species['V'], dx_/2:, 0:dx_/4] = kf*omega/math.sqrt(k1*kf)/2.

if (False):
    np_state_new[species['U'],:,:] = 250.
    np_state_new[species['V'],:,:] = 0.

    np_state_new[species['U'],dx_/2-5:dx_/2, 0:dx_/4] = 0.
    np_state_new[species['V'],dx_/2-5:dx_/2, 0:dx_/4] = 500.

    #refractory
    np_state_new[species['U'],dx_/2:, :] = 75.
    np_state_new[species['V'],dx_/2:, :] = 2.
