# This python code snippet serves to initialize the species array
import numpy
import sys
import numpy.random

# print out information about the model
print 'Init script for majority vote model called!'
print 'Runtime information is: ', runtimeInformation
print 'Optional parameters are: ', parameters
print "Species list is:", species
sys.stdout.flush()

# extract some parameters we need
nx = runtimeInformation['nx']
ny = runtimeInformation['ny']
nSpecies = runtimeInformation['nspecies']
p = parameters['p']

# make a numpy array out of it
if (runtimeInformation['solver']=='stochastic_inhomogeneous'):
    # we have a stochastic simulation going
    print "Detected stochastic inhomogeneous simulation!"
    numpy_state = numpy.array(state, copy=False)
elif (runtimeInformation['solver']=='stochastic_homogeneous'):
    # we have a stochastic simulation going
    print "Detected stochastic homogeneous simulation!"
    numpy_state = numpy.array(state, copy=False)
else:
    # this is for the deterministic simulation
    print "Detected deterministic simulation!"
    numpy_state = numpy.array(deterministicState, copy=False)

# reshape the array
np_state_new = numpy_state.reshape(nSpecies, nx, ny)
np_state_new[:,:,:]=0

# distribute particles
nmols = int(parameters['numMolecules'])

print "nmols:{0:d}, nx:{1:d}, ny:{2:d} p:{3:g}.".format(nmols, nx, ny, p)

# create 2D random array
ngreens = numpy.random.binomial(nmols, p, (nx, ny))
nreds = numpy.ones((nx, ny))*nmols-ngreens

# and distribute them
if ('R00' in species) and ('G00' in species):
    print "Found 'R00' and 'G00' in species list."
    np_state_new[species['R00'], :, :] = nreds
    np_state_new[species['G00'], :, :] = ngreens
elif ('R0' in species) and ('G0' in species):
    print "Found 'R0' and 'G0' in species list."
    np_state_new[species['R0'], :, :] = nreds
    np_state_new[species['G0'], :, :] = ngreens
else:
    print "ERROR: Did not find suitable species in species list."
