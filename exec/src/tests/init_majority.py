# This python code snippet serves to initialize the species array
import numpy
import sys
import random

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
random.seed()
nmols = int(parameters['numMolecules'])

print "nmols:{0:d}, nx:{1:d}, ny:{2:d}.".format(nmols, nx, ny)
nreds = 0
ngreens=0

for i in range(nmols):
    # get position
    posx=random.randint(0, nx-1)
    posy=random.randint(0, ny-1)

    # red or green?
    if (random.random()<0.5):
        nreds = nreds + 1
        # red
        np_state_new[species['R00'], posx, posy] = np_state_new[species['R00'], posx, posy] + 1
        print "Adding red particle at {0:d},{1:d}.".format(posx,posy)
    else:
        ngreens=ngreens+1
        np_state_new[species['G00'], posx, posy] = np_state_new[species['G00'], posx, posy] + 1
        print "Adding green particle at {0:d},{1:d}.".format(posx,posy)

print "Distributed {0:d} red and {1:d} green molecules (total={2:d})".format(nreds, ngreens, nreds+ngreens)
