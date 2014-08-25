# This python code snippet serves to initialize the species array
import numpy
import sys

print 'Init script for Grey-Scott model called!'
sys.stdout.flush()
print 'Grid topology is (dx,dy)=({0:d}, {1:d}). #species = {2:d}.'.format(dx_, dy_, nSpecies_)

#print state

# the species are available as a dictionary
print species

# make a numpy array out of it

#numpy_state = numpy.array(state, copy=False)
# -- this is for the deterministic simulation
numpy_state = numpy.array(deterministicState, copy=False)

# display shape
print numpy.shape(numpy_state)

# reshape it
np_state_new = numpy_state.reshape(nSpecies_, dx_, dy_)

print numpy.shape(np_state_new)

# zero it
np_state_new[:,:,:] = 0.
np_state_new[:,:,:] = 448.685

# convert to number counts per subvolume and set it
#scale = 5e-9
#length = 250.
#subvolume = (length*scale/float(dx_))**3*1e3;
#na = 6.02214e23;
#M = na*subvolume;

# number of channels per subvolume
#nchannelssv = 200.*(1./2.4*length/float(dx_))**2

#np_state_new[species['S'],:,:] = 2.25

#np_state_new[:,:,:] = 0.01*4.48685e9*scale
#np_state_new[species['c'],:,:] = 0.1*M
#np_state_new[species['y'],:,:] = 0.1*nchannelssv
