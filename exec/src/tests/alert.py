# This python code snippet serves to initialize the species array
import numpy

print 'Alert called at time={0:g}'.format(time_)
print 'Grid topology is (dx,dy)=({0:d}, {1:d}). #species = {2:d}.'.format(dx_, dy_, nSpecies_)
print species

# make a numpy array out of it
numpy_state = numpy.array(state, copy=False)

# display shape
print numpy.shape(numpy_state)

# reshape it
np_state_new = numpy_state.reshape(nSpecies_, dx_, dy_)

print numpy.shape(np_state_new)

# trigger another wave this time coming from the bottom
#np_state_new[species['Z'],:, 0:2]=147

# and set it to steady-state solution
#scale = 4e-4

#np_state_new[:,:, 0:40] = 0
#np_state_new[:,:, 80:] = 0

#np_state_new[0, 0:40, : ] = 2.26364e11*scale**3
#np_state_new[1, 0:40, : ] = 4.57189e13*scale**3
#np_state_new[2, 0:40, : ] = 1.13182e13*scale**3
#np_state_new[3, 0:40, : ] = 0
#np_state_new[4, 0:40, : ] = 2.28039e11*scale**3

#np_state_new[0, 80:, : ] = 2.26364e11*scale**3
#np_state_new[1, 0:40, : ] = 4.57189e13*scale**3
#np_state_new[2, 0:40, : ] = 1.13182e13*scale**3
#np_state_new[3, 0:40, : ] = 0
#np_state_new[4, 0:40, : ] = 2.28039e11*scale**3


#print np_state_new[0, :, :]
#print np_state_new[1, :, :]
#print np_state_new[2, :, :]

# and change it 
#np_state_new[2, :, :] = 13
#np_state_new[2, 3, :] = 15

#print np_state_new[2, :, :]

# manipulate numpy state
#numpy_state[0] = 1.
#numpy_state[1] = 2.
#numpy_state[2] = 3.

#print "Numpy state:"
#print numpy_state

#print "State:"
#print state

#import _gpgmp
#val = _gpgmp.get_state()

#print val


#oarr = numpy.array([1,2,3])

#print "sqrt of array is:"
#print numpy.sqrt(arr)
