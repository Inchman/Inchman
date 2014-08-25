# This python code snippet serves to initialize the species array
import numpy
import sys
import math

# print out information about the model
print 'Init script for Slit (after Cai et al., 2006) model called!'
print 'Runtime information is: ', runtimeInformation
print 'Optional parameters are: ', parameters
print "Species list is:", species
sys.stdout.flush()

# extract some parameters we need
nx = runtimeInformation['nx']
ny = runtimeInformation['ny']
length = runtimeInformation['length']
nSpecies = runtimeInformation['nspecies']

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

if (False):
    # this is for the test problems

    # set background
    np_state_new[species['Cells'],:,:] = 0

    # set source 
    # note that index1 is y and index2 is x
    h = length/float(nx)
    halfwidth = int(5./h)
    print "h={0:g}, dh={1:g}.".format(h, halfwidth)  
    np_state_new[species['Cells'],:,nx/2-halfwidth:nx/2+halfwidth] = int(2000.*h)

if (True):
    # set background
    ntemp = np_state_new[species['Cells'],:,:]
    ntemp[:,:] = 0
    
    # concentration (from Cai)
    hcai = 5. # 5 mum is the cell size
    h = length/float(nx)
    u0 = 14./hcai**2
    
    # compute coordinate system
    # create coordinate system
    xc=numpy.mgrid[-nx/2:nx/2]/float(nx)*length
    yc=numpy.mgrid[-nx/2:nx/2]/float(nx)*length
    gxc=numpy.outer(xc, numpy.ones_like(xc))
    gyc=numpy.outer(numpy.ones_like(yc), yc)
    rs = gxc**2+gyc**2
    
    # and set explant value in the middle
    r0 = 100.
    ntemp[numpy.where(rs<=r0**2)] = u0*h**2
