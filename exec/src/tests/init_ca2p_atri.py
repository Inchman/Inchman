# This python code snippet serves to initialize the species array
import numpy
import sys

# print out information about the model
print 'Init script for Ca2+-waves (after Atri, 1993) model called!'
print 'Runtime information is: ', runtimeInformation
print 'Optional parameters are: ', parameters
print "Species list is:", species
sys.stdout.flush()

# extract some parameters we need
nx = runtimeInformation['nx']
ny = runtimeInformation['ny']
nSpecies = runtimeInformation['nspecies']
scale = parameters['omega']
nmax = parameters['nmax']
muM = parameters['muM']

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

# convert to number counts per subvolume and set it
nM = 1e-3*muM

if (True):
    # this one works for spirals, yay!
    # set background
    np_state_new[species['Ca'],:,:] = 0.129*muM
    np_state_new[species['n'],:,:] = 0.96627*nmax
    np_state_new[species['p'],:,:] = 0.095*muM

    # now set initial region to start wave
    #np_state_new[species['Ca'], 0:10, :]=0.920*muM
    np_state_new[species['Ca'], dx_/2-5:dx_/2, 0:dx_/4]=0.920*muM
    #np_state_new[species['Ca'], dx_/2:, :] = 0.020*muM
    np_state_new[species['Ca'], dx_/2:, :] = 0.19*muM
    np_state_new[species['n'], dx_/2:, :] = 0.71*nmax
    
if (False):
    # set background
    np_state_new[species['Ca'],:,:] = 0.328*muM
    np_state_new[species['n'],:,:] = 0.820*nmax
    np_state_new[species['p'],:,:] = 0.095*muM

    # now set initial region to start wave
    np_state_new[species['Ca'],:,0:10] = 0.01*muM
    #np_state_new[species['Ca'], dx_/2-5:dx_/2, 0:dx_/2]=0.01*muM
    #np_state_new[species['Ca'], dx_/2:dx_/2+5, 0:dx_/2]=0.15*muM


if (False):
    # this one generates a circular wave by IP3 increase
    # set background
    np_state_new[species['Ca'],:,:] = 130.*nM
    np_state_new[species['n'],:,:] = 0.96585*nmax
    np_state_new[species['p'],:,:] = 0.

    # now set initial region to start wave
    np_state_new[species['p'], dx_/2-5:dx_/2+5, dx_/2-5:dx_/2+5]=2.5*muM
    #np_state_new[species['p'], :, 0:10]=2.5*muM
    
if (False):
    # this one generates a plane wave
    # set background
    np_state_new[species['Ca'],:,:] = 130.*nM
    np_state_new[species['n'],:,:] = 0.96585*nmax
    np_state_new[species['p'],:,:] = 0

    # now set initial region to start wave
    np_state_new[species['Ca'], :, 0:10]=0.9*muM

if (False):
    # this one generates a circular wave by Ca2+ increase
    # set background
    np_state_new[species['Ca'],:,:] = 130.*nM
    np_state_new[species['n'],:,:] = 0.96585*nmax
    np_state_new[species['p'],:,:] = 0

    # now set initial region to start wave
    np_state_new[species['Ca'], dx_/2-5:dx_/2+5, dx_/2-5:dx_/2+5]=0.9*muM

    # we try to get a spiral wave
    #np_state_new[species['Ca'], dx_/2-20:dx_/2, dx_/2-20:dx_/2] = 10*nM
    
if (False):
    print "Generating spiral waves.."
    # this one generates a spirsl wave by Ca2+ increase
    # set background
    np_state_new[species['Ca'],:,:] = 130.*nM
    np_state_new[species['n'],:,:] = 0.96585*nmax
    np_state_new[species['p'],:,:] = 95.*nM

    # and initiate circular wave
    np_state_new[species['Ca'], dx_/2-5:dx_/2+5, dx_/2-5:dx_/2+5] = 920*nM

    # and this is the refractory region
    #np_state_new[species['Ca'], dx_/2-20:dx_/2, dx_/2-20:dx_/2] = 10*nM
    #np_state_new[species['p'],:,:] = 0
