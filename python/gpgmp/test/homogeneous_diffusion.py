# import matplotlib
import matplotlib
import matplotlib.pyplot as plt

# import system libs and numpy
import gpgmp.io
import numpy

import gpgmp.common.jobs
import pickle
import importlib

import os.path
import sys

# gives the analytic solution for the homogeneous diffusion test problem
def diffusion(x, t, y=0, z=0, diffConst=1, dim=1.):
    return 1./(4.*diffConst*numpy.pi*t)**(0.5*dim)*numpy.exp(-(x*x+y*y+z*z)/(4.*diffConst*t))

# compute rmse given the data
def rmse(xb, n, t, y=0, z=0, diffusionConstant=1, dim=1):
    return numpy.sqrt(numpy.sum((n-diffusion(xb,t,y, z, diffusionConstant, dim))**2)/numpy.shape(xb)[0])

def run(options, datadir):

    # get executable name
    executable = options['inchmanexec']

    # parameters
    params = "{0:s}/Homogeneous\ Diffusion.xml --cl-source {1:s} --output test_homogeneous_diffusion.h5".format(options['inchmanexampledir'], options['inchmancldir'])

    # set output files
    newoptions = options
    options['outlog'] = 'homogeneous_diffusion.out'
    options['errlog'] = 'homogeneous_diffusion.err'

    # run the job
    gpgmp.common.jobs.run(options, executable, datadir, params)

# analyze
# for backwards compability..
def analyze(options, datadir, filename='test_homogeneous_diffusion.h5'):
    analyze(datadir, filename)

def analyze(datadir, filename):
    # parameters
    diffusionConstant = 1.
    numMolecules = 1e5
    length = 100.
    nx = 128
    runtime = 40.

    # read out run time
    timefile = open(datadir+'runtime.pkl', 'r')
    runTime = pickle.load(timefile)
    timefile.close()

    # read out the results
    n, times, species, nruns = gpgmp.io.read_gmp_hdf5(datadir+filename)

    # sort it after time
    timesarg=numpy.argsort(times)
    sampleTime=times[timesarg[-1]]

    # compute coordinate system
    dx=length/nx
    xgmp=(numpy.mgrid[0:nx])*dx-length/2.

    # average over runs
    ns = n[:,timesarg[-1], 0, :,:]
    nav = numpy.mean(ns,0)

    # compute rmse for slices in x and y
    rmsex = rmse(xgmp, nav[nx/2,:]/(numMolecules*dx**2), sampleTime,
                 y=xgmp[nx/2], dim=2, diffusionConstant=diffusionConstant)
    rmsey = rmse(xgmp, nav[:,nx/2]/(numMolecules*dx**2), sampleTime,
                 y=xgmp[nx/2], dim=2, diffusionConstant=diffusionConstant)

    if (rmsex<=5e-4 and rmsey<=5e-4):
        passed = True
    else:
        passed = False

    # plot it if needed
    hasmatplotlib = False        
    try:
        matplotlib = importlib.import_module('matplotlib')
        #matplotlib.use('PDF')
        plt = importlib.import_module('matplotlib.pyplot')
        hasmatplotlib = True        
    except:
        print("Could not find matplotlib. Skipping all plots..")

    if hasmatplotlib:
        plt.clf()
        plt.plot(xgmp, nav[:,nx/2], "r^", label="GPGMP")
        plt.plot(xgmp, diffusion(xgmp, sampleTime, y=xgmp[nx/2], dim=2, diffConst=diffusionConstant)*numMolecules*dx**2,"b--")
        plt.savefig(datadir+"homogeneous_diffusion_xslice.eps", format="eps")
        plt.clf()
        plt.plot(xgmp, nav[nx/2,:], "r^", label="GPGMP")
        plt.plot(xgmp, diffusion(xgmp, sampleTime, y=xgmp[nx/2], dim=2, diffConst=diffusionConstant)*numMolecules*dx**2,"b--")
        plt.savefig(datadir+"homogeneous_diffusion_yslice.eps", format="eps")
        
    return passed, runTime

# analyze it
if __name__ == "__main__":
    # print command line arguments
    dirname, filename = os.path.split(sys.argv[1])

    if dirname=='':
        dirname = './'
            
    analyze(dirname, filename)
