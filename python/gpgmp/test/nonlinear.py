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


# returns the analytic solution for the nonlinear model (Zhang et al., 1997)
def analytic(xc, t, nldiff=1., omega=1., theta=1., x0=0.):
    mx = x0*numpy.exp(-(omega+theta)*t)
    st = nldiff/omega*(1.-numpy.exp(-2.*omega*t))
    return 1./numpy.sqrt(2.*numpy.pi*st)*numpy.exp(-(xc-mx)**2/(2.*st))

# compute rmse
def compute_rmse(state,theo):
    numMolecules = numpy.sum(state)
    # this is only 1d
    return numpy.sqrt(numpy.sum((state-theo)**2)/(numpy.shape(state)[0]))/numMolecules

def run(options, datadir):

    # get executable name
    executable = options['inchmanexec']

    # set output files
    newoptions = options
    newoptions['outlog'] = 'nonlinear.out'
    newoptions['errlog'] = 'nonlinear.err'

    # parameters
    params = "{0:s}/Nonlinear.xml --cl-source {1:s} --output nonlinear.h5".format(options['inchmanexampledir'], options['inchmancldir'])

    # run the job
    gpgmp.common.jobs.run(newoptions, executable, datadir, params)

def analyze(options, datadir, filename='nonlinear.h5'):
    analyze(datadir, filename)

def analyze(datadir, filename):

    # Parameters
    length = 32.
    nx = 64
    tmax = 5
    diffx= 1.
    diffy= 1.
    theta = 1.
    omega = 1.
    numMolecules = 10000
    numruns = 10

    # read out the results
    n, times, species, nruns = gpgmp.io.read_gmp_hdf5(datadir+filename)

    # read out run time
    timefile = open(datadir+'runtime.pkl', 'r')
    runtime = pickle.load(timefile)
    timefile.close()

    # sort it after time
    timesarg=numpy.argsort(times)
    sampleTime=times[timesarg[-1]]

    # average over runs
    nstate=numpy.transpose(numpy.mean(n[:,timesarg[-1],0,:,:], axis=0))

    # compute coordinate system
    xc, dx = (float(length)/float(nx)*(numpy.mgrid[0:nx]-nx/2.)+length/8., float(length)/float(nx))
    yc = xc
    x0 = xc[nx/2]
    y0 = yc[nx/2]

    # compute analytic solution (1d)
    theo = analytic(xc, sampleTime, diffx, omega, theta, x0)*dx*numMolecules

    # compute rmse
    rmse = compute_rmse(nstate[:,nx/2], theo)

    # and plot it
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
        plt.plot(xc, nstate[:,nx/2], "r^", label="GPGMP")
        plt.plot(xc, theo,"b--")
        plt.savefig(datadir+"/nonlinear.eps", format="eps")

    if rmse<=1e-3:
        passed=True
    else:
        passed = False

    return passed, runtime

# analyze it
if __name__ == "__main__":
    # print command line arguments
    dirname, filename = os.path.split(sys.argv[1])

    if dirname=='':
        dirname = './'

    analyze(dirname, filename)
