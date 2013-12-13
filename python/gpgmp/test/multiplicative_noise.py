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

# analytic solution
def analytic_2d(x, x0, y, y0, sigmax, sigmay,  t, mux=0, muy=0):
    # make grid from x and y
    gx=numpy.outer(x, numpy.ones_like(x))
    gy=numpy.outer(numpy.ones_like(y), y)

    return 1/(2.*numpy.exp((((t*(-2.*mux + sigmax**2))/2. + numpy.log(gx/x0))**2./sigmax**2 + ((t*(-2.*muy + sigmay**2))/2. + numpy.log(gy/y0))**2/sigmay**2)/(2.*t))*numpy.pi*t*gx*gy*sigmax*sigmay)

def compute_rmse(state, theo):
    numMolecules = numpy.sum(state)
    return numpy.sqrt(numpy.sum((state-theo)**2)/(numpy.shape(state)[0]*numpy.shape(state)[1]))/numMolecules

# run the job
#def run(testdir, plot=False, executable='gpgmp', deterministic=False, aqss=False):
def run(options, datadir):

    # get executable name
    executable = options['inchmanexec']

    # set output files
    newoptions = options
    newoptions['outlog'] = 'multiplicative_noise.out'
    newoptions['errlog'] = 'multiplicative_noise.err'

    # parameters
    params = "{0:s}/Multiplicative\ Noise.xml --cl-source {1:s} --output multiplicative_noise.h5".format(options['inchmanexampledir'], options['inchmancldir'])

    # run the job
    gpgmp.common.jobs.run(newoptions, executable, datadir, params)

def analyze(options, datadir, filename='multiplicative_noise.h5'):
    analyze(datadir, filename)

def analyze(datadir, filename):

    # Parameters
    length = 32.
    nx = 64
    tmax = 0.0006
    cx = 5
    cy = 4
    mux = 0.1
    muy = 0.2
    numMolecules = 100000
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
    xc, dx = (float(length)/float(nx)*(numpy.mgrid[0:nx]-nx/2.)+length, float(length)/float(nx))
    yc = xc

    # compute analytic solution
    x0 = xc[nx/4]
    y0 = yc[nx/4]
    theo = analytic_2d(xc, x0, yc, y0, cx, cy, sampleTime, mux, muy)*dx**2*numMolecules

    # compute RMSE
    rmse = compute_rmse(nstate,theo)

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
        plt.plot(xc, nstate[:,nx/4], "r^", label="GPGMP")
        plt.plot(xc, theo[:,nx/4],"b--")
        plt.savefig(datadir+"multiplicative_noise_xslice.eps", format="eps")
        plt.clf()
        plt.plot(xc, nstate[nx/4,:], "r^", label="GPGMP")
        plt.plot(xc, theo[nx/4,:],"b--")
        plt.savefig(datadir+"multiplicative_noise_yslice.eps", format="eps")

    if rmse<1e-4:
        passed=True
    else:
        passed=False

    return passed, runtime

# analyze it
if __name__ == "__main__":
    # print command line arguments
    dirname, filename = os.path.split(sys.argv[1])

    if dirname=='':
        dirname = './'

    analyze(dirname, filename)
