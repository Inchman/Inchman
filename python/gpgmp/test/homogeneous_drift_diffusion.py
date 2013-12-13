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

# compute theoretical solution for 2D biased diffusion
def biased_diffusion_2d(xc, yc, t, diffx, diffy, qx, qy):
    # make grid from x and y
    gxc=numpy.outer(xc, numpy.ones_like(xc))
    gyc=numpy.outer(numpy.ones_like(yc), yc)

    px = 1./(4.*diffx*numpy.pi*t)**(0.5)*numpy.exp(-((gxc-t*qx)**2)/(4.*diffx*t))
    py = 1./(4.*diffy*numpy.pi*t)**(0.5)*numpy.exp(-((gyc-t*qy)**2)/(4.*diffy*t))
    return px*py

# compute rmse
def compute_rmse(theo, data):
    numMolecules = numpy.sum(data)

    return numpy.sqrt(numpy.sum((data/numMolecules-theo)**2)/(numpy.shape(data)[0]*numpy.shape(data)[1]))

def run(options, datadir, individual=False):

    # get executable name
    executable = options['inchmanexec']

    # parameters
    if individual:
        params = "{0:s}/Homogeneous\ Drift\ Individual.xml --cl-source {1:s} --output test_homogeneous_drift_individual.h5".format(options['inchmanexampledir'], options['inchmancldir'])
    else:
        params = "{0:s}/Homogeneous\ Drift.xml --cl-source {1:s} --output test_homogeneous_drift.h5".format(options['inchmanexampledir'], options['inchmancldir'])

    # set output files
    newoptions = options
    if individual:
        options['outlog'] = 'homogeneous_drift_individual.out'
        options['errlog'] = 'homogeneous_drift_individual.err'
    else:
        options['outlog'] = 'homogeneous_drift.out'
        options['errlog'] = 'homogeneous_drift.err'

    # run the job
    gpgmp.common.jobs.run(options, executable, datadir, params)

# run the job
#def _deprec(testdir, initscriptdir, plot=False, norun=False, executable='gpgmp', deterministic=False, aqss=False):
# analyze
# for backwards compability..
def analyze(options, datadir, filename='test_homogeneous_drift_diffusion.h5', individual=False, field=False):
    analyze(datadir, filename)

def analyze(datadir, filename):

    # parameters
    length = 32.
    nx = 512
    tmax = 5.
    diffusionConstantX = 0.6
    diffusionConstantY = 0.8
    rx = 0.8
    ry = 0.4
    numMolecules = 100000

    numruns = 10

#    # run it and get the timing
#    if deterministic:
#        numruns = 1
#        if aqss:
#            solver = 6
#        else:
#            solver = 5
#    else:
#        solver = 0

    # run it and get the timing
#    if norun==False:
#        runTime = gpgmp.job.run_job(testdir, executable, parameters='{0:g} {1:g} {2:g} {3:g} {4:g}'.format(diffusionConstantX, diffusionConstantY, rx, ry, numMolecules),
 #                                   problem=7, length=length, nx=nx, runtime=tmax, dtout=1., solver=solver, numRuns=numruns,
 #                                   datafile='test_homogeneous_drift.h5', outlog='homogeneous_drift.out', errlog='homogeneous_drift.err')

    # read out the results
    n, times, species, nruns = gpgmp.io.read_gmp_hdf5(datadir+filename, deterministic=False)

    # read out run time
    timefile = open(datadir+'runtime.pkl', 'r')
    runTime = pickle.load(timefile)
    timefile.close()

    # sort it after time
    timesarg=numpy.argsort(times)
    sampleTime=times[timesarg[-1]]

    # compute coordinate system
    dx = length/float(nx)
    xc = (numpy.mgrid[0:nx])*dx-length/2.
    yc = (numpy.mgrid[0:nx])*dx-length/2.

    # average over runs
    nstate=numpy.mean(n[:,timesarg[-1],0,:,:], axis=0)

    # get analytic solution
    #theo = numpy.transpose(biased_diffusion_2d(xc, yc, sampleTime, diffusionConstantY, diffusionConstantX, ry, rx)*dx**2)
    theo = (biased_diffusion_2d(xc, yc, sampleTime, diffusionConstantY, diffusionConstantX, ry, rx)*dx**2)

    # compute rmse 
    rmse = compute_rmse(theo, nstate)

    # and plot if needed
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
        plt.plot(xc, theo[:,nx/2]*numMolecules,"b--")
        plt.savefig(datadir+"homogeneous_drift_xslice.eps", format="eps")
        plt.clf()
        plt.plot(xc, nstate[nx/2,:], "r^", label="GPGMP")
        plt.plot(xc, theo[nx/2,:]*numMolecules,"b--")
        plt.savefig(datadir+"homogeneous_drift_yslice.eps", format="eps")

    if rmse<=1e-4:
        passed = True
    else:
        passed = False

    return passed, runTime

def run_field(options, datadir):

    # get executable name
    executable = options['inchmanexec']

    params = "{0:s}/Homogeneous\ Drift\ Field.xml --cl-source {1:s} --output test_homogeneous_drift_field.h5".format(options['inchmanexampledir'], options['inchmancldir'])

    # set output files
    newoptions = options
    options['outlog'] = 'homogeneous_drift_field.out'
    options['errlog'] = 'homogeneous_drift_field.err'

    # run the job
    gpgmp.common.jobs.run(options, executable, datadir, params)

# analyze it
if __name__ == "__main__":
    # print command line arguments
    dirname, filename = os.path.split(sys.argv[1])

    if dirname=='':
        dirname = './'

    analyze(dirname, filename)
