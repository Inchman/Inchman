# import matplotlib
import matplotlib
import matplotlib.pyplot as plt

# import system libs and numpy
import gpgmp.io
import gpgmp.job
import numpy

# analytic solution
def source_diffusion_x(xb,yy,phi,ll):
    numCompartments=numpy.shape(xb)[0]

    retval = numpy.zeros((numCompartments), dtype=float)

    for i in range(0, numCompartments):
        xx = xb[i]
        for k in range(0,200) :
            retval[i] = retval[i] + 4.*phi/numpy.pi*numpy.exp(-(2.*k+1.)*numpy.pi*yy/ll)*numpy.sin((2.*k+1.)*numpy.pi*xx/ll)/(2.*k+1.)
        
    return retval

# run the job
def run(testdir, initscriptdir, plot=False, norun=False, executable='gpgmp'):

    # parameters
    length=40.
    nx=64
    tmax=1000
    diffusionConstant=1.
    sourceStrength = 10.
    dx = length/nx
    numRuns = 100

    # run it and get the timing
    runTime = -1
    if (norun==False):
        runTime = gpgmp.job.run_job(testdir, executable, parameters='{0:g} {1:g}'.format(diffusionConstant, sourceStrength),
                                    problem=4, length=length, nx=nx, runtime=tmax, dtout=20., solver=2, numRuns = numRuns,
                                    datafile='test_infinite_slab.h5', outlog='infinite_slab.out', errlog='infinite_slab.err')

    # read out the results
    n, times, species, nruns = gpgmp.io.read_gmp_hdf5('test_infinite_slab')

    # sort it after time
    timesarg=numpy.argsort(times)
    sampleTime=times[timesarg[-1]]

    # average over runs
    na = numpy.mean(n[:, timesarg[-1], 0, :, :], 0)

    # compute geometry
    xb = numpy.mgrid[0:nx]*(length/nx)

    # compute analytic solution and rmse
    analytic = source_diffusion_x(xb+dx,xb[1], sourceStrength, length+dx)
    rmse = numpy.sqrt(numpy.sum((na[:,0]-analytic)**2)/numpy.shape(xb)[0])

    if rmse<=1:
        passed = True
    else:
        passed = False

    # plot it if necessary
    if plot:
        plt.clf()
        plt.plot(xb, na[:,0], "ro")
        plt.plot(xb, analytic)
        plt.savefig("infinite_slab.eps", format="eps")
                
    return passed, runTime

