# matplotlib
import matplotlib.pyplot as plt

# for io
import gpgmp.io
import os

# for constant drift/diff. solution
import gpgmp.test.homogeneous_drift_diffusion

# numpy
import numpy

# need error function
from scipy.special import erf

import gpgmp.common.jobs
import pickle
import importlib

import os.path
import sys

# analytic solution
def analytic(xc, d, t, mu, sigma):
    ana = 1./(numpy.exp((xc-t*mu)**2/(2.*t*(2*d+t*sigma**2)))*numpy.sqrt(2.*numpy.pi)*numpy.sqrt(t*(2.*d + t*sigma**2)))

    return ana

# analytic solution for the discrete contribution:
def analytic_discrete(xc, d, t, mu, sigma, k):
    # standard values for the cutoffs
    # todo: this will depend on the variance
    # todo: maybe make these parameters as well? to guarantee consistency
    # between the analytic solution and the code
    bl = -10.
    br = 10.
    db = (br-bl)/k
    
    xcr = numpy.zeros_like(xc)

    for i in range(k):
        xcr = xcr + (-erf((bl + ((-bl + br)*i)/k - mu)/(numpy.sqrt(2)*sigma))
                      + erf((br + br*i + bl*(-1 - i + k) - k*mu)/(numpy.sqrt(2)*k*sigma))) \
                      /(4.*numpy.exp((((bl - br)*(1 + 2*i) - 2*bl*k)*t + 2*k*xc)**2/(16.*d*k**2*t))
                        *numpy.sqrt(numpy.pi)*numpy.sqrt(d*t))
    
    return xcr

def run(options, datadir):

    # get executable name
    executable = options['inchmanexec']

    # set output files
    newoptions = options
    newoptions['outlog'] = 'randomdrift.out'
    newoptions['errlog'] = 'randomdrift.err'

    # parameters
    params = "{0:s}/Random\ Drift.xml --cl-source {1:s} --output random_drift.h5".format(options['inchmanexampledir'], options['inchmancldir'])

    # run the job
    gpgmp.common.jobs.run(newoptions, executable, datadir, params)


def analyze(options, datadir, filename='random_drift.h5'):
    analyze(datadir, filename)

def analyze(datadir, filename):
    #parameters
    diffx = 1.
    diffy = 0.5
    tmax = 2.
    mux = 0.
    sigmax = 2.
    muy = 0.
    sigmay = 1.
    discK = 5
    numMols = 10000.
    numRuns = 1
    length = 40.
    nx = 512
    runTime = 0.

    # read out run time
    timefile = open(datadir+'runtime.pkl', 'r')
    runtime = pickle.load(timefile)
    timefile.close()

    # get non-invidual species
    nua, times, species, nruns = gpgmp.io.read_gmp_hdf5(datadir+filename, deterministic=False)
    
    # take latest time
    targ = (numpy.argsort(times))[-1]

    # read in results
    cellpos, ids = gpgmp.io.read_gmp_individuals(datadir+filename, 2, 'A')

    # now compute 2d cell pos
    x = numpy.int_(cellpos)/512
    y = cellpos % 512

    # make histogram
    minx = -length/2.
    maxx = length/2.
    dx = (maxx-minx)/float(nx)
    edges=numpy.mgrid[0:nx+1]*dx+minx
    xc = edges[:-1]+dx
    
    # and convert cell indices into positions
    xpos = x * dx - length/2.
    ypos = y * dx - length/2.

    # average
    n = numpy.mean(nua,axis=0)
    statea=numpy.squeeze(n[2,numpy.where(species == 'A'),:,:])
    stateb=numpy.zeros((nx,nx))

    for ix in range(discK):
        for iy in range(discK):
            name = "B{0:d}{1:d}".format(ix, iy)            
            stb = numpy.squeeze(n[2,numpy.where(species == name),:,:])
            stateb = stateb + stb

    numMolsB = numpy.sum(stateb)

    # get analytic solutions
    anay = analytic(xc, diffx, tmax, mux, sigmax)
    anax = analytic(xc, diffy, tmax, muy, sigmay)
    anadiscy = analytic_discrete(xc, diffx, tmax, mux, sigmax, discK)
    anadiscx = analytic_discrete(xc, diffy, tmax, muy, sigmay, discK)

    # and plot it
    hasmatplotlib = False        
    try:
        matplotlib = importlib.import_module('matplotlib')
        #matplotlib.use('PDF')
        plt = importlib.import_module('matplotlib.pyplot')
        hasmatplotlib = True        
    except:
        print("Could not find matplotlib. Skipping all plots..")

    if hasmatplotlib:
        width1c = 2.*3.27
        fig=plt.figure(1, figsize=(2.*width1c, 2.*width1c), dpi=75)

        plt.subplot(221)
        pdf, bins, patches = plt.hist(xpos, bins=edges, normed=True, edgecolor='none')
        plt.plot(xc, numpy.sum(statea, axis=1)/numMols/dx, "m")
        plt.plot(xc, anax, "r")

        plt.subplot(222)
        pdf, bins, patches = plt.hist(ypos, bins=edges, normed=True, edgecolor='none')
        plt.plot(xc, numpy.sum(statea, axis=0)/numMols/dx, "m")
        plt.plot(xc, anay, "r")
        #plt.plot(xc, numpy.sum(anahom, axis=0), "g")

        plt.subplot(223)
        plt.plot(xc, numpy.sum(stateb, axis=1)/numMolsB/dx)
        plt.plot(xc, anax, "r")
        plt.plot(xc, anadiscx, "g")
        #plt.plot(xc, numpy.sum(anahom, axis=1), "g")

        plt.subplot(224)
        plt.plot(xc, numpy.sum(stateb, axis=0)/numMolsB/dx)
        plt.plot(xc, anay, "r")
        plt.plot(xc, anadiscy, "g")
        #plt.plot(xc, numpy.sum(anahom, axis=0), "g")

        plt.savefig(datadir+"/random_drift.eps", format="eps")

    # compute RMSE
    rmsex = numpy.sqrt(numpy.sum((numpy.sum(statea, axis=1)/numMols/dx-anax)**2)/(numpy.shape(statea)[0]*numpy.shape(statea)[1]))
    rmsey = numpy.sqrt(numpy.sum((numpy.sum(statea, axis=0)/numMols/dx-anay)**2)/(numpy.shape(statea)[0]*numpy.shape(statea)[1]))

    passed = False
    
    if (rmsex < 5e-4) and (rmsey < 5e-4):
        passed = True

    return passed, runtime

# analyze it
if __name__ == "__main__":
    # print command line arguments
    dirname, filename = os.path.split(sys.argv[1])

    if dirname=='':
        dirname = './'

    analyze(dirname, filename)
