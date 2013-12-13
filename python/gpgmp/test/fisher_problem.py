# import matplotlib
import matplotlib
import matplotlib.pyplot as plt

# import analysis helpers
import gpgmp.constants

# import system libs and numpy
import gpgmp.io
import numpy
import shutil

import common.jobs
import pickle
import importlib

import os.path
import sys

# analytic solutions
def fisher(x, t, c=0.5):
      z = x + 5./numpy.sqrt(6.)*t
      return 1./(1.+c*numpy.exp(-z/numpy.sqrt(6.)))**2

def rmse_fisher(xb, species, t, c=1):
    return numpy.sqrt(sum((species-fisher(xb,t,c))**2)/numpy.shape(xb)[0])

# for debug - gives the reaction rate etc in terms of simulation parameters
def getRates(baseConcentration=1e-6, length=100., nx = 256.):
      # subvolume in litres
      subvolume = (length/nx*1e-6)**3*1e3
      
      # avogadro constant
      avogadro = 6.0221415e23

      # reaction rate in physical units
      ratemol = 1./baseConcentration

      # reaction rate in 1/s
      rate = ratemol/(avogadro*subvolume)
# run
def run(options, datadir, individual=False):
      # get executable name
      executable = options['inchmanexec']

      # set output files
      newoptions = options
      if individual:
            newoptions['outlog'] = 'fisher_individual.out'
            newoptions['errlog'] = 'fisher_individual.err'
      else:
            newoptions['outlog'] = 'fisher.out'
            newoptions['errlog'] = 'fisher.err'

      # parameters
      if individual:
            params = "{0:s}/Fisher\ Problem\ Individual.xml --cl-source {1:s} --output test_fisher_individual.h5".format(options['inchmanexampledir'], options['inchmancldir'])
      else:
            params = "{0:s}/Fisher\ Problem.xml --cl-source {1:s} --output test_fisher.h5".format(options['inchmanexampledir'], options['inchmancldir'])

      # run the job
      gpgmp.common.jobs.run(options, executable, datadir, params)

def analyze(options, datadir, filename='test_fisher.h5', individual=False):
    analyze(datadir, filename, individual)

def analyze(datadir, filename, individual):

    nx = 256
    length = 100.
    #baseConc = 6.4e-5
    if individual:
          baseConc=0.0020971537215416936*5e-5
    else:
          baseConc=0.0020971537215416936

    # read out the results
    n, times, species, nruns = gpgmp.io.read_gmp_hdf5(datadir+filename, deterministic=False)

    # read out run time
    timefile = open(datadir+'runtime.pkl', 'r')
    runtime = pickle.load(timefile)
    timefile.close()

    # sort it after time
    timesarg=numpy.argsort(times)

    # get time index closest to desired time
    timeIndex = (numpy.where(times >= 5.))[0][0]

    # average over number of runs
    nt = (numpy.mean(n,0))[timeIndex, 0, :, :]

    # average over y values and number of runs
    nav=numpy.mean(nt,0)

    # for individuals we also plot the collective species
    if individual:
          ntc = (numpy.mean(n,0))[timeIndex, 1, :, :]
          navc=numpy.mean(ntc,0)

    # compute coordinates values
    dx=length/nx
    xgmp=(numpy.mgrid[0:nx]+0.5)*dx-3.*length/4.

    # number of particles
    volumeLitres=(dx*1e-6)**3*1e3
    #nparticles=(nx/64)**3*baseConc*gpgmp.constants.na*volumeLitres
    #nparticles=75277
    nparticles=baseConc*gpgmp.constants.na*volumeLitres

    # plot it
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
          plt.plot(xgmp, nav/nparticles, "r^")
          plt.plot(xgmp, fisher(xgmp, 5., c=1))
          plt.plot(xgmp, fisher(xgmp, 0., c=1), "g-")
          if individual:
                plt.plot(xgmp, navc/nparticles, "go")
                plt.savefig(datadir+"fisher_problem_individual.eps", format="eps")
          else:
                plt.savefig(datadir+"fisher_problem.eps", format="eps")

    # compute rms
    rmse = rmse_fisher(xgmp, nav/nparticles, 5.)

    if (rmse<=0.0001):
        passed = True
    else:
        passed = False

    return True, runtime


# analyze it
if __name__ == "__main__":
    # print command line arguments
    dirname, filename = os.path.split(sys.argv[1])

    if dirname=='':
        dirname = './'
            
    individual = False

    if 'individual' in sys.argv:
        individual = True
        
    analyze(dirname, filename, individual)
