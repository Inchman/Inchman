# import matplotlib
import matplotlib
import matplotlib.pyplot as plt

# import system libs and numpy
import gpgmp.io
#import gpgmp.job
import gpgmp.constants
import numpy
import common.jobs
import pickle
import importlib

import os.path
import sys

# run the job
def run(options, datadir):
    
    # get executable name
    executable = options['inchmanexec']

    # set output files
    newoptions = options
    newoptions['outlog'] = 'annihilation.out'
    newoptions['errlog'] = 'annihilation.err'

    # parameters
    params = "{0:s}/AB\ Annihilation.xml --cl-source {1:s} --output abannihilation.h5".format(options['inchmanexampledir'], options['inchmancldir'])

    # run the job
    common.jobs.run(options, executable, datadir, params)

# for backwards compability..
def analyze(options, datadir, filenam='abannihilation'):
    analyze(datadir, filename)

def analyze(datadir, filename):
    # read out the results
    n, times, species, nruns = gpgmp.io.read_gmp_hdf5(datadir+filename, deterministic=False)

    # read out run time
    timefile = open(datadir+'runtime.pkl', 'r')
    runtime = pickle.load(timefile)
    timefile.close()

    # sort it after time
    timesarg=numpy.argsort(times)
    sampleTime=times[timesarg[-1]]

    # average over runs
    nstate=numpy.mean(n[:,timesarg[-1],:,:,:], axis=0)

    # make plots if needed
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
        fig=plt.figure(1, figsize=(gpgmp.constants.width2c,gpgmp.constants.height1c), dpi=300)
        plt.clf()

        fig.subplots_adjust(hspace=0.1, bottom=0.15, right=0.95, left=0.1)
        ax1 = fig.add_subplot(1, 2, 1)
        ax2 = fig.add_subplot(1, 2, 2, sharex=ax1, sharey=ax1)

        ax1.imshow(nstate[0,:,:])
        ax2.imshow(nstate[1,:,:])
        plt.savefig(datadir+"annihilation_2d.eps", format="eps")

    return runtime

# analyze it
if __name__ == "__main__":
    # print command line arguments
    dirname, filename = os.path.split(sys.argv[1])

    if dirname=='':
        dirname = './'
            
    
    analyze(dirname, filename)
