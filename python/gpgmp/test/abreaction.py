# matplotlib
import matplotlib.pyplot as plt

# for io
import gpgmp.io

# numpy
import numpy

import gpgmp.common.jobs
import pickle
import importlib

import os.path
import sys

# run the job
#def run(testdir, plot=False, executable='gpgmp', deterministic=False, aqss=False):
def run(options, datadir, individual=False, new=False):

    # get executable name
    executable = options['inchmanexec']

    # set output files
    newoptions = options
    if individual:
        if new:
            newoptions['outlog'] = 'abreaction_individual_new.out'
            newoptions['errlog'] = 'abreaction_individual_new.err'
        else:
            newoptions['outlog'] = 'abreaction_individual.out'
            newoptions['errlog'] = 'abreaction_individual.err'
    else:
        newoptions['outlog'] = 'abreaction.out'
        newoptions['errlog'] = 'abreaction.err'

    # parameters
    if individual:
        if new:
            params = "{0:s}/New\ Individuals\ Test.xml --cl-source {1:s} --output abreaction_individual_new.h5".format(options['inchmanexampledir'], options['inchmancldir'])
        else:
            params = "{0:s}/AB\ Reaction\ Individual.xml --cl-source {1:s} --output abreaction_individual.h5".format(options['inchmanexampledir'], options['inchmancldir'])

    else:
        params = "{0:s}/AB\ Reaction.xml --cl-source {1:s} --output abreaction.h5".format(options['inchmanexampledir'], options['inchmancldir'])

    # run the job
    gpgmp.common.jobs.run(options, executable, datadir, params)

# analyze
# for backwards compability..
def analyze(options, datadir, filename='ab_reaction', individual=False):
    analyze(datadir, filename, individual)

def analyze(datadir, filename, individual):

    n, times, species, nruns = gpgmp.io.read_gmp_hdf5(datadir+filename, deterministic=False)

    # read out run time
    timefile = open(datadir+'runtime.pkl', 'r')
    runtime = pickle.load(timefile)
    timefile.close()

    # sort it after time
    timesarg=numpy.argsort(times)

    # compute mean
    na = n[0,timesarg,0,:,:]
    nb = n[0,timesarg,1,:,:]
    namean = numpy.mean(numpy.mean(na,1),1)
    nbmean = numpy.mean(numpy.mean(nb,1),1)

    if individual:
        nac = n[0,timesarg,2,:,:]
        nbc = n[0,timesarg,3,:,:]
        nacmean = numpy.mean(numpy.mean(nac,1),1)
        nbcmean = numpy.mean(numpy.mean(nbc,1),1)


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
        plt.plot(times[timesarg],namean, "ro", label="Species A")
        plt.plot(times[timesarg],nbmean, "bo", label="Species B")

        if individual:
            plt.plot(times[timesarg],nacmean, "r-", label="Species A")
            plt.plot(times[timesarg],nbcmean, "b-", label="Species B")

        plt.title('A+B reaction problem without diffusion')
        plt.xlabel('t [s]')
        plt.ylabel('N')
        plt.legend(loc='lower right')

        plt.savefig(datadir+"aplusb_reaction.eps", format="eps")

    # now do we pass the test
    #if deterministic:
    #    acca = numpy.abs(namean[-1]-10.)/10.
    #    accb = numpy.abs(nbmean[-1]-10.)/10.
    #else:
    acca = numpy.abs(namean[-1]-9.6)/9.6
    accb = numpy.abs(nbmean[-1]-12.1)/12.1
    print("Mean: a={0:g}, b={1:g}".format(namean[-1], nbmean[-1]))

    if (acca<=0.01 and accb<=0.01):
        passed = True
    else:
        passed = False

    return passed, runtime

# analyze it
if __name__ == "__main__":
    # print command line arguments
    dirname, filename = os.path.split(sys.argv[1])

    if dirname=='':
        dirname = './'
            
    individual = False
    new = False

    if 'individual' in sys.argv:
        individual = True
        
    analyze(dirname, filename, individual)
