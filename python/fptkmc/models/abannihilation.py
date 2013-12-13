'''
Created on 18/10/2012

@author: matthias
'''
import common.jobs
import os
import pickle
import numpy
import math
import matplotlib
matplotlib.use('PS')

import matplotlib.pyplot as plt

def poisson(n, k1=0.2, k2=0.02, vol=8, b0=1):
    # need stupid work around to avoid using scipy.misc.factorial .. which
    # doesn't work on poorly configured cluster!
    retarr = numpy.zeros_like(n, dtype=numpy.float)
    len = numpy.shape(n)[0]
    
    for t in range(len):
        retarr[t] = 1./float(math.factorial(n[t]))*numpy.power((k2*vol**2/(k1*b0)),n[t])*numpy.exp(-k2*vol**2/(k1*b0))

    return retarr

def average(k1, k2, vol, b0):
    return k2*vol**2/(k1*b0)

def getVfRange():
    return [1., 5., 10.]

def getDtRange():
    return [0.1, 0.01, 0.001, 0.0001, 0.00001]

def run(options, datadir):
    # get executable name
    executable = options['fptkmcexec']
    
    # loop over volfactor
    vfrange=getVfRange()
    
    for volfactor in vfrange:
        # get output directory name
        ddirloop = datadir+'vf_{0:g}/'.format(volfactor)
        
        # make directory
        os.makedirs(ddirloop)
        
        # compute parameters and run        
        parameters = '{0:g}'.format(volfactor)
        
        # write parameter information 
        paramfile = open(ddirloop+'parameters.pkl', 'w')
        pickle.dump(volfactor, paramfile)
        paramfile.close()
        
        common.jobs.run(options, executable, ddirloop, parameters)

def createSmoldynInputFile(datadir, volfactor, dt):
    # compute parameters using the given volfactos
    size = 2.*math.pow(volfactor, 1./3.)
    k2 = 0.02/volfactor**2
    timeend = 60000.*volfactor
    dtout = 300.*volfactor
    
    # this is the PBS script template
    template = '''# 3D a+b annihilation

dim 3
species speciesA speciesB

difc speciesA 10
difc speciesB 10

time_start 0
time_stop {2:g}
time_step {4:g}

boundaries 0 0 {0:g} r
boundaries 1 0 {0:g} r
boundaries 2 0 {0:g} r

reaction anni speciesA + speciesB -> speciesB 0.2
reaction creation 0 -> speciesA {1:g}

mol 1 speciesB 1 1 1

output_files outfile.txt

cmd i 0 {2:g} {3:g} molcount outfile.txt

end_file
    '''.format(size, k2, timeend, dtout, dt)
    
    # write template to file    
    f = open(datadir+'/aplusb.input', 'w')
    f.write(template)
    f.close()
    
def runSmoldyn(options, datadir):
    # get executable name
    executable = options['smoldynexec']
    
    # loop over volfactor
    vfrange=getVfRange()
    
    # use standard dt
    dt = 0.01
    
    for volfactor in vfrange:
        # get output directory name
        ddirloop = datadir+'vf_{0:g}/'.format(volfactor)
        
        # make directory
        os.makedirs(ddirloop)
        
        # parameters for smoldyn is just the input file name        
        parameters = ddirloop+'/aplusb.input'
        
        # write parameter information 
        paramfile = open(ddirloop+'parameters.pkl', 'w')
        pickle.dump(volfactor, paramfile)
        pickle.dump(dt, paramfile)
        paramfile.close()
        
        # create input file
        createSmoldynInputFile(ddirloop, volfactor, dt)
        
        # run
        common.jobs.run(options, executable, ddirloop, parameters)

def runSmoldynAccuracy(options, datadir):
    # get executable name
    executable = options['smoldynexec']
    
    # loop over accuracy
    arange = getDtRange()
    
    # use normal volfactor
    volfactor = 1.
    
    for dt in arange:
        # get output directory name
        ddirloop = datadir+'dt_{0:g}/'.format(dt)
        
        # make directory
        os.makedirs(ddirloop)
        
        # parameters for smoldyn is just the input file name        
        parameters = ddirloop+'/aplusb.input'
        
        # write parameter information 
        paramfile = open(ddirloop+'parameters.pkl', 'w')
        pickle.dump(volfactor, paramfile)
        pickle.dump(dt, paramfile)
        paramfile.close()
        
        # create input file
        createSmoldynInputFile(ddirloop, volfactor, dt)
        
        # run
        common.jobs.run(options, executable, ddirloop, parameters)

def plotSmoldynAccuracy(optons, datadir):
     # loop over accuracy
    arange = getDtRange()

    # use standard volfactor
    volfactor = 1.
    
    # get data from directories
    for dt in arange:
        # get output directory name
        ddirloop = datadir+'dt_{0:g}/'.format(dt)
        
        # read in data
        dat=numpy.genfromtxt(ddirloop+'/outfile.txt')
        dumptime=dat[:,0]
        nb = dat[1:, 1]
        
        # extract all dump times which are a multiple of the autocorrelation time
        # (and strip first element as well)
        #dumpindices = numpy.where(dumptime.astype(int) % 300 == 0)
        #nb = ((dat[dumpindices, 1])[0])[1:]

        # clear plot
        plt.clf()
        
        # plot histogram
        hist,bins = numpy.histogram(nb, bins=range(17), normed=True)
        width = 0.7*(bins[1]-bins[0])
        center=(bins[:-1]+bins[1:])/2
        plt.bar(center,hist,align='center',width=width)
        
        # plot analytic solution
        size = 1.240700981798
        vol = 4./3.*numpy.pi*size**3*volfactor
        ana = poisson(center, k1=0.2, k2=0.02/volfactor**2, vol=vol)
        plt.plot(center, ana, "r")
        plt.savefig("./aplusb_dt_{0:g}.eps".format(dt), format="eps")
        
        # compute rms and KS
        nelements = float(numpy.size(hist))
        rms = numpy.sqrt(numpy.sum((hist-ana)**2)/nelements)
        print "For dt = {0:g}  we find: RMS={1:g}.".format(dt, rms)

def plotTimeStep(options, datadir, logfile='log'):
    
    # input name
    inputfile = datadir+'/'+logfile
    binaryname = inputfile+'.npy'
    
    # check if numpy binary file exists..
    if (os.path.isfile(binaryname)):
        # read in from binary
        print "Binary found. Reading in from {0:s}.".format(binaryname)
        data = numpy.load(binaryname)
    else:    
        print "No binary found. Reading in from ASCII file {0:s}.".format(inputfile)
        
        # read in data
        #data = numpy.genfromtxt(inputfile, dtype='float',  skip_header=5, skip_footer=1)
        
        # create empty numpy array
        data = numpy.array([0.])
        lineno=0
        
        file = open(inputfile, 'r')
        for line in file:
            line = line.strip()
            columns = line.split()
            # check if it's a comment and append if not
            try:
                dz = [float(x) for x in columns]
                
                # we only save dt ..
                oldsize = numpy.shape(data)[0]
                data.resize(oldsize+1)
                data[oldsize] = dz[3]

                if ((lineno % 5000) == 0 ):
                    print "Reading in line {0:d}".format(lineno)

            except ValueError:
                print "Skipping line {0:d}.".format(lineno)
            

            lineno = lineno + 1
            #if lineno>20:
            #    break
            
        file.close()
        
        # remove first entry
        data = data[1:]
        
        # save as binary
        numpy.save(inputfile+'.npy', data)
    
    # extract time step
    print data
    ts = data
    print numpy.min(ts), numpy.max(ts)
    
    # limits
    minx = -8
    maxx = 0
        
    # plot histogram
    plt.clf()
    n, bins, patches = plt.hist(ts, bins=numpy.logspace(minx, maxx, 50))
    plt.gca().set_xscale("log")
    plt.gca().set_xlim([10**minx, 10**maxx])
    plt.gca().set_xlabel("timestep")
    plt.gca().set_ylabel("# occurrences")
    
    plt.show()
    plt.savefig(datadir+"/timestep.eps", format="eps")
    
def plotSingle(options, datadir):
    # read in data
    dat=numpy.genfromtxt(datadir+'/output.dat')[50:,2]

    # plot histogram
    plt.hist(dat, bins=range(17), normed=True)

    # plot analytic solution
    volfactor = 100.
    size = 1.240700981798
    vol = 4./3.*numpy.pi*size**3*volfactor
    plt.plot(range(16), poisson(range(16), k1=0.2, k2=0.02/volfactor**2, vol=vol), "r")
    
    plt.savefig(datadir+"/aplusb.eps", format="eps")
    
def plotSingleSmoldyn(options, datadir):
    # read in data
    dat=numpy.genfromtxt(datadir+'/outfile.txt')[50:,1]

    # plot histogram
    plt.hist(dat, bins=range(17), normed=True)

    # plot analytic solution
    volfactor = 1.
    size = 1.240700981798
    vol = 4./3.*numpy.pi*size**3*volfactor
    plt.plot(range(16), poisson(range(16), k1=0.2, k2=0.02/volfactor**2, vol=vol), "r")
    
    plt.savefig(datadir+"/aplusb_smoldyn.eps", format="eps")
