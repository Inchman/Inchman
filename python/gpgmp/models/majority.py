import numpy
import gpgmp.io

import common.jobs
import common.options
import os
import pickle
import gc
import h5py

import random
import math
import shutil

import subprocess

import importlib

import numpy.random

def readFlamePosition(datadir, run, iteration):
    # reads in the flame position output and returns the red and green positions
    pos = numpy.genfromtxt(datadir+'/position_{0:d}.dat'.format(run))
    return pos[numpy.where(numpy.logical_and(pos[:,0]==iteration, pos[:,-1]==1)), 2:4][0,:,:], pos[numpy.where(numpy.logical_and(pos[:,0]==iteration, pos[:,-1]==0)), 2:4][0,:,:] 


def plotCorrelationComparison(dirs, resultsdir, iteration):
        # and plot them
    try:
        matplotlib = importlib.import_module('matplotlib')
        #matplotlib.use('PDF')
        plt = importlib.import_module('matplotlib.pyplot')
        hasmatplotlib = True        
    except:
        print("Could not find matplotlib. Skipping all plots..")

    if hasmatplotlib:
        # Make SDE comparison
        matplotlib.rc('text', usetex=True)
        plt.figure(figsize=(4*1.5, 3*1.5), dpi=600)

        # go through all directories and plot the correlation function
        colorcodes = ['b', 'r', 'g', 'k']
        for i, datadir in enumerate(dirs):

            # read in correlation function
            df = numpy.load(datadir+'/correlation.npz')
            g = df['g']
            
            # and plot it
            plt.plot(g[iteration, :, 0], g[iteration, :, 1], '-'+colorcodes[i])
            
        # set limits and labels
        plt.xlim(g[iteration, -1, 0], g[iteration, 0, 0])
        plt.ylim([0, 3.5])
        plt.xlabel(r'$r$')
        plt.ylabel(r'$g(r)$')
        
        # save figure
        plt.savefig(resultsdir+'/correlation_comparison.pdf')

        # delete plots
        plt.clf()

        # turn of latex support
        matplotlib.rc('text', usetex=False)

def reduceFlamePositionData(datadir):
    # computes the two-point correlation function for all experiments from the stored position data

    # read in parameters
    paramfile = open(datadir+'parameters.pkl', 'r')
    nruns = pickle.load(paramfile)
    tot = pickle.load(paramfile)
    paramfile.close()
    
    # these are the iteration numbers we store the position data for.. they're hard-wired into flame code
    iterations = [1, 100, 500, 1000, 2000]

    # output data will - for each stored iteration - contain the grid and the two-point correlation function
    output = numpy.zeros((len(iterations), 101, 2))
    count = numpy.zeros((len(iterations)))

    # collect position data from all experiments
    for run in range(nruns):
        # read in position data
        pos = numpy.genfromtxt(datadir+'/position_{0:d}.dat'.format(run))

        # go through iterations
        for i, iteration in enumerate(iterations):
            # get all reds from current iteration
            index = numpy.where(numpy.logical_and(pos[:,0]==iteration, pos[:,-1]==1))

            # extract the positions
            posred = (pos[index, 2:4])[0,:,:]

            # and add them to the output array
            if len(posred)>0:
                r, g, gamma = gpgmp.models.majority.computeTwoPointCorrelationRed(posred, h=0.01)
                output[i,:,0] = r
                output[i,:,1] += g
                count[i] += 1
            #else:
            #    print "Warning: No red robots found for run {0:d} (iteration {1:d}).".format(run, iteration)

    # do average
    for i in range(len(iterations)):
        output[i,:,1] /= count[i]

    # save the average correlation function
    numpy.savez(datadir+'correlation.npz', g=output)

    return output, count

def createPointSets(num, setType):
    # creates a number of point sets that can be used for testing purpose
    if setType=='regular':
        # creates a regular grid
        xx = numpy.tile(numpy.mgrid[0:1:complex(0,num)], num)
        yy = (numpy.transpose(numpy.tile(numpy.mgrid[0:1:complex(0,num)], (num,1)))).flatten()
        return numpy.transpose(numpy.array([xx,yy]))

    if setType=='poisson':
        # creates a Poisson process (actually a binomial process as the total number of particles is fixed..)
        return numpy.random.rand(num, 2)

def boxKernel(x, h=0.1):
    return numpy.asarray(numpy.logical_and(x <= h, x >=-h), dtype=float)/(2.*h)

def indicatorFunction(x, r):
    # returns 1 for each element in x that is inside disk around origin
    return numpy.asarray(x <= r, dtype=float)

def computeK(pos):
    # this gives an estimate of Ripley's K-function

    # check if we only got an x/y array
    if (numpy.shape(pos)[1]==2):
        posred = pos
    else:
        # get red index
        indred=numpy.where(pos[:,-1]==1)

        # extract the positions
        posred = (pos[indred, 2:4])[0,:,:]

    # get number of particles inside the zone
    tot = len(posred)
    #tot = len(numpy.where(numpy.logical_and(posred[:,0] <=a, posred[:,1] <=a))[0])

    # compute intensity
    a = 1.
    lam = tot / a**2

    # define the grid
    npoints=100
    ra = 0.01
    rb = numpy.sqrt(1./numpy.pi)

    rgrid = numpy.mgrid[0:npoints+1]/float(npoints)*(rb-ra)+ra
    k = numpy.zeros((npoints+1))

    # compute the pair distance and remove the self-distance indices
    try:
        scipy = importlib.import_module('scipy')
        hasscipy = True        
    except:
        print("Could not find scipy. Skipping all integrating..")

    if hasscipy:
        distr=scipy.spatial.distance.squareform(scipy.spatial.distance.pdist(posred))
        dist = distr[numpy.where(distr > 0)]

    for i,r in enumerate(rgrid):
        # this gives the total number of points which are closer than r
        k[i] = numpy.sum(gpgmp.models.majority.indicatorFunction(dist, r))/tot
    # so this one gives the *average* number of points inside distance r
    #k = tpc/tot

    return rgrid, k/lam

def computeTwoPointCorrelationRed(pos, h=None):

    # check if we only got an x/y array
    if (numpy.shape(pos)[1]==2):
        posred = pos
    else:
        # get red index
        indred=numpy.where(pos[:,-1]==1)

        # extract the positions
        posred = (pos[indred, 2:4])[0,:,:]

    # get length
    tot = len(posred)

    # we normalize the distances to the maximum packing distance
    # we assume it's a unit area so we take
    #rmax = numpy.sqrt(1./(2.*numpy.sqrt(3)*tot))
    #rmax = 1.

    # compute the  intensity
    a = 1.
    lam = tot / a**2
  
    # kernel width estimate
    if h==None:
        h = 0.1/numpy.sqrt(lam)

    # define the  grid
    npoints=100
    rb = h
    ra = a
    rgrid = numpy.mgrid[0:npoints+1]/float(npoints)*(rb-ra)+ra
    g = numpy.zeros((npoints+1))

    # compute the normalized pair distance and remove the self-distance indices
    try:
        scipy = importlib.import_module('scipy')
        hasscipy = True        
    except:
        print("Could not find scipy. Skipping all integrating..")

    if hasscipy:
        distr=scipy.spatial.distance.squareform(scipy.spatial.distance.pdist(posred))
        dist = distr[numpy.where(distr > 0)]

    # compute the isotropized set covariance
    gamma = a**2 - 2.*rgrid/numpy.pi*2.*a+rgrid**2/numpy.pi

    # define the width of the Gaussian kernel
    #sigma = 0.25

    #norm = tot*(tot-1.)/2.
    norm = tot**2

    for i,r in enumerate(rgrid):
        g[i] = numpy.sum(gpgmp.models.majority.boxKernel(r-dist, h))/(2.*numpy.pi*r*gamma[i])/lam**2
        #g[i] = 1./rmax**2/(2.*numpy.pi*r*norm)*numpy.sum(numpy.exp(-(r-dist)**2/sigma**2)/numpy.sqrt(numpy.pi)/sigma)

    # use scipy to get the two-point correlation
    #kdt = scipy.spatial.KDTree(posred)
    #tpc = kdt.count_neighbors(kdt, r)

    return rgrid, g, gamma

def computeExitTimeFPE(datadir, test=False, delta=0.1):
    # computes the exit time by using the homogeneous process equation
    # (5.5.21) from Gardiner
    if test:
        # set standard values for diffusivity and drift so we can compare it with
        # analytic solution
        tot = 150
        v = 1.
        diff = 0.1
        alpha = numpy.ones((tot+1))*v
        beta2 = numpy.ones((tot+1))*2.*diff # factor two cause definition in Gardiner is B/2 d2p/dx2..
        # set up the coordinate system
        x=numpy.mgrid[0:tot+1]/float(tot)

        # theroetical solution [from Farkas (2001)]
        nr = (1.-numpy.exp(-x*v/diff))/(1.-numpy.exp(-v/diff))
        tau = (nr-x)/v
    else:
        # read in coefficients
        alphar = numpy.genfromtxt(datadir+'/alpha.dat')
        beta2r = numpy.genfromtxt(datadir+'/beta2.dat')
        totr = numpy.shape(alphar)[0]-1

        # set up the coordinate system
        xx=numpy.mgrid[0:totr+1]/float(totr)

        # now we need to find the indices which are closest to delta
        id1 = (numpy.abs(xx-delta)).argmin()
        id2 = (numpy.abs(xx-(1.-delta))).argmin()

        # and we just use these values to compute the exit time
        x = xx[id1:id2+1]
        alpha = alphar[id1:id2+1]
        beta2 = beta2r[id1:id2+1]

    # try to import scipy
    try:
        scipy = importlib.import_module('scipy')
        hasscipy = True        
    except:
        print("Could not find scipy. Skipping all integrating..")

    if hasscipy:
        # compute the variuos integrals
        eps = scipy.integrate.cumtrapz(2.*alpha/beta2, x, initial=0.)
        iy = scipy.integrate.cumtrapz(numpy.exp(eps)/beta2, x, initial=0.)
        iax = scipy.integrate.cumtrapz(numpy.exp(-eps)*iy, x, initial=0.)
        ixb = -scipy.integrate.cumtrapz((numpy.exp(-eps))[::-1], x[::-1], initial=0.)[::-1]
        ixb1 = -scipy.integrate.cumtrapz((numpy.exp(-eps)*iy)[::-1], x[::-1], initial=0.)[::-1]
        iax1 = scipy.integrate.cumtrapz(numpy.exp(-eps), x, initial=0.)
        t = 2.*(iax1*ixb1 - ixb*iax)/iax1[-1]

        # for splitting probability we need
        psim1 = scipy.integrate.cumtrapz(numpy.exp(-eps), x, initial=0.)
        pi = psim1/psim1[-1]

    # and save it as *csv (compability with mathematica)
    if test:
        return t, tau
    else:
        numpy.savetxt(datadir+"decision_time.csv", numpy.column_stack([x, t]), delimiter=',')
        numpy.savetxt(datadir+"splitting_probability.csv", numpy.column_stack([x, pi]), delimiter=',')

def computeBirthDeath(datadir):
    # computes the splitting probability and exit time from the
    # transition probabilities under the assumption that the underlying
    # process is a birth-death process

    # read in transition probabilities - these should be prob per time
    # W+/W- in the Gillespie (1992) dictus
    pplus  = numpy.genfromtxt(datadir+'/pplus.dat')
    pminus = numpy.genfromtxt(datadir+'/pminus.dat')
    nagents = numpy.shape(pminus)[0]-1

    # we set the probability to a small value wherever it's zero
    if (numpy.shape(numpy.where(pminus == 0))[1] > 0):
        pminus[numpy.where(pminus == 0)]=1e-5
    if (numpy.shape(numpy.where(pplus == 0))[1] > 0):
        pplus[numpy.where(pplus == 0)]=1e-5
    
    # define helper arrays - see Gillespie 6.7.12./6.7.13
    t1 = numpy.zeros((nagents+1))
    t2 = numpy.zeros((nagents+1))
    p1 = numpy.zeros((nagents+1))
    p2 = numpy.zeros((nagents+1))

    # compute t1,t2
    t1[1]         = 1./pminus[1]
    t2[nagents-1] = 1./pplus[nagents-1]
    for i in range(2, nagents):
        t1[i]         = (1.+pplus[i-1]*t1[i-1])/pminus[i]
        t2[nagents-i] = (1.+pminus[nagents-i+1]*t2[nagents-i+1])/pplus[nagents-i]

    # compute splitting probabilities
    p1 = t2/(t1+t2)
    p2 = t1/(t1+t2)
    p1[0] = 1.
    p2[0] = 0.
    p1[nagents] = 0.
    p2[nagents] = 1.

    # finally compute first-exit time
    dectime = numpy.zeros((nagents+1))
    for i in range(1, nagents):
        dectime[i] = p1[i]*numpy.sum(t1[1:i+1]) + p2[i]*numpy.sum(t2[i+1:nagents])

    if False:
        # -- below is the way via the potential function --
        # we set the probability to a small value wherever it's zero
        if (numpy.shape(numpy.where(pminus == 0))[1] > 0):
            pminus[numpy.where(pminus == 0)]=1e-5
        if (numpy.shape(numpy.where(pplus == 0))[1] > 0):
            pplus[numpy.where(pplus == 0)]=1e-5

        # compute potential function 6.4.12 in Gillespie (1992)
        phi=numpy.cumsum(numpy.log(pminus/pplus)[1:])
        phi=numpy.insert(phi, 0, 0.)
        norm = numpy.sum(numpy.exp(phi[:-1]))

        # compute splitting probs. (6.7.18a) in Gillespie (1992)
        p1 = numpy.cumsum(numpy.exp(phi[::-1][1:]))[::-1]/norm # goes from 0 to N-1 (N values)
        p2 = numpy.cumsum(numpy.exp(phi))[:-1]/norm # goes from 1 to N (N values)
        pi = numpy.insert(p2, 0, 0.) # actual spliting probability .. has 0 at the beginning

        # compute decision time (6.7.19a) in Gillespie (1992)
        # yeargh.. use recursive function might be easier!
        #dectime1 = numpy.cumsum((p2*norm*numpy.exp(-phi)/pplus)[1:])*p1[1:] # goes from 1 to N-1 (N-1 values)
        #dectime2 = numpy.cumsum((p1*norm*numpy.exp(-phi)/pplus)[::-1][1:])[::-1][1:]*p2[:-2] 
        #dectime=dectime1[:-2]+dectime2[1:]
        dc1con = numpy.zeros((9)) # goes from 1 to N-1 (N-1 values)
        dc2con = numpy.zeros((8)) # goes from 1 to N-2 (N-2 values)
        for i in range(1,nagents):
            dc1con[i-1] = numpy.sum(((numpy.exp(-phi)/pplus)[1:-1]*(p2*norm)[:-1]*p1[1:])[0:i])
        for i in range(1,nagents-1):
            dc2con[i-1] = numpy.sum(((numpy.exp(-phi)/pplus)[1:-2]*(p1*norm)[1:-1]*p2[:-2])[i+1:])

        dc = numpy.zeros((nagents+1))
        dc[1:nagents] = dc1con
        dc[1:nagents-1] += dc2con

    # and save it
    numpy.savetxt(datadir+'/splitting_probability.dat', p2)
    numpy.savetxt(datadir+'/decision_time.dat', dectime)

def makeComparisonPlot(dirs, resultsdir, useFPE=True):
        # and plot them
    try:
        matplotlib = importlib.import_module('matplotlib')
        #matplotlib.use('PDF')
        plt = importlib.import_module('matplotlib.pyplot')
        hasmatplotlib = True        
    except:
        print("Could not find matplotlib. Skipping all plots..")

    if hasmatplotlib:
        # Make SDE comparison
        matplotlib.rc('text', usetex=True)
        plt.figure(figsize=(4*2.*1.5, 3*2*1.5), dpi=600)
        plt.subplot(221)
        plt.xlabel(r'$s$')
        if useFPE:
            plt.ylabel(r'$\alpha(s)/10^{-4}$')
        else:
            plt.ylabel(r'$W_+(s)$')
        plt.subplot(222)
        plt.xlabel(r'$s$')
        if useFPE:
            plt.ylabel(r'$\beta^2(s)/10^{-4}$')
        else:
            plt.ylabel(r'$W_-(s)$')
        plt.subplot(223)
        plt.xlabel(r'$s$')
        plt.ylabel(r'$\pi(s)$')
        plt.subplot(224)
        plt.xlabel(r'$s$')
        plt.ylabel(r'$T(s)$')

        # go through datadirs
        colorcodes = ['b', 'r', 'g', 'k']
        for i, datadir in enumerate(dirs):

            # read in parameters
            alpha = numpy.genfromtxt(datadir+'/alpha.dat')
            beta2 = numpy.genfromtxt(datadir+'/beta2.dat')
    
            # bin position for the symmetry parameters
            tot = numpy.shape(alpha)[0]-1
            si = numpy.asfarray(numpy.mgrid[0:tot+1])/tot

            plt.subplot(221)
            if useFPE:
                plt.plot(si, alpha/1e-4, colorcodes[i]+'-')
            else:
                pplus  = numpy.genfromtxt(datadir+'/pplus.dat')
                plt.plot(si, pplus, colorcodes[i]+'-')
            plt.xlim([0.,1.])

            plt.subplot(222)
            if useFPE:
                plt.plot(si, beta2/1e-4, colorcodes[i]+'-')
            else:
                pminus = numpy.genfromtxt(datadir+'/pminus.dat')
                plt.plot(si, pminus, colorcodes[i]+'-')
            plt.xlim([0.,1.])

            # read in splitting probability
            sp = numpy.genfromtxt(datadir+'/splitting_probability_experiment.dat')
    
            plt.subplot(223)
            if useFPE:
                se = numpy.genfromtxt(datadir+'/splitting_probability.csv', delimiter=',')
                plt.plot(se[:,0], se[:,1], colorcodes[i]+'-')
                plt.plot(si, sp, colorcodes[i]+'+')
            else:
                se = numpy.genfromtxt(datadir+'/splitting_probability.dat')
                plt.plot(si, se, colorcodes[i]+'-')
                plt.plot(si, sp, colorcodes[i]+'o')                

            # read in decision time
            decTimeData = numpy.genfromtxt(datadir+'/decision_time_experiment.dat')

            plt.subplot(224)
            if useFPE:
                decTimeInt = numpy.genfromtxt(datadir+'/decision_time.csv', delimiter=',')
                plt.plot(decTimeInt[:,0], decTimeInt[:,1], colorcodes[i]+'-')
                plt.plot(si, decTimeData, colorcodes[i]+'+')
            else:
                decTimeInt = numpy.genfromtxt(datadir+'/decision_time.dat')
                plt.plot(si, decTimeInt, colorcodes[i]+'-')
                plt.plot(si, decTimeData, colorcodes[i]+'o')

        plt.savefig(resultsdir+'/basic_comparison.pdf')

        # delete plots
        plt.clf()

        # turn of latex support
        matplotlib.rc('text', usetex=False)

def analyzeGeoff(options, datadir):
    # this contains all directories we reduce
    dirs=[]

    # go through all directories and reduce the files in there
    for dirpath, dirnames, files in os.walk(datadir):
        for name in dirnames:
            rdir = datadir+'/'+name+'/'
            dirs.append(rdir)
            print("Reducing dir {0:s}".format(rdir))
            analyzeGeoffSingle(rdir, options)

    # now we do the comparison plot   
    makeComparisonPlot(sorted(dirs), datadir, False)

def analyzeGeoffSingle(datadir, options=None):

    # get standard options if needed
    if options==None:
        options, op = common.options.getSystemOptions(argv=None, parser=None)

    # this will contain the change in symmetry parameter for each channel
    nrobots = 10
    deltas  = numpy.zeros((nrobots+1))
    deltass = numpy.zeros((nrobots+1))
    pp  = numpy.zeros((nrobots+1))
    pm = numpy.zeros((nrobots+1))
    cs  = numpy.zeros((nrobots+1)) # counts the number of values we have in each bin so we can compute the average later
    ds  = numpy.zeros((nrobots+1)) # counts the number of values we have in each bin when a decision was made
    tt = numpy.zeros((nrobots+1)) # counts the sum of the deltat for averaging

    # this one will contain the decision time
    timedec = numpy.zeros((nrobots+1))

    # and this one is for the splitting probability
    splitting = numpy.zeros((nrobots+1))

    # go through directory and reduce all 'txt' files
    for dirpath, dirnames, files in os.walk(datadir):
        for file in files:
            if file.endswith('.txt'):
                f=open(dirpath+file,'r')
                #print "Reducing file {0:s}".format(f)

                start = True
                # read in each line
                for line in f:
                    #print "Reducing line {0:s}.".format(line)
                    if start:
                        start = False;
                        exited = False;
                        rho0=int(line.split(' ')[0])
                        rho = rho0
                        mins=int(line.split(' ')[1].split(':')[0])
                        secs=int(line.split(' ')[1].split(':')[1].strip())
                        currTime=0
                        currTime=60*mins+secs
                        t0 = currTime
                        #print "Experiment started with {0:d} red robots at t={1:g}".format(rho0, currTime)
                    else:
                        mins=int(line.split(':')[0])
                        secs=int(line.split(':')[1].split(  ' ')[0])
                        op=line.split(' ')[1].strip()
                        time=60*mins+secs
                        deltat = time-currTime
                        if deltat==0:
                            deltat = 1
                        currTime = time
                        if op == '+':
                            deltas[rho]  += 0.1/deltat
                            deltass[rho] += 0.1**2/deltat
                            pp[rho] += 1.
                            cs[rho] += 1
                            tt[rho] += deltat
                            rho += 1
                        elif op == '-':
                            deltas[rho]  -= 0.1/deltat
                            deltass[rho] += 0.1**2/deltat
                            pm[rho] += 1.
                            cs[rho] += 1
                            tt[rho] += deltat
                            rho -=1
                        elif op == 's':
                            # not sure what this is for
                            start = True
                        if ((rho == 10 or rho==0) and (not exited)):
                            #print "Exit at {0:g} (for s0={1:d})".format(currTime, rho0)
                            timedec[rho0] += currTime 
                            splitting[rho0] += float(rho)/10.
                            ds[rho0]+=1
                            exited = True
                        #print "We currently have {0:d} robots at t={1:g}.".format(rho, currTime)
                f.close()

    # save the arrays
    numpy.savetxt(datadir+'/alpha.dat', deltas/cs)
    numpy.savetxt(datadir+'/beta2.dat', deltass/cs)
    numpy.savetxt(datadir+'/decision_time_experiment.dat', timedec/ds)
    numpy.savetxt(datadir+'/splitting_probability_experiment.dat', splitting/ds)
    numpy.savetxt(datadir+'/pplus.dat', pp/tt)
    numpy.savetxt(datadir+'/pminus.dat', pm/tt)

    # compute splitting probability and decision time from birth-death master equation
    computeBirthDeath(datadir)

    # now we need to run the mathematica script
    if False:
        mathparam = "-script {0:s}/mk_split_probabilities_picture.m {1:s} {2:s}".format(options['majoritymathscripts'], datadir+'/split_probabilities.pdf', datadir+'/')

        # run the mathematica script
        mathcommand = options['mathematicabin']+' '+mathparam
        print("Running mathematica script using {0:s}".format(mathcommand))
        
        # create files to capture output
        outfile = open(datadir+'/mathematica_out.txt', 'w')
        errfile = open(datadir+'/mathematica_err.txt', 'w')
    
        # and run it
        subprocess.call(mathcommand, stdout=outfile, stderr=errfile, shell=True)

def computeDecisionTime(rho, t, nagents, delta=0.1):
    # extract number of runs and dumps
    nruns = numpy.shape(rho)[0]
    ndumps = numpy.shape(t)[0]
    
    # duplicate times array
    tm = numpy.reshape(numpy.tile(t, nruns), (nruns, ndumps))

    # now pick the times where no decision was made
    di = numpy.where(numpy.logical_and(rho > delta, rho < 1-delta))

    # these will be set to infinity in the array
    tm[di] = numpy.inf

    # and now we pick the minimum for each experiment
    mtime = numpy.min(tm,1)

    # mtime now contains the decision time for each experiment
    # we still have to sort them according to the experiment start condition
    # create edges
    nbins = nagents+1
    edges = getEdges(nbins)
    
    # create output array
    ret = numpy.zeros((nbins))

    # for each bin sum up the decision time and compute the average
    init = rho[:,0]
    for i in range(nbins):
        index = numpy.where(numpy.logical_and(numpy.logical_and(init>edges[i], init<edges[i+1]), mtime < numpy.inf))[0]
        si = numpy.shape(index)[0]
        if (si>0):
            ret[i] = numpy.sum(mtime[index])/float(si)
            if (ret[i] == numpy.inf):
                ret[i] = numpy.nan
        else:
            ret[i] = numpy.nan

    return ret

def plotComparison(datadir):
    # and plot them
    try:
        matplotlib = importlib.import_module('matplotlib')
        #matplotlib.use('PDF')
        plt = importlib.import_module('matplotlib.pyplot')
        hasmatplotlib = True        
    except:
        print("Could not find matplotlib. Skipping all plots..")

    if hasmatplotlib:
        # prepare plot
        matplotlib.rc('text', usetex=True)
        plt.figure(figsize=(4*2., 3*1.), dpi=600)

        # read in pstar
        pstarexp = numpy.genfromtxt(datadir+'/pstar.dat')
        pstarnum = numpy.genfromtxt(datadir+'/pstar.csv', delimiter=',')

        # and plot it
        plt.subplot(121)
        plt.plot(pstarexp[:,0], pstarexp[:,1], 'r+')
        plt.plot(pstarnum[:,0], pstarnum[:,1], 'b-')
        plt.ylim([0,1])
        plt.xlabel(r'$s$')
        plt.ylabel(r'$p^\ast(s)$')

        # read in experimental splitting probability
        sp = numpy.genfromtxt(datadir+'/splitting_probability_experiment.dat')

        # create grid for experimental
        tot = numpy.shape(sp)[0] - 1
        si = numpy.asfarray(numpy.mgrid[0:tot+1])/tot

        # read in theoretical splitting probability
        se = numpy.genfromtxt(datadir+'/splitting_probability.csv', delimiter=',')

        # and plot splitting probability
        plt.subplot(122)
        plt.plot(se[:,0], se[:,1], "b-")
        plt.plot(si, sp, "r+")
        plt.xlabel(r'$s$')
        plt.ylabel(r'$\pi(s)$')

        # and save it
        plt.savefig(datadir+'/pstar_pi_comparison.pdf')

        # delete plots
        plt.clf()

        # turn of latex support
        matplotlib.rc('text', usetex=False)

def computeSplittingProbability(init, end, nagents):
    # create edges
    nbins = nagents+1
    edges = getEdges(nbins)
    
    # create output array
    ret = numpy.zeros((nbins))

    # for each bin sum up the end results and comptue the average
    for i in range(nbins):
        index = numpy.where(numpy.logical_and(init>edges[i], init<edges[i+1]))[0]
        si = numpy.shape(index)[0]
        if (si>0):
            ret[i] = numpy.sum(end[index])/float(si)
        else:
            ret[i] = numpy.nan

    # and return it
    return ret

def makePDFEvolutionPlot(datadir, times):
    # read in data for PDF evolution from experiment
    pdf = numpy.genfromtxt(datadir+'/evolutionPDF.dat')
    x = numpy.mgrid[0:151]/150.
    t = numpy.genfromtxt(datadir+'/times.dat')

    # read in PDF evolution data from numerical integration
    intpdf = numpy.genfromtxt(datadir+'/fpe_solution.dat', delimiter=',')

    # find the indices that match the required times
    tindex = []
    for tt in times:
        tindex.append(numpy.where(t == tt)[0][0])

    # dynamically import matplotlib (if possible)
    hasmatplotlib = False

    try:
        matplotlib = importlib.import_module('matplotlib')
        #matplotlib.use('PDF')
        plt = importlib.import_module('matplotlib.pyplot')
        hasmatplotlib = True
    except:
        print("Could not find matplotlib. Skipping all plots..")

    if hasmatplotlib:
        # plot simulation results and overplot integrated values
        plt.plot(x, numpy.transpose(pdf[tindex]), "r-")
        plt.plot(x, numpy.transpose(intpdf[tindex]), "b-")
        plt.ylim([0, 15])

        # save ut
        plt.savefig(datadir+'/pdf_comparison.pdf')
        plt.clf()

def mkDecisionErrorPlot(errorfiles, recError):
    # dynamically import matplotlib (if possible)
    hasmatplotlib = False

    try:
        matplotlib = importlib.import_module('matplotlib')
        #matplotlib.use('PDF')
        plt = importlib.import_module('matplotlib.pyplot')
        hasmatplotlib = True
    except:
        print("Could not find matplotlib. Skipping all plots..")

    if hasmatplotlib:
        # prepare plot
        matplotlib.rc('text', usetex=True)
        plt.figure(figsize=(4*2, 3*2))

        # read in all files and plot them
        for ef in errorfiles:
            e = numpy.genfromtxt(ef)
            # check if all data points are there .. if not we assume it's the first x ones
            np = numpy.shape(e)[0]
            plt.plot(recError[0:np], e/1e-2, 'o')

        plt.xlabel(r'$r$')
        plt.ylabel(r'$\epsilon/10^{-2}$')

        # and save it
        plt.savefig('recognitionError.pdf')
        plt.clf()

        # turn of latex support
        matplotlib.rc('text', usetex=False)

def getEdges(nbins):
    dx = 1./float(nbins)
    edges = numpy.mgrid[0:nbins+1]*dx
    return edges

def getEvolutionPDFRho(rho, times, nagents):

    # create edges
    nbins = nagents+1
    edges = getEdges(nbins)
    
    # create output array
    ndumps = numpy.shape(times)[0]
    ret = numpy.zeros((ndumps, nbins))
    
    # create histogram for all
    for i in range(ndumps):
        # check if sum is zero
        if (numpy.sum(rho[:,i])==0):
            print("Warning: getEvolutionPDFRho finds zero count for i={0:d}.".format(i))
        
        hist, tedges = numpy.histogram(rho[:,i], bins=edges, density=True)
        ret[i,:] = hist
        
    return ret, edges

def getStationaryRho(rho, times, nagents=150):
    # get time indices for threshold time
    threshold = 1000.
    indices = numpy.where(times>threshold)

    # create edges for histogram
    nbins = nagents+1
    edges = getEdges(nbins)

    # create histogram for indices
    nindices = numpy.shape(indices)[0]
    h = numpy.zeros((nindices, nbins))
    for i in range(nindices):
        # check if sum is zero
        if (numpy.sum(rho[:,indices[i]])==0):
            print("Warning: getStationaryRho finds zero count for i={0:d}.".format(i))
        hist, tedges = numpy.histogram(rho[:,indices[i]], bins=edges, density=True)
        h[i,:] = hist

    # and average over i
    ret = numpy.mean(h, axis=0)

    # get middle values of edges
    em = (edges[1:]-edges[:-1])/2. + edges[:-1]

    return ret, em
    
def computeFptHomogeneous(rho, deltat, n, t0=0):
    # we assume time-homogeneity as in Geoffrey's thesis
    # we compute rho and delta rho where we omit the first time value    
    r = (rho[:,t0:-1]).flatten()
    dr = (rho[:,t0+1:]-rho[:,t0:-1]).flatten()
    
    # now get the relative and absolute changes
    sr, sa, sr2, pplus, pminus, p0 = computeChange(r, dr, deltat, n)
    
    # and compute alpha and beta
    alpha = sr
    beta = sa - numpy.abs(sr)
    beta2 = sr2
    
    return alpha, beta, beta2, pplus, pminus, p0
    
def computeChange(rho, drho, dt, n):
    # we assume that rho and drho are flat arrays
    # this might be because it's one time slice or because
    # we assume it's time-homogeneous (as in Geoffreys's thesis)
    
    # compute relative and absolute change
    srel=drho/dt
    srel2=drho**2/dt
    sabs=numpy.abs(srel)

    # now for each value of rho el (0, 1/N, 2/N, .., 1) we find all corresponding values
    # of sabs and srel and average them (see ref [13] from Bernd's paper)
    srelnew = numpy.zeros((n+1))
    srel2new = numpy.zeros((n+1))
    sabsnew = numpy.zeros((n+1))
    pplus = numpy.zeros((n+1))
    pminus = numpy.zeros((n+1))
    p0 = numpy.zeros((n+1))
    
    # we use the same edges as the histogram
    edges = getEdges(n+1)
        
    for i in range(n+1):
        index = numpy.where(numpy.logical_and(rho>=edges[i], rho<edges[i+1]))
        
        if (numpy.shape(index)[1]>0):
            srelnew[i] = numpy.mean(srel[index])
            srel2new[i] = numpy.mean(srel2[index])
            sabsnew[i] = numpy.mean(sabs[index])

            # now we compute the change probability
            srt = srel[index]
            pplus[i] = float(numpy.sum(srt[numpy.where(srt>0)]*n))/float(numpy.shape(index)[1])
            pminus[i] = -float(numpy.sum(srt[numpy.where(srt<0)]*n))/float(numpy.shape(index)[1])
            p0[i] = float(numpy.sum(srt[numpy.where(srt==0)]*n))/float(numpy.shape(index)[1])
            #pplus[i] = float(numpy.shape()[1])/float(numpy.shape(index)[1])
            #pminus[i] = float(numpy.shape(numpy.where(srel[index]<0))[1])/float(numpy.shape(index)[1])
            #p0[i] = float(numpy.shape(numpy.where(srel[index]==0))[1])/float(numpy.shape(index)[1])
            
    return srelnew, sabsnew, srel2new, pplus, pminus, p0

def reducePBS(resultsdir, dirs=None):

    if dirs==None:
        # we first read in the directory list
        pkl_file = open(resultsdir+'parameters.pkl', 'rb')
        dirs = pickle.load(pkl_file)
        pkl_file.close()

    # browse through directory
    for cdir in dirs:
        # we change into this dir (so we can use tar)
        olddir = os.getcwd()
        currentDir = resultsdir+'/'+cdir
        os.chdir(currentDir)

        # get parameters for this run
        paramfile = open(currentDir+'parameters.pkl', 'r')
        nruns = pickle.load(paramfile)
        tot = pickle.load(paramfile)
        paramfile.close()
        
        # loop through runs
        print("Cleaning up dir {0:s} with {1:d} runs.".format(currentDir, nruns))
        for i in range(nruns):
            if os.path.isfile(currentDir+'/rho_{0:d}.dat'.format(i)):
                # read in rho
                rho = numpy.loadtxt(currentDir+'rho_{0:d}.dat'.format(i))
                nmax = 5000 # TODO: This should be read it from the parameter file!
                time = numpy.mgrid[0:nmax]
                numpy.savez(currentDir+'rho_{0:d}.npz'.format(i), rho=rho, time=time)
                os.remove(currentDir+"rho_{0:d}.dat".format(i))

        # create tars for log files and delete originals
        if (not os.path.isfile(currentDir+'/logs.tar.gz')):
            subprocess.call("tar czf {0:s} {1:s}".format("logs.tar.gz", "run_pbs_*"), shell=True)
            subprocess.call("rm -rf run_pbs_*", shell=True)
        if (not os.path.isfile(currentDir+'/models.tar.gz')):
            subprocess.call("tar czf {0:s} {1:s}".format("models.tar.gz", "model_*"), shell=True)
            subprocess.call("rm -rf model_*", shell=True)
    
    # change into previous dir
    os.chdir(olddir)

def reduceHeiko(resultsdir):
    # read in results file
    resultsRaw = numpy.genfromtxt('ratioTimeDump.dat')

    # and reshape it
    results=numpy.reshape(resultsRaw, ((100, 10001, 5)))

    # compute rho and t
    ndumps = numpy.shape(results)[1]
    t = numpy.mgrid[0:ndumps]
    rho = results[:,:,4]

    # and save it
    numpy.savez(resultsdir+'/rho.npz', rho=rho, time=t)

    # read in change history file
    nbins = 101
    changeRaw = numpy.genfromtxt(resultsdir+'/changeHisto')
    ntimedumps = numpy.shape(changeRaw)[0]/nbins
    change=numpy.reshape(changeRaw, (-1, nbins, 6))
    print("Found {0:d} time dumps.".format(ntimedumps))

    # members of change:
    # - first index: time (accumulated over dt=fr_chg_delta_t)
    # - second index: s-bin (0, 0.01, .., 1.)
    # - third index: [t, t+dt, bin, abs, rel, count]
    # .. abs/rel are averaged over samples and time bins => [abs]=[rel]=1/s
    # .. count is not .. count is ntimesteps per time bin * samples

    # compute mean .. for consistency only use the numbers where count > 0
    crelm = numpy.zeros((nbins))
    cabsm = numpy.zeros((nbins))
    for i in range(nbins):
        nzindex = numpy.where(change[:,i,5]>0)
        crelm[i] = numpy.mean(change[nzindex,i,4])
        cabsm[i] = numpy.mean(change[nzindex,i,3])

    return rho, t, change, crelm, cabsm

def analyze(options, resultsdir):

    # check if the totals file exists:
    if os.path.isfile(resultsdir+'/rho.npz'):
        print("Found total file {0:s} and will use it.".format(resultsdir+'/rho.npz'))
        
        # read in file
        df = numpy.load(resultsdir+'/rho.npz')
        #nr = df['nreds']
        #ng = df['ngreens']
        rho = df['rho']
        t = df['time']
    else:
        # we need to compute the totals first
        print("No total file found.. computing.")
        
        # read in totals
        nru, ngu, tu = gpgmp.models.majority.read_totals(resultsdir+'output.h5')
    
        # we need to time-order them
        targs = numpy.argsort(tu[0,:])
        nr = nru[:, targs]
        ng = ngu[:, targs]
        t = tu[0, targs]
        
        # now get the total number of agents and compute rho
        tot = int(nr[0,0]+ng[0,0])
        rho = nr/tot

        # save totals
        numpy.savez(resultsdir+'/rho.npz', rho=rho, time=t)

    # read in parameters
    paramfile = open(resultsdir+'parameters.pkl', 'r')
    nruns = pickle.load(paramfile)
    tot = pickle.load(paramfile)
    paramfile.close()

    # check if we should use FPE or Master equation
    useFPE = False
    if (tot>100):
        print("Found {0:d} agents. Using FPE for integration.".format(tot))
        useFPE = True
    else:
        print("Found {0:d} agents. Using Master equation for integration.".format(tot))
    # dynamically import matplotlib (if possible)
    hasmatplotlib = False

    try:
        matplotlib = importlib.import_module('matplotlib')
        #matplotlib.use('PDF')
        plt = importlib.import_module('matplotlib.pyplot')
        hasmatplotlib = True
    except:
        print("Could not find matplotlib. Skipping all plots..")

    #import matplotlib.pyplot as plt
    #import matplotlib.cm as cm
    #import matplotlib

    # compute stationary PDF
    rs, em = gpgmp.models.majority.getStationaryRho(rho, t, tot)
    rhostar =  numpy.column_stack((em, rs))
    
    # and save it
    numpy.savetxt(resultsdir+'/pstar.dat', rhostar)

    # export PDF
    ev, edges =gpgmp.models.majority.getEvolutionPDFRho(rho, t, tot)
    numpy.savetxt(resultsdir+'/evolutionPDF.dat', ev)
    numpy.savetxt(resultsdir+'/times.dat', t)

    # compute splitting probability and save it
    init=rho[:,0]
    end=rho[:,-1]
    sp = gpgmp.models.majority.computeSplittingProbability(init, end, tot)
    numpy.savetxt(resultsdir+'/splitting_probability_experiment.dat', sp)

    if hasmatplotlib:
        # prepare plot
        matplotlib.rc('text', usetex=True)
        plt.figure(figsize=(4, 3), dpi=300)
    
        # plot evolution of PDF
        plt.imshow(numpy.transpose(ev), extent=[0., numpy.max(t), edges[0], edges[-1]], aspect='auto')
        plt.xlabel(r'$t\ \mathrm{[s]}$')
        plt.ylabel(r'$s$')
        plt.savefig(resultsdir+'/pdfEvolution.pdf')
        plt.clf()

    # compute dt
    dt = t[1]-t[0]
    
    # now get alpha and beta assuming time-homogenity
    alpha, beta, beta2, pplus, pminus, p0 = gpgmp.models.majority.computeFptHomogeneous(rho, dt, tot)
    
    # we export them
    numpy.savetxt(resultsdir+'/alpha.dat', alpha)
    numpy.savetxt(resultsdir+'/beta.dat', beta)
    numpy.savetxt(resultsdir+'/beta2.dat', beta2)
    numpy.savetxt(resultsdir+'/pplus.dat', pplus)
    numpy.savetxt(resultsdir+'/pminus.dat', pminus)
    numpy.savetxt(resultsdir+'/p0.dat', p0)

    # start the MATLAB/Mathematica scripts
    if useFPE:
        # run matlab script to integrate FPE
        matcom = '{0:s} -nodesktop -nosplash -r "cd \'{1:s}\'; solveFPE(\'{2:s}\'); exit "'.format(options['matlabbin'], options['majoritymathscripts'], resultsdir)
        print("Running matlab script using the command {0:s}".format(matcom))
        
        # create files to capture output
        outfile = open(resultsdir+'/matlab_out.txt', 'w')
        errfile = open(resultsdir+'/matlab_err.txt', 'w')
    
        # and run it
        subprocess.call(matcom, stdout=outfile, stderr=errfile, shell=True)

        # plot PDF evolution comparison
        makePDFEvolutionPlot(resultsdir, [0, 300, 500])

        # and compute the exit time
        computeExitTimeFPE(resultsdir, False)

        # this is the parameter string for the mathematica script
        #mathparam = "-script {0:s}/mk_split_probabilities_picture.m {1:s} {2:s}".format(options['majoritymathscripts'], resultsdir+'/split_probabilities.pdf', resultsdir+'/')

        # Run the mathematica script
        #mathcommand = options['mathematicabin']+' '+mathparam
        #print "Running mathematica script using {0:s}".format(mathcommand)
    
        # create files to capture output
        #outfile = open(resultsdir+'/mathematica_out.txt', 'w')
        #errfile = open(resultsdir+'/mathematica_err.txt', 'w')
    
        # and run it
        #subprocess.call(mathcommand, stdout=outfile, stderr=errfile, shell=True)
        
    else:
        computeBirthDeath(resultsdir)

    # bin position for the symmetry parameters
    si = numpy.asfarray(numpy.mgrid[0:tot+1])/tot
    
    # compute decision time from data
    if useFPE:
        decTimeData = computeDecisionTime(rho, numpy.float_(t), tot)
    else:
        # for this we set delta to zero ..
        decTimeData = computeDecisionTime(rho, numpy.float_(t), tot, 0)

    # and save it
    numpy.savetxt(resultsdir+'/decision_time_experiment.dat', decTimeData)

    # and read it in from mathematica output
    if useFPE:
        decTimeInt = numpy.genfromtxt(resultsdir+'/decision_time.csv', delimiter=',')

        if hasmatplotlib:
            # plot coefficients
            plt.figure(figsize=(4*4,3*2), dpi=300)
    
            plt.subplot(121)
            plt.plot(si, alpha/1e-4)
            plt.xlabel(r'$s$')
            plt.ylabel(r'$\alpha(s)/10^{-4}$')
        
            plt.subplot(122)
            plt.plot(si, beta2/1e-4)
            plt.xlabel(r'$s$')
            plt.ylabel(r'$\beta^2(s)/10^{-4}$')
            plt.savefig(resultsdir+'/sde.pdf')

            # plot decision time
            plt.figure(figsize=(4*1.5,3*1.5), dpi=600)
            plt.plot(si, decTimeData, "r+")
            plt.plot(decTimeInt[:,0], decTimeInt[:,1], "b-")
            plt.xlabel(r'$s$')
            plt.ylabel(r'$T(s)$')
            plt.savefig(resultsdir+'/decision_time.pdf')

            # turn of latex support
            matplotlib.rc('text', usetex=False)

def writeFlameInputFile(filename, tot, sigma=0, recError=0, outputid=-1, speed=0.01, avoidance=0.01, length=1., track_positions=False):

    # seed RNG
    random.seed();
    
    # open input file
    inputFile = open(filename, 'w')
    
    # write header
    inputFile.write("<states>\n")
    inputFile.write("  <itno>0</itno>\n")
    inputFile.write("  <environment>\n")
    #inputFile.write("    <diffusivity>2.5e-3</diffusivity>\n")
    inputFile.write("    <speed>{0:g}</speed>\n".format(speed))
    inputFile.write("    <size_x>{0:g}</size_x>\n".format(length))
    inputFile.write("    <size_y>{0:g}</size_y>\n".format(length))
    inputFile.write("    <avoidance_radius>{0:g}</avoidance_radius>\n".format(avoidance))
    inputFile.write("    <sigma>{0:g}</sigma>\n".format(sigma))
    inputFile.write("    <recognition_error>{0:g}</recognition_error>\n".format(recError))
    inputFile.write("    <output_id>{0:d}</output_id>\n".format(outputid))
    if track_positions:
        inputFile.write("    <track_positions>1</track_positions>\n")
    else:
        inputFile.write("    <track_positions>0</track_positions>\n")

    inputFile.write("  </environment>\n")

    
    # write agent data
    inputFile.write("  <agents>\n")
    for i in range(tot):
        inputFile.write("    <xagent>\n")
        inputFile.write("      <name>agent</name>\n")
        inputFile.write("      <my_id>{0:d}</my_id>\n".format(i))
        inputFile.write("      <my_position>{{ {0:g}, {1:g} }}</my_position>\n".format(random.uniform(0.001, length - 0.001 ), random.uniform(0.001, length-0.001)))
        inputFile.write("      <my_velocity>{0:g}</my_velocity>\n".format(random.uniform(0, 2.*math.pi)))
        color = 0
        if (random.random() < 0.5):
             color = 1
        inputFile.write("      <my_color>{0:d}</my_color>\n".format(color))
        inputFile.write("      <my_green_encounters>0</my_green_encounters>\n")
        inputFile.write("      <my_red_encounters>0</my_red_encounters>\n")
        inputFile.write("      <my_encounter_ids>{-1, -1, -1, -1, -1}</my_encounter_ids>\n")
        inputFile.write("    </xagent>\n")
        
    # close all groups
    inputFile.write("</agents>\n")
    inputFile.write("</states>\n")
    inputFile.close()

def runFlameAvoidance(options, datadir, diffusive, track_positions):
    # create list of dirs for pickling
    dirs=[]
    legend = []

    # loop over avoidance radius
    for radius in  [0.01, 0.05, 0.075, 0.1]:
        # compute actual avoidance radius
        #radius = 0.01*float(avoidance)/10.+0.01

        # create dir
        curdir = datadir+'/radius_{0:g}/'.format(radius)
        dirs.append('/radius_{0:g}/'.format(radius))
        os.makedirs(curdir)
        runFlameSingle(options, curdir, diffusive, 150, 5000, 0, 0, avoidance=radius, track_positions=track_positions)
        #runFlameSingle(options, curdir, diffusive, 150, 1, 0, 0, avoidance=radius, track_positions=track_positions)
        legend.append(r'$r_a={0:g}$'.format(radius))

    # and pickle the directory list
    paramfile = open(datadir+'parameters.pkl', 'w')
    pickle.dump(dirs, paramfile)
    pickle.dump(legend, paramfile)
    paramfile.close()
    
def runFlameVelocity(options, datadir, diffusive):
    # create list of dirs for pickling
    dirs=[]
    legend = []

    # loop over sigma
    for sigma in range(10):
        # create dir
        curdir = datadir+'/sigma_{0:g}/'.format(sigma)
        dirs.append('/sigma_{0:g}/'.format(sigma))
        os.makedirs(curdir)
        runFlameSingle(options, curdir, diffusive, 150, 1000, sigma, 0)
        legend.append(r'$\sigma={0:g}$'.format(sigma))

    # and pickle the directory list
    paramfile = open(datadir+'parameters.pkl', 'w')
    pickle.dump(dirs, paramfile)
    pickle.dump(legend, paramfile)
    paramfile.close()

def runFlameRecognition(options, datadir, diffusive):
    # create list of dirs for pickling
    dirs=[]
    legend = []
    
    # loop over sigma
    for error in range(10):
        recError = float(error)/10.*0.2

        # create dir
        curdir = datadir+'/error_{0:g}/'.format(recError)
        dirs.append(curdir)
        os.makedirs(curdir)
        runFlameSingle(options, curdir, diffusive, 150, 1000, 0, recError)
        legend.append(r'$r={0:g}$'.format(recError))

    # and pickle the directory list
    paramfile = open(datadir+'parameters.pkl', 'w')
    pickle.dump(dirs, paramfile)
    pickle.dump(legend, paramfile)
    paramfile.close()

def runFlameExperiment(options, datadir, diffusive):
    # create list of dirs and legend for pickling
    dirs=[]
    legend = []

    # non-looping parameters
    tot = 10
    nruns = 1000
    sigma = 0
    recError = 0
    speed = 1. # in cm s^-1 // should be 1 cm s^-1 .. but we change the time base by this
    length = 60. # in cm
    nmax = 2000
    radiuslist = [8., 6., 4.]  # in cm. 8 cm corresponds to 127

    # loop over avoidance radius
    for radius in range(3):
        avoidance = radiuslist[radius]

        # create dir
        curdir = datadir+'/radius_{0:g}/'.format(avoidance)
        dirs.append('/radius_{0:g}/'.format(avoidance))
        legend.append(r'$d_a={0:g}$'.format(avoidance))
        os.makedirs(curdir)
        
        runFlameSingle(options, curdir, diffusive, tot, nruns, sigma, recError, speed, avoidance, length, nmax)

    # and pickle the directory list
    paramfile = open(datadir+'parameters.pkl', 'w')
    pickle.dump(dirs, paramfile)
    pickle.dump(legend, paramfile)
    paramfile.close()

def analyzeToyModel(options, datadir):
    # get data
    pkl_file = open(datadir+'parameters.pkl', 'rb')
    tot = pickle.load(pkl_file)
    nred = pickle.load(pkl_file)
    k = pickle.load(pkl_file)
    pkl_file.close()

    # reduce the data
    analyzeGpgmpNondiffusiveSingle(options, datadir, toyModel=True)

    # read in plus and minus probabilities
    pplus  = numpy.genfromtxt(datadir+'/pplus.dat')
    pminus  = numpy.genfromtxt(datadir+'/pminus.dat')

    # generate xgrid
    #tot=numpy.shape(pplus)[0] - 1
    s = numpy.mgrid[0:tot+1]/float(tot)

    # compute analytic solution from the toy model
    gamma = numpy.sqrt(1.-3.*(s-2.)*s)
    rho = numpy.sqrt(4.-3.*s*s)
    pmt = k*(tot*(-(tot*(-1. + 2.*s)*(-3. + s + gamma)*(-2. + s + rho)) + 4.*(-2. + rho + s*(4. - s*gamma + (-2. + s)*rho))))/(-4. + 8.*s)
    ppt = k*(tot*(-(tot*(-1. + 2.*s)*(1. + s - gamma)*(2. + s - rho)) + 4.*(-2. + rho + s*(4. - s*gamma + (-2. + s)*rho))))/(-4. + 8.*s)
    
    # set to limit solution at s=0.5
    s05 = numpy.where(s==0.5)[0][0]
    pmt[s05] = k*(0.0527756 *tot* (-0.204417 + 1. *tot))
    ppt[s05] = k*(0.0527756 *tot* (-0.204417 + 1. *tot))

    # read in mesoscopic solution
    ppmeso = numpy.genfromtxt("/home/matthias/data/models/majority/gpgmp/virtual/avoidance_radius/results/radius_0.01/pplus.dat")
    pmmeso = numpy.genfromtxt("/home/matthias/data/models/majority/gpgmp/virtual/avoidance_radius/results/radius_0.01/pminus.dat")

    # and plot the whole shebang
    try:
        matplotlib = importlib.import_module('matplotlib')
        #matplotlib.use('PDF')
        plt = importlib.import_module('matplotlib.pyplot')
        hasmatplotlib = True
        
    except:
        print("Could not find matplotlib. Skipping all plots..")

    if hasmatplotlib:
        # prepare plot
        matplotlib.rc('text', usetex=True)
        plt.figure(figsize=(4*2, 3*2), dpi=600)
        
        plt.plot(s, pplus, "ro")
        plt.plot(s, pminus, "bo")
        plt.plot(s, ppt, "r-")
        plt.plot(s, pmt, "b-")
        plt.plot(s, ppmeso, "r+")
        plt.plot(s, pmmeso, "b+")

        plt.xlabel(r'$s$')
        plt.ylabel(r'$W_\pm(s)$')
        plt.savefig(datadir+'/analytic_comparison.pdf')
        plt.clf()
        
        # turn of latex support
        matplotlib.rc('text', usetex=False)

def analyzeFlame(options, datadir, gpgmp=False):
    # read in directories
    print("Reading in directories from {0:s}.".format(datadir+'parameters.pkl'))
    pkl_file = open(datadir+'parameters.pkl', 'rb')
    dirs = pickle.load(pkl_file)
    legend = pickle.load(pkl_file)
    pkl_file.close()

    # these will hold the arrays for the SDE coefficients
    a=[]
    b=[]
    useFPE = False

    # loop over dirs
    for draw in dirs:
        d=datadir+'/'+draw+'/'
        if gpgmp:
            analyzeGpgmpNondiffusiveSingle(options, d)
        else:
            print("Reducing dir {0:s}.".format(d))
            reduce(d)
            print("Analyzing dir {0:s}.".format(d))
            analyze(options, d)

        # read in data files
        alphacurr = numpy.genfromtxt(d+'/alpha.dat')
        b2curr = numpy.genfromtxt(d+'/beta2.dat')
        nagents = numpy.shape(b2curr)[0]-1
        if (nagents>100):
            print("Found {0:d} agents. Using FPE for integration.".format(nagents))
            useFPE = True
        else:
            print("Found {0:d} agents. Using Master equation for integration.".format(nagents))

        a.append(alphacurr)
        b.append(b2curr)


    # convert the python arrays into numpy
    alphatot = numpy.array(a)
    beta2tot = numpy.array(b)

    # and plot them
    try:
        matplotlib = importlib.import_module('matplotlib')
        #matplotlib.use('PDF')
        plt = importlib.import_module('matplotlib.pyplot')
        hasmatplotlib = True
        
    except:
        print("Could not find matplotlib. Skipping all plots..")

    if useFPE:
        if hasmatplotlib:
            # prepare plot
            matplotlib.rc('text', usetex=True)
            plt.figure(figsize=(4*4, 3*2), dpi=600)

            ss = numpy.shape(alphatot)[1]
            sgrid = numpy.mgrid[0:ss]/float(ss-1)
            styles = ['bo-', 'ro-', 'go-']

            # and plot evolution
            plt.subplot(121)
            for i in range(numpy.shape(alphatot)[0]):
                plt.plot(sgrid, alphatot[i,:]/1e-4, styles[i % len(styles)], label=legend[i])

            plt.xlabel(r'$s$')
            plt.ylabel(r'$\alpha(s)/10^{-4}$')
            plt.legend(loc='lower right')

            # and beta^2
            plt.subplot(122)
            for i in range(numpy.shape(beta2tot)[0]):
                plt.plot(sgrid, beta2tot[i,:]/1e-4, styles[i % len(styles)], label=legend[i])

            plt.xlabel(r'$s$')
            plt.ylabel(r'$\beta^2(s)/10^{-4}$')
            plt.legend(loc='lower center')

            plt.savefig(datadir+'/total_sde.pdf')
            plt.clf()

            # turn of latex support
            matplotlib.rc('text', usetex=False)

    # make the comparison plot
    dd = [datadir+'/'+d for d in dirs]
    makeComparisonPlot(dd, datadir, useFPE)


def runFlameSingle(options, datadir, diffusive, tot, nruns, sigma, recError,
                   speed=0.01, avoidance=0.01, length=1, nmax=5000, track_positions=False):
    # get exec name
    if diffusive:
        executable = '/home/matthias/Software/flame/models/majority_diff/main'
    else:
        executable = options['flameexec']
    print("Using flame executable {0:s}.".format(executable))

    # write parameter information 
    paramfile = open(datadir+'parameters.pkl', 'w')
    pickle.dump(nruns, paramfile)
    pickle.dump(tot, paramfile)
    paramfile.close()

    # we change into this dir (so we can use tar)
    olddir = os.getcwd()
    os.chdir(datadir)

    # we need to convert walltime into hr::mn::sc
    min, sec = divmod(options['pbsWalltime'], 60)
    h, min = divmod(min, 60)
    walltime = "{0:d}:{1:02d}:{2:02d}".format(h, min, sec)
    
    # if we run it using gps we first need to create all input files
    # and add them to the script before submitting
    if options['pbs']:
        # number of experiments per batch
        nexp = 250

        # we need to create our custom PBS script
        template = '''#!/bin/bash
#$ -l h_rt={1:s}
#$ -o {0:s}
#$ -e {0:s}

cd {0:s}
'''.format(datadir, walltime)
        for j in range(0, nruns/nexp):
            # write script header for this batch
            scriptname = datadir+'/run_pbs_{0:d}'.format(j)
            f = open(scriptname, 'w')
            f.write(template)

            # now create the input files for batch and append them to script
            for i in range(nexp):
                # create model input file
                inputFileName = 'model_{0:d}_0.xml'.format(j*nexp+i);
                writeFlameInputFile(datadir+'/'+inputFileName, tot, sigma, recError, j*nexp+i, speed, avoidance, length, track_positions=track_positions)

                # create the flame command
                params = "{0:d} {1:s} -f {0:d}".format(nmax, datadir+'/'+inputFileName)

                # append the current job to the PBS script
                command = '''{0:s} {1:s}\n'''.format(executable, params)
                f.write(command)
                
            # and run it
            f.close()
            subprocess.call(['qsub', scriptname])
    
    # this is for single node jobs
    else:
        for i in range(0,nruns):
        
            # create model input file
            inputFileName = 'model_{0:d}_0.xml'.format(i);
            writeFlameInputFile(datadir+'/'+inputFileName, tot, sigma, recError, i, speed, avoidance, length, track_positions=track_positions)

            # run it
            params = "{0:d} {1:s} -f {0:d}".format(nmax, datadir+'/'+inputFileName)

            print("Running job {0:d} in {1:s}.".format(i, datadir))
            common.jobs.run(options, executable, datadir, params)
        
            # read in rho
            rho = numpy.loadtxt(datadir+'rho_{0:d}.dat'.format(i))
            time = numpy.mgrid[0:nmax]
            numpy.savez(datadir+'rho_{0:d}.npz'.format(i), rho=rho, time=time)
        
            # delete rho and rename files
            shutil.move(datadir+"stderr.txt", datadir+"stderr_{0:d}.txt".format(i))
            shutil.move(datadir+"stdout.txt", datadir+"stdout_{0:d}.txt".format(i))
            os.remove(datadir+"rho_{0:d}.dat".format(i))
        
            # add model and log files to archive
            subprocess.call("tar rf {0:s} {1:s} {2:s}".format("logs.tar", "stderr_{0:d}.txt".format(i), "stdout_{0:d}.txt".format(i)), shell=True)
            subprocess.call("tar rf {0:s} {1:s}".format("models.tar", inputFileName), shell=True)
            os.remove(datadir+"stderr_{0:d}.txt".format(i))
            os.remove(datadir+"stdout_{0:d}.txt".format(i))
            os.remove(datadir+'/'+inputFileName)


    # zip the tar fles
    #print "Calling {0:s}.".format("gzip {0:s}".format(datadir+"logs.tar"))
    #subprocess.call("gzip {0:s}".format(datadir+"logs.tar"))
    #subprocess.call("gzip {0:s}".format(datadir+"models.tar"))

    # change into previous dir
    os.chdir(olddir)

    # copy over model spec and agent_functions.c
    shutil.copy(options['flamedir']+'/majority.xml', datadir)
    shutil.copy(options['flamedir']+'/agent_functions.c', datadir)

def plotSDEComparison(datadir1, datadir2, resultsdir):
    # read in parameters
    pkl_file = open(datadir1+'parameters.pkl', 'rb')
    tot = pickle.load(pkl_file)
    pkl_file.close()

    # read in data from both dirs
    alpha1 = numpy.loadtxt(datadir1+'/alpha.dat')
    beta1 = numpy.loadtxt(datadir1+'/beta.dat')
    beta21 = numpy.loadtxt(datadir1+'/beta2.dat')
    alpha2 = numpy.loadtxt(datadir2+'/alpha.dat')
    beta2 = numpy.loadtxt(datadir2+'/beta.dat')
    beta22 = numpy.loadtxt(datadir2+'/beta2.dat')

    # dynamically import matplotlib (if possible)
    hasmatplotlib = False

    try:
        matplotlib = importlib.import_module('matplotlib')
        #matplotlib.use('PDF')
        plt = importlib.import_module('matplotlib.pyplot')
        hasmatplotlib = True
    except:
        print("Could not find matplotlib. Skipping all plots..")

    # and plot them
    si = numpy.asfarray(numpy.mgrid[0:tot+1])/tot
    
    if hasmatplotlib:
        plt.figure(figsize=(4*4,3*2), dpi=300)
    
        plt.subplot(121)
        plt.plot(si, alpha1/1e-4)
        plt.plot(si, alpha2/1e-4)
        plt.xlabel(r'$s$')
        plt.ylabel(r'$\alpha(s)/10^{-4}$')
        
        plt.subplot(122)
        plt.plot(si, beta21/1e-4)
        plt.plot(si, beta22/1e-4)
        plt.xlabel(r'$s$')
        plt.ylabel(r'$\beta^2(s)/10^{-4}$')
    
        plt.savefig(resultsdir+'/sde_comparison.pdf')

        # turn of latex support
        matplotlib.rc('text', usetex=False)


def analyzeGpgmpNondiffusiveSingle(options, datadir, toyModel=False):

    # read in parameters
    pkl_file = open(datadir+'parameters.pkl', 'rb')
    if toyModel:
        tot = pickle.load(pkl_file)
        nred = pickle.load(pkl_file)
        k0 = pickle.load(pkl_file)
    else:
        tot = pickle.load(pkl_file)
        nred = pickle.load(pkl_file)
        depth = pickle.load(pkl_file)
        r0 = pickle.load(pkl_file)
        v0 = pickle.load(pkl_file)
        omega = pickle.load(pkl_file)
        k0 = pickle.load(pkl_file)
    pkl_file.close()

    # check if we should use FPE or Master equation
    useFPE = False
    if (tot>100):
        print("Found {0:d} agents. Using FPE for integration.".format(tot))
        useFPE = True
    else:
        print("Found {0:d} agents. Using Master equation for integration.".format(tot))

    # read in data
    nru, ngu, tu = gpgmp.models.majority.read_totals_no_diffusion(datadir+'output.h5', numRuns=1)

    # get times - we assume it was only one run
    times=tu[0,:]

    # sort them
    args=numpy.argsort(times)

    # compute density
    rho=numpy.transpose(nru[0,args,:])/tot

    # get evolution histograms and store them
    ev, edges =gpgmp.models.majority.getEvolutionPDFRho(rho, times, tot)
    numpy.savetxt(datadir+'/evolutionPDF.dat', ev)
    numpy.savetxt(datadir+'/times.dat', times[args])

    # dynamically import matplotlib (if possible)
    hasmatplotlib = False

    try:
        matplotlib = importlib.import_module('matplotlib')
        #matplotlib.use('PDF')
        plt = importlib.import_module('matplotlib.pyplot')
        hasmatplotlib = True
    except:
        print("Could not find matplotlib. Skipping all plots..")

    # and plot it
    if hasmatplotlib:
        # prepare plot
        matplotlib.rc('text', usetex=True)
        plt.figure(figsize=(4, 3), dpi=300)

        plt.imshow(numpy.transpose(ev), extent=[0., numpy.max(times), edges[0], edges[-1]], aspect='auto')
        plt.xlabel(r'$t\ \mathrm{[s]}$')
        plt.ylabel(r'$s$')
        plt.savefig(datadir+'/pdfEvolution.pdf')
        plt.clf()

    # compute dt
    dt = times[1]-times[0]
    
    # now get alpha and beta assuming time-homogenity
    alpha, beta, beta2, pplus, pminus, p0 = gpgmp.models.majority.computeFptHomogeneous(rho, dt, tot)
    
    # we export them
    numpy.savetxt(datadir+'/alpha.dat', alpha)
    numpy.savetxt(datadir+'/beta.dat', beta)
    numpy.savetxt(datadir+'/beta2.dat', beta2)
    numpy.savetxt(datadir+'/pplus.dat', pplus)
    numpy.savetxt(datadir+'/pminus.dat', pminus)
    numpy.savetxt(datadir+'/p0.dat', p0)

    # compute splitting probability and save it
    init=rho[:,0]
    end=rho[:,-1]
    sp = gpgmp.models.majority.computeSplittingProbability(init, end, tot)
    numpy.savetxt(datadir+'/splitting_probability_experiment.dat', sp)

    if useFPE:
        # run matlab script to integrate FPE
        matcom = '{0:s} -nodesktop -nosplash -r "cd \'{1:s}\'; solveFPE(\'{2:s}\'); exit "'.format(options['matlabbin'], options['majoritymathscripts'], datadir)

        # create files to capture output
        outfile = open(datadir+'/matlab_out.txt', 'w')
        errfile = open(datadir+'/matlab_err.txt', 'w')
    
        # and run it
        subprocess.call(matcom, stdout=outfile, stderr=errfile, shell=True)

        # we also need to compute the decision time
        computeExitTimeFPE(datadir, False)

    else:
        computeBirthDeath(resultsdir)

    # and plot them
    si = numpy.asfarray(numpy.mgrid[0:tot+1])/tot
    t = times[args]

    if useFPE:
        decTimeData = computeDecisionTime(rho, t, int(tot))
    else:
        # for this we set delta to zero ..
        decTimeData = computeDecisionTime(rho, t, int(tot), 0)

    # and save it
    numpy.savetxt(datadir+'/decision_time_experiment.dat', decTimeData)
    
    if hasmatplotlib:
        plt.figure(figsize=(4*4,3*2), dpi=300)
    
        plt.subplot(121)
        plt.plot(si, alpha/1e-4)
        plt.xlabel(r'$s$')
        plt.ylabel(r'$\alpha(s)/10^{-4}$')
        
        plt.subplot(122)
        plt.plot(si, beta2/1e-4)
        plt.xlabel(r'$s$')
        plt.ylabel(r'$\beta^2(s)/10^{-4}$')
    
        plt.savefig(datadir+'/sde.pdf')

        # turn of latex support
        matplotlib.rc('text', usetex=False)

def runGpgmpNondiffusiveAvoidance(options, datadir):
    # create list of dirs for pickling
    dirs=[]
    legend = []
    
    # loop over avoidance radius
    for radius in   [0.01, 0.05, 0.075, 0.1]:
        #radius = 0.01*float(avoidance)/10.+0.01

        # create dir
        reldir = 'radius_{0:g}/'.format(radius)
        curdir = datadir+'/'+reldir
        dirs.append(reldir)
        os.makedirs(curdir)
        runGpgmpNondiffusiveSingle(options, curdir, 150, 0, avoidanceRadius=radius)
        legend.append(r'$r_a={0:g}$'.format(radius))

    # and pickle the directory list
    paramfile = open(datadir+'parameters.pkl', 'w')
    pickle.dump(dirs, paramfile)
    pickle.dump(legend, paramfile)
    paramfile.close()
    

def runGpgmpNondiffusive(options, datadir):
    # create list of dirs for pickling
    dirs=[]
    legend = []
    
    # loop over sigma
    for error in range(10):
        recError = float(error)/10.*0.2

        # create dir
        curdir = datadir+'/error_{0:g}/'.format(recError)
        dirs.append(curdir)
        os.makedirs(curdir)
        runGpgmpNondiffusiveSingle(options, curdir, 150, recError)
        legend.append(r'$r={0:g}$'.format(recError))

    # and pickle the directory list
    paramfile = open(datadir+'parameters.pkl', 'w')
    pickle.dump(dirs, paramfile)
    pickle.dump(legend, paramfile)
    paramfile.close()

def runGpgmpToyModel(options, datadir):
    # get executable name
    executable = options['gpgmpexec']
    
    # compute parameters
    tot = 150
    nred = 75
    k0 = 0.0002

    # write parameter information 
    paramfile = open(datadir+'parameters.pkl', 'w')
    pickle.dump(tot, paramfile)
    pickle.dump(nred, paramfile)
    pickle.dump(k0, paramfile)
    paramfile.close()

    # create parameters
    params = "-p 29 -x 128 -y 128 -l 1 -t 4000 -n 1000000 -r 1 -o 10 -s 0 -c {0:s} -i {1:s} {2:d} {3:d} {4:g}".format(options['gpgmpcldir'], options['gpgmpinitdir'],
                                                                                                                          tot, nred, k0)
    # run the job
    common.jobs.run(options, executable, datadir, params)

def runGpgmpNondiffusiveSingle(options, datadir, tot, recError, avoidanceRadius=0.01):
    # get executable name
    executable = options['gpgmpexec']
    
    # compute parameters
    tot = 150
    nred = 75
    depth = 5

    # equation here comes from rate_constants_gillespie.nb
    r0 = avoidanceRadius # avoidance radius
    v0 = 0.01 # robot speed
    omega = 1 # domain area
    #k0 = 0.0002
    k0 = 2*r0*v0/omega

    # write parameter information 
    paramfile = open(datadir+'parameters.pkl', 'w')
    pickle.dump(tot, paramfile)
    pickle.dump(nred, paramfile)
    pickle.dump(depth, paramfile)
    pickle.dump(r0, paramfile)
    pickle.dump(v0, paramfile)
    pickle.dump(omega, paramfile)
    pickle.dump(k0, paramfile)
    paramfile.close()

    # create parameters
    params = "-p 28 -x 128 -y 128 -l 1 -t 4000 -n 1000000 -r 1 -o 10 -s 0 -c {0:s} -i {1:s} {2:d} {3:d} {4:d} {5:g} {6:g}".format(options['gpgmpcldir'], options['gpgmpinitdir'],
                                                                                                                          tot, nred, depth, k0, recError)
    # run the job
    common.jobs.run(options, executable, datadir, params)

def run(options, datadir):
    # get executable name
    executable = options['gpgmpexec']
    
    # create parameters
    tot = 150
    params = "-p 27 -x 64 -y 64 -l 1 -t 2000 -n 10000000 -r 1 -o 10 -s 0 -c {0:s} -i {1:s} 3e-03 3e-03 0 0. 150 50 5 0.485675".format(options['gpgmpcldir'], options['gpgmpinitdir'] )
    
    # do number of runs and process the output accordingly
    nruns = 500

    # write parameter information 
    paramfile = open(datadir+'parameters.pkl', 'w')
    pickle.dump(nruns, paramfile)
    pickle.dump(tot, paramfile)
    paramfile.close()

    for i in range(0,nruns):
        print("Running job {0:d}.".format(i))
        common.jobs.run(options, executable, datadir, params)
        print("done.")
        
        #os.rename('output.h5', 'majority_output_{0:d}.h5'.format(i))

        # reduce the single file
        #nreds, ngreens, times = reduce_single(datadir+'output')
        nr, ng, t = read_totals(datadir+'output.h5')
        nreds = nr[0,:]
        ngreens = ng[0,:]
        times = t[0,:]
        
        # do garbage collection .. to free up memory
        gc.collect()

        # get total count
        tot = nreds[0]+ngreens[0]
    
        # time ordering
        args=numpy.argsort(times)
    
        # compute rho and ordered time
        rho=nreds[args]/tot
        time = times[args]

        # and save it
        numpy.savez(datadir+'rho_{0:d}.npz'.format(i), rho=rho, time=time)
        
        # and delete the HDF file
        os.remove(datadir+'output.h5')

def read_totals_no_diffusion(filename, numRuns=0):
        # read in file
        print("Processing file {0:s}.".format(filename))

        # Open file
        f=h5py.File(filename, 'r')

        # get model group
        modelGroup = f['Model']

        # get dumps
        runGroup = f['Runs']

        # get number of dumps
        crg = runGroup['0']
        numDumps = len(crg)

        iRun = 0

        # get shape
        currentRunGroup = runGroup['{0:d}'.format(iRun)]
        dataGroup = currentRunGroup['Dump_0']
        numCells = numpy.shape((dataGroup['Red'][:,:]).flatten())[0]
        print("Found {0:d} runs, {1:d} cells and {2:d} time dumps.".format(numRuns, numCells, numDumps))

        # output arrays
        nreds   = numpy.zeros((numRuns, numDumps, numCells))
        ngreens = numpy.zeros((numRuns, numDumps, numCells))
        times   = numpy.zeros((numRuns, numDumps))
                
        # loop over dumps - there should be only one run!
        iDump = 0
        for name in currentRunGroup:

            # get time
            dataGroup = currentRunGroup[name]
            time = float(dataGroup.attrs['time'][0])
            times[iRun, iDump] = time

            # get reds
            dataGroup = currentRunGroup[name]
            reds = (dataGroup['Red'][:,:]).flatten()
            nreds[iRun, iDump, :] = reds

            # get greens
            dataGroup = currentRunGroup[name]
            greens = (dataGroup['Green'][:,:]).flatten()
            ngreens[iRun, iDump, :] = greens

            # and increase
            iDump = iDump + 1
        
        f.close()
        
        return nreds, ngreens, times

def read_totals(filename, numRuns=0):
        # read in file
        print("Processing file {0:s}.".format(filename))

        # Open file
        f=h5py.File(filename, 'r')

        # get model group
        modelGroup = f['Model']

        # get number of runs
        if numRuns==0:
            numRuns=int(modelGroup.attrs['numRuns'])

        # get dumps
        runGroup = f['Runs']

        # get number of dumps
        crg = runGroup['0']
        numDumps = len(crg)

        print("Found {0:d} runs and {1:d} time dumps.".format(numRuns, numDumps))

        # output arrays
        nreds   = numpy.zeros((numRuns, numDumps))
        ngreens = numpy.zeros((numRuns, numDumps))
        times   = numpy.zeros((numRuns, numDumps))

        # loop over runs
        for iRun in range(numRuns):
            currentRunGroup = runGroup['{0:d}'.format(iRun)]
                
            # loop over dumps - there should be only one run!
            iDump = 0
            for name in currentRunGroup:
            
                # get time
                dataGroup = currentRunGroup[name]
                time = float(dataGroup.attrs['time'][0])
                times[iRun, iDump] = time

                # get reds
                dataGroup = currentRunGroup[name]
                reds = float(dataGroup.attrs['Red'][0])
                nreds[iRun, iDump] = reds

                # get greens
                dataGroup = currentRunGroup[name]
                greens = float(dataGroup.attrs['Green'][0])
                ngreens[iRun, iDump] = greens

                # and increase
                iDump = iDump + 1
        
            # increase run number
            iRun = iRun + 1

        # close file
        f.close()
        
        return nreds, ngreens, times

def reduce(datadir, nruns=0):
    # check if rho.npz exists .. in which case we don't reduce
    if os.path.isfile(datadir+'/rho.npz'):
        print("Found {0:s}.. skipping reduce step.".format(datadir+'/rho.npz'))
        return

    # read in number of runs if not given
    if nruns==0:
        pkl_file = open(datadir+'parameters.pkl', 'rb')
        nruns = pickle.load(pkl_file)
        pkl_file.close()
    
    print("Reading in {0:d} runs.".format(nruns))
    
    for i in range(nruns):
        
        # reduce the single file
        #nreds, ngreens, times = reduce_single(datadir+'majority_output_{0:d}'.format(i))
        # do garbage collection .. to free up memory
        #gc.collect()

        # the HDF5 files should already have been reduced to numpy arrays .. read those in!
        if (i%1)==0:
            print("Reading in run {0:s}.".format(datadir+'rho_{0:d}.npz'.format(i)))
            
        df = numpy.load(datadir+'rho_{0:d}.npz'.format(i) )
        rhoc  = df['rho']
        times = df['time']
        
        if i==0:
            # create output array
            nt = numpy.shape(times)[0]
            rho = numpy.zeros((nruns, nt))
            print("Found {0:d} time dumps.".format(numpy.shape(times)[0]))
            time = times
            rho[i,:]=rhoc    

        # check if time agrees ..
        if i>0:
            if ((numpy.shape(times)[0] != numpy.shape(time)[0]) or (numpy.shape(rhoc)[0]!=numpy.shape(rho)[1])):
                print("WARNING: Dump {0:d} has incorrect shape ({1:d} instead of {2:d}).. dismissing.".format(i,numpy.shape(times)[0],numpy.shape(time)[0] ))
            else:
                # compute rho and ordered time
                rho[i,:]=rhoc    
                
    # and save as binary
    numpy.savez(datadir+'rho.npz', rho=rho, time=time)

    # finally delete all intermediate density files
    subprocess.call("rm -rf rho_*", shell=True)

# -------------------------------------------------------------------------------------------------
# all stuff below should be DEPRECATED!!

def computeFPT(rho, time, n, ii):
    # first compute dt
    dt = time[1]-time[0]
    
    # compute relative and absolute change where we omit the first time value
    srel=(rho[:,1:]-rho[:,:-1])/dt
    sabs=numpy.abs(srel)

    # now srel and sabs are given in the center of the cell .. we need to compute
    # the corresponding rho (which we assume to be the same as at the start of the time step..
    # sort of like an Ito definition)
    rhonew=rho[:,:-1]

    # pick out the right time index
    sr = srel[:,ii]
    sa = sabs[:,ii]
    rn = rhonew[:, ii]
    
    # now for each value of rho el (0, 1/N, 2/N, .., 1) we find all corresponding values at a certain time
    # of sabs and srel and average them (see ref [13] from Bernd's paper)
    srelnew = numpy.zeros((n+1))
    sabsnew = numpy.zeros((n+1))
    for i in range(n+1):
        if (i==0):
            index = numpy.where(rn==0)
        elif (i==n):
            index = numpy.where(rn==1)
        else:
            index = numpy.where(numpy.logical_and(rn<(float(i)+0.25)/float(n), rn>(float(i)-0.25)/float(n)))
        
        if (numpy.shape(index)[1]>0):
            srelnew[i] = numpy.mean(sr[index])
            sabsnew[i] = numpy.mean(sa[index])

    return srelnew, sabsnew

def count(n, species):
    # output arrays
    nindices = numpy.shape(n)[1]
    nreds = numpy.zeros((nindices))
    ngreens = numpy.zeros((nindices))

    # get all red and green species
    redspecies = [numpy.where(species == x)[0][0] for x in species if x[0]=='R']
    greenspecies = [numpy.where(species == x)[0][0] for x in species if x[0]=='G']

    # count all reds
    for red in redspecies:
        dreds = numpy.sum(numpy.sum(n[:,:,red, :, :], axis=2), axis=2)
        nreds = nreds + dreds
        #print "For {0:s}:".format(species[red]), dreds

    # count all greens
    for green in greenspecies:
        dgreens = numpy.sum(numpy.sum(n[:,:,green, :, :], axis=2), axis=2)
        ngreens = ngreens + dgreens
        #print "For {0:s}:".format(species[green]), dgreens

    # and return the results
    return (nreds, ngreens)

def reduce_single(filename):
        # read in file
        print("Processing file {0:s}.".format(filename))
        #n, times, species, nruns = gpgmp.io.read_gmp_hdf5(filename)
        # count
        #nreds, ngreens = count(n, species)

        # we need to do it by hand .. so we don't need to load massive files into memory!
        f=h5py.File('output.h5', 'r')
        
        # get model group
        modelGroup = f['Model']

        # get species list
        numSpecies=modelGroup.attrs['numSpecies']
        species=modelGroup['Species'][:]
        
        # get all red and green species
        redspecies = [x for x in species if x[0]=='R']
        greenspecies = [x for x in species if x[0]=='G']

        # get dumps
        runGroup = f['Runs']
        currentRunGroup = runGroup['0']
        numDumps = len(currentRunGroup)
        
        # output arrays
        nreds   = numpy.zeros((numDumps))
        ngreens = numpy.zeros((numDumps))
        times   = numpy.zeros((numDumps))
        
        # loop over dumps - there should be only one run!
        iDump = 0
        for name in currentRunGroup:
                dataGroup = currentRunGroup[name]
                time = float(dataGroup.attrs['time'][0])
                #print "Found dataset at t={0:g}.".format(time)
                
                times[iDump] = time
                
                # count all reds
                nredscur = 0
                for red in redspecies:
                    dreds = numpy.sum(dataGroup[red])
                    nredscur = nredscur + dreds

                # count all greens
                ngreenscur = 0
                for green in greenspecies:
                    dgreens =  numpy.sum(dataGroup[green])
                    ngreenscur = ngreenscur + dgreens

                # and put them into output array
                #print "Found {0:d} reds and {1:g} greens (total={2:d})".format(nredscur, ngreenscur, nredscur+ngreenscur)
                nreds[iDump] = nredscur
                ngreens[iDump] = ngreenscur
                
                # and increase
                iDump = iDump + 1
        
        # close file
        f.close()
        
        return nreds, ngreens, times

    
def plot(fn='output'):
    # read in file
    n, times, species, nruns = gpgmp.io.read_gmp_hdf5(fn)

    # count
    nreds, ngreens = count(n, species)

    # get total count .. hopefully conserved!
    tot = nreds[0]+ngreens[0]

    # time ordering
    args=numpy.argsort(times)

    # get number of runs
    ind=numpy.shape(n)[0]

    # and plot

    for i in range(ind):
        plt.plot(times[args], nreds[i,args]/tot, color=cm.prism(i))

    #plt.plot(times[args], ngreens[args], "g")
    #plt.plot(times[args], nreds[args]+ngreens[args], "b")
    plt.show()

def computeRho(nr, ng, times):
    # we first need to get the ordered times .. we assume the dump times are the same for all
    targs = numpy.argsort(times[0,:])
    
    # total number of particles . should be the same for all!
    tot = nr[0,targs[0]] + ng[0, targs[0]]
    
    return nr[:,targs]/float(tot), times[0,targs]


def getEvolutionPDF(nr, ng, times, nbins=50):
    # we first need to get the ordered times .. we assume the dump times are the same for all
    targs = numpy.argsort(times)
    
    # total number of particles . should be the same for all!
    tot = nr[0,targs[0]] + ng[0, targs[0]]
    
    # create edges
    edges = numpy.mgrid[0:1:complex(nbins+1)]
    
    # create output array
    ndumps = numpy.shape(targs)[0]
    ret = numpy.zeros((ndumps, nbins))
    
    # create histogram for all
    for i in range(ndumps):
        hist, tedges = numpy.histogram(nr[:,targs[i]]/float(tot), bins=edges, density=True)
        ret[i,:] = hist
        
    return ret, edges

