# This runs the random drift test for the FPKMC/RW algorithm
# for gpgmp job control
import gpgmp.job

# matplotlib
import matplotlib.pyplot as plt

# for io
import gpgmp.io
import os

# numpy
import numpy

# need to add problem specific elements
import xml.etree.ElementTree as ET

# need error function
from scipy.special import erf

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

# run the job
def run(testdir, cldir, plot=False, executable='rw', norun=False, discrete=False):

    # run it and get the timing
    if norun==False:            
        disel = ET.Element("discrete")
        if discrete==True:
            disel.text='5'
        else:
            disel.text='0'
        runtime = gpgmp.job.run_rw(testdir, executable, cldir, 'RANDOM_DRIFT', problemElements=[disel])

        # rename output data files
        os.rename('position.dat', "random_drift.pos.dat")
    else:
        runtime = 0.


    # read in results
    pos=numpy.genfromtxt("random_drift.pos.dat")[:8000]

    #parameters
    diff = 1.
    time = 2.
    mu = 0.
    sigma = 2.

    if plot:
        # make histogram
        minx = -20.
        maxx = 20.
        ncells = 128.
        dx = (maxx-minx)/ncells
        edges=numpy.mgrid[0:ncells+1]*dx+minx
        xc = edges[:-1]+dx

        # compute analytic solution
        if discrete == True:
            ana = analytic_discrete(xc, diff, time, mu, sigma, 5)
        else:
            ana = analytic(xc, diff, time, mu, sigma)

        # plot 
        width1c = 2.*3.27
        fig=plt.figure(1, figsize=(2.*width1c,width1c), dpi=300)
        plt.clf()
        plt.title('Random drift problem')

        plt.subplot(121)
        pdf, bins, patches = plt.hist(pos[:,0], bins=edges, normed=True, edgecolor='none')
        plt.plot(xc, ana, "r")

        plt.subplot(122)
        pdf, bins, patches = plt.hist(pos[:,1], bins=edges, normed=True, edgecolor='none')
        plt.plot(xc, ana, "r")

        #plt.xlabel('x')
        #plt.ylabel('y')
        #plt.imshow(H, extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]])

        if discrete==True:
            plt.savefig("random_drift_discrete.pdf", format="pdf")
        else:
            plt.savefig("random_drift.pdf", format="pdf")

    # check if passed
    meanxx = numpy.sqrt(numpy.mean((pos[:,0]**2+pos[:,1]**2+pos[:,2]**2)))
    mxxtheo = time*(2*diff + time*(mu**2 + sigma**2))

    if ((meanxx-mxxtheo)/mxxtheo <= 1e-3):
        passed = True
    else:
        passed = False

    return passed, runtime
