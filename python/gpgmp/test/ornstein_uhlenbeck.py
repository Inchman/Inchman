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

# compute mean
def mean(diffx, diffy, g11, g12, g21, g22, t, x0, y0):
    # for ease of notation
    detg = numpy.sqrt(4.*g12*g21 + (g11 - g22)**2)

    mx = (detg*(1 + numpy.exp(detg*t))*x0 - (-1 + numpy.exp(detg*t))*(g11*x0 - g22*x0 + 2*g12*y0))/(2.*detg*numpy.exp(((detg + g11 + g22)*t)/2.))
    my = (-2*(-1 + numpy.exp(detg*t))*g21*x0 + (detg - g11 + numpy.exp(detg*t)*(detg + g11 - g22) + g22)*y0)/(2.*detg*numpy.exp(((detg + g11 + g22)*t)/2.))

    return (mx, my)

# compute variance
def sigma(diffx, diffy, g11, g12, g21, g22, t0):

    # for ease of notation
    detg = numpy.sqrt(4.*g12*g21 + (g11 - g22)**2)

    s11 = ((4*diffy*g12**2 + diffx*(detg + g11 - g22)**2)/(detg + g11 + g22) + \
               (-8*diffy*g12**2 + 2*diffx*(detg + g11 - g22)*(detg - g11 + g22))/(g11 + g22) - \
               (4*diffy*g12**2 + diffx*(detg - g11 + g22)**2)/(detg - g11 - g22) +  \
               ((2*numpy.exp(detg*t0)*(-(detg**2*diffx) + 4*diffy*g12**2 + diffx*(g11 - g22)**2))/(g11 + g22) - \
                    (4*diffy*g12**2 + diffx*(detg + g11 - g22)**2)/(detg + g11 + g22) + \
                    (numpy.exp(2*detg*t0)*(4*diffy*g12**2 + diffx*(detg - g11 + g22)**2))/(detg - g11 - g22))/numpy.exp((detg + g11 + g22)*t0))/(2.*detg**2)

    s12 = -((-2*diffy*g12 + 2*diffx*g21 - (2*(diffy*g12 - diffx*g21)*(g11 - g22))/(g11 + g22) + \
                  (2*(diffy*g11*g12 + diffx*(detg - g11)*g21))/(-detg + g11 + g22) + \
                  (2*diffy*g11*g12 - 2*diffx*(detg + g11)*g21)/(detg + g11 + g22) + \
                  (diffy*g12 - diffx*g21 + (numpy.exp(2*detg*t0)*(detg*(diffy*g12 + diffx*g21) + (diffy*g12 - diffx*g21)*(g11 - g22)))/(detg - g11 - g22) + \
                       (2*numpy.exp(detg*t0)*(diffy*g12 - diffx*g21)*(g11 - g22))/(g11 + g22) + \
                       (-2*diffy*g11*g12 + 2*diffx*(detg + g11)*g21)/(detg + g11 + g22))/numpy.exp((detg + g11 + g22)*t0))/detg**2)
    # should be s12=s21
    s21 = -((-2*diffy*g12 + 2*diffx*g21 - (2*(diffy*g12 - diffx*g21)*(g11 - g22))/(g11 + g22) + \
                  (2*(diffy*g11*g12 + diffx*(detg - g11)*g21))/(-detg + g11 + g22) + \
                  (2*diffy*g11*g12 - 2*diffx*(detg + g11)*g21)/(detg + g11 + g22) + \
                  (diffy*g12 - diffx*g21 + (numpy.exp(2*detg*t0)*(detg*(diffy*g12 + diffx*g21) + (diffy*g12 - diffx*g21)*(g11 - g22)))/(detg - g11 - g22) + \
                       (2*numpy.exp(detg*t0)*(diffy*g12 - diffx*g21)*(g11 - g22))/(g11 + g22) + \
                       (-2*diffy*g11*g12 + 2*diffx*(detg + g11)*g21)/(detg + g11 + g22))/numpy.exp((detg + g11 + g22)*t0))/detg**2)

    s22 = (-((4*diffx*g21**2 + diffy*(detg + g11 - g22)**2)/(detg - g11 - g22)) + \
                (-8*diffx*g21**2 + 2*diffy*(detg + g11 - g22)*(detg - g11 + g22))/(g11 + g22) + \
                (4*diffx*g21**2 + diffy*(detg - g11 + g22)**2)/(detg + g11 + g22) + \
                ((numpy.exp(2*detg*t0)*(4*diffx*g21**2 + diffy*(detg + g11 - g22)**2))/(detg - g11 - g22) + \
                     (2*numpy.exp(detg*t0)*(-(detg**2*diffy) + 4*diffx*g21**2 + diffy*(g11 - g22)**2))/(g11 + g22) - \
                     (4*diffx*g21**2 + diffy*(detg - g11 + g22)**2)/(detg + g11 + g22))/numpy.exp((detg + g11 + g22)*t0))/(2.*detg**2)

    return (s11, s12, s21, s22)


# analytic solution
def ornstein_uhlenbeck_2d(x, x0, y, y0, t, diffx, diffy, g11, g12, g21, g22):
    
    # create coordinate matrix
    gx=numpy.outer(x, numpy.ones_like(x))
    gy=numpy.outer(numpy.ones_like(y), y)

    # mean
    mux, muy = mean(diffx, diffy, g11, g12, g21, g22, t, x0, y0)

    # covariance matrix
    s11, s12, s21, s22 = sigma(diffx, diffy, g11, g12, g21, g22, t)

    # sigmas
    sigmax = numpy.sqrt(s11)
    sigmay = numpy.sqrt(s22)
    rho = s12/(sigmax*sigmay)

    # and compute gaussian
    gauss = 1./(2.*numpy.pi*sigmax*sigmay*numpy.sqrt(1.-rho**2)) \
        * numpy.exp(-1./(2.*(1-rho**2))*((gx-mux)**2/sigmax**2+(gy-muy)**2/sigmay**2-2.*rho*(gx-mux)*(gy-muy)/(sigmax*sigmay)))

    return gauss

def compute_rmse(state, theo):
    numMolecules = numpy.sum(state)
    return numpy.sqrt(numpy.sum((state-theo)**2)/(numpy.shape(state)[0]*numpy.shape(state)[1]))/numMolecules

#def run(testdir, plot=False, executable='gpgmp', deterministic=False, aqss=False):
def run(options, datadir):

    # get executable name
    executable = options['inchmanexec']

    # set output files
    newoptions = options
    newoptions['outlog'] = 'ornstein_uhlenbeck.out'
    newoptions['errlog'] = 'ornstein_uhlenbeck.err'

    # parameters
    params = "{0:s}/Ornstein\ Uhlenbeck.xml --cl-source {1:s} --output ornstein_uhlenbeck.h5".format(options['inchmanexampledir'], options['inchmancldir'])

    # run the job
    gpgmp.common.jobs.run(newoptions, executable, datadir, params)


def analyze(options, datadir, filename='ornstein_uhlenbeck.h5'):
    analyze(datadir, filename)

def analyze(datadir, filename):

    # Parameters
    length = 32.
    nx = 128
    tmax = 10
    diffx = 1.
    diffy = 1.
    gxx = 0.5
    gxy = 0.6
    gyx = 0.7
    gyy = 0.8
    numMolecules = 10000
    numruns = 50

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
    nstate=numpy.transpose(numpy.mean(n[:,timesarg[-1],0,:,:], axis=0))

    # compute coordinate system
    xc, dx = (float(length)/float(nx)*(numpy.mgrid[0:nx]-nx/2.), float(length)/float(nx))
    yc, dy = xc, dx
    x0 = xc[nx/2]
    y0 = yc[nx/2]

    # compute analytic solution
    theo = ornstein_uhlenbeck_2d(xc, x0, yc, y0, sampleTime, diffx, diffy, gxx, gxy, gyx, gyy)*dx*dy*numMolecules

    # compute rmse
    rmse = compute_rmse(nstate, theo)

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
        plt.clf()
        plt.plot(xc, nstate[:,nx/2], "r^", label="GPGMP")
        plt.plot(xc, theo[:,nx/2],"b--")
        plt.savefig(datadir+"/ornstein_uhlenbeck_xslice.eps", format="eps")
        plt.clf()
        plt.plot(xc, nstate[nx/2,:], "r^", label="GPGMP")
        plt.plot(xc, theo[nx/2,:],"b--")
        plt.savefig(datadir+"/ornstein_uhlenbeck_yslice.eps", format="eps")

    if rmse<=1e-4:
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
