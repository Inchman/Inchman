import gpgmp.job
import gpgmp.io
import gpgmp.helpers
import numpy
import matplotlib.pyplot as plt

# plot a single problem
def plot_single(hdffile, deterministic=False, sampletime=1000., speciesid=1., length=89.4):

    # read in problem
    n, times, species, nspecies = gpgmp.io.read_gmp_hdf5(hdffile, deterministic=deterministic)
    im=numpy.squeeze(n[0,numpy.where(times==sampletime),speciesid,:,:])

    # get plot specs and create figure
    specs = gpgmp.helpers.getPlotSpecs('plos')
    fig=plt.figure(1, figsize=(specs['width2c'], specs['width2c']*3./4.), dpi=300)
    plt.clf()

    # plot it
    plt.imshow(im, extent=[0,length,0,length])

    # remove ticklabels
    ax = plt.gca()
    for tick in ax.xaxis.get_major_ticks():
        tick.label1On = False
        tick.label2On = False
    for tick in ax.yaxis.get_major_ticks():
        tick.label1On = False
        tick.label2On = False
 
# run Gray-Scott problem
def run_single(problemdir, clSourceDir, initDir,
               executable='test',
               deterministic=False, aqss=False,
               length=89.4, nx=256, tmax=1000, numRuns=1, dtOut=50, U=0., V=0.005,
               k = 0.0225, F = 0.0025, Omega = 250):

    # Parameters
    #./gpgmp -x 256 -y 256 -t 1000 -l 89.4 -r 1 -o 50 -s 2 -p 20 0 0.005 0.0225 0.0025 250
    
    # run it and get the timing
    # we assume that the init script and the executable are already in the current directory..
    if deterministic:
        numruns = 1
        if aqss:
            solver = 4
        else:
            solver = 3
    else:
        solver = 2

    time = gpgmp.job.run_job(problemdir, executable, clSourceDir = clSourceDir, initDir=initDir,
                             parameters='{0:g} {1:g} {2:g} {3:g} {4:g}'.format(U, V, k, F, Omega),
                             problem=20, length=length, nx=nx, runtime=tmax, dtout=dtOut, solver=solver,
                             datafile='gray_scott.h5', outlog='gray_scott.out', errlog='gray_scott.err')

    return time

def run_range(datadir, clSourceDir, initDir, executable='test', deterministic=False):

    # problem parameters
    diffU = 0.
    diffV = 0.005 
    k = 0.0225
    F = 0.0025
    omegabase = 250.

    # geometry
    length = 44.7*2

    # runtime
    tmax = 5000.
    tout = 50.

    # loop over scale
    for scale in range(5):
        omega = omegabase*10**(scale-2)

        # create directory
        fulldir = basedir+'/{0:d}/run_{1:g}/'.format(nx, 10**scale)
        os.makedirs(fulldir)

        time = gpgmp.job.run_job(datadir, executable, clSourceDir=clSourceDir, initDir=initDir,
                                 parameters='{0:g} {1:g} {2:g} {3:g} {4:g}'.format(U, V, k, F, Omega),
                                 problem=20, length=length, nx=nx, runtime=tmax, dtout=dtOut, solver=solver,
                                 datafile='gray_scott.h5', outlog='gray_scott.out', errlog='gray_scott.err')
