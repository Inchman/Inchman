import gpgmp.job
import gpgmp.io
import gpgmp.helpers
import numpy
import matplotlib.pyplot as plt

# plot a single problem
def plot_single(resultdir, deterministic=False, sampletime=50., speciesid=0, length=0.54):
    
    # read in problem
    n, times, species, nspecies = gpgmp.io.read_gmp_hdf5(resultdir+'/oregonator', deterministic=deterministic)
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
 

# run single oregonator problem
def run_single(problemdir, executable='test', deterministic=False, aqss=False,
               length=0.54, nx=64, tmax=50, numRuns=1, dtOut=1,
               Omega=0.00005):
    # ./gpgmp -x 256 -y 256 -t 50 -l 0.54 -r 1 -o 5 -s 2 -p 19 0.00005

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

    time = gpgmp.job.run_job(problemdir, executable, parameters='{0:g}'.format(Omega),
                             problem=19, length=length, nx=nx, runtime=tmax, dtout=dtOut, solver=solver,
                             datafile='oregonator.h5', outlog='oregonator.out', errlog='oregonator.err')

    return time


    
