# import matplotlib
import matplotlib
import matplotlib.pyplot as plt

# use latex for rendering
from matplotlib import rc, mpl

# import system libs and numpy
import gpgmp.io
import gpgmp.job
import numpy
import shutil

# use for grid layout of subplots
from mpl_toolkits.axes_grid1 import ImageGrid

# run the job
def run(testdir, initscriptdir, plot=False, norun=False, executable='gpgmp', deterministic=False, aqss=False):
    
    # Parameters
    nx = 256
    tmax = 172800
    diffusionType = 1
    numruns = 1
    length = 1200.

    runTime = 0.

    # run it and get the timing
    if deterministic:
        numruns = 1
        if aqss:
            solver = 6
        else:
            solver = 5
    else:
        solver = 0

    # run it and get the timing
    if norun==False:
        # we need to get the init script
        shutil.copy(initscriptdir+'/init_slit.py', testdir)

        runTime = gpgmp.job.run_job(testdir, executable, parameters='{0:d}'.format(diffusionType),
                                    problem=8, length=length, nx=nx, runtime=tmax, dtout=5000, solver=solver, numRuns=numruns,
                                    datafile='problem_slit.h5', outlog='slit.out', errlog='slit.err')

    # read in results
    n, times, species, nspecies = gpgmp.io.read_gmp_hdf5('problem_slit')

    # and plot it
    width1c = 3.27
    height1c = width1c
    width2c = 6.83
    fs = 10

    if plot:
        # Create figure
        fig=plt.figure(1, figsize=(width2c,height1c), dpi=300)

        grid = ImageGrid(fig, 111, nrows_ncols = (1, 2), axes_pad=0.5, cbar_mode='single')

        # find index for t=24h and 48h
        idx24=(numpy.abs(times-3600.*24)).argmin()
        idx48=(numpy.abs(times-3600.*48)).argmin()

        maxcount = max(numpy.max(n[0,idx24,0,:,:]), numpy.max(n[0,idx48,0,:,:]))
        n24 = numpy.array(n[0,idx24,0,:,:], dtype=float)
        n48 = numpy.array(n[0,idx48,0,:,:], dtype=float)
        norm = mpl.colors.BoundaryNorm([0, 1,2,3,4,5,6,7], plt.cm.spectral_r.N)
        im = grid[0].imshow(n24, extent=[-600, 600, -600, 600], cmap=plt.cm.spectral_r, norm=norm)
        grid[1].imshow(n48, extent=[-600, 600, -600, 600], cmap=plt.cm.spectral_r, norm=norm)

        # prepare colorbar
        cb = grid[0].cax.colorbar(im)
        grid[0].cax.set_yticks([2, 4, 6])
        grid[0].cax.set_ylabel('counts')

        # set font size of ticklabels
        ax = plt.gca()
        for tick in ax.xaxis.get_major_ticks():
            tick.label1.set_fontsize(fs)
        for tick in ax.yaxis.get_major_ticks():
            tick.label1.set_fontsize(fs)

        # add labels
        grid[0].set_xlabel(r'$x\,[\mu\mathrm{m}]$', fontsize=fs)
        grid[1].set_xlabel(r'$x\,[\mu\mathrm{m}]$', fontsize=fs)
        grid[0].set_ylabel(r'$y\,[\mu\mathrm{m}]$', fontsize=fs)

        # set ticks
        ticks = [-600, -400, -200, 0, 200, 400, 600]
        ticklabels = [r'-600', '', '', r'0', '', '', '600']
        grid[0].set_xticks(ticks)
        grid[0].set_xticklabels(ticklabels)
        grid[1].set_xticks(ticks)
        grid[1].set_xticklabels(ticklabels)
        grid[0].set_yticks(ticks)
        grid[0].set_yticklabels(ticklabels)

        # title subplots
        grid[0].set_title(r'24 h')
        grid[1].set_title(r'48 h')

        plt.savefig("slit.pdf", format="pdf")
    
    return True, runTime
