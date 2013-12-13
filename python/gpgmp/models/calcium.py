import gpgmp.io
import gpgmp.helpers
import numpy
import importlib

import os.path
import sys

# plot a single problem
def plot_single(resultdir, filename, deterministic=False, sampletime=100., speciesid=0, length=250.):
    
    # read in problem
    n, times, species, nspecies = gpgmp.io.read_gmp_hdf5(resultdir+filename, deterministic=deterministic)
    im=numpy.squeeze(n[0,numpy.where(times==sampletime),speciesid,:,:])

    # get plot specs and create figure
    specs = gpgmp.helpers.getPlotSpecs('plos')

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

        plt.title('Intracellular Calcium')
        plt.legend(loc='lower right')

        plt.savefig(resultdir+"calcium.eps", format="eps")


# analyze it
if __name__ == "__main__":
    # print command line arguments
    dirname, filename = os.path.split(sys.argv[1])

    if dirname=='':
        dirname = './'
            
    plot_single(dirname, filename)
