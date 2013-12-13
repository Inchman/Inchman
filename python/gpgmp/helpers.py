import gpgmp.constants
import numpy

# Returns a concentration (in mol/l) from a given
# particle number and cell width (in mum)
def numberToConcentration(nmols, dx):
    return nmols/(gpgmp.constants.na*subvolumeInLiters(dx))

# Returns the subvolume size in liters from a given
# cell width (in mum)
def subvolumeInLiters(dx):
    return dx**3*1e-15

# Returns standard plot specifications for plots
def getPlotSpecs(journal='plos'):

    if journal=='plos':
        return {'width': 3.27, 'height': 3.27, 'width2c': 6.83}

# finds index and value of element that closest to the value given
def find_nearest(array,value):
    idx=(numpy.abs(array-value)).argmin()
    return idx, array[idx]

# tests the segmented scan implementation
def test_segmentedScan():
    # read in values
    ssinitial=numpy.genfromtxt("segscan_initial.dat")
    ssfinal=numpy.genfromtxt("segscan_final.dat")

    # compute number of values
    n = numpy.shape(ssinitial)[0]

    # find flag positions
    flags = numpy.where(ssinitial[:,2]==1)[0]
    nflags = numpy.shape(flags)[0]
    print("Finding {0:d} elements and {1:d} flags.".format(n, nflags))

    # compute proper partial sum
    partialSum = numpy.zeros((n))

    # start array parts
    il = 0
    ir = 0

    if (nflags > 0):
        ir = flags[0]

        for i in range(1,nflags):
            #print "Pass {0:d} working in range ({1:d},{2:d})".format(i, il, ir)
            # compute cumulative sum (inclusive)
            # in range (il, ir-1) = (il, flags-1)
            cumsum = numpy.cumsum(ssinitial[il:ir,1])

            # and put it into destination array
            # (shifted one to the right to get exclusive scan)
            partialSum[il+1:ir] = cumsum[:-1]

            # which array part are we dealing with next
            il = ir
            ir = flags[i]

        # and do last one
        #print "Pass {0:d} working in range ({1:d},{2:d})".format(nflags, il, ir)

        # compute cumulative sum (inclusive)
        # in range (il, ir-1) = (il, flags-1)
        cumsum = numpy.cumsum(ssinitial[il:ir,1])

        # and put it into destination array
        # (shifted one to the right to get exclusive scan)
        partialSum[il+1:ir] = cumsum[:-1]

    # and do really last one
    il = ir
    ir = n
    #print "Pass {0:d} working in range ({1:d},{2:d})".format(nflags+1, il, ir)

    # compute cumulative sum (inclusive)
    # in range (il, ir-1) = (il, flags-1)
    cumsum = numpy.cumsum(ssinitial[il:ir,1])
    
    # and put it into destination array
    # (shifted one to the right to get exclusive scan)
    partialSum[il+1:ir] = cumsum[:-1]

    # now check if it's all good
    diff = (numpy.where(partialSum != ssfinal[:,1]))[0]
    if (numpy.shape(diff)[0]) > 0 :
        print("FAILED !")
    else:
        print("PASSED.")

    return (ssinitial[:,1], flags, ssfinal[:,1], partialSum, diff)
