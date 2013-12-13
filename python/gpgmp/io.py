import csv
from numpy import *
# for hdf5
import h5py
import gpgmp.individuals
import sys
import os.path

# deprecated!
#def read_individual(filename):
#    kvp = genfromtxt(filename)
#    
#    particleindex = kvp[:,1]
#    cellpos = kvp[:,2]

#    # now compute 2d cell pos
#    cpy = int_(cellpos)/512
#    cpx = cellpos % 512
#    
#    return cpx, cpy

def read_cpugmp(filename, sizex=64, sizey=64, sizez=1, numSpecies=1, numRuns=1):

    # open datafile
    datafile=open(filename, 'r')

    # create initial output array
    n=[]

    for run in range(0, numRuns):
    #initial time index is 0
        time = 0

        ntempouter = zeros((1, sizex,sizey,sizez, numSpecies), dtype=int)

        # loop over times
        while True:
            # read in line
            line = datafile.readline()

            # is end of run or end of simulation?
            if (line=="# END OF SIMULATION\n"):
                #print "Detected new simulation.."
                break

            if  (len(line)==0):
                #print "Detected EOF."
                break

            # throw away time line (and print it)
            #print(line)

            # create temp array for this step
            ntemp = zeros((1, sizex, sizey, sizez, numSpecies), dtype=int)

            # loop over coordinates
            for i in range(0, sizex):
                for j in range(0, sizey):
                    for k in range(0, sizez):

                        # read in line
                        dline = datafile.readline()

                        # split it
                        linelist = dline.split('\t')

                        # get coordinates
                        ic = int(linelist[0])
                        jc = int(linelist[1])
                        kc = int(linelist[2])

                        # get species count
                        for s in range(0, numSpecies):
                            ntemp[0, ic, jc, kc, s] = int(linelist[3+s])

            # append temporary array to target array
            if (time == 0):
                ntempouter[0,:,:,:,:]=ntemp[0,:,:,:,:]
            elif (time > 0):
                ntempouter=append(ntempouter, ntemp, axis=0)

            # read next line to skip empty break and increment time index
            datafile.readline()
            time = time + 1

        # append to list
        n.append(ntempouter)

        # throw next line
        #print datafile.readline()
        
    # close file
    datafile.close()

    # now make an array out of it
    nret = zeros((numRuns, shape(n[0])[0], sizex,sizey,sizez, numSpecies), dtype=int)

    for i in range(0, numRuns):
        nret[i,:,:,:,:,:]=n[i]

    # return the array
    return nret

def gmp_hdf5_get_compartment_list(filename):
    # open file
    f=h5py.File(filename, 'r')

    # model group
    modelGroup = f['Model']

    # compartment group
    compartmentGroup = modelGroup['Compartments']

    # create return array
    retlist = []
    compnames = []

    # loop through all compartments
    for comp in compartmentGroup:
        item = compartmentGroup[comp][:]
        retlist.append(item)
        compnames.append(comp)

    # and return it as a map
    return (compnames, retlist)

# return value is given by
# n[#run, #dump, #species, #x, #y]

def read_field(filename, fieldname):
    # determine shape

    # read in file
    f=h5py.File(filename, 'r')

    # model group
    modelGroup = f['Model']

    # get number of runs
    numRuns = modelGroup.attrs['numRuns'][0]

    # get array size
    runGroup = f['Runs']
    firstRun = runGroup['0']
    numDumps = len(firstRun)
    timeStamps = zeros((numDumps), dtype=float)

    # get shape
    dataGroup = firstRun['Dump_0']
    fieldsGroup = dataGroup['Fields']
    field = fieldsGroup[fieldname]

    # size
    dx = field.shape[0]
    dy = field.shape[1]

    # create output array
    n = zeros((numRuns, numDumps, dx, dy), dtype=float)

    # loop over number of runs
    for run in range(0, numRuns):
        # get run group
        currentRunGroup = runGroup['{0:d}'.format(run)]

        # get all dumps
        iDump = 0

        for name in currentRunGroup:
            dataGroup = currentRunGroup[name]
            time = float(dataGroup.attrs['time'][0])
            if (iDump < numDumps):
                timeStamps[iDump]=time

            # get field
            fieldsGroup = dataGroup['Fields']
            field = fieldsGroup[fieldname]

            # and store it
            n[run, iDump, :, :] = field

            # increase dump counter
            iDump = iDump + 1

    # close file
    f.close()

    return n

def read_gmp_hdf5(filenamebase, deterministic=False, verbose=False):

    # determine shape

    # check if filename has ".h5" extension and add if necessary
    # (this is only for backward compability)
    filename = os.path.splitext(filenamebase)[0]+".h5"

    # read in file
    if verbose:
        print("Reading in {0:s}.".format(filename))

    f=h5py.File(filename, 'r')

    # model group
    modelGroup = f['Model']

    # find out if we're using old version or new version
    if ('numRuns' in modelGroup.attrs):
        # new file version
        if verbose:
            print("Detected new HDF5 structure. Reading in.")

        # get number of runs
        numRuns = modelGroup.attrs['numRuns'][0]

        # read in list of species
        numSpecies=modelGroup.attrs['numSpecies'][0]
        specieslist=modelGroup['Species'][:]

        # get array size
        ft = f
        if 'Jobs' in ft:
            if verbose:
                print("Found Jobs key!")
            f = (ft['Jobs'])['1']
            if verbose:
                print("Done")
        else:
            f = ft

        runGroup = f['Runs']

        firstRun = runGroup['0']
        numDumps = len(firstRun)
        timeStamps = zeros((numDumps), dtype=float)

        # get shape
        # read in list of species
        specieslist=modelGroup['Species'][:]

        # for python version >= 3.0, the species names are given as bytes objects,
        # we need to decode first
        if sys.version_info >= (3, 0):
            specieslist = array([x.decode() for x in specieslist])

        # get number of species
        numSpecies=specieslist.shape[0]
        dataGroup = firstRun['Dump_0']

        # use first available species to determine grid size
        species=dataGroup[str(specieslist[0])]
        dx = species.shape[0]
        dy = species.shape[1]

        # create output array
        if deterministic:
            n = zeros((numRuns, numDumps, numSpecies, dx, dy), dtype=float)
        else:
            n = zeros((numRuns, numDumps, numSpecies, dx, dy), dtype=int)

        # loop over number of runs
        for run in range(0, numRuns):

            if verbose:
                print("Reading run {0:d}..".format(run))

            # get run group
            currentRunGroup = runGroup['{0:d}'.format(run)]

            # get all dumps
            iDump = 0

            for name in currentRunGroup:
                dataGroup = currentRunGroup[name]
                time = float(dataGroup.attrs['time'][0])
                #print "Found dataset at t={0:g}.".format(time)
                if (iDump < numDumps):
                    timeStamps[iDump]=time

                    # loop over species
                    for iSpecies in range(0, numSpecies):
                        data=dataGroup[str(specieslist[iSpecies])][:,:]
                        # it seems to be tranposed so we'll get it back - no it's not anymore
                        #n[run, iDump, iSpecies, :, :] = transpose(data)
                        n[run, iDump, iSpecies, :, :] = data
                else:
                    print("(gillespieio.read_gmp_hdf5) WARNING: Found unexpected dataset at t={0:g}.. dropping!".format(time))

                iDump = iDump+1

        f.close()
        return (n, timeStamps, specieslist, numRuns)

    else:
        #old file version
        if verbose:
            print("Detected old HDF5 structure. Reading in.")

        # get number of runs
        numRuns = modelGroup['numRuns'][0]

        # get array size
        runGroup = f['Run_0']
        numDumps = len(runGroup)
        timeStamps = zeros((numDumps), dtype=float)
        dataGroup = runGroup['Data_0']

        # read in list of species
        specieslist=modelGroup['Species'][:]
        numSpecies=specieslist.shape[0]
        species=dataGroup[str(specieslist[0])]

        # size
        dx = species.shape[0]
        dy = species.shape[1]

        # create output array
        n = zeros((numRuns, numDumps, numSpecies, dx, dy), dtype=int)

        # close file
        f.close()

        # read in file
        filename = filenamebase+".h5"
        f=h5py.File(filename, 'r')

        # loop over number of runs
        for run in range(0, numRuns):

            #print "Reading run {0:d}..".format(run)

            # get run group
            runGroup = f['Run_{0:d}'.format(run)]

            # get all dumps
            iDump = 0
            for name in runGroup:
                dataGroup = runGroup[name]
                time = dataGroup['time'][0]
                #print "Found dataset at t={0:g}.".format(time)
                timeStamps[iDump]=time

                # loop over species
                for iSpecies in range(0, numSpecies):
                    data=dataGroup[str(specieslist[iSpecies])][:,:]
                    # it seems to be tranposed so we'll get it back
                    n[run, iDump, iSpecies, :, :] = transpose(data)

                iDump = iDump+1

        # close file
        f.close()

        return (n, timeStamps, specieslist, numRuns)

def read_gmp(filenamebase, numRuns=1, readInitial=False):

    retlist=[]

    for i in range(0, numRuns):
        if (readInitial):
            filename=filenamebase+"_initial"
        else:
            filename=filenamebase+"_{0:d}".format(i)

        datafile=open(filename, 'r')

        #open reader
        reader=csv.reader(datafile, delimiter=" ", quoting=csv.QUOTE_NONNUMERIC, skipinitialspace=True)

        # read the stuff in
        el = reader.next()
        nlist=[el[:len(el)-1]]
        for row in reader:
            # throw out last element (it's a quoted whitespace)
            nlist.append(row[:len(row)-1])

        #close it
        datafile.close()

        retlist.append(nlist)

    return array(retlist)

def read_ssc(filename, dx=32, dim=1):
    datafile=open(filename, 'r')

    #throw away header line
    datafile.readline()
    reader=csv.reader(datafile, delimiter="\t", quoting=csv.QUOTE_NONNUMERIC)

    # read the stuff in (throw away time argument)
    nlist=[reader.next()[1:]]
    for row in reader:
        nlist.append(row[1:])

    #close it
    datafile.close()

    # remove empty lines
    nlist.remove([])
    n = array(nlist)

    # reshape if necessary
    if dim==2:
        ndim = shape(n)
        ret = n.reshape((ndim[0], dx, dx))
    elif dim==1:
        ret = n

    return ret

def read_smartcell(dirname, timestep, numSpecies=1, numRuns=1):
    #dirname = "/home/matthias/.wine/drive_c/Program Files/Smartcell/windows/Output/sd/TimeSeries_0/"
    #timestep=50

    # read in slot file
    dirnameTime=dirname+"TimeSeries_0/"
    datafile=open(dirnameTime+"SlotMap.txt", 'r')

        #throw away header line
    datafile.readline()
    reader=csv.reader(datafile, delimiter="\t", quoting=csv.QUOTE_MINIMAL)

    fileNames=[]
    compNames=[]
    xc=[]
    yc=[]
    zc=[]
    cellSize=[]
    for row in reader:
        if (len(row)==8):
            fileNames.append(row[0])
            compNames.append(row[1])
            xc.append(row[3])
            yc.append(row[4])
            zc.append(row[5])
            cellSize.append(row[6])

        #close file
    datafile.close()

        # create coordinate system
    xc = map(int, xc)
    yc = map(int, yc)
    zc = map(int, zc)

    sizex=len(list(set(xc)))
    sizey=len(list(set(yc)))
    sizez=len(list(set(zc)))

    # create target array
    n = zeros((sizex,sizey,sizez, numSpecies, numRuns), dtype=int)

    # read stuff in and average over time series
    for ts in range(0, numRuns):
        print("Reading in "+dirname+"/TimeSeries_{0:d}/".format(ts)+"\n")
        for i in range(0, len(fileNames)):
            voxelFile=open(dirname+"/TimeSeries_{0:d}/".format(ts)+compNames[i]+"/"+fileNames[i]+".voxel", 'r')

            # throw first line
            voxelFile.readline()

            # read all stuff
            reader=csv.reader(voxelFile, delimiter="\t", quoting=csv.QUOTE_NONE)

            vn=[]

            for row in reader:
                vn.append(row)

            # close file
            voxelFile.close()

            for s in range(1, numSpecies+1):
                # pick only the right step
                voxelNum=vn[timestep][s]

                # and put it in the right position
                n[xc[i]][yc[i]][zc[i]][s-1][ts]=int(voxelNum)

    # and return the array
    return n

def read_multiple_mesord(filenameshape, sizex, numRuns, geometryname='geometry.txt'):

    # read in first to suss it out
    state = read_mesord_full(filenameshape.format(0), geometryname, sizex)

    # create array
    stateShape = shape(state)
    n = zeros((numRuns, stateShape[0], stateShape[1], stateShape[2]), dtype=float)
    n[0,:,:,:] = state

    for i in range(1, numRuns):
        state = read_mesord_full(filenameshape.format(i), geometryname, sizex)
        n[i, :, :, :] = state

    return n

def read_mesord_full(filename, geometryname, sizex):
    # read in state
    state = read_mesord(filename, sizex)
    
    # read in geometry file
    datafile=open(geometryname, 'r')

    #open reader
    reader=csv.reader(datafile, delimiter=" ", quoting=csv.QUOTE_ALL, skipinitialspace=True)

    nlist=[]

    for row in reader:
        # throw out last element (it's a quoted whitespace)
        nlist.append([int(float(x)-0.5) for x in row[:len(row)-1]])

    #close it
    datafile.close()

    # only get those where the z coordinate is 0
    nnew = [x for x in nlist if x[2]==0.0]

    # find minimum
    nmin = min(nnew)

    # and add it
    nlist = [[x[0]+abs(nmin[0]), x[1]+abs(nmin[1])] for x in nnew]

    # create array
    tmax=shape(state)[0]
    n = zeros((tmax, sizex,sizex), dtype=float)
    
    # reshape state
    s = reshape(state, (tmax, sizex*sizex))

    # and save in the right order
    for i in range(0, sizex*sizex):
        n[:,nlist[i][0], nlist[i][1]] = s[:, i]
                
    return n

def read_mesord(filename, sizex):
    datafile=open(filename, 'r')

    #open reader
    reader=csv.reader(datafile, delimiter=" ", quoting=csv.QUOTE_NONNUMERIC, skipinitialspace=True)

    nlist=[]

    for row in reader:
        # throw out last element (it's a quoted whitespace)
        nlist.append(row[:len(row)-1])

        #close it
    datafile.close()

    # create array
    tmax=shape(nlist)[0]
    n = zeros((tmax, sizex,sizex), dtype=float)

    for t in range(0, tmax):
        for i in range(0, sizex):
            for j in range(0, sizex):
                # need to sum over z
                n[t,i,j] = (nlist[t][i*sizex*2+2*j]+nlist[t][i*sizex*2+2*j+1])
    return n

def read_crank_nicholson(filename, sizex, split=' '):

    # open datafile
    datafile=open(filename, 'r')

    # create initial output array
    n=zeros((1,sizex), dtype=float)
    step = 0

    # loop over times
    while True:
        # read in line
        line = datafile.readline()

        # EOF?
        if len(line)==0:
            break

        # create temp array for this step
        ntemp = zeros((1, sizex), dtype=float)

        # split it
        linelist = line.split(split)

        # to array
        for i in range(0, sizex):
            ntemp[0,i] = float(linelist[i])

        # append temporary array to target array
        if (step == 0):
            n[0,:]=ntemp[0,:]
        elif (step > 0):
            n=append(n, ntemp, axis=0)

        # increase step counter
        step = step +1

    # return array
    return n

def read_gmp_individuals_all(filenamebase, species, runNo=0,  deterministic=False, verbose=False):
    # read in file
        # check if filename has ".h5" extension and add if necessary
    # (this is only for backward compability)
    filename = os.path.splitext(filenamebase)[0]+".h5"

    if verbose:
        print("Reading in {0:s}.".format(filename))

    f=h5py.File(filename, 'r')

    # model group
    modelGroup = f['Model']

    # find out if we're using old version or new version
    if ('numRuns' in modelGroup.attrs):
        # new file version
        if verbose:
            print("Detected new HDF5 structure. Reading in.")

        # get number of runs
        numRuns = modelGroup.attrs['numRuns'][0]

        # read in list of species
        numSpecies=modelGroup.attrs['numSpecies'][0]
        specieslist=modelGroup['Species'][:]
        if sys.version_info >= (3, 0):
            specieslist = array([x.decode() for x in specieslist])

        # get array size
        ft = f
        if 'Jobs' in ft:
            if verbose:
                print("Found Jobs key!")
            f = (ft['Jobs'])['1']
            if verbose:
                print("Done")
        else:
            f = ft

        runGroup = f['Runs']

        firstRun = runGroup['0']
        numDumps = len(firstRun)
        timeStamps = zeros((numDumps), dtype=float)

        # get run group
        currentRunGroup = runGroup['{0:d}'.format(runNo)]

        # result dictionary
        result = {}

        # go through all datasets and save in dictionary
        for name in currentRunGroup:
            dataGroup = currentRunGroup[name]
            igroup = dataGroup['Individuals']
            if (len(igroup)>0):
                idata = igroup[species]            
                result[name] = idata[:]

        return result

def read_gmp_individuals(filenamebase, dumpNo, species, runNo=0,  deterministic=False, verbose=False):

    # read in file
    # read in file
        # check if filename has ".h5" extension and add if necessary
    # (this is only for backward compability)
    filename = os.path.splitext(filenamebase)[0]+".h5"

    if verbose:
        print("Reading in {0:s}.".format(filename))

    f=h5py.File(filename, 'r')

    # model group
    modelGroup = f['Model']

    # find out if we're using old version or new version
    if ('numRuns' in modelGroup.attrs):
        # new file version
        if verbose:
            print("Detected new HDF5 structure. Reading in.")

        # get number of runs
        numRuns = modelGroup.attrs['numRuns'][0]

        # read in list of species
        numSpecies=modelGroup.attrs['numSpecies'][0]
        specieslist=modelGroup['Species'][:]
        if sys.version_info >= (3, 0):
            specieslist = array([x.decode() for x in specieslist])

        # get array size
        ft = f
        if 'Jobs' in ft:
            if verbose:
                print("Found Jobs key!")
            f = (ft['Jobs'])['1']
            if verbose:
                print("Done")
        else:
            f = ft

        runGroup = f['Runs']

        firstRun = runGroup['0']
        numDumps = len(firstRun)
        timeStamps = zeros((numDumps), dtype=float)

        # get run group
        currentRunGroup = runGroup['{0:d}'.format(runNo)]

        # now iterate through until we get the desired dump no ..
        # this is a bit cumbersome
        iDump = 0
        result = 0
        found = False

        for name in currentRunGroup:
            if verbose:
                print("Reading in group {0:s} for dump nr {1:d}.".format(name, iDump))
            dataGroup = currentRunGroup[name]

            if (iDump == dumpNo):
                if verbose:
                    print("Dump found. Reading in species {0:s}.".format(species))
                igroup = dataGroup['Individuals']
                if (len(igroup)>0):
                    idata = igroup[species]
                    result = idata[:]
                    found = True
                else:
                    print("WARNING: Empty dataset at dump {0:s} ({1:d})".format(name, iDump))

            iDump = iDump +1

        if (not found):
            print("ERROR: Can't find dataset!")

        return gpgmp.individuals.split_keys(result)
