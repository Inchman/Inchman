import numpy
import gpgmp.io

def split_keys(keys):
    # create bit mask array
    mask=numpy.ones_like(keys)
    mask.fill(0xffffffff)

    # first is cell index, second is id
    return keys & mask, keys >> 32

def read_trajectory(filenamebase, species, particleId, listInactive=False, runNo=0):
    
    # read in all particles
    poslist = gpgmp.io.read_gmp_individuals_all('output', species, runNo)

    # create empty result array
    rtime = numpy.zeros((0))
    result = numpy.zeros((0))

    # go through all dumps
    for name in poslist:
        # get time
        time =  float(name.partition('_')[2])
        
        # get positions and ids
        pos, ids = split_keys(poslist[name])

        # find my particle
        index = numpy.where(ids == particleId)[0]

        # if it's there, we append it to results
        if (len(index) > 0):
            rtime = numpy.append(rtime, time)
            result = numpy.append(result, pos[index])
        else:
            if listInactive:
                # append max float
                rtime = numpy.append(rtime, time)
                result = numpy.append(result, 0xffffffff)

    # and return it
    return rtime, result
