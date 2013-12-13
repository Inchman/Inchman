# for gpgmp job control
import gpgmp.job
# matplotlib
import matplotlib.pyplot as plt

# for io
import gpgmp.io

import gpgmp.helpers

# numpy
import numpy

# run the job
def run(testdir, plot=False, norun=False, executable='gpgmp', deterministic=False, aqss=False):

    # parameters we use
    length = 1000.
    nx=64
    tmax=1800.
    k1=9.4096e+09
    k2=9.4096e+10
    k3=1.27529e-13
    k4=1.06274e-13
    diffusionConstant=100.
    xminA = 0.
    xmaxA = 900.
    xminB = 400.
    xmaxB = 1000.

    if deterministic:
          if aqss:
                solver = 4
          else:
                solver =3
    else:
          solver = 2 
    
    # run it and get the timing
    runTime=0
    if norun==False:
        runTime = gpgmp.job.run_job(testdir, executable,
                                    parameters='{0:g} {1:g} {2:g} {3:g} {4:g} {5:g} {6:g} {7:g} {8:g}'.format(diffusionConstant,
                                                                                                              k1,k2,k3,k4,
                                                                                                              xminA,xmaxA,xminB,xmaxB),
                                    problem=3, length=length, nx=nx, runtime=tmax, dtout=50., solver=solver, datafile='test_localised_ab.h5')

    # read out the results
    n, times, species, nruns = gpgmp.io.read_gmp_hdf5('test_localised_ab', deterministic=deterministic)

    # sort it after time
    timesarg=numpy.argsort(times)
    sampleTime=times[timesarg[-1]]

    # average over runs
    nmean=numpy.mean(n[:,timesarg[-1],:,:,:], axis=0)
    nstate=numpy.mean(nmean, axis=2)

    # read in deterministic results (from the Crank-Nicholson solver)
    concA = numpy.array([2.60316e-12, 2.60316e-12, 2.60314e-12, 2.60309e-12, 2.60295e-12, 2.60261e-12, 2.60178e-12, 2.5997e-12, 2.5946e-12, 2.58221e-12, 2.55264e-12, 2.48443e-12, 2.33641e-12, 2.05087e-12, 1.61921e-12, 1.29519e-12, 1.05499e-12, 8.61224e-13, 6.94822e-13, 5.51069e-13, 4.31271e-13, 3.36188e-13, 2.63618e-13, 2.09328e-13, 1.68807e-13, 1.38287e-13, 1.14964e-13, 9.68452e-14, 8.23855e-14, 6.79906e-14, 3.25562e-15, 1.44636e-16, 6.46373e-18, 6.46373e-18])
    concB = numpy.array([1.39556e-20, 1.39556e-20, 4.53566e-20, 1.78812e-19, 7.14592e-19, 2.85812e-18, 1.14312e-17, 4.57115e-17, 1.82709e-16, 7.29459e-16, 2.90432e-15, 1.14873e-14, 4.47386e-14, 1.68341e-13, 5.90424e-13, 8.62843e-13, 1.12557e-12, 1.43963e-12, 1.85122e-12, 2.40178e-12, 3.12518e-12, 4.04426e-12, 5.17124e-12, 6.51118e-12, 8.0657e-12, 9.83526e-12, 1.18201e-11, 1.40205e-11, 1.64367e-11, 1.90663e-11, 2.186e-11, 2.3762e-11, 2.47143e-11, 2.47143e-11])

    # coordinate systems
    dx = length/nx
    xgmp=(numpy.mgrid[0:nx])*dx
    dxcr = length/(32+1)
    xcg = (numpy.mgrid[0:32+2])*dxcr

    # and plot it if needed
    if plot:
        plt.clf()
        plt.plot(xgmp, gpgmp.helpers.numberToConcentration(nstate[0], dx), "ro")
        plt.plot(xgmp, gpgmp.helpers.numberToConcentration(nstate[1], dx), "bo")
        plt.plot(xcg, concA, "r")
        plt.plot(xcg, concB, "b")
        plt.savefig("localized_aplusb_reaction.eps", format="eps")

    return True, runTime
