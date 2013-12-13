# This runs the simple A+B -> A reaction test for the FPKMC/RW algorithm
# for gpgmp job control
import gpgmp.job

# matplotlib
import matplotlib.pyplot as plt

# for io
import gpgmp.io

# numpy
import numpy

# run the job
def run(testdir, cldir, plot=False, executable='rw', norun=False):
    
    # run it and get the timing
    if norun==False:            
        runtime = gpgmp.job.run_rw(testdir, executable, cldir, 'SIMPLE_AB')
    else:
        runtime = 0.

    # read in results
    tot=numpy.genfromtxt("total.dat")
        
    # take time points from output
    time = tot[:,0]
        
    # compute reaction rate from reaction radius
    k1 = 8e-3*4*numpy.pi*(1+1)
    volsub = 2**3

    # compute variance
    var=numpy.exp(-k1/volsub*500.*time)*500.*(1.-numpy.exp(-k1/volsub*500.*time))

    if plot:
        # plot it
        plt.clf()
        plt.title('A+B reaction problem without diffusion')
        plt.xlabel('t [s]')
        plt.ylabel('N')


        # and overplot it
        plt.plot(time, numpy.exp(-k1/volsub*500.*time)*500, "b-", label='analytic')

        # plot variance lines
        plt.plot(time, numpy.exp(-k1/volsub*500.*time)*500+numpy.sqrt(var), "b--", label='variance')
        plt.plot(time, numpy.exp(-k1/volsub*500.*time)*500-numpy.sqrt(var), "b--")
        plt.plot(time, numpy.exp(-k1/volsub*500.*time)*500+2.*numpy.sqrt(var), "b:")
        plt.plot(time, numpy.exp(-k1/volsub*500.*time)*500-2.*numpy.sqrt(var), "b:")

        # plot rw solution
        plt.plot(tot[:,0], tot[:,2]-500, "r-", linewidth=2, label='rw')

        plt.legend(loc='upper right')
        plt.savefig("simple_ab.eps", format="eps")

    # check if passed
    varlower = numpy.exp(-k1/volsub*500.*time[-1])*500-2.*numpy.sqrt(var[-1])
    varupper = numpy.exp(-k1/volsub*500.*time[-1])*500+2.*numpy.sqrt(var[-1])
    if ((tot[-1,2]-500 > varlower) and (tot[-1,2]-500 < varupper)):
        passed = True
    else:
        passed = False

    return passed, runtime
