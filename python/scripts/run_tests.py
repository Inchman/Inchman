#!/usr/bin/python

# interface to os
import os
import sys

# parse options
from optparse import OptionParser

# need math
import math

import shutil

import time

# for gpgmp job control
import gpgmp.job

# test modules
import gpgmp.test.abreaction
import gpgmp.test.fisher_problem
import gpgmp.test.homogeneous_diffusion
import gpgmp.test.localised_abreaction
import gpgmp.test.infinite_slab
import gpgmp.test.homogeneous_drift_diffusion
import gpgmp.test.random_drift
import gpgmp.test.multiplicative_noise
import gpgmp.test.ornstein_uhlenbeck
import gpgmp.test.nonlinear
import gpgmp.test.annihilation_2d
import gpgmp.test.annihilation_2d_drift
import gpgmp.test.slit

#TODO: support all opencl paths (-c command line argument..)
#TODO: support all init script paths (-i command line argument..)

def main(argv=None):
    # get arguments from command line if none are given (if it's started as a script)
    if argv is None:
        argv = sys.argv[1:]

    # get current directory as reference for standard tree
    cwd = os.getcwd()

    # check if running under windows
    windows = False
    if os.name=='nt':
        windows=True

    # add all possible options
    parser = OptionParser()
    parser.add_option("-p", "--path", action="store", type="string",
                      dest="gpgmppath", default=cwd+'/')
    parser.add_option("--datadir", action="store", type="string",
                      dest="datadir", default=cwd+'/')
    parser.add_option("--use-pbs", action="store_true", dest="use_pbs")
    parser.add_option("--opencl", action="store_true", dest="opencl")
    parser.add_option("--norun", action="store_true", dest="norun", default=False)
    parser.add_option("--all", action="store_true", dest="all", default=False)
    parser.add_option("--deterministic", action="store_true", dest="deterministic", default=False)
    parser.add_option("--aqss", action="store_true", dest="aqss", default=False)

    # problems
    parser.add_option("--aplusb", action="store_true", dest="aplusb")
    parser.add_option("--individual-aplusb", action="store_true", dest="individual_aplusb")
    parser.add_option("--fisher", action="store_true", dest="fisher")
    parser.add_option("--homogeneous-diffusion", action="store_true", dest="homogeneous_diffusion")
    parser.add_option("--homogeneous-drift", action="store_true", dest="homogeneous_drift")
    parser.add_option("--random-drift", action="store_true", dest="random_drift")
    parser.add_option("--localised-ab", action="store_true", dest="localised_aplusb")
    parser.add_option("--infinite-slab", action="store_true", dest="infinite_slab")
    parser.add_option("--multiplicative-noise", action="store_true", dest="multiplicative_noise")
    parser.add_option("--ornstein-uhlenbeck", action="store_true", dest="ornstein_uhlenbeck")
    parser.add_option("--nonlinear", action="store_true", dest="nonlinear")
    parser.add_option("--all-homogeneous", action="store_true", dest="all_homogeneous")
    parser.add_option("--all-inhomogeneous", action="store_true", dest="all_inhomogeneous")
    parser.add_option("--annihilation-2d", action="store_true", dest="annihilation_2d")
    parser.add_option("--annihilation-2d-drift", action="store_true", dest="annihilation_2d_drift")
    parser.add_option("--slit", action="store_true", dest="slit")

    # query options and set them
    (options, args) = parser.parse_args(argv)

    # set problem sets
    if options.all:
        options.aplusb = True
        options.fisher = True
        options.homogeneous_diffusion = True
        options.homogeneous_drift = True
        options.localised_aplusb = True
        #options.infinite_slab = True
        options.multiplicative_noise = True
        options.ornstein_uhlenbeck = True
        options.nonlinear = True
        options.random_drift = True

    if options.all_homogeneous:
        options.aplusb = True
        options.fisher = True
        options.homogeneous_diffusion = True
        options.localised_aplusb = True
        #options.infinite_slab = True

    if options.all_inhomogeneous:
        options.homogeneous_drift = True
        options.multiplicative_noise = True
        options.ornstein_uhlenbeck = True
        options.nonlinear = True

    # set paths
    path = options.gpgmppath
    datadir = options.datadir
    pbs = options.use_pbs

    print "Welcome to the automated gpgmp test suite."
    print "Path to gpgmp is {0:s}".format(path)
    print "Path to data is {0:s}".format(datadir)
    if pbs:
        print "We use the PBS queuing system to execute jobs."

    # time it
    tstart = time.time()

    # make test directory if runs are performed
    if options.norun:
        testdir = datadir+'results/'
        print "Not running tests. Results directory is {0:s}".format(testdir)
        os.chdir(testdir)
        if windows:
            binname = 'tests.exe'
        else:
            binname = 'tests'

    else:
        # copy opencl files if needed
        if options.opencl:
            if windows:
                binname = 'tests.exe'
            else:
                binname = 'tests'
            testdir = gpgmp.job.prepare_directories(path, datadir, 'results/', executable=binname)
        else:
            binname = 'gpgmp'
            testdir = gpgmp.job.prepare_directories(path, datadir, 'results/', executable=binname)


    # finally add the whole path to the executable so it can be found
    executable = path+'/'+binname

    # ---------------- homogeneous problems

    # A + B test
    if options.aplusb:
        print "Running A+B reaction test ...\t\t\t",
        sys.stdout.flush()
        passed, runtime = gpgmp.test.abreaction.run(testdir, plot=True, executable=executable,
                                                    deterministic=options.deterministic, aqss=options.aqss)
        if passed:
            print "\x1b[32m \x1b[1m passed.\x1b[0m",
        else:
            print "\x1b[31m \x1b[1m failed.\x1b[0m",

        print "Runtime: {0:g} seconds.".format(runtime)

    # 2D annihilation without drift
    if options.annihilation_2d:
        print "Running 2D annihilation without drift ...\t\t\t",
        sys.stdout.flush()
        passed, runtime = gpgmp.test.annihilation_2d.run(testdir, plot=True, executable=executable,
                                                         deterministic=options.deterministic, aqss=options.aqss)
        print "\x1b[33m \x1b[1m (no comparison) \x1b[0m"

        print "Runtime: {0:g} seconds.".format(runtime)

    # localised A + B test
    if options.localised_aplusb:
        print "Running localised A+B reaction test ...\t\t",
        sys.stdout.flush()
        passed, runtime = gpgmp.test.localised_abreaction.run(testdir, plot=True, executable=executable,
                                                              deterministic=options.deterministic, aqss=options.aqss)
        if passed:
            print "\x1b[32m \x1b[1m passed.\x1b[0m",
        else:
            print "\x1b[31m \x1b[1m failed.\x1b[0m",
            
        print "Runtime: {0:g} seconds.".format(runtime)

    # Fisher problem
    if options.fisher:
        # copy init file
        #shutil.copy(path+"init_fisher.py", datadir+'results/')
        print "Running Fisher problem test ...\t\t\t",
        sys.stdout.flush()
        passed, runtime = gpgmp.test.fisher_problem.run(testdir, path, plot=True, executable=executable,
                                                        deterministic=options.deterministic, aqss=options.aqss)
        if passed:
            print "\x1b[32m \x1b[1m passed.\x1b[0m",
        else:
            print "\x1b[31m \x1b[1m failed.\x1b[0m",

        print "Runtime: {0:g} seconds.".format(runtime)

    # homogeneous diffusion
    if options.homogeneous_diffusion:
        print "Running homogeneous diffusion test ...\t\t",
        sys.stdout.flush()
        passed, runtime = gpgmp.test.homogeneous_diffusion.run(testdir, path, plot=True, executable=executable,
                                                               deterministic=options.deterministic)
        if passed:
            print "\x1b[32m \x1b[1m passed.\x1b[0m",
        else:
            print "\x1b[31m \x1b[1m failed.\x1b[0m",

        print "Runtime: {0:g} seconds.".format(runtime)

    # semi-infinite slab problem
    if options.infinite_slab:
        if options.deterministic:
            print "\x1b[33m\x1b[1mWarning: deterministic solver does not support source BCs.. skipping semi-infinite slab test.\x1b[0m"
        else:
            print "Running semi-infinite slab test ...\t\t",
            sys.stdout.flush()
            passed, runtime = gpgmp.test.infinite_slab.run(testdir, path, plot=True, norun=options.norun, executable=executable)
            if passed:
                print "\x1b[32m \x1b[1m passed.\x1b[0m",
            else:
                print "\x1b[31m \x1b[1m failed.\x1b[0m",

            print "Runtime: {0:g} seconds.".format(runtime)

    # ---------------- inhomogeneous problems
                                               
    # 2D annihilation with drift
    if options.annihilation_2d_drift:
        print "Running 2D annihilation with drift ...\t\t\t",
        sys.stdout.flush()
        passed, runtime = gpgmp.test.annihilation_2d_drift.run(testdir, plot=True, executable=executable,
                                                              deterministic=options.deterministic, aqss=options.aqss)

        print "\x1b[33m \x1b[1m (no comparison) \x1b[0m"

        print "Runtime: {0:g} seconds.".format(runtime)


    # homogeneous drift - diffusion
    if options.homogeneous_drift:
        print "Running inhomogeneous diffusion test ...\t",
        sys.stdout.flush()
        passed, runtime = gpgmp.test.homogeneous_drift_diffusion.run(testdir, path, plot=True, norun=options.norun, executable=executable,
                                                                     deterministic=options.deterministic, aqss=options.aqss)
        if passed:
            print "\x1b[32m \x1b[1m passed.\x1b[0m",
        else:
            print "\x1b[31m \x1b[1m failed.\x1b[0m",

        print "Runtime: {0:g} seconds.".format(runtime)

    # random drift test (individual)
    if options.random_drift:
        print "Running random drift test ...\t",
        sys.stdout.flush()
        passed, runtime = gpgmp.test.random_drift.run(testdir, path, plot=True, norun=options.norun, executable=executable,
                                                      deterministic=options.deterministic, aqss=options.aqss)
        if passed:
            print "\x1b[32m \x1b[1m passed.\x1b[0m",
        else:
            print "\x1b[31m \x1b[1m failed.\x1b[0m",

        print "Runtime: {0:g} seconds.".format(runtime)

    # individual A + B test
    if options.individual_aplusb:
        print "Running individual A+B reaction test ...\t\t\t",
        sys.stdout.flush()
        passed, runtime = gpgmp.test.abreaction.run_individual(testdir, plot=True, executable=executable, norun=options.norun)
        if passed:
            print "\x1b[32m \x1b[1m passed.\x1b[0m",
        else:
            print "\x1b[31m \x1b[1m failed.\x1b[0m",

        print "Runtime: {0:g} seconds.".format(runtime)

    # multiplicative noise
    if options.multiplicative_noise:
        print "Running multiplicative noise test ...\t\t",
        sys.stdout.flush()
        passed, runtime = gpgmp.test.multiplicative_noise.run(testdir, path, plot=True, norun=options.norun, executable=executable,
                                                              deterministic=options.deterministic, aqss=options.aqss)
        if passed:
            print "\x1b[32m \x1b[1m passed.\x1b[0m",
        else:
            print "\x1b[31m \x1b[1m failed.\x1b[0m",

        print "Runtime: {0:g} seconds.".format(runtime)

    # Ornstein-Uhlenbeck
    if options.ornstein_uhlenbeck:
        print "Running Ornstein-Uhlenbeck test ...\t\t",
        sys.stdout.flush()
        passed, runtime = gpgmp.test.ornstein_uhlenbeck.run(testdir, path, plot=True, norun=options.norun, executable=executable,
                                                            deterministic=options.deterministic, aqss=options.aqss)
        if passed:
            print "\x1b[32m \x1b[1m passed.\x1b[0m",
        else:
            print "\x1b[31m \x1b[1m failed.\x1b[0m",

        print "Runtime: {0:g} seconds.".format(runtime)

    # non linear
    if options.nonlinear:
        print "Running nonlinear test ...\t\t\t",
        sys.stdout.flush()
        passed, runtime = gpgmp.test.nonlinear.run(testdir, path, plot=True, norun=options.norun, executable=executable,
                                                   deterministic=options.deterministic, aqss=options.aqss)
        if passed:
            print "\x1b[32m \x1b[1m passed.\x1b[0m",
        else:
            print "\x1b[31m \x1b[1m failed.\x1b[0m",

        print "Runtime: {0:g} seconds.".format(runtime)

    # ------ model problems
    # non linear
    if options.slit:
        print "Running SLIT problem ...\t\t\t",
        sys.stdout.flush()
        passed, runtime = gpgmp.test.slit.run(testdir, path, plot=True, norun=options.norun, executable=executable,
                                                   deterministic=options.deterministic, aqss=options.aqss)

        print "\x1b[33m \x1b[1m (no comparison) \x1b[0m",

        print "Runtime: {0:g} seconds.".format(runtime)

    totaltime = time.time()-tstart
    print "Done. Total runtime {0:g} seconds.".format(totaltime)

# only run when scripted
if __name__ == "__main__":
    # run main
    main()
