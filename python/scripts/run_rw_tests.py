#!/usr/bin/python

# interface to os
import os
import sys

# parse options
from optparse import OptionParser

import time

# for gpgmp job control
import gpgmp.job

# problems
import fpkmc.test.simpleab
import fpkmc.test.random_drift

def main(argv=None):
    # get arguments from command line if none are given
    # (if it's started as a script)
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
                      dest="path", default=cwd+'/')
    parser.add_option("--datadir", action="store", type="string",
                      dest="datadir", default=cwd+'/')
    parser.add_option("--use-pbs", action="store_true", dest="use_pbs")
    parser.add_option("--norun", action="store_true", dest="norun", default=False)
    parser.add_option("--all", action="store_true", dest="all", default=False)

    # problems
    parser.add_option("--simple-ab", action="store_true", dest="simple_ab")
    parser.add_option("--random-drift", action="store_true", dest="random_drift")
    parser.add_option("--random-drift-discrete", action="store_true", dest="random_drift_discrete")

    # query options and set them
    (options, args) = parser.parse_args(argv)

    # set problem sets
    if options.all:
        options.simple_ab = True

    # set paths
    path = options.path
    datadir = options.datadir
    pbs = options.use_pbs

    print "Welcome to the automated RW test suite."
    print "Path to rw is {0:s}".format(path)
    print "Path to data is {0:s}".format(datadir)
    if pbs:
        print "We use the PBS queuing system to execute jobs."

    # make test directory if runs are performed
    if options.norun:
        testdir = datadir+'results/'
        print "Not running tests. Results directory is {0:s}".format(testdir)
        os.chdir(testdir)
        if windows:
            binname = 'rw.exe'
        else:
            binname = 'rw'

    else:
        binname = 'rw'
        testdir = gpgmp.job.prepare_directories(path, datadir, 'results/', executable=binname)

    # finally add the whole path to the executable so it can be found
    executable = path+'/'+binname

    # simple A + B -> A
    if options.simple_ab:
        print "Running simple A+B reaction test ...\t\t\t",
        passed, runtime = fpkmc.test.simpleab.run(testdir, plot=True, executable=executable, cldir=path, norun=options.norun)
        if passed:
            print "\x1b[32m \x1b[1m passed.\x1b[0m",
        else:
            print "\x1b[31m \x1b[1m failed.\x1b[0m",
            
        if (options.norun):
            print " (did not run test)"
        else:
            print " (runtime: {0:g} seconds)".format(runtime)

    # random drift problems
    if options.random_drift:
        print "Running random drift test ...\t\t\t",
        passed, runtime = fpkmc.test.random_drift.run(testdir, plot=True, executable=executable, cldir=path, norun=options.norun)
        if passed:
            print "\x1b[32m \x1b[1m passed.\x1b[0m",
        else:
            print "\x1b[31m \x1b[1m failed.\x1b[0m",
            
        if (options.norun):
            print " (did not run test)"
        else:
            print " (runtime: {0:g} seconds)".format(runtime)

    if options.random_drift_discrete:
        print "Running discrete random drift test ...\t\t",
        passed, runtime = fpkmc.test.random_drift.run(testdir, plot=True, executable=executable, cldir=path, norun=options.norun, discrete=True)
        if passed:
            print "\x1b[32m \x1b[1m passed.\x1b[0m",
        else:
            print "\x1b[31m \x1b[1m failed.\x1b[0m",
            
        if (options.norun):
            print " (did not run test)"
        else:
            print " (runtime: {0:g} seconds)".format(runtime)

# only run when scripted
if __name__ == "__main__":
    # run main
    main()
