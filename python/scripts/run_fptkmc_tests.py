#!/usr/bin/python

'''
Created on 18/10/2012

@author: matthias
'''
import sys

import common.options
import common.jobs

import fptkmc.models.abannihilation

# parse options
from optparse import OptionParser

def main(argv=None):
    # get arguments from command line if none are given (if it's started as a script)
    if argv is None:
        argv = sys.argv[1:]        

    # parse output for additional parameters
    parser = OptionParser()
    parser.add_option("--smoldyn", action="store_true", dest="smoldyn")
    parser.add_option("--smoldyn-accuracy", action="store_true", dest="smoldynaccuracy")
    parser.add_option("--timestep", action="store_true", dest="timestep")
    
    ops, options = common.options.getSystemOptions(argv, parser)

    
    # query options and set them
    #(options, args) = parser.parse_args(argv)

    # check if we need to run the sims
    if ops['runSimulations']:
        # run them
        print "Running FPTKMC/abannihilation.."
        
        # prepare directories
        datadir = common.jobs.prepare_directories(ops, 'results/', ops['fptkmcdir'])
        
        # and run a+b annihilation
        if options.smoldyn:
            fptkmc.models.abannihilation.runSmoldyn(ops, datadir)
        elif options.smoldynaccuracy:  
            fptkmc.models.abannihilation.runSmoldynAccuracy(ops, datadir)      
        else:
            fptkmc.models.abannihilation.run(ops, datadir)

    # check if we need to analyze
    if ops['analyzeSimulations']:
        # do it
        print "Analyzing FPTKMC/abannihilation."
        
        # get results dir
        resultsdir = ops['datadir']+'results/'

        # and plot
        if options.smoldyn:
            fptkmc.models.abannihilation.plotSingleSmoldyn(ops, ".")
        elif options.smoldynaccuracy:
            fptkmc.models.abannihilation.plotSmoldynAccuracy(ops, ops['datadir']+"/aplusb/smoldyn/accuracy/")
        elif options.timestep:
            fptkmc.models.abannihilation.plotTimeStep(ops, ops['datadir'])
        else:
            fptkmc.models.abannihilation.plotSingle(ops, resultsdir)

# only run when scripted
if __name__ == "__main__":
    # run main
    main()
