'''
Created on 18/10/2012

@author: matthias
'''

import platform
import os

# parse options
from optparse import OptionParser

def getSystemOptions(argv=None, parser=None):
    # set all option values to default
    # default means cluster (since its hard to
    # identify the cluster nodes..)
    sysoptions = {'datadirBase': '.',
                  'fptkmcdir': '/nfs/monash/home/mvigeliu/checkout/gpgmp/fptkmc/',
                  'fptkmcexec': '/nfs/monash/home/mvigeliu/checkout/gpgmp/fptkmc/build/main',
                  'gpgmpdir': '/nfs/monash/home/mvigeliu/checkout/gpgmp/gpgmp_cl/',
                  'gpgmpexec': '/nfs/monash/home/mvigeliu/checkout/gpgmp/gpgmp_cl/build/tests',
                  'gpgmpcldir': '/nfs/monash/home/mvigeliu/checkout/gpgmp/gpgmp_cl/src/',
                  'gpgmpinitdir': '/nfs/monash/home/mvigeliu/checkout/gpgmp/gpgmp_cl/tests/',
                  'inchmanexec': '/usr/local/bin/inchman-exec',
                  'inchmancldir': '/usr/local/share/inchman/cl_source/',
                  'inchmanexampledir': '/nfs/monash/home/mvigeliu/checkout/gpgmp/inchman/examples/',
                  'pbs': False,
                  'outlog': 'stdout.txt',
                  'errlog': 'stderr.txt',
                  'timingFile': 'runtime.pkl',
                  'runSimulations': False,
                  'analyzeSimulations': False,
                  'fakeRun': False,
                  'smoldynDir': '/nfs/monash/home/mvigeliu/software/smoldyn-2.28/build/',
                  'smoldynexec': '/nfs/monash/home/mvigeliu/software/smoldyn-2.28/build/smoldyn',
                  'flamedir': '/nfs/monash/home/mvigeliu/checkout/gpgmp/documentation/majority/data/flame_models/majority/',
                  'flameexec': '/nfs/monash/home/mvigeliu/checkout/gpgmp/documentation/majority/data/flame_models/majority/main',
                  'mathematicabin': '/opt/sw/Wolfram/Mathematica-8.0.1/Executables/math',
                  'majoritymathscripts': '/nfs/monash/home/mvigeliu/checkout/gpgmp/documentation/majority/notebooks/'
                  }

    sysoptions['datadir'] = os.getcwd()+'/'

    # check first if we recognize the system
    #hostname = platform.node()
    #if (hostname == "XibalbaII"):
    #    #print("Host XibalbaII found.")
    #    sysoptions['datadir'] = os.getcwd()+'/'
    #    sysoptions['fptkmcdir'] = '/home/matthias/checkout/fptkmc/'
    #    sysoptions['fptkmcexec'] = sysoptions['fptkmcdir']+'build/main'
    #    sysoptions['gpgmpdir'] = '/home/matthias/checkout/gpgmp_cl/'
    #    sysoptions['gpgmpexec'] = sysoptions['gpgmpdir']+'build/tests'
    #    sysoptions['mathematicabin'] = '/usr/local/bin/math'
    #    sysoptions['matlabbin'] = '/usr/local/MATLAB/R2013a/bin/matlab'
    #    sysoptions['majoritymathscripts'] = '/home/matthias/checkout/majority/notebooks/'

    #if (hostname == "Xibalba"):
    #    #print("Host Xibalba found.")
    #    sysoptions['datadir'] = os.getcwd()+'/'
    #    sysoptions['fptkmcdir'] = '/home/matthias/checkout/fptkmc/'
    #    sysoptions['fptkmcexec'] = sysoptions['fptkmcdir']+'build/main'
    #    sysoptions['gpgmpdir'] = '/home/matthias/checkout/gpgmp_cl/'
    #    sysoptions['gpgmpexec'] = sysoptions['gpgmpdir']+'build/tests'
    #    sysoptions['gpgmpcldir'] = sysoptions['gpgmpdir']+'src/'
    #    sysoptions['gpgmpinitdir'] = sysoptions['gpgmpdir']+'tests/'
    #    sysoptions['flameexec'] = '/home/matthias/checkout/majority/data/flame_models/majority/main'
    #    sysoptions['flamedir'] = '/home/matthias/checkout/majority/data/flame_models/majority/'
    #    sysoptions['flameexec'] = '/home/matthias/checkout/majority/data/flame_models/majority/main'
    #    sysoptions['inchmanexec'] = '/home/matthias/checkout/inchman/build/src/inchman/inchman-exec'
    #    sysoptions['inchmancldir'] = '/home/matthias/checkout/inchman/exec/src/gpgmp/src/'
    #    sysoptions['inchmanexampledir'] = '/home/matthias/checkout/inchman/examples/'
    #    sysoptions['mathematicabin'] = '/usr/local/bin/math'
    #    sysoptions['matlabbin'] = '/usr/local/MATLAB/R2013a/bin/matlab'
    #    sysoptions['majoritymathscripts'] = '/home/matthias/checkout/majority/notebooks/'
        
    # parse output for standard parameters
    if not parser:
        parser = OptionParser()
        
    parser.add_option("--run", action="store_true", dest="runSimulations")
    parser.add_option("--analyze", action="store_true", dest="analyzeSimulations")
    #parser.add_option("--fake", action="store_true", dest="fakeRun")
    parser.add_option("--pbs", action="store_true", dest="pbs")
    parser.add_option("--pbs-walltime", type="int", dest="walltime", default=600)
    parser.add_option("--datadir", type="string", dest="datadir")
    
    # query options and set them
    (options, args) = parser.parse_args(argv)
    if options.runSimulations:
        sysoptions['runSimulations'] = True
    if options.analyzeSimulations:
        sysoptions['analyzeSimulations'] = True
    #if options.fakeRun:
    #    sysoptions['fakeRun'] = True
    if options.pbs:
        sysoptions['pbs'] = True
    sysoptions['datadir'] = options.datadir
    sysoptions['pbsWalltime'] = options.walltime
    # print info
    #print("Using datadir {0:s}.".format(sysoptions['datadir']))

    # return options and arguments
    return sysoptions, options, args
