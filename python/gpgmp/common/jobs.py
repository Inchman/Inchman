'''
Created on 18/10/2012

@author: matthias
'''
import os
import errno
import uuid
import glob
import shutil
import sys
import subprocess
import time
import pickle
import common.pbs

def prepare_directories(options, extension, subversiondir=None):
    # extract datadir from options
    datadir = options['datadir']
    
    print("Creating directory {0:s}.".format(datadir+extension))

    # recursively create directory
    try:
        os.makedirs(datadir+extension)
    except OSError as err:
        # check if the error is because the dir exists
        if (err.errno == errno.EEXIST):
            print("Directory \"{0:s}\" exists. Moving contents to backup dir.".format(datadir+extension))

            # in that case we just create a unique directory and copy all the old stuff there
            olddir = datadir+'/'+str(uuid.uuid4())+'.backup/'
            os.makedirs(olddir)
            # need to expand wildcards first
            for file in glob.glob(datadir+extension+'*'):
                shutil.move(file, olddir)
        else:
            print("Error \"{0:s}\" while creating directory.".format(err.strerror))
            sys.exit(1)
    
    # change into dir
    os.chdir(datadir+extension)

    # put in subversion information (if requested)
    if subversiondir:
        # open info file
        infofile = open("svninfo", 'w')
        
        # run svn info in subversiondir
        subprocess.call(["svn", "info", subversiondir], stdout=infofile)
        
        # close infofile
        infofile.close()
        
    # and return a path to the full dir
    return datadir+extension

def runNode(options, executable, datadir, parameters):
    """Runs a job on the current machine."""
    command = executable + ' ' + parameters;

    print("Command is {0:s}".format(command))

    #change to full directory
    os.chdir(datadir)
    
    # run only if it's not set to fake mode
    if not options['fakeRun']:
        # create files to capture output
        outfile = open(options['outlog'], 'w')
        errfile = open(options['errlog'], 'w')
    
        # and time it
        ts = time.time()
        subprocess.call(command, stdout=outfile, stderr=errfile, shell=True)
        t = time.time()-ts
    
        # close log files
        outfile.close()
        errfile.close()
    else:
        print("Running (fakemode): {0:s} in directory {1:s}.".format(command, datadir))
        t = 1.
        
    # write timing information 
    timingfile = open(options['timingFile'], 'w')
    pickle.dump(t, timingfile)
    timingfile.close()
    
def run(options, executable, datadir, parameters):
    """Runs a job specified by executable. The name of the executable must contain the complete path.
    If options['pbs'] is set to True, the jobs will be submitted to the PBS scheduler."""
    if options['pbs']:
        common.pbs.run(options, executable, datadir, parameters)
    else:
        runNode(options, executable, datadir, parameters)
