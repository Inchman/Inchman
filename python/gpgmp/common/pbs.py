'''
Created on 22/10/2012

@author: matthias
'''
import sys
import time
import pickle
import subprocess
import os
import os.path
import datetime

# submit a job to the PBS queue
def run(options, executable, datadir, parameters):
    # create the command to execute
    command = executable + ' ' + parameters;
        
    #change to full directory
    os.chdir(datadir)
    
    # we need to create the pbs script
    scriptname='pbs_run_gmp'
    f = open(datadir+scriptname, 'w')

    # we need to convert walltime into hr::mn::sc
    min, sec = divmod(options['pbsWalltime'], 60)
    h, min = divmod(min, 60)
    walltime = "{0:d}:{1:02d}:{2:02d}".format(h, min, sec)

    # this is the PBS script template
    template = '''#!/bin/bash
# To give your job a name, replace "MyJob" with an appropriate name
#$ -N {2:s}

# Redirect stdout/stderr
#$ -e {0:s}stderr.txt
#$ -o {0:s}stdout.txt

# set your minimum acceptable walltime=hours:minutes:seconds
#$ -l h_rt={4:s}

# Changes directory to your execution directory
cd {0:s}

# set python search path
setenv PYTHONPATH /nfs/monash/home/mvigeliu/checkout/gpgmp/python/

# load required modules
module load python/2.7.2

# Command to run a job, either mpi or serial :
python -mcommon.pbs "{1:s}" "{3:s}"

# do any postprocessing needed
'''.format(datadir, command, os.path.basename(executable), options['timingFile'], walltime)

    f.write(template)

    # and close the file
    f.close()

    # and submit to PBS
    subprocess.call(['qsub', scriptname])
    
    return 0

# When scripted, just execute the command given as the first argument, time it
# and save output time
if __name__ == "__main__":
    # print command line arguments
    argv = sys.argv[1:]
            
    # and time it
    ts = time.time()
    subprocess.call(argv[0], shell=True)
    t = time.time()-ts
    
    # write timing information 
    timingfile = open(argv[1], 'w')
    pickle.dump(t, timingfile)
    timingfile.close()
