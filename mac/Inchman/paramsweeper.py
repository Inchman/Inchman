#!/usr/bin/python
#
# Inchman Parameter Sweeper
#
# Created by Aidan Lane, 6 Feb 2013
# Copyright (c) 2013 Monash University
#
# This script creates and executes a cross-product of all parameters within an
# Inchman annotated SBML file. We call this a full parameter sweep.
#
# This has been designed to work in the same was as Nimrod/G, as detailed at:
# http://messagelab.monash.edu.au/NimrodG/ExperimentManagement
# The following technical description is based on that document.
#
# domain - the following paremeter domains are available:
#
#   single (value)
#     Parameter has only a single value. This is the default domain,
#     i.e. it is used if no other domain has been specified.
#     The syntax for this domain is default <value>
#
#   range
#     This parameter has a number of values (points) that are within a specified
#     range, i.e. values are generated within the lower and upper bounds
#     (inclusive) for this parameter.
#     For this type of parameter we also need to specify how to generate
#     parameter values and there are two ways to do this.
#     The first one is to specify a number of uniformly distributed points where
#     the second one is to use a step value to generate points.
#     The range domain is:
#          range from <value> to <value> points <value> 
#          range from <value> to <value> step <value>
#
#   random
#     This domain of parameter has a number of random points generated between
#     specified lower and upper bounds from a uniform distribution.
#     Important type information:
#          integer : Discrete uniform distribution over the closed interval [from, to]
#            float : Floats uniformly distributed over the half-open [from, to)
#     The random domain is:
#          random from <value> to <value> [points <value>]
#     If the 'points' option is specified, 'points' random numbers are generated
#     as a dimension/degree-of-freedom for your parameter set, so each random
#     number in the set will be paired against each combination of all other
#     parameters. Without the 'points' option a single distinct random value
#     will be generated for every job.
#     Examples:
#     1. random from 1.0 to 15.0 points 6
#        This will yield 6 randomly generated values between 1.0 and 15.0, and
#        increase the parameter-set cross product size 6-fold.
#     2. random from 1.0 to 15.0
#        In this case there will be a distinct random value generated for every
#        job of the experiment at the time when your Inchman file is parsed.
#        For the float datatype it is unlikely any jobs will share the same
#        value. The parameter-set size remains unaffected.
#
# FUTURE: allow for different random distribution types
#

import sys, collections, copy, subprocess, numpy
import xml.etree.ElementTree as ET

class Param:
    pass
    
    
def get_attrib_as_type(element, value_name, type_name):
    v = element.get(value_name, None)
    if v:
        v = float(v) if (type_name == 'float') else int(v)
    return v


def find_params(filename):
    param_list = []
    tree = ET.parse(filename)
    root = tree.getroot()
    sbml_ns = root.tag[1:].split('}')[0] # SBML namespace (allows for different versions - e.g. http://www.sbml.org/sbml/level2/version3)
    for sbml_p in root.findall('{%s}model/{%s}listOfParameters/{%s}parameter' % (sbml_ns, sbml_ns, sbml_ns)):
        im_p = sbml_p.find('{%s}annotation/{http://www.csse.monash.edu.au/~berndm/inchman}parameter' % sbml_ns)
        p = Param()
        p.id           = sbml_p.get('id')
        p.value        = sbml_p.get('value')
        p.domain       = im_p.get('domain')
        p.type_name    = im_p.get('type')
        p.value_from   = get_attrib_as_type(im_p, 'from',   p.type_name)
        p.value_to     = get_attrib_as_type(im_p, 'to',     p.type_name)
        p.value_step   = get_attrib_as_type(im_p, 'step',   p.type_name)
        p.value_points = get_attrib_as_type(im_p, 'points', 'integer') # force always as integer
        #print "Found parameter: %s" % p.__dict__
        param_list.append(p)
        
    return param_list


def sweep_params(params):
    # ordering isn't critical, but it's nice for the user :)
    return recursive_sweep_params(collections.OrderedDict(), params)

# recursive function to build job list by sweeping over the parameter list
def recursive_sweep_params(current_job, remaining_params):
    if remaining_params == []:
        return [current_job] # end of branch -> job definition is complete -> return it
    
    popped_remaining_params = copy.deepcopy(remaining_params)
    p = popped_remaining_params.pop(0) # pop front / index 0
    
    values = []
    
    if p.domain == 'single' or not p.domain: # 'single' and if domain not specified
        values = [p.value];
    
    elif p.domain == 'range':
        if p.value_step and p.value_step > 0 and p.value_from > p.value_to:
            raise Exception("Invalid step value (%s), as from value (%s) > to value (%s) but step is positive" % (p.value_step, p.value_from, p.value_to))
        num_points = p.value_points if p.value_points else ((p.value_to-p.value_from)/p.value_step)+1 # |lin_points| includes the endpoint
        values = numpy.linspace(p.value_from, p.value_to, num_points, endpoint=True) # |lin_points| includes the endpoint
        
    elif p.domain == 'random':
        num_points = p.value_points if p.value_points else 1
        if p.type_name == 'integer':
            # Discrete uniform distribution over the closed interval [low, high]
            values = numpy.random.random_integers(p.value_from, p.value_to, num_points)
        elif p.type_name == 'float':
            # Floats uniformly distributed over the half-open [from, to)
            values = (p.value_to - p.value_from) * numpy.random.random_sample(num_points) + p.value_from
        else:
            raise Exception("Invalid parameter type: %s" % p.type_name)
    
    else:
        raise Exception("Invalid parameter domain: %s" % p.domain)
    
    new_jobs = []
    for v in values:
        j = copy.deepcopy(current_job)
        j[p.id] = float(v) if p.type_name=='float' else int(v) # add/update current param and force correct type
        new_jobs += recursive_sweep_params(j, popped_remaining_params)

    return new_jobs
    

def print_jobs(jobs):
    print "Queued Jobs:"
    for idx, job in enumerate(jobs):
        print "  %3d: " % idx,
        for k, v in job.iteritems():
            print "%6s=%-6s" % (k, str(v)),
        print "" # force newline
        
        
def exec_jobs(inchman_exec, filename, jobs):
    for idx, job in enumerate(jobs):
        params = []
        for k, v in job.iteritems():
            params += ['--p', '%s=%s' % (k, str(v))]
        print "Executing Inchman job %d with: %s" % (idx, str(params))
        subprocess.call([inchman_exec, '--job', str(idx)] + params + [filename])


def exec_sweep(inchman_exec, filename):
    params = find_params(filename)
    jobs = sweep_params(params)
    print_jobs(jobs)
    exec_jobs(inchman_exec, filename, jobs)


def main(argv):
    exec_sweep(argv[1], argv[2])


if __name__ == '__main__':
    main(sys.argv)
