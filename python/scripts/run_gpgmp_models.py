import sys

# parse options
from optparse import OptionParser

# for job processing and options
import common.options
import common.jobs

import gpgmp.models.majority

# tests
import gpgmp.test.abreaction
import gpgmp.test.annihilation_2d
import gpgmp.test.annihilation_2d_drift
import gpgmp.test.fisher_problem
import gpgmp.test.homogeneous_diffusion
import gpgmp.test.homogeneous_drift_diffusion
import gpgmp.test.multiplicative_noise
import gpgmp.test.ornstein_uhlenbeck
import gpgmp.test.nonlinear
import gpgmp.test.random_drift

def main(argv=None):
    # get arguments from command line if none are given (if it's started as a script)
    if argv is None:
        argv = sys.argv[1:]        

    # parse output for additional parameters
    parser = OptionParser()

    # for the tests
    parser.add_option("--test-aplusb", action="store_true", dest="test_aplusb")
    parser.add_option("--test-aplusb-individual", action="store_true", dest="test_aplusb_individual")
    parser.add_option("--test-aplusb-individual-new", action="store_true", dest="test_aplusb_individual_new")
    parser.add_option("--test-annihilation-2d", action="store_true", dest="test_annihilation_2d")
    parser.add_option("--test-annihilation-2d-drift", action="store_true", dest="test_annihilation_2d_drift")
    parser.add_option("--test-fisher", action="store_true", dest="test_fisher")
    parser.add_option("--test-fisher-individual", action="store_true", dest="test_fisher_individual")
    parser.add_option("--test-homogeneous-diffusion", action="store_true", dest="test_homogeneous_diffusion")
    parser.add_option("--test-homogeneous-drift", action="store_true", dest="test_homogeneous_drift")
    parser.add_option("--test-homogeneous-drift-individual", action="store_true", dest="test_homogeneous_drift_individual")
    parser.add_option("--test-homogeneous-drift-field", action="store_true", dest="test_homogeneous_drift_field")
    parser.add_option("--test-multiplicative-noise", action="store_true", dest="test_multiplicative_noise")
    parser.add_option("--test-ornstein-uhlenbeck", action="store_true", dest="test_ornstein_uhlenbeck")
    parser.add_option("--test-nonlinear", action="store_true", dest="test_nonlinear")
    parser.add_option("--test-random-drift", action="store_true", dest="test_random_drift")

    # for the majority vote model
    parser.add_option("--majority", action="store_true", dest="majority")
    parser.add_option("--majority-flame-avoidance", action="store_true", dest="majority_flame_avoidance")
    parser.add_option("--majority-flame-velocity", action="store_true", dest="majority_flame_velocity")
    parser.add_option("--majority-flame-recognition", action="store_true", dest="majority_flame_recognition")
    parser.add_option("--majority-flame-diffusive", action="store_true", dest="majority_flame_diffusive")
    parser.add_option("--majority-flame-experiment", action="store_true", dest="majority_flame_experiment")
    parser.add_option("--majority-gpgmp-nondiffusive", action="store_true", dest="majority_gpgmp_nondiffusive")
    parser.add_option("--majority-gpgmp-avoidance", action="store_true", dest="majority_gpgmp_avoidance")
    parser.add_option("--majority-geoff", action="store_true", dest="majority_geoff")
    parser.add_option("--majority-toy-model", action="store_true", dest="majority_toy_model")
    parser.add_option("--majority-track-positions", action="store_true", dest="majority_track_positions")

    ops, options = common.options.getSystemOptions(argv, parser)

    # check if we need to run the sims
    if ops['runSimulations']:
        # prepare directories
        datadir = common.jobs.prepare_directories(ops, 'results/')

        # tests
        if options.test_aplusb:
            gpgmp.test.abreaction.run(ops, datadir)
        if options.test_aplusb_individual:
            gpgmp.test.abreaction.run(ops, datadir, individual=True)
        if options.test_aplusb_individual_new:
            gpgmp.test.abreaction.run(ops, datadir, individual=True, new=True)
        if options.test_annihilation_2d:
            gpgmp.test.annihilation_2d.run(ops, datadir)
        if options.test_annihilation_2d_drift:
            gpgmp.test.annihilation_2d_drift.run(ops, datadir)
        if options.test_fisher:
            gpgmp.test.fisher_problem.run(ops, datadir)
        if options.test_fisher_individual:
            gpgmp.test.fisher_problem.run(ops, datadir, individual=True)
        if options.test_homogeneous_diffusion:
            gpgmp.test.homogeneous_diffusion.run(ops, datadir)
        if options.test_homogeneous_drift:
            gpgmp.test.homogeneous_drift_diffusion.run(ops, datadir)
        if options.test_homogeneous_drift_individual:
            gpgmp.test.homogeneous_drift_diffusion.run(ops, datadir, individual=True)
        if options.test_homogeneous_drift_field:
            gpgmp.test.homogeneous_drift_diffusion.run_field(ops, datadir)
        if options.test_multiplicative_noise:
            gpgmp.test.multiplicative_noise.run(ops, datadir)
        if options.test_ornstein_uhlenbeck:
            gpgmp.test.ornstein_uhlenbeck.run(ops, datadir)
        if options.test_nonlinear:
            gpgmp.test.nonlinear.run(ops, datadir)
        if options.test_random_drift:
            gpgmp.test.random_drift.run(ops, datadir)

        # majority vote
        if options.majority:
            gpgmp.models.majority.run(ops, datadir)
        if options.majority_flame_avoidance:
            gpgmp.models.majority.runFlameAvoidance(ops, datadir, False, track_positions = options.majority_track_positions)
        if options.majority_flame_velocity:
            gpgmp.models.majority.runFlameVelocity(ops, datadir, False)
        if options.majority_flame_recognition:
            gpgmp.models.majority.runFlameRecognition(ops, datadir, False)
        if options.majority_flame_experiment:
            gpgmp.models.majority.runFlameExperiment(ops, datadir, False)

        if options.majority_flame_diffusive:
            gpgmp.models.majority.runFlame(ops, datadir, True)
            
        if options.majority_gpgmp_nondiffusive:
            gpgmp.models.majority.runGpgmpNondiffusive(ops, datadir)

        if options.majority_gpgmp_avoidance:
            gpgmp.models.majority.runGpgmpNondiffusiveAvoidance(ops, datadir)
            
        if options.majority_toy_model:
            gpgmp.models.majority.runGpgmpToyModel(ops, datadir)
            
    # check if we need to analyze
    if ops['analyzeSimulations']:

        # get results dir
        #resultsdir = ops['datadir']+'results/'
        datadir = ops['datadir']

        # tests
        if options.test_aplusb:
            passed, runtime = gpgmp.test.abreaction.analyze(ops, datadir)
            if passed:
                print("A+B reaction test \t\t\t \x1b[32m \x1b[1m passed.\x1b[0m Runtime: {0:g} seconds".format(runtime))
            else:
                print("A+B reaction test \t\t\t \x1b[31m \x1b[1m failed.\x1b[0m Runtime: {0:g} seconds".format(runtime))

        if options.test_aplusb_individual:
            passed, runtime = gpgmp.test.abreaction.analyze(ops, datadir, individual=True)
            if passed:
                print("A+B reaction test (individual) \t\t\t \x1b[32m \x1b[1m passed.\x1b[0m Runtime: {0:g} seconds".format(runtime))
            else:
                print("A+B reaction test (individual) \t\t\t \x1b[31m \x1b[1m failed.\x1b[0m Runtime: {0:g} seconds".format(runtime))

        if options.test_aplusb_individual_new:
            passed, runtime = gpgmp.test.abreaction.analyze(ops, datadir, individual=True, new=True)
            if passed:
                print("A+B reaction test (individual new) \t\t\t \x1b[32m \x1b[1m passed.\x1b[0m Runtime: {0:g} seconds".format(runtime))
            else:
                print("A+B reaction test (individual new) \t\t\t \x1b[31m \x1b[1m failed.\x1b[0m Runtime: {0:g} seconds".format(runtime))

        if options.test_annihilation_2d:
            runtime = gpgmp.test.annihilation_2d.analyze(ops, datadir)
            print("A+B annihilation test \t\t\t \x1b[33m \x1b[1m (no comparison) \x1b[0m Runtime: {0:g} seconds.".format(runtime))

        if options.test_annihilation_2d_drift:
            runtime = gpgmp.test.annihilation_2d_drift.analyze(ops, datadir)
            print("A+B annihilation test with drift \t\t\t \x1b[33m \x1b[1m (no comparison) \x1b[0m Runtime: {0:g} seconds.".format(runtime))

        if options.test_fisher:
            passed, runtime = gpgmp.test.fisher_problem.analyze(ops, datadir)
            if passed:
                print("Fisher test \t\t\t \x1b[32m \x1b[1m passed.\x1b[0m Runtime: {0:g} seconds".format(runtime))
            else:
                print("Fisher test \t\t\t \x1b[31m \x1b[1m failed.\x1b[0m Runtime: {0:g} seconds".format(runtime))

        if options.test_fisher_individual:
            passed, runtime = gpgmp.test.fisher_problem.analyze(ops, datadir, individual=True)
            if passed:
                print("Fisher test (individual) \t\t\t \x1b[32m \x1b[1m passed.\x1b[0m Runtime: {0:g} seconds".format(runtime))
            else:
                print("Fisher test (individual) \t\t\t \x1b[31m \x1b[1m failed.\x1b[0m Runtime: {0:g} seconds".format(runtime))


        if options.test_homogeneous_diffusion:
            passed, runtime = gpgmp.test.homogeneous_diffusion.analyze(ops, datadir)
            if passed:
                print("Homogeneous diffusion test \t\t\t \x1b[32m \x1b[1m passed.\x1b[0m Runtime: {0:g} seconds".format(runtime))
            else:
                print("Homogeneous diffusion test \t\t\t \x1b[31m \x1b[1m failed.\x1b[0m Runtime: {0:g} seconds".format(runtime))

        if options.test_homogeneous_drift:
            passed, runtime = gpgmp.test.homogeneous_drift_diffusion.analyze(ops, datadir)
            if passed:
                print("Homogeneous drift-diffusion test \t\t\t \x1b[32m \x1b[1m passed.\x1b[0m Runtime: {0:g} seconds".format(runtime))
            else:
                print("Homogeneous drift-diffusion test \t\t\t \x1b[31m \x1b[1m failed.\x1b[0m Runtime: {0:g} seconds".format(runtime))

        if options.test_homogeneous_drift_field:
            passed, runtime = gpgmp.test.homogeneous_drift_diffusion.analyze(ops, datadir, field=True)
            if passed:
                print("Homogeneous drift-diffusion test with field \t\t\t \x1b[32m \x1b[1m passed.\x1b[0m Runtime: {0:g} seconds".format(runtime))
            else:
                print("Homogeneous drift-diffusion test with field \t\t\t \x1b[31m \x1b[1m failed.\x1b[0m Runtime: {0:g} seconds".format(runtime))

        if options.test_homogeneous_drift_individual:
            passed, runtime = gpgmp.test.homogeneous_drift_diffusion.analyze(ops, datadir, individual=True)
            if passed:
                print("Homogeneous drift-diffusion test (individual) \t\t\t \x1b[32m \x1b[1m passed.\x1b[0m Runtime: {0:g} seconds".format(runtime))
            else:
                print("Homogeneous drift-diffusion test (individual) \t\t\t \x1b[31m \x1b[1m failed.\x1b[0m Runtime: {0:g} seconds".format(runtime))

        if options.test_multiplicative_noise:
            passed, runtime = gpgmp.test.multiplicative_noise.analyze(ops, datadir)
            if passed:
                print("Multiplicative noise test \t\t\t \x1b[32m \x1b[1m passed.\x1b[0m Runtime: {0:g} seconds".format(runtime))
            else:
                print("Multiplicative noise test \t\t\t \x1b[31m \x1b[1m failed.\x1b[0m Runtime: {0:g} seconds".format(runtime))

        if options.test_ornstein_uhlenbeck:
            passed, runtime = gpgmp.test.ornstein_uhlenbeck.analyze(ops, datadir)
            if passed:
                print("Ornstein Uhlenbeck test \t\t\t \x1b[32m \x1b[1m passed.\x1b[0m Runtime: {0:g} seconds".format(runtime))
            else:
                print("Ornstein Uhlenbeck test \t\t\t \x1b[31m \x1b[1m failed.\x1b[0m Runtime: {0:g} seconds".format(runtime))

        if options.test_nonlinear:
            passed, runtime = gpgmp.test.nonlinear.analyze(ops, datadir)
            if passed:
                print("Non-linear test \t\t\t \x1b[32m \x1b[1m passed.\x1b[0m Runtime: {0:g} seconds".format(runtime))
            else:
                print("Non-linear test \t\t\t \x1b[31m \x1b[1m failed.\x1b[0m Runtime: {0:g} seconds".format(runtime))

        if options.test_random_drift:
            passed, runtime = gpgmp.test.random_drift.analyze(ops, datadir)
            if passed:
                print("Random drift test \t\t\t \x1b[32m \x1b[1m passed.\x1b[0m Runtime: {0:g} seconds".format(runtime))
            else:
                print("Random drift test \t\t\t \x1b[31m \x1b[1m failed.\x1b[0m Runtime: {0:g} seconds".format(runtime))

        # and plot
        if options.majority:
            gpgmp.models.majority.reduce(datadir)
        if options.majority_flame_recognition or options.majority_flame_velocity or options.majority_flame_experiment or options.majority_flame_avoidance:
            if ops['pbs']:
                # we first need to reduce the single runs
                gpgmp.models.majority.reducePBS(datadir)

            gpgmp.models.majority.analyzeFlame(ops, datadir)

        if options.majority_gpgmp_nondiffusive or options.majority_gpgmp_avoidance:
            gpgmp.models.majority.analyzeFlame(ops, datadir, True)
#            gpgmp.models.majority.analyzeGpgmpNondiffusive(ops, datadir)
        if  options.majority_toy_model:
            gpgmp.models.majority.analyzeToyModel(ops, datadir)
        if options.majority_geoff:
            gpgmp.models.majority.analyzeGeoff(ops, datadir)

# only run when scripted
if __name__ == "__main__":
    # run main
    main()
