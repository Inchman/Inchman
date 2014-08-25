#include "AnnihilationProblem.h"
#include "SemiInfiniteSlabProblem.h"
#include "DiffusionProblem.h"
#include "LocalizedAnnihilationProblem.h"
#include "FisherProblem.h"
#include "ReactionProblem.h"
#include "DiffusionModel.h"
#include "MultiplicativeNoiseProblem.h"
#include "OrnsteinUhlenbeckProblem.h"
#include "NonlinearProblem.h"
#include "SlitModel.h"
#include "GrayScottModel.h"
#include "OregonatorModel.h"
#include "CalciumModel.h"
#include "RandomDriftProblem.h"
#include "IndividualReactionProblem.h"
#include "MajorityVoteProblem.h"
#include "MajorityFPEProblem.h"

#include "Compartment.h"
#include "Species.h"
#include "Reaction.h"
//#include "Version.h"

#include <cstdio>
#include <iostream>

#include <boost/math/constants/constants.hpp>
#include <boost/python.hpp>
#include <boost/detail/lightweight_test.hpp>
#include <boost/filesystem.hpp>

using namespace std;
using namespace gpgmp;

namespace python = boost::python;

#define INIT_SCRIPT "init.py"
#define EVENTS_SCRIPT "alert.py"

#if defined(_WIN32) || defined(_WIN64)
//int _tmain(int argc, TCHAR** argv)
int main(int argc, char **argv)
{
#else
int main(int argc, char **argv)
{
#endif // defined(_WIN32) || defined(_WIN64)

	boost::filesystem::path programPath = boost::filesystem::absolute(boost::filesystem::path(argv[0]).parent_path());

#if defined(_WIN32) || defined(_WIN64)
    // Tell the Python runtime where it 'home' is :-)
	//boost::filesystem::path pythonPath = programPath / "python";
	//boost::filesystem::path pythonPath = "C:\\Users\\Matthias\\Documents\\Projects\\Inchman\\python";
	boost::filesystem::path pythonPath = "C:\\Users\\Matthias\\Anaconda";
	if (!boost::filesystem::exists(pythonPath)) {
		cerr << "Error: failed to locate Python" << endl;
		exit(1);
	}
	string pythonPathStr = pythonPath.string();
	std::vector<char> pythonPathChars;
	pythonPathChars.assign(pythonPathStr.begin(), pythonPathStr.end());
	pythonPathChars.push_back('\0');
	Py_SetPythonHome(pythonPathChars.data());
#endif // defined(_WIN32) || defined(_WIN64)

    // Set default OpenCL source and Python init scripts paths
    string oclSourcePath = programPath.string();
    string initScriptsPath = programPath.string();

    // print out the command line
    std::cout <<"Welcome. GPGMP" /*<< Version::toString() <<*/ " was invoked with the following command line:\n";
    for (int i=0; i<argc; i++)
        cout <<" "<<argv[i];
    std::cout <<"\n";

    // parameters
    int dx=32, dy=32;
    int numRuns=1;
    int numSteps=500000;
    float time=100.;
    Real length=40.;

    //opterr = 0;
    int c;

    int oclDeviceIndexNo=0;
    int problem=0;
    float p0 = 0.1f;
    float dt = 0.1f;

    OutputFormat outputFormat = OutputHdf5;
    //InhomogeneousDiffusionParameters::AlgorithmType algType = InhomogeneousDiffusionParameters::AT_CA;

    SolverType solverType = gpgmp::stochastic_homogeneous;

	/*
    while ((c = getopt(argc, argv, "x:y:t:l:n:r:c:i:d:p:w:o:f:s:")) != -1)
	{
        switch (c)
        {
        case 'x':
            dx = atoi(optarg);
            break;
        case 'y':
            dy = atoi(optarg);
            break;
        case 't':
            time = atof(optarg);
            break;
        case 'l':
            length = atof(optarg);
            break;
        case 'n':
            numSteps = atoi(optarg);
            break;
        case 'r':
            numRuns = atoi(optarg);
            break;
        case 'c':
            // OpenCL source path
            oclSourcePath = optarg;
            break;
        case 'i':
            // Python init scripts path
            initScriptsPath = optarg;
            break;
        case 'd':
            // OpenCL device to run on
            oclDeviceIndexNo = atoi(optarg);
            break;
        case 'p':
            problem = atoi(optarg);
            break;
        case 'o':
            dt = atof(optarg);
            break;
        case 'w':
            p0 = atof(optarg);
            break;
        case 'f':
            switch(atoi(optarg)) {
            case 0:
                cout <<"Setting output format to ASCII.\n";
                outputFormat = OutputAscii;
                break;
            case 1:
                cout <<"Setting output format to HDF5.\n";
                outputFormat = OutputHdf5;
                break;
            }
            break;
        case 's':
            switch(atoi(optarg)) {
            case 0:
                cout <<"Setting solver type to stochastic inhomogeneous (CA).\n";
                solverType = gpgmp::stochastic_inhomogeneous;
                break;
            case 1:
                cout <<"Setting solver type to stochastic inhomogeneous (FPE).\n";
                solverType = gpgmp::stochastic_inhomogeneous_fpe;
                break;
            case 2:
                cout <<"Setting solver type to stochastic homogeneous.\n";
                solverType = gpgmp::stochastic_homogeneous;
                break;
            case 3:
                cout <<"Setting solver type to deterministic homogeneous (RK4).\n";
                solverType = gpgmp::deterministic_homogeneous_RK4;
                break;
            case 4:
                cout <<"Setting solver type to deterministic homogeneous (alpha-QSS).\n";
                solverType = gpgmp::deterministic_homogeneous_aqss;
                break;
            case 5:
                cout <<"Setting solver type to deterministic inhomogeneous RK4.\n";
                solverType = gpgmp::deterministic_inhomogeneous_RK4;
                break;
            case 6:
                cout <<"Setting solver type to deterministic inhomogeneous alpha-QSS.\n";
                solverType = gpgmp::deterministic_inhomogeneous_aqss;
                break;
            }
            break;

        case '?':
            if (optopt == 'x' || optopt =='y' || optopt == 't'
                    || optopt =='l' || optopt == 'n' || optopt=='r'
                    || optopt=='c'  || optopt=='i' || optopt=='d'
                    || optopt=='p' || optopt == 'w' || optopt=='s')
                fprintf (stderr, "Option -%c requires an argument.\n", optopt);
            else if (isprint (optopt))
                fprintf (stderr, "Unknown option `-%c'.\n", optopt);
            else
                fprintf (stderr,
                         "Unknown option character `\\x%x'.\n",
                         optopt);
            return 1;
        default:
            abort ();
        }
    }
	*/

    int index=0;


    // This is the model we will use.
    DiffusionModel *model = nullptr;
	problem = 30;

    if (problem == 0) {
        // -------------------------------- A plus B reaction model
        float diffusionConstant;
        float k1, k2, k3, k4;
        if (index == argc-5) {
            diffusionConstant=atof(argv[index++]);
            k1=atof(argv[index++]);
            k2=atof(argv[index++]);
            k3=atof(argv[index++]);
            k4=atof(argv[index++]);
        } else if (index==argc) {
            cout <<"No parameters given, using default.\n";
            length=100.;
            dx=128;
            dy=128;
            p0=0.1;
            time=500.;
            k1=1e-3;
            k2=1e-2;
            k3=1.2;
            k4=1.;
            diffusionConstant=0.;
        } else {
            cout <<"Mismatching parameters. Usage:\n";
            cout <<"gpgmp -p 0 [diffusionConstant k1 k2 k3 k4].\n";
            exit(1);
        }
        cout <<"Setting up A plus B reaction problem with diffusion Constant :"
            <<diffusionConstant<<"\n";
        cout <<"k1 = "<<k1<<", k2 = "<<k2<<", k3 = "<<k3<<", k4 = "<<k4<<"\n";

        // create reaction problem
        model = new ReactionProblem(length, dx, dy,
                                    diffusionConstant,
                                    k1, k2, k3, k4);

    } // Simple reaction problem

    else if (problem == 1) {
        // FISHER PROBLEM
        float baseConcentration = (index>=argc ? 1.0e-6 : atof(argv[index++])); // provides sane default
        cout <<"Setting up Fisher Problem with concentration "<<baseConcentration<<".\n";

        model = new FisherProblem(length, dx, dy, baseConcentration, initScriptsPath);
    } // Fisher problem

    else if (problem == 2) {
            // SIMPLE DIFFUSION
            float diffusionConstant;
            int numMolecules;

            if (index == argc-2) {
                diffusionConstant=atof(argv[index++]);
                numMolecules = atoi(argv[index++]);
            } else if (index==argc) {
                cout <<"No parameters given, using default.\n";
                length=100.;
                dx=128;
                dy=128;
                p0=0.1;
                dt = 5.;
                time=20.;
                diffusionConstant=1.;
                numMolecules = 1e5;
            } else {
                cout <<"Mismatching parameters. Usage:\n";
                cout <<"gpgmp -p 2 [diffusionConstant numMolecules].\n";
                exit(1);
            }

            cout <<"Setting up simple diffusion problem with diffusion Constant :"<<diffusionConstant
                    <<", nMol:"<<numMolecules<<" \n";

            model = new DiffusionProblem(length, dx, dy, diffusionConstant, numMolecules);
    } // simple diffusion problem

    else if (problem ==3) {
        // LOCALIZED A+B ANNIHILATION
        cout <<"Setting up localized A plus B reaction problem.\n";
        float diffusionConstant;
        float k1, k2, k3, k4;
        float xminA, xmaxA,xminB, xmaxB;

        if (index == argc-9) {
            diffusionConstant=atof(argv[index++]);
            k1=atof(argv[index++]);
            k2=atof(argv[index++]);
            k3=atof(argv[index++]);
            k4=atof(argv[index++]);
            xminA=atof(argv[index++]);
            xmaxA=atof(argv[index++]);
            xminB=atof(argv[index++]);
            xmaxB=atof(argv[index++]);
        } else if (index==argc) {
            cout <<"No parameters given, using default.\n";
            length=1000.;
            dx=32;
            dy=32;
            p0=0.1;
            time=1800.;
            // annihilation rates are in M^-1 s^-1
            k1=9.4096e+09;
            k2=9.4096e+10;
            // creation rates are in M s^-1
            k3=1.27529e-13;
            k4=1.06274e-13;
            diffusionConstant=100.;
            xminA = 0.;
            xmaxA = 900.;
            xminB = 400.;
            xmaxB = 1000.;
        } else {
            cout <<"Mismatching parameters. Usage:\n";
            cout <<"gpgmp -p 3 [diffusionConstant k1 k2 k3 k4 xminA xmaxA xminB xmaxB].\n";
            exit(1);
        }
        cout <<"Parameters : diffusion Constant :"<<diffusionConstant<<"\n";
        cout <<"k1 = "<<k1<<", k2 = "<<k2<<", k3 = "<<k3<<", k4 = "<<k4<<"\n";
        cout <<"xminA = "<<xminA<<", xmaxA = "<<xmaxA<<", xminB = "<<xminB<<", xmaxB = "<<xmaxB<<"\n";

        model = new LocalizedAnnihilationProblem(length, dx, dy,
                                                 diffusionConstant, k1, k2, k3, k4,
                                                 xminA, xmaxA, xminB, xmaxB);
    } // localised annihilation problem

    else if (problem == 4)
    {
        // SEMI-INFINITE SLAB problem
        cout <<"Setting up semi-infinite slab problem.\n";
        float diffusionConstant;
        float sourceNumber;

        if (index == argc-2) {
            diffusionConstant=atof(argv[index++]);
            sourceNumber=atof(argv[index++]);
        } else if (index==argc) {
            cout <<"No parameters given, using default.\n";
            length=40.;
            dx=64;
            dy=64;
            p0=0.1;
            time=1000;
            diffusionConstant=1.;
            sourceNumber = 10.;
            dt=500.;
        } else {
            cout <<"Mismatching parameters. Usage:\n";
            cout <<"gpgmp -p 4 [diffusionConstant sourceNumber].\n";
            exit(1);
        }
        cout <<"Setting up semi-infinite slab problem.\n";
        cout <<"Source number is:"<<sourceNumber<<", diffusion Constant:"
                <<diffusionConstant<<"\n";

        model = new SemiInfiniteSlabProblem(length, dx, dy, sourceNumber, diffusionConstant);
    } // semi-infinite slab problem

    else if (problem==6)
    {
        // 2D A+B annihilation problem without drift
        std::cout <<"Setting up 2D annihilation problem without drift\n";

        // default parameters
        Real diffusionConstant = 1.;
        Real rate;
        int numMolecules;

        // check if more parameters are given - if not take defaults
        if (index == argc-3) {
            diffusionConstant=atof(argv[index++]);
            rate=atof(argv[index++]);
            numMolecules = atoi(argv[index++]);
        } else if (index == argc) {
            // use defaults
            cout <<"No parameters given ... using defaults.\n";
            length=12.8;
            dx=64;
            dy=64;
            p0=0.1;
            time=4.;
            rate = 1e8;
            numMolecules = 100000.;
        } else {
            cout <<"Mismatching parameters. Usage:\n";
            cout <<"gpgmp -p 24 diffusionConstant drift rate numMolecules.\n";
            exit(1);
        }
        std::cout <<"Annihilation rate is:"<<rate<<", diffusion Constant:"
                <<diffusionConstant<<", number of Molecules:"
		  <<numMolecules<<".\n" <<std::flush;

        // set up the model
	std::cout <<"About to construct ..\n"<<std::flush;
        model = new AnnihilationProblem(length, dx, dy,
                                        numMolecules, diffusionConstant, rate);
	std::cout <<"About to construct finished ..\n"<<std::flush;
    } // Annihilation problem (without drift)

    else if (problem == 7)
    {
        // SIMPLE INHOMOGENEOUS DIFFUSION
        float diffusionConstantX, diffusionConstantY;
        float rx, ry;
        int numMolecules;

        if (index == argc-5) {
            diffusionConstantX=atof(argv[index++]);
            diffusionConstantY=atof(argv[index++]);
            rx = atof(argv[index++]);
            ry = atof(argv[index++]);
            numMolecules = atoi(argv[index++]);
        } else {
            cout <<"No parameters given ... using defaults.\n";
            length = 32.;
            dx = 64;
            dy = 64;
            time = 5.;
            diffusionConstantX = 0.3;
            diffusionConstantY = 0.4;
            rx = 0.8;
            ry = 0.4;
            numMolecules = 100000;
        }


        cout <<"Setting up simple inhomogeneous diffusion problem with diffusion Constant :"
                <<diffusionConstantX <<", and "<<diffusionConstantY
                <<", nMol:"<<numMolecules<<" and drift field "<<
                "("<<rx<<","<<ry<<").\n";

        model = new DiffusionProblem(length, dx, dy, diffusionConstantX,
                                     numMolecules, diffusionConstantY,
                                     rx, ry);

    } // homogeneous drift-diffusion problem
    else if (problem == 8)
    {
        // Slit Problem
        model = new SlitModel(length, dx, dy, initScriptsPath);
    } // slit problem
    else if (problem == 13)
    {

        // multiplicative noise
        Real cx, cy, mux, muy;
        int numMolecules;

        if (index == argc-5) {
            cx=atof(argv[index++]);
            cy=atof(argv[index++]);
            mux=atof(argv[index++]);
            muy=atof(argv[index++]);
            numMolecules = atoi(argv[index++]);
        } else {
            cout <<"No parameters given ... using defaults.\n";
            length = 32.;
            dx = 64;
            dy = 64;
            time = 6e-4;
	    dt = 1e-4;
            cx = 5.;
            cy = 4.;
            mux = 0.1;
            muy = 0.2;
            numMolecules = 100000;
	    numRuns = 10;
	    solverType = gpgmp::stochastic_inhomogeneous;
        }

        cout <<"Setting up multiplicative noise problem with cx="
                <<cx <<", cy="<< cy <<", mux = "<<mux<<", muy="<<muy<<" and "
                <<"nMol="<<numMolecules<<".\n";

        // simple model
        model = new MultiplicativeNoiseProblem(length, dx, dy,
                                               numMolecules,
                                               cx, cy, mux, muy);
    } // multiplicative-noise

    else if (problem == 14)
    {
        // Ornstein-Uhlenbeck process
        Real gxx, gxy, gyx, gyy;
        Real diffx, diffy;
        int numMolecules;

        if (index == argc-7) {
            diffx=atof(argv[index++]);
            diffy=atof(argv[index++]);
            gxx=atof(argv[index++]);
            gxy=atof(argv[index++]);
            gyx=atof(argv[index++]);
            gyy=atof(argv[index++]);
            numMolecules = atoi(argv[index++]);
        } else {
            cout <<"No parameters given ... using defaults.\n";
            length = 32.;
            dx = 64;
            dy = 64;
            time = 0.001;
            diffx=1.; diffy=1.;
            gxx = 0.5; gxy=0.6;
            gyx = 0.7; gyy=0.8;
            numMolecules = 10000;
        }

        cout <<"Setting OU problem with D=("<<diffx<<","<<diffy<<"), "
                << "G = ("<<gxx<<","<<gxy<<","<<gyx<<","<<gyy<<"),"
                <<"nMol="<<numMolecules<<".\n";

        model = new OrnsteinUhlenbeckProblem(length, dx, dy,
                                             diffx, diffy,
                                             gxx, gxy, gyx, gyy,
                                             numMolecules);
    } // Ornstein-Uhlenbeck

    else if (problem == 15)
    {
        // Non-linear
        Real theta, omega;
        Real diffx, diffy;
        int numMolecules;

        if (index == argc-5) {
            diffx=atof(argv[index++]);
            diffy=atof(argv[index++]);
            theta=atof(argv[index++]);
            omega=atof(argv[index++]);
            numMolecules = atoi(argv[index++]);
        } else {
            cout <<"No parameters given ... using defaults.\n";
            length = 32.;
            dx = 64;
            dy = 64;
            time = 0.001;
            diffx=1.; diffy=1.;
            theta = 1.;
            omega = 1.;
            numMolecules = 10000;
        }

        cout <<"Setting nonlinear problem with D=("<<diffx<<","<<diffy<<"), "
                << "theta = "<<theta<<", omega="<<omega<<", "
                <<"nMol="<<numMolecules<<".\n";


        model = new NonlinearProblem(length, dx, dy,
                                     diffx, diffy,
                                     theta, omega,
                                     numMolecules);
    } // non-linear problem
   else if (problem == 19) {

            Real scale = 3e-4f;

            if (index == argc-1) {
                scale=atof(argv[index++]);
            } else {
                cout <<"No parameters given ... using defaults.\n";
                length = 0.27*2.f;
                dx = 128*2;
                dy = 128*2;
            }

            std::cout <<"Running dimensional three-state oregonator problem..";

            model = new OregonatorModel(length, dx, dy, scale, initScriptsPath);

        } // Oregonator Model
    else if (problem == 20) {
            // Gray Scott model
            float diffusionConstantU, diffusionConstantV;
            float k1;
            float k2;
            float kf;
            float u0;
            float omega;

            Real F, k;

            if (index == argc-5) {
                diffusionConstantU=atof(argv[index++]);
                diffusionConstantV=atof(argv[index++]);

                k=atof(argv[index++]);
                F=atof(argv[index++]);
                omega=atof(argv[index++]);
            } else if (index==argc) {
                cout <<"No parameters given, using default.\n";

                // this one is Muratov
                length = 44.7f; // l*100
                dx=256;
                dy=256;

                omega = 250.f;
                p0=0.1f;

                // this is Wang et al. for spirals
                k = 0.0225f;
                F = 0.0025f;
                diffusionConstantU=0;
                diffusionConstantV=0.005;
            } else {
                cout <<"Mismatching parameters. Usage:\n";
                cout <<"gpgmp -p 20 [diffusionConstantu diffusionConstantv k F omega].\n";
                exit(1);
            }

            model = new GrayScottModel(length, dx, dy,
                                       diffusionConstantU, diffusionConstantV,
                                       k, F, omega,
                                       initScriptsPath);
        } // Gray-Scott model

    else if (problem == 22) {
        float scale;
        float nmax;
        float kp, betac, diffc, diffp;

        std::cout <<"Setting up Ca2+ problem (after Atri et al., 1993).\n"<<std::flush;

        if (index == argc-6) {
            scale = atof(argv[index++]);
            nmax = atof(argv[index++]);
            kp = atof(argv[index++]);
            betac = atof(argv[index++]);
            diffc=atof(argv[index++]);
            diffp=atof(argv[index++]);
        } else if (index == argc) {
            std::cout <<"No parameters given - using default.\n";
            length=250.; // from Pearson (1993)
            dx=128;
            dy=128;
            p0=0.1;
            diffc =20.;
            diffp = 0.;

            scale = 1e-6;
            nmax = 1000;
            betac = 0.02;
            kp = 0.;
        } else {
            std::cout <<"Mismatching parameters.\n";
            std::cout <<"USAGE: gpgmp -p 22 [scale nmax kp beta diffc diffp]\n";
        }

        std::cout <<"Setting up Ca2+ problem (Atri) with "
                <<"scale = "<<scale<<", nmax="<<nmax<<", kp="<<kp<<", betac="
                <<betac<<", diffc = "<<diffc<<", diffp="<<diffp<<"\n"<<std::flush;


    // create simple Ca2+ model
    model = new CalciumModel(length, dx, dy,
                             diffc, diffp,
                             nmax, kp, betac, scale,
                             initScriptsPath);
    } // calcium model
    else if (problem == 24) {
        // 2D A+B annihilation problem with drift
        std::cout <<"Setting up 2D annihilation problem with drift\n";

        // default parameters
        Real diffusionConstant = 1., drift=0.5;
        Real rate;
        int numMolecules;

        // check if more parameters are given - if not take defaults
        if (index == argc-4) {
            diffusionConstant=atof(argv[index++]);
            drift=atof(argv[index++]);
            rate=atof(argv[index++]);
            numMolecules = atoi(argv[index++]);
        } else if (index == argc) {
            // use defaults
            cout <<"No parameters given ... using defaults.\n";
            length=12.8;
            dx=64;
            dy=64;
            p0=0.1;
            time=4.;
            rate = 1e8;
            numMolecules = 100000;
        } else {
            cout <<"Mismatching parameters. Usage:\n";
            cout <<"gpgmp -p 24 diffusionConstant drift rate numMolecules.\n";
            exit(1);
        }
        std::cout <<"Annihilation rate is:"<<rate<<", diffusion Constant:"
                <<diffusionConstant<<", number of Molecules:"
                <<numMolecules<<"and drift rate is "<<drift<<".\n";

        model = new AnnihilationProblem(length, dx, dy,
                                        numMolecules, diffusionConstant, rate, drift);
    } // A+B annihilation with drift
    else if (problem == 25)
    {
        // RANDOM DRIFT MODEL
        float diffusionConstantX, diffusionConstantY;
        float mux, sigmax, muy, sigmay;

        int numMolecules, discreteK;

        if (index == argc-8) {
            diffusionConstantX=atof(argv[index++]);
            diffusionConstantY=atof(argv[index++]);
            mux = atof(argv[index++]);
            sigmax = atof(argv[index++]);
            muy = atof(argv[index++]);
            sigmay = atof(argv[index++]);
            numMolecules = atoi(argv[index++]);
            discreteK = atoi(argv[index++]);
        } else {
            cout <<"No parameters given ... using defaults.\n";
            length = 40.;
            dx = 512;
            dy = 512;
            time = 2.;
            diffusionConstantX = 1.;
            diffusionConstantY = 0.5;
            mux = 0.; sigmax=2.;
            muy = 0.; sigmay=1.;
            numMolecules = 16384;
            discreteK = 5;
        }


        std::cout <<"Setting up random drift problem with parameters"
                    "diffX = " <<diffusionConstantX << ", diffY = " << diffusionConstantY
                 << "sigmax = " <<sigmax << ", sigmay = " << sigmay
                 << "numMols = " <<numMolecules <<", discreteK = " << discreteK <<".\n";

        model = new RandomDriftProblem(length, dx, dy,
                                       diffusionConstantX, diffusionConstantY,
                                       mux, sigmax, muy, sigmay,
                                       numMolecules, discreteK);

        model->setHasContinuousParameterSpecies(true);

    } // random drift model
    else if (problem == 26)
    {
        // Individual reaction model
        float diffusionConstantX, diffusionConstantY, k0;
        int numMolecules;
        float k1, k2, k3, k4;

        if (index == argc-6) {
            diffusionConstantX=atof(argv[index++]);
            diffusionConstantY=atof(argv[index++]);
            k1=atof(argv[index++]);
            k2=atof(argv[index++]);
            k3=atof(argv[index++]);
            k4=atof(argv[index++]);
        } else {
            cout <<"No parameters given ... using defaults.\n";
            length=100.;
            dx=128;
            dy=128;
            p0=0.1;
            time=500.;
            k1=1e-3;
            k2=1e-2;
            k3=1.2;
            k4=1.;
            diffusionConstantX = 1.;
            diffusionConstantY = 1.;
        }


        std::cout <<"Setting up individual reaction problem with parameters"
                    "diffX = " <<diffusionConstantX << ", diffY = " << diffusionConstantY
                  <<"k1 = "<<k1<<", k2 = "<<k2<<", k3 = "<<k3<<", k4 = "<<k4<<"\n";

        model = new IndividualReactionProblem(length, dx, dy,
                                              diffusionConstantX, diffusionConstantY,
                                              k1, k2, k3, k4);

        model->setHasContinuousParameterSpecies(true);

    } // individual reaction problem

    else if (problem == 27)
    {
        // Majority vote model
        Real diffx, diffy, mux, muy, reactionProb;
        uint numMolecules, nReds, numEncounters;

        if (index == argc-8) {
            diffx = atof(argv[index++]);
            diffy = atof(argv[index++]);
            mux = atof(argv[index++]);
            muy = atof(argv[index++]);
            numMolecules  = atoi(argv[index++]);
            nReds = atoi(argv[index++]);
            numEncounters = atoi(argv[index++]);
            reactionProb = atof(argv[index++]);
        } else {
            cout <<"No parameters given ... using defaults.\n";
            /*
             these are for the honours theses:
            diffx = 0.25; diffy = 0.25; // (4 D dt)^1/2 = speed dt =>  D = speed^2 dt^2/4 dt = speed^2 dt/4 = 1/4 (speed=dt=1)
            mux = 0.; muy = 0.;
            numMolecules = 10;
            nReds = 6;
            numEncounters = 5   ;
            reactionProb = 2.41085; // this roughly corresponds to d=4cm ..
            length = 80.; // 80x80 square domain
            dx = 32; // todo: area of grid square should (roughly) correspond to 4./3. pi r^2 with r=comm. dist
            dy = 32;
            */
            // these are for the paper:
            diffx = 3e-3; diffy = 3e-3; // (4 D dt)^1/2 = speed dt =>  D = speed^2 dt^2/4 dt = speed^2 dt/4 = 1/4 (speed=dt=1)
            mux = 0.; muy = 0.;
            numMolecules = 150;
            nReds = 50.; // this is not actually being respected!
            numEncounters = 5   ;
            reactionProb = 0.485675; // this roughly corresponds to d=0.01
            length = 1.; // 1x1 square domain
            dx = 64; // todo: area of grid square should (roughly) correspond to 4./3. pi r^2 with r=comm. dist
            dy = 64;
        }


        std::cout <<"Setting up majority vote problem with parameters"
                 << "D=("<<diffx<<","<<diffy<<"), mu=("<<mux<<","<<muy<<"), "
                 << "n="<<numMolecules<<" (reds="<<nReds<<"), nencounters="<<numEncounters
                 << ", k="<<reactionProb<<"\n";

        model = new MajorityVoteProblem(length, dx, dy,
                                        diffx, diffy, mux, muy,
                                        numMolecules, nReds, numEncounters,
                                        reactionProb, initScriptsPath);
        model->setHasContinuousParameterSpecies(false);

    } // majority vote model
    else if (problem == 28)
    {
        // sets up a majority vote model without diffusion
        // Majority vote model
        Real reactionProb, recognitionError;
        uint numMolecules, nReds, numEncounters;

        if (index == argc-5) {
            numMolecules  = atoi(argv[index++]);
            nReds = atoi(argv[index++]);
            numEncounters = atoi(argv[index++]);
            reactionProb = atof(argv[index++]);
            recognitionError = atof(argv[index++]);
        } else {
            cout <<"No parameters given ... using defaults.\n";
            // these are for the paper:
            numMolecules = 150;
            nReds = 75.;
            numEncounters = 5   ;
            reactionProb = 0.485675; // this roughly corresponds to d=0.01
            recognitionError = 0.;
            length = 1.; // 1x1 square domain
            dx = 64; // todo: area of grid square should (roughly) correspond to 4./3. pi r^2 with r=comm. dist
            dy = 64;
        }


        std::cout <<"Setting up majority vote problem (no diffusion) with parameters"
                 << "n="<<numMolecules<<" (reds="<<nReds<<"), nencounters="<<numEncounters
                 << ", k="<<reactionProb<<", recError="<<recognitionError<<"\n";

        model = new MajorityVoteProblem(length, dx, dy,
                                        numMolecules, nReds, numEncounters,
                                        reactionProb, recognitionError, initScriptsPath);
        model->setHasContinuousParameterSpecies(false);
    } // majority vote model

    else if (problem == 29)
    {
        // sets up the majority vote toy model without diffusion
        // Majority vote model
        Real reactionProb;
        uint numMolecules, nReds;

        if (index == argc-3) {
            numMolecules  = atoi(argv[index++]);
            nReds = atoi(argv[index++]);
            reactionProb = atof(argv[index++]);
        } else {
            cout <<"No parameters given ... using defaults.\n";
            // these are for the paper:
            numMolecules = 150;
            nReds = 75.;
            reactionProb = 0.0002;
            length = 1.; // 1x1 square domain
            dx = 64;
            dy = 64;
        }


        std::cout <<"Setting up majority vote problem (no diffusion) with parameters"
                 << "n="<<numMolecules<<" (reds="<<nReds<<"), "
                 << "k="<<reactionProb<<"\n";

        model = new MajorityVoteProblem(length, dx, dy,
                                        numMolecules, nReds,
                                        reactionProb, initScriptsPath);
        model->setHasContinuousParameterSpecies(false);
    } // majority vote model
	else if (problem == 30)
	{
		// Solving the FPE for the majority papter
		length = 1.;
		dx = 256;
		dy = 256;
		time = 2000.;
		unsigned int numMolecules = 100000;
		solverType = gpgmp::stochastic_inhomogeneous;
		numRuns = 1;
		dt = 50;
		oclSourcePath = "C:\\Users\\Matthias\\Documents\\Projects\\Inchman\\exec\\src\\gpgmp\\src";
		initScriptsPath + "C:\\Users\\Matthias\\Documents\\Projects\\Inchman\\exec\\src\\tests";
		std::cout << "Setting up majority vote FPE problem... ";
		std::cout << std::flush;

		model = new MajorityFPEProblem(length, dx, dy,
			numMolecules);
		std::cout << "Done. \n";
		std::cout << std::flush;

	}
    // run the models with the user-provided parameters
    if (model != 0) {
	std::cout <<"About to solve..\n"<<std::flush;
        model->setOutputFormat(outputFormat);
        model->setP0(p0);
        model->setNumRuns(numRuns);
        model->setMaxSteps(numSteps);
        model->setMaxTime(time);
        model->setOutputInterval(dt);
        model->setGridDims(dx, dy);
        model->setPhysicalLength(length);
        model->setSolver(solverType);
        model->setOclSourcePath(oclSourcePath);
        model->setOclDeviceIndex(oclDeviceIndexNo);
	model->logModel();
	std::cout << std::flush;
        model->solve();
	std::cout <<"Finished..\n"<<std::flush;
    } else {
        std::cerr <<"ERROR: No test problem provided!\n";
    }
} // main
