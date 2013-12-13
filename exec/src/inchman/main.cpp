/*
 * main.cpp
 *
 *  Created on: Sep 1, 2010
 *      Author: aidan
 */

#include "InchmanDiffusionModel.h"

#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/filesystem.hpp>
#include <boost/foreach.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/program_options.hpp>
#include <boost/log/trivial.hpp>

#include <iostream>

#if defined(_WIN32) || defined(_WIN64)
# include <conio.h>
#else
# include <termios.h>
#endif


using namespace std;
namespace po = boost::program_options;


void waitForAnyKey();


int main(int argc, char **argv)
{
    boost::filesystem::path programPath = boost::filesystem::absolute(boost::filesystem::path(argv[0]).parent_path());
    
#if defined(_WIN32) || defined(_WIN64)
    // Tell the Python runtime where it 'home' is :-)
	boost::filesystem::path pythonPath = programPath / "python";
	if (!boost::filesystem::exists(pythonPath)) {
        BOOST_LOG_TRIVIAL(fatal) << "Error: failed to locate Python";
		exit(1);
	}
	string pythonPathStr = pythonPath.string();
	std::vector<char> pythonPathChars;
	pythonPathChars.assign(pythonPathStr.begin(), pythonPathStr.end());
	pythonPathChars.push_back('\0');
	Py_SetPythonHome(pythonPathChars.data());
#endif // defined(_WIN32) || defined(_WIN64)
    
    // Set default OpenCL source path
#ifdef __APPLE__
    string oclSourcePath = programPath.string() + "/../Frameworks/cl_source/";
#else
    string oclSourcePath = programPath.string() + "/../share/inchman/cl_source/";
#endif
    
     BOOST_LOG_TRIVIAL(info) << "Welcome to Inchman.";
     std::cout << "Invoked using \"";
     for (int i=0; i<argc; i++)
         cout <<" "<<argv[i];
     std::cout <<"\"\n";

    // Declare a group of options that will be 
    // allowed only on command line
    po::options_description generic;//("Generic options");
    generic.add_options()
    ("help", "Print usage information and exit.")
    ;
    
    // Declare a group of options that will be 
    // allowed both on command line and in
    // config file
    po::options_description config;//("Configuration");
    generic.add_options()
    // TODO: get the default output filename from a CDiffusionModel static
    ("output", po::value<string>(), "Output filename (default: \"output.h5\").")
    ("p", po::value< vector<string> >(), "<name>=<value> Override parameter value.")
    ("wait-on-exit", "Wait for user to press any key on exit - helpful when launching from a graphical desktop.")
    ("cl-source", po::value<string>(), "Path to Open-CL source files (default:\".\").")
    ;
    
    // Hidden options, will be allowed both on command line and
    // in config file, but will not be shown to the user.
    po::options_description hidden;//("Hidden options");
    hidden.add_options()
    ("model-file", po::value<string>(), "Model file.")
    ;
    
    po::options_description cmdline_options;
    cmdline_options.add(generic).add(config).add(hidden);
    
    /*
     po::options_description config_file_options;
     config_file_options.add(config).add(hidden);
     */
    po::options_description visible("Usage:\n\n"
                                    "  inchman [options] <model-file>\n\n"
                                    "Options");
    visible.add(generic).add(config);
    
    try
    {
        po::positional_options_description p;
        p.add("model-file", -1);
        
        po::variables_map vm;
        po::store(po::command_line_parser(argc, argv).
                  options(cmdline_options).positional(p).run(), vm);
        po::notify(vm);
        
        if (vm.count("wait-on-exit"))
            atexit(waitForAnyKey);
        
        if (vm.count("help")) {
            cout << visible << endl;
            return 1;
        }
        
        // Get the input filename
        if (vm.count("model-file") != 1) {
            cerr << "Error: an single model file must be provided." << endl << endl << visible << endl;
            return 1;
        }
        string filename = vm["model-file"].as<string>();
        
        // Get the parameter overrides
        map<string, double> paramOverrides;
        if (vm.count("p") > 0) { // check before trying to cast it to vec
            BOOST_FOREACH(const string &var, vm["p"].as< vector<string> >()) {
                // TODO: error handling
                vector<string> parts;
                boost::split(parts, var, boost::is_any_of("="));
                paramOverrides.insert(make_pair(parts[0], boost::lexical_cast<double>(parts[1])));
            }
        }
        
        // Attempt to load and process the model
        gpgmp::InchmanDiffusionModel *imModel = gpgmp::InchmanDiffusionModel::fromFile(filename, paramOverrides);
        if (imModel)
        {
            // Get the output filename override, if it exists
            if (vm.count("output"))
                imModel->setOutputFilename(vm["output"].as<string>());            

            if (vm.count("cl-source"))
                oclSourcePath = vm["cl-source"].as<string>();

            // set source path for ocl files
            BOOST_LOG_TRIVIAL(debug) <<"Using CL source path "<<oclSourcePath<<".";
            imModel->setOclSourcePath(oclSourcePath);
        
            // output configuration
            imModel->logModel();

            // Solve!
            imModel->solve();
            
            delete imModel;
        }
        
        BOOST_LOG_TRIVIAL(info) << "Inchman job finished successfully."; // use cout, so that it can be easily logged
    }
    catch (exception &e) {
        BOOST_LOG_TRIVIAL(error) << "Error: " << e.what() << endl << endl << visible << endl;
        return 1;
    }
    
    return 0;
}


void waitForAnyKey()
{
    cerr << "Press any key to continue..." << endl; // use cerr, so that it user knows the reason for the pause, even when cout is redirected to a file
    
#if defined(_WIN32) || defined(_WIN64)
    _getch(); // conio.h: read a single character from the console without echoing the character. (enter/return does NOT need to be pressed)
#else
    /*
     This code attempts to mimic the getch function within DOS.
     It sets the terminal into non-canonical mode, thus disabling line buffering,
     reads a character from stdin and then restores the old terminal status. 
     
     Note: There's also a getch function in the curses library, but it is not
     equivalent to the DOS getch and may only be used within real curses
     applications (ie: it only works in curses "WINDOW"s).
     */
    termios oldt, newt;
    tcgetattr(STDIN_FILENO, &oldt);
    newt = oldt;
    newt.c_lflag &= ~(ICANON | ECHO);
    tcsetattr(STDIN_FILENO, TCSANOW, &newt);
    getchar();
    tcsetattr( STDIN_FILENO, TCSANOW, &oldt);
#endif
}
