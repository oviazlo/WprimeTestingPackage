#include "xAODRootAccess/Init.h"
#include "SampleHandler/SampleHandler.h"
#include "SampleHandler/ToolsDiscovery.h"
#include "EventLoop/Job.h"
#include "EventLoop/ProofDriver.h"
#include "SampleHandler/DiskListLocal.h"
#include <TSystem.h>

#include "MyAnalysis/MyxAODAnalysis.h"
#include <EventLoopAlgs/NTupleSvc.h>
#include <EventLoop/OutputStream.h>

/// boost
#include "boost/program_options.hpp"
#include <boost/algorithm/string.hpp>

namespace 
{ 
  const size_t ERROR_IN_COMMAND_LINE = 1; 
  const size_t SUCCESS = 0; 
  const size_t ERROR_UNHANDLED_EXCEPTION = 2; 
  const size_t HELP_CALLED = 3;
  
  namespace po = boost::program_options;
 
} /// namespace 

int parseOptionsWithBoost(po::variables_map &vm, int argc, char* argv[]);

int main( int argc, char* argv[] ) {

	/// get global input arguments:
	po::variables_map vm; 
	const size_t returnedMessage = parseOptionsWithBoost(vm,argc,argv);
	if (returnedMessage!=SUCCESS) std::exit(returnedMessage);

	/// Take the submit directory from the input if provided:
	std::string submitDir = "submitDir";
	if ( vm.count("folder") ){
		submitDir = vm["folder"].as<std::string>();  
	}

	int nEvents = -1;
	if ( vm.count("nEvents") ){
		nEvents = vm["nEvents"].as<unsigned int>();  
	}
	
	// Set up the job for xAOD access:
	xAOD::Init().ignore();

	// Construct the samples to run on:
	SH::SampleHandler sh;

	const char* inputFilePath = gSystem->ExpandPathName ("/afs/cern.ch/work/o/oviazlo/Wprime/datasets");
	//const char* inputFilePath = gSystem->ExpandPathName ("$ALRB_TutorialData/r5591"); // using $ALRB_TutorialData previously defined
	SH::DiskListLocal list (inputFilePath);
	SH::scanDir (sh, list, "DAOD_EXOT9*root*");
	//SH::scanDir (sh, list, "DAOD_EXOT9.04655450._000002.pool.root.1");
	//SH::scanDir (sh, list, "AOD.01494882._113691.pool.root.1"); // specifying one particular file for testing

	/// Set the name of the input TTree. It's always "CollectionTree"
	/// for xAOD files.
	sh.setMetaString( "nc_tree", "CollectionTree" );

	/// Print what we found:
	sh.print();

	/// Create an EventLoop job:
	EL::Job job;
	job.sampleHandler( sh );

	/// Specify that we only want to run on 1k events
	if (nEvents!=-1)
		job.options()->setDouble(EL::Job::optMaxEvents, nEvents);

	/// define an output and an ntuple associated to that output 
	//    EL::OutputStream output  ("myOutput");
	//    job.outputAdd (output);
	//    EL::NTupleSvc *ntuple = new EL::NTupleSvc ("myOutput");
	//    job.algsAdd (ntuple);

	/// Add our analysis to the job:
	MyxAODAnalysis* alg = new MyxAODAnalysis();
	job.algsAdd( alg );
	/// all configuration to the algorithm should be passed after attachment to the job (practically after this comment)

	//    alg->outputName = "myOutput"; // give the name of the output to our algorithm
	//alg->m_useHistObjectDumper = false;
	alg->m_useMuonCalibrationAndSmearingTool = false;

	if ( vm.count("proof") ){/// Run the job using the local/direct driver:
		EL::ProofDriver driver;
		if ( vm.count("nWorkers") ){
			driver.numWorkers = vm["nWorkers"].as<unsigned int>();  
		}
		driver.submit( job, submitDir );
	}
	else{/// Run the job using the local/direct driver:
		EL::DirectDriver driver;
		driver.submit( job, submitDir );
	}
	
	return 0;
}



/// parse input arguments
int parseOptionsWithBoost(po::variables_map &vm, int argc, char* argv[]){
  
  try 
  { 
    /** Define and parse the program options 
     */ 
    po::options_description desc("Options"); 
    desc.add_options()
      ("help,h", "Print help messages") 
      ("folder,f", po::value<string>(), "output working-folder name")
      ("nWorkers,w", po::value<unsigned int>(), "number of workers")
      ("proof,p", "enable PROOF-Lite mode") 
      ("nEvents,n", po::value<unsigned int>(), "number of events to proceed")
      ;
    try 
    { 
      po::store(po::parse_command_line(argc, argv, desc),  
                vm); /// can throw 
 
      /** --help option 
       */ 
      if ( vm.count("help")  ) 
      { 
        std::cout << "Basic Command Line Parameter App" << std::endl 
                  << desc << std::endl; 
        return HELP_CALLED; 
      } 
      
      po::notify(vm); /// throws on error, so do after help in case 
                      /// there are any problems 
    } 
    catch(po::error& e) 
    { 
      std::cerr << "ERROR: " << e.what() << std::endl << std::endl; 
      std::cerr << desc << std::endl; 
      return ERROR_IN_COMMAND_LINE; 
    } 
 
    /// application code here // 
    return SUCCESS;
 
  } 
  catch(std::exception& e) 
  { 
    std::cerr << "Unhandled Exception reached the top of main: " 
              << e.what() << ", application will now exit" << std::endl; 
    return ERROR_UNHANDLED_EXCEPTION; 
 
  } 
  
}

