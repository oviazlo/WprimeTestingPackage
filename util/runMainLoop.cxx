#include "xAODRootAccess/Init.h"
#include "SampleHandler/SampleHandler.h"
#include "SampleHandler/ToolsDiscovery.h"
#include "EventLoop/Job.h"
#include "EventLoop/ProofDriver.h"
#include "EventLoop/DirectDriver.h"
#include "EventLoop/LSFDriver.h"
#include "SampleHandler/DiskListLocal.h"
#include <TSystem.h>

#include "MyAnalysis/MyxAODAnalysis.h"
#include <EventLoopAlgs/NTupleSvc.h>
#include <EventLoop/OutputStream.h>

/// boost
#include "boost/program_options.hpp"
#include <boost/algorithm/string.hpp>
#include "boost/filesystem.hpp"

#include <stdlib.h>     /* getenv */

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
  
  /// Set up the job for xAOD access:
  xAOD::Init().ignore();

  /// Construct the samples to run on:
  SH::SampleHandler sh;

  const char* inputFilePath;
  
  /// get hostname of the running machine
  /// check if it is cern or lunarc machines
  char* hostNameChArr;
  hostNameChArr = getenv("HOSTNAME");
  string hostName(hostNameChArr);
  
  /// define which input-path to use
  /// look for HOSTNAME env. variable
  std::size_t found = hostName.find("cern");
  if (found!=std::string::npos){
    /// cern machines
    inputFilePath = gSystem->ExpandPathName
    ("/afs/cern.ch/work/o/oviazlo/Wprime/AnalysisFramework/rel20/data");
  }
  else{
    std::size_t found = hostName.find("alarik");
    if (found!=std::string::npos){
      /// alarik cluster machines
      inputFilePath = gSystem->ExpandPathName
      ("/lunarc/nobackup/users/oviazlo/xAOD/p2425");
    }
    else{
      /// iridium cluster
      inputFilePath = gSystem->ExpandPathName
      ("/nfs/shared/pp/oviazlo/xAOD/cutFlow"); 
    }
  }

  SH::DiskListLocal list (inputFilePath);
  SH::scanDir (sh, list, "DAOD_EXOT9.*root*");
  
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
  
  /// all configuration to the algorithm should be passed after attachment to 
  /// the job (practically after this comment)
  ///
  /// syntax of passing configurations:
  ///   alg->m_myMegaFlag = true;
  ///
  /// [List of possble flags to use]
  /// bool m_useHistObjectDumper; - do not make default plots by HistObjectDumper
  /// bool m_useBitsetCutflow; - do not save cutflow
  /// bool m_useCalibrationAndSmearingTool; - don't do muon calibration and 
  ///                                         smearing TODO not in the code now?
  /// bool m_runElectronChannel; - run electron cycle instead of muon one
  /// bool m_doWprimeTruthMatching; - do truth matching to identify Wprime decay 
  ///                                 to muon/electron channel
  /// bool m_doNotApplyTriggerCuts; - do not apply triggers in MC
  /// string outputName; - name of output tree TODO not implemented yet
  
  alg->m_useCalibrationAndSmearingTool = true;
  if (vm.count("noSmearing"))
    alg->m_useCalibrationAndSmearingTool = false;
  
  alg->m_doWprimeTruthMatching = false;

  if (vm.count("electronChannel"))
    alg->m_runElectronChannel = true;
  
  /// overwrite output folder
  if (vm.count("overwrite")){
    const boost::filesystem::path path( submitDir); 
    if ( boost::filesystem::exists( path ) ) 
    {
      boost::filesystem::remove_all( path );
    }
  }
  
  if ( vm.count("proof") ){/// Run the job using the local/direct driver:
    EL::ProofDriver driver;
    if ( vm.count("nWorkers") ){
      driver.numWorkers = vm["nWorkers"].as<unsigned int>();  
    }
    driver.submit( job, submitDir );
  }
  else if(vm.count("slurm")){
    system("mkdir -p ~/bin/; ln -s /usr/bin/sbatch ~/bin/bsub;"
    " export PATH=$PATH:~/bin");
    std::string slurmOptions = "-n 1 --cpus-per-task 1 --mem=2000"
    " -p short -t 2:00:00";
    EL::Driver* driver = new EL::LSFDriver;
    job.options()->setBool(EL::Job::optResetShell, false);
    job.options()->setString(EL::Job::optSubmitFlags, slurmOptions);
    driver->submit(job, submitDir);
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
      ("noSmearing", "don't do lepton calibration and smearing") 
      ("electronChannel,e", "run electron selection") 
      ("overwrite,o", "overwrite output folder") 
      ("slurm,s", "run on batch system SLURM") 
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

