/* Copyright 2016 Oleksandr Viazlo */

///*****************************************************************************
///
/// Main script. Runs over samples, makes analysis selection and produce 
/// histograms. Can be run in parallel (EventLoop) or standalone.
///
///*****************************************************************************

/// WARNING HOWTO run the script for MC/data
///
/// runMainLoop -p data15_13TeV.0027* -t p2495 -f finalData_0027 -o
/// runMainLoop -p mc15_13TeV.301*EXOT9* -t p2495_r7 -f finalMC_massBinnedDY -o
/// 
/// -p <pattern> 
/// run over samples which satisfy this pattern.
/// inclusive MC pattern:
/// diboson: 3610*
/// top: 4100*
/// mass-binned DY: 301*
/// inclusive W: 3611*L1_W*
/// inclusive Z: 3611*L1_Z*
///
/// -t <dir>
/// specify directory where to look for data.
/// prefix is hardcoded in defaultConfigurations.h file
/// /nfs/shared/pp/oviazlo/xAOD/<dir>
/// 
/// -f <outDir>
/// specify output directory where EventLoop will write histograms
///
/// -o
/// allow overwritting of output directory.
///
/// --directDriver
/// run script in standalone mode (w/o paralleling). useful for debugging.
///
/// description for other flags is in parseOptionsWithBoost() function.

/// TODO
/// 1) default patterns... and call them with separate flags
/// 2) cross check electron channel

/// EventLoop/xAOD headers
#include "xAODRootAccess/Init.h"
#include "EventLoop/Job.h"
#include "EventLoop/ProofDriver.h"
#include "EventLoop/DirectDriver.h"
#include "EventLoop/LSFDriver.h"
#include "SampleHandler/DiskListLocal.h"
#include "SampleHandler/ToolsSplit.h"
#include "SampleHandler/SampleHandler.h"
#include "SampleHandler/ToolsDiscovery.h"
#include "SampleHandler/ScanDir.h"
#include "SampleHandler/Sample.h"
#include "SampleHandler/ToolsJoin.h"
#include <EventLoopAlgs/NTupleSvc.h>
#include <EventLoop/OutputStream.h>

/// ROOT
#include <TSystem.h>

/// private
#include "MyAnalysis/RecoAnalysis.h"
#include "defaultConfigurations.h"

/// boost
#include "boost/program_options.hpp"
#include <boost/algorithm/string.hpp>
#include "boost/filesystem.hpp"
#include "boost/assign.hpp" /// map_list_of

/// C/C++ 
#include <stdlib.h>     /* getenv */
#include <map>

namespace 
{ 
  const size_t ERROR_IN_COMMAND_LINE = 1; 
  const size_t SUCCESS = 0; 
  const size_t ERROR_UNHANDLED_EXCEPTION = 2; 
  const size_t HELP_CALLED = 3;
  
  namespace po = boost::program_options;
 
  const size_t CERN=0;
  const size_t IRIDIUM=1;
  const size_t AURORA_HEP=2;
  const size_t AURORA_SNIC=3;
  
} /// namespace 

map<size_t,string> systemMap = boost::assign::map_list_of (CERN,"CERN") 
(IRIDIUM,"IRIDIUM") (AURORA_HEP,"AURORA_HEP") (AURORA_SNIC,"AURORA_SNIC");

int parseOptionsWithBoost(po::variables_map &vm, int argc, char* argv[]);
size_t systemType;

int main( int argc, char* argv[] ) {
 
  /// get global input arguments:
  po::variables_map vm; 
  const size_t returnedMessage = parseOptionsWithBoost(vm,argc,argv);
  if (returnedMessage!=SUCCESS) std::exit(returnedMessage);

  /// Take the submit directory from the input if provided:
  
  std::string submitDir = configMap["submitDir"];
  if ( vm.count("folder") ){
    submitDir = vm["folder"].as<std::string>();
    if (submitDir.find("submitDirs/")==std::string::npos)
      submitDir = "submitDirs/" + submitDir;
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
  
  
  std::string strSamplePattern = configMap["strSamplePattern"];
  if ( vm.count("samplePattern") ){
    strSamplePattern = vm["samplePattern"].as<std::string>();
  }
  
  std::size_t found = hostName.find("cern");
  if (found!=std::string::npos)
    systemType = CERN;
  else{
    /// WARNING before was aurora
    std::size_t found = hostName.find("au");
    if (found!=std::string::npos){
      found = hostName.find("aurora"); /// front-end node --> use SNIC nodes (data storage is visible only from HEP nodes) 
      if (found!=std::string::npos)
        systemType = AURORA_SNIC;
      else{
        string tmpStrAurora (hostName.begin()+2,hostName.end());
        std::string::size_type sz;   // alias of size_t
        int intFromStr = std::stoi (tmpStrAurora,&sz);
        if ((intFromStr>=193) && (intFromStr<=216)) /// HEP nodes
          systemType = AURORA_HEP;
        else
          systemType = AURORA_SNIC;
      }
    }
    else
      systemType = IRIDIUM;
  }
  
  cout << "[JobSetup]\tCode is running on system: " << systemMap[systemType] 
  << endl;
  
  string pathToExtend = "";
  
  if (systemType == CERN){
    pathToExtend = configMap["pathToExtend_CERN"];
  }
  else if (systemType == AURORA_HEP){
    pathToExtend = configMap["pathToExtend_AURORA_HEP"];
  }
  else if (systemType == AURORA_SNIC){
    pathToExtend = configMap["pathToExtend_AURORA_SNIC"];
  }
  else if (systemType == IRIDIUM){
    pathToExtend = configMap["pathToExtend_IRIDIUM"];
  }

  if(vm.count("sampleTag"))
    pathToExtend += vm["sampleTag"].as<std::string>();
  else
    pathToExtend += "cutFlow";
  
  inputFilePath = gSystem->ExpandPathName(pathToExtend.c_str());

  cout << "[JobSetup]\tLooking for a sample pattern: " << strSamplePattern 
  << endl << endl;
  
  SH::ScanDir()
  .samplePattern (strSamplePattern)
  .scan (sh, inputFilePath);

  if ( vm.count("mergeSamples") ){
    string sampleMergePattern;
    if (strSamplePattern.find("data")!=std::string::npos)
      sampleMergePattern = "data15_13TeV.*";
    else
      sampleMergePattern = "mc15_13TeV.*";
    cout << "[JobSetup]\tMake attampt of merging sample with pattern: " 
    << sampleMergePattern << " to one sample with name: " << 
    vm["mergeSamples"].as<std::string>() << endl;
    SH::mergeSamples (sh, vm["mergeSamples"].as<std::string>(), 
                      sampleMergePattern); 
  }
  
  /// Print what we found:
  if (vm.count("debug"))
    sh.print();
  
  /// Set the name of the input TTree. It's always "CollectionTree"
  /// for xAOD files.
  sh.setMetaString( "nc_tree", "CollectionTree" );

  
  /// scan the number of events in each root file 
  if ( vm.count("nEventsPerJob") )
    SH::scanNEvents (sh);
  
  /// Create an EventLoop job:
  EL::Job job;
  job.sampleHandler( sh );

  /// specify nEventsPerJob
  if ( vm.count("nEventsPerJob") )
    sh.setMetaDouble (EL::Job::optEventsPerWorker, 
                      vm["nEventsPerJob"].as<unsigned int>());
  
  /// Specify nEvents to run on
  if (nEvents!=-1)
    job.options()->setDouble(EL::Job::optMaxEvents, nEvents);

  if ( vm.count("nFilesPerJob") )
    job.options()->setDouble (EL::Job::optFilesPerWorker, 
                              vm["nFilesPerJob"].as<unsigned int>());

  /// TODO to find optimal values to use
  job.options()->setDouble (EL::Job::optCacheSize, 10*1024*1024);
  job.options()->setDouble (EL::Job::optCacheLearnEntries, 50);
    
  /// Add our analysis to the job:
  RecoAnalysis* alg = new RecoAnalysis();
  job.algsAdd( alg );
  
  /// all configuration to the algorithm should be passed after attachment to 
  /// the job (practically after this comment)
  ///
  /// syntax of passing configurations:
  ///   alg->m_myMegaFlag = true;
  ///
  /// [List of possble flags to use]
  /// bool m_useHistObjectDumper; - do not make default plots by 
  ///                                           HistObjectDumper
  /// bool m_useBitsetCutflow; - do not save cutflow
  /// bool m_useCalibrationAndSmearingTool; - don't do muon calibration and 
  ///                                         smearing TODO not in the code 
  ///                                         now?
  /// bool m_runElectronChannel; - run electron cycle instead of muon one
  /// bool m_doWprimeTruthMatching; - do truth matching to identify Wprime 
  ///                                 decay to muon/electron channel
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
  
  if (vm.count("info"))
    alg->setMsgLevel (MSG::INFO);
    
  if ( vm.count("proofDriver") ){/// Run the job using the local/direct driver:
    EL::ProofDriver driver;
    if ( vm.count("nWorkers") ){
      driver.numWorkers = vm["nWorkers"].as<unsigned int>();  
    }
    driver.submit( job, submitDir );
  }
  else if(vm.count("directDriver")){
    EL::DirectDriver driver;
    driver.submit( job, submitDir );
  }
  else{/// Local/Direct Driver:
    EL::LSFDriver* driver = new EL::LSFDriver;
    std::string slurmSystemDependentOptions;
    
    if (systemType == AURORA_HEP){
      system("mkdir -p ~/bin/; ln -s /usr/bin/sbatch"
      " ~/bin/bsub; export PATH=$PATH:~/bin");
      slurmSystemDependentOptions = "-n 1 --cpus-per-task 1"
    " -p hep -t 48:00:00 -A hep2016-1-1 ";
    }
    else if (systemType == AURORA_SNIC){
      system("mkdir -p ~/bin/; ln -s /usr/bin/sbatch"
      " ~/bin/bsub; export PATH=$PATH:~/bin");
      slurmSystemDependentOptions = "-n 1 --cpus-per-task 1"
    " -p snic -t 48:00:00 -A snic2015-2-63";
    }
    else if (systemType == IRIDIUM){
      system("mkdir -p ~/bin/; ln -s /usr/bin/sbatch ~/bin/bsub;"
      " export PATH=$PATH:~/bin");
      slurmSystemDependentOptions = "-n 1 --cpus-per-task 1"
    " -p long -t 2:00:00";
    }
    else if (systemType == CERN){
      slurmSystemDependentOptions = "-L /bin/bash -q 1nh";
      driver->shellInit = "export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/"
      "repo/ATLASLocalRootBase && source ${ATLAS_LOCAL_ROOT_BASE}/user/"
      "atlasLocalSetup.sh";
    }
    
    job.options()->setBool(EL::Job::optResetShell, false);
    job.options()->setString(EL::Job::optSubmitFlags, 
                             slurmSystemDependentOptions);
    driver->submit(job, submitDir);
    
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
      ("nFilesPerJob", po::value<unsigned int>(), "number of files per job")
      ("nEventsPerJob", po::value<unsigned int>(), "number of events per job")
      ("proofDriver", "run with ProofDriver - PROOF-Lite mode") 
      ("noSmearing", "don't do lepton calibration and smearing") 
      ("electronChannel,e", "run electron selection") 
      ("overwrite,o", "overwrite output folder") 
      ("directDriver", "run with DirectDriver") 
      ("info", "set message level to INFO") 
      ("debug", "print debug printouts from code") 
      ("mergeSamples", po::value<string>(),"merge everything in one sample; "
      "specify final sample name")
      ("samplePattern,p", po::value<string>(),"specify Sample Pattern")
      ("sampleTag,t", po::value<string>(),"specify Sample tag to use")
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

