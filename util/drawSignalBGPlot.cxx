#include "xAODRootAccess/Init.h"
#include "SampleHandler/SampleHandler.h"
#include "SampleHandler/Sample.h"
#include "SampleHandler/ToolsDiscovery.h"
#include "EventLoop/Job.h"
#include "EventLoop/ProofDriver.h"
#include "EventLoop/DirectDriver.h"
#include "SampleHandler/DiskListLocal.h"
#include <SampleHandler/ScanDir.h>
#include <SampleHandler/ToolsJoin.h>
#include <TSystem.h>
#include <TH1F.h>
#include <TCanvas.h>

//~ #include "AtlasStyle.C"
//~ #include "AtlasUtils.C"
//~ #include "AtlasUtils.h"

#include "MyAnalysis/MyxAODAnalysis.h"
#include <EventLoopAlgs/NTupleSvc.h>
#include <EventLoop/OutputStream.h>

/// boost
#include "boost/program_options.hpp"
#include <boost/algorithm/string.hpp>

#include <HelpFunctions.cxx>

map<string,string> sampleMap;

namespace 
{ 
  const size_t ERROR_IN_COMMAND_LINE = 1; 
  const size_t SUCCESS = 0; 
  const size_t ERROR_UNHANDLED_EXCEPTION = 2; 
  const size_t HELP_CALLED = 3;
  
  namespace po = boost::program_options; 
 
} /// namespace  

int parseOptionsWithBoost(po::variables_map &vm, po::options_description desc, 
                          int argc, char* argv[]);
TStyle* AtlasStyle();
void SetAtlasStyle();

int main( int argc, char* argv[] ) {

  po::options_description desc("Options"); 
  desc.add_options()
    ("help,h", "Print help messages") 
    ("folder,f", po::value<string>(), "name of folder to read")
    ("drawHists,d", "draw histograms")
    ;
  
  /// get global input arguments:
  po::variables_map vm; 
  const size_t returnedMessage = parseOptionsWithBoost(vm, desc,argc,argv);
  if (returnedMessage!=SUCCESS) std::exit(returnedMessage);

  /// Take the submit directory from the input if provided:
  std::string folder = "submitDir";
  if (vm.count("folder"))
    folder = vm["folder"].as<std::string>();
  else
    cout << "[INFO]\tread cutflow from default directory: " << folder << endl;

  /// Construct the samples to run on:
  SH::SampleHandler sh;

  sh.load ((folder + "/hist").c_str());

  cout << "[INFO]\tSize of read sample: " << sh.size() << endl;
  if (sh.size()<1){
    cout << "[ERROR]\tSize of sample is zero... Aborting execution!" << endl;
    return 0;
  }
  
  SH::Sample* mySample = *(sh.begin());
  
  std::size_t found = mySample->name().find("mc15");
  string mergePattern;
  if (found!=std::string::npos)
    mergePattern = "mc15_13TeV.*";
  else
    mergePattern = "data15_13TeV.*";

  cout << "[INFO]\tNumber of samples in directory: " << sh.size() << endl;
  cout << "[INFO]\tUse merge pattern: " << mergePattern << endl;
  SH::mergeSamples (sh, "final", mergePattern);  
  cout << "[INFO]\tNumber of samples after merging: " << sh.size() << endl;

  if (sh.size()!=1){
    cout << "[ERROR]\tNumber of samples after merging are not equal 1. "
    "Print all samples after merging:" << endl;
    sh.print();
    return -1;
  }
  

  if (vm.count("drawHists")){
    SetAtlasStyle();
    
    /// read pt, MET, Mt histograms, make them pretty and save them
    TH1D* h_pt = (TH1D*)mySample->readHist ("h_pt");
    TH1D* h_MET = (TH1D*)mySample->readHist ("h_met");
    TH1D* h_Mt = (TH1D*)mySample->readHist ("h_mt");
    
    TCanvas *can = new TCanvas("can","can",800,800);
    gPad->SetLogy();
    
    h_pt->Draw();
    h_pt->GetXaxis()->SetRangeUser(0,1000.0);
    can->SaveAs("pT.png");
    
    h_MET->GetXaxis()->SetRangeUser(0,1000.0);
    h_MET->Draw();
    can->SaveAs("MET.png");
    
    h_Mt->GetXaxis()->SetRangeUser(0,1000.0);
    h_Mt->Draw();
    can->SaveAs("mT.png");
  }
  
  return 0;
}



/// parse input arguments
int parseOptionsWithBoost(po::variables_map &vm, po::options_description desc, 
                          int argc, char* argv[]){

  try 
  { 
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


