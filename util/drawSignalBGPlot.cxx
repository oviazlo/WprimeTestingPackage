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

#include "HelpFunctions.h"

#include <fstream>      // std::ifstream

map<string,string> sampleMap;

int main( int argc, char* argv[] ) {

  po::options_description desc("Options"); 
  desc.add_options()
    ("help,h", "Print help messages") 
    ("sampleList,l", po::value<string>(), "file with list of samples and tags")
    ;
  
  /// get global input arguments:
  po::variables_map vm; 
  const size_t returnedMessage = parseOptionsWithBoost(vm, desc,argc,argv);
  if (returnedMessage!=SUCCESS) std::exit(returnedMessage);

  /// Take the submit directory from the input if provided:
  std::string sampleList = "sampleList.txt";
  if (vm.count("sampleList"))
    sampleList = vm["sampleList"].as<std::string>();
  else
    cout << "[INFO]\tread list of samples from file: " << sampleList << endl;

  std::ifstream sampleListStream;
  sampleListStream.open (sampleList.c_str(), std::ifstream::in);

  while(sampleListStream.good()){
    std::string sample, tag;
    sampleListStream >> sample >> tag; 
  }
  
  /*
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
  
  */
  
  return 0;
}


