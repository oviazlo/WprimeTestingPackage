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

#include "TStyle.h"

#include "TROOT.h"

//~ #include "AtlasStyle.C"
//~ #include "AtlasUtils.C"
//~ #include "AtlasUtils.h"

#include "MyAnalysis/MyxAODAnalysis.h"
#include <EventLoopAlgs/NTupleSvc.h>
#include <EventLoop/OutputStream.h>

/// boost
#include "boost/program_options.hpp"
#include <boost/algorithm/string.hpp>

#include "HelpFunctions.h"

map<string,string> sampleMap;


int main( int argc, char* argv[] ) {

  po::options_description desc("Options"); 
  desc.add_options()
    ("help,h", "Print help messages") 
    ("folder,f", po::value<string>(), "name of folder to read")
    ("drawHists,d", "draw histograms")
    ;
  
  /// get global input arguments:
  po::variables_map vm; 
  const size_t returnedMessage = parseOptionsWithBoost(vm,desc,argc,argv);
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
  
  /// get histogram
  TH1I* cutflowHist = (TH1I*)mySample->readHist ("cutflow_hist");
  cout << "[DEBUG]\tcutflowHist = " << cutflowHist << endl;
  
  for (int i=1; i<=cutflowHist->GetNbinsX(); i++){
    int binContent = cutflowHist->GetBinContent(i);
    if (binContent<=0) break;
    string binLabel = cutflowHist->GetXaxis()->GetBinLabel(i);
    cout << binLabel << ":\t" << binContent << endl;
  }

  if (vm.count("drawHists")){
    SetAtlasStyle();
    
    /// read pt, MET, Mt histograms, make them pretty and save them
    TH1D* h_pt = (TH1D*)mySample->readHist ("h_pt");
    TH1D* h_MET = (TH1D*)mySample->readHist ("h_met");
    TH1D* h_Mt = (TH1D*)mySample->readHist ("h_mt");
    
//     h_pt->GetXaxis()->SetTitle("p_{T}^{#mu} [GeV]");
//     h_MET->GetXaxis()->SetTitle("E_{T}^{miss} [GeV]");
//     h_Mt->GetXaxis()->SetTitle("m_{T} [GeV]");
//     
//     h_pt->GetYaxis()->SetTitle("Counts");
//     h_MET->GetYaxis()->SetTitle("Counts");
//     h_Mt->GetYaxis()->SetTitle("Counts");
//     
//     h_pt->SetTitle("");
//     h_MET->SetTitle("");
//     h_Mt->SetTitle("");
    
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

