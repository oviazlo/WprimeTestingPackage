/// EventLoop/xAOD headers
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
#include <EventLoopAlgs/NTupleSvc.h>
#include <EventLoop/OutputStream.h>

/// ROOT headers
#include <TSystem.h>
#include <TH1F.h>
#include <TCanvas.h>
#include "TStyle.h"
#include "TROOT.h"

/// std C/C++ headers
#include <stdlib.h>

/// boost
#include "boost/program_options.hpp"
#include <boost/algorithm/string.hpp>

/// private headers
#include "MyAnalysis/RecoAnalysis.h"
#include "HelpFunctions.h"

map<string,string> sampleMap;

int main( int argc, char* argv[] ) {

  po::options_description desc("Options"); 
  desc.add_options()
    ("help,h", "Print help messages") 
    ("folder,f", po::value<string>(), "name of folder to read")
    ("drawHists,d", "draw histograms")
    ("useOnlyOneFile", "use only one file from folder")
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

  cout << "[INFO]\tRead samples from dir: " << folder << endl;
  sh.load (folder + "/hist");
 
  if (sh.size()==0){
    cout << "sh.size() = " << sh.size() << endl;
    cout << "Terminate!!!"<< endl;
    return 0;
  }
  
  SH::Sample* mySample = *(sh.begin());
  TH1I* cutflowHist = (TH1I*)mySample->readHist ("cutflow_hist");
  cout << "[INFO]\tsh.size() = " << sh.size() << endl;
  
  if (sh.size()>1){
  
    if (vm.count("useOnlyOneFile")){
      int counter = 0;
      cout << "[WARNING]\tThere are more than one sample!" << endl;
      for (SH::SampleHandler::iterator iter = sh.begin(); iter != sh.end(); 
        ++ iter){
        cout << counter << ": " << (*iter)->name() << endl;
        counter++;
      }
      int sampleToUse = 0;
      cout << "[WARNING]\tPlease choose which one to use from list above: ";
      cin >> sampleToUse;
      if (sampleToUse<0)
        sampleToUse = 0;
      if (sampleToUse>=sh.size())
        sampleToUse = sh.size()-1;
      
      mySample = sh.at(sampleToUse);
      cutflowHist = (TH1I*)mySample->readHist ("cutflow_hist");
      
      cout << "[INFO]\tUsing samle " << sampleToUse << ": " 
      << mySample->name() << endl;
    }
    else{ /// get cutflow_hist from each file add sum them up
      for (SH::SampleHandler::iterator iter = sh.begin()+1; iter != sh.end(); 
        ++ iter){
        TH1I* tmpHist = (TH1I*)(*iter)->readHist ("cutflow_hist");
        cutflowHist->Add(tmpHist);
      }
    }
  }
    
  cout << "[DEBUG]\tcutflowHist = " << cutflowHist << endl;
  
  for (int i=1; i<=cutflowHist->GetNbinsX(); i++){
    int binContent = cutflowHist->GetBinContent(i);
    if (binContent<=0) break;
    string binLabel = cutflowHist->GetXaxis()->GetBinLabel(i);
    cout << binLabel << ":\t" << binContent << endl;
  }

  /// FIXME deprecaged functionality; remove it
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
    can->SaveAs("pictures/pT.png");
    
    h_MET->GetXaxis()->SetRangeUser(0,1000.0);
    h_MET->Draw();
    can->SaveAs("pictures/MET.png");
    
    h_Mt->GetXaxis()->SetRangeUser(0,1000.0);
    h_Mt->Draw();
    can->SaveAs("pictures/mT.png");
  }
  
  return 0;
}

