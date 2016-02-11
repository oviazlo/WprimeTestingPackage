/// EventLoop/xAOD
#include "xAODRootAccess/Init.h"
#include "SampleHandler/SampleHandler.h"
#include "SampleHandler/Sample.h"
#include "SampleHandler/ToolsDiscovery.h"
#include "SampleHandler/DiskListLocal.h"
#include <SampleHandler/ScanDir.h>
#include <SampleHandler/ToolsJoin.h>
#include "EventLoop/Job.h"
#include "EventLoop/ProofDriver.h"
#include "EventLoop/DirectDriver.h"
#include <EventLoop/OutputStream.h>
#include <EventLoopAlgs/NTupleSvc.h>

/// ROOT
#include <TSystem.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TColor.h>
#include <Rtypes.h>
#include <TLegend.h>
#include <THStack.h>

/// std C/C++
#include <fstream>
#include <sstream>

/// private
#include "MyAnalysis/RecoAnalysis.h"
#include "HelpFunctions.h"
#include "histPlotter/WprimeMergedSample.h"

map<string,string> sampleMap;
Color_t colorArr[] = {kBlack, kBlue, kGreen, kViolet, kOrange, kRed, kYellow, 
  kGreen, kRed, kGreen, kRed, kGreen, kRed, kGreen, kRed, kGreen, kRed, 
  kGreen, kRed, kGreen, kRed, kGreen, kRed, kGreen, kRed, kGreen, kRed,
  kGreen, kRed, kGreen, kRed, kGreen, kRed, kGreen, kRed, kGreen, kRed,
  kGreen, kRed, kGreen, kRed, kGreen, kRed, kGreen, kRed, kGreen, kRed,
};

void setHistStyle(TH1D* inHist, Color_t kColor);

int main( int argc, char* argv[] ) {

  po::options_description desc("Options"); 
  desc.add_options()
    ("help,h", "Print help messages") 
    ("sampleList,l", po::value<string>(), "file with list of samples and tags")
    ("zoom,z","zoom x-axis")
    ("yAxis,y",po::value<string>(),"set y-range; syntax: y1-y2")
    ("z2","zoom x-axis to lower values")
    ("noW","do not make summing up of first sample to others")
    ("drawSeparately,s", "draw each sample separately")
    ;
  
  /// get global input arguments:
  po::variables_map vm; 
  const size_t returnedMessage = parseOptionsWithBoost(vm, desc,argc,argv);
  if (returnedMessage!=SUCCESS) std::exit(returnedMessage);

  /// Take the submit directory from the input if provided:
  std::string sampleList = "sampleList.txt";
  if (vm.count("sampleList"))
    sampleList = vm["sampleList"].as<std::string>();
  
  cout << "[INFO]\tread list of samples from file: " << sampleList << endl;

  std::ifstream sampleListStream;
  sampleListStream.open (sampleList.c_str(), std::ifstream::in);

  vector<string> samples, tags;
  
  while(sampleListStream.good()){
    std::string sample = "";
    std::string tag = "";
    sampleListStream >> sample >> tag; 
    if (sample == "" || tag == "")
      continue;
    samples.push_back(sample);
    tags.push_back(tag);
    cout << "[INFO]\tread sample-tag pair: {" << sample << "," << tag << "}" 
         << endl;
  }
  
  if (samples.size()!=tags.size() || samples.size()==0){
    cout << "[ERROR]\tsamples.size() = " << samples.size() << 
            ", tags.size() = " << tags.size() << endl;
            return -1;
  }
  
  WprimeMergedSample *mergedSample = new WprimeMergedSample();
  
  for (int i=0; i<samples.size(); i++){
    SH::SampleHandler sh;
    cout << "[INFO]\tRead samples from dir: " << samples[i] << endl;
    sh.load (samples[i] + "/hist");
    mergedSample->AddSampleHandler(sh,samples[i]);
  }

  vector<string> plotsToDraw = {"h_mgen","h_mt","h_pt","h_met"};
  
  SetAtlasStyle();
//   gStyle->SetHistTopMargin(0.);
  TCanvas* tmpCan = new TCanvas("c","c",3508*150/300.,2480*150/300.);
  
  
  for (int i=0; i<plotsToDraw.size(); i++){
      
    THStack *hs = new THStack("hs","Stacked 1D histograms");
    TH1D* testHist = mergedSample->GetMergedHist("wmunu",plotsToDraw[i]);
    TH1D* testHist2 = mergedSample->GetMergedHist("test2",plotsToDraw[i]);
    
    setHistStyle(testHist,kOrange);
    setHistStyle(testHist2,kRed);
    
    hs->Add(testHist);
    hs->Add(testHist2);
    
    hs->Draw("HIST");
    hs->SetMinimum(10E-10);
    hs->SetMaximum(1000);
    hs->GetXaxis()->SetTitle(testHist->GetXaxis()->GetTitle());
    hs->GetYaxis()->SetTitle(testHist->GetYaxis()->GetTitle());
    
    gPad->Update();
    gPad->SetLogy();
    gPad->SetLogx();
    
    string outFileName = "pictures/testHistPlotter.ps";
    if (i==0&&plotsToDraw.size()>1)
        outFileName += "(";
    if (i==(plotsToDraw.size()-1)&&plotsToDraw.size()>1)
        outFileName += ")";
    tmpCan->SaveAs(outFileName.c_str());
  }
    
  return 0;
}

void setHistStyle(TH1D* inHist, Color_t kColor){
  if (inHist==NULL)
    return;
  inHist->SetFillColor(kColor);
  inHist->SetMarkerStyle(21);
  inHist->SetMarkerColor(kColor);

}
