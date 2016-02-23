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
Color_t colorArr[] = {kOrange, kRed, kGreen, kYellow, kBlue, kViolet,
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
    ("samplesToDraw,s",po::value<string>()," list of samples to draw, separated by coma")
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
  
  vector<string> samplesToDraw;
  if (vm.count("samplesToDraw")){
    samplesToDraw = GetWords(vm["samplesToDraw"].as<std::string>(),',');
    cout << "[INFO]\tRunning over saples: ";
    vector<string> supportedSamples = mergedSample->GetAllSupportedGlobalSampleTags();
    std::vector<string>::iterator it;
    for (int i=0; i<samplesToDraw.size(); i++){
      cout << samplesToDraw[i] << " ";
      it = find (supportedSamples.begin(), supportedSamples.end(), samplesToDraw[i]);
      if (it == supportedSamples.end()){
        cout << "Sample *" << samplesToDraw[i] << "* is not supported by code!!!" << endl;
        cout << "Available options:" << endl;
        for (int k=0; k<supportedSamples.size(); k++)
          cout << supportedSamples[k] << " ";
        cout << endl;
        cout << "Terminate execution!!!" << endl;
        return 0;
      }
    }
    cout << endl;
  }
  
  for (int i=0; i<samples.size(); i++){
    SH::SampleHandler sh;
    cout << "[INFO]\tRead samples from dir: " << samples[i] << endl;
    sh.load (samples[i] + "/hist");
    mergedSample->AddSampleHandler(sh,samples[i]);
  }

  string prefix = "muon/stage_final/hObjDump_";
  vector<string> plotsToDraw = {"pt","met","mt","eta","phi"};
  
  SetAtlasStyle();
//   gStyle->SetHistTopMargin(0.);
  TCanvas* tmpCan = new TCanvas("c","c",3508*150/300.,2480*150/300.);
  
  
  for (int i=0; i<plotsToDraw.size(); i++){
      
    THStack *hs = new THStack("hs","Stacked 1D histograms");
    
    if (!vm.count("samplesToDraw"))
      samplesToDraw = {"diboson","z","top","w"};
    
    TH1D* testHist;
    
    for (int k=0; k<samplesToDraw.size(); k++){
      testHist = mergedSample->GetMergedHist(samplesToDraw[k],prefix+plotsToDraw[i]);
      if (testHist!=NULL){
        setHistStyle(testHist,colorArr[k]);
        if (samplesToDraw[k]=="wmunu_massbinned")
          testHist->Scale(10); /// FIXME hardcoded scaling
        hs->Add(testHist);
      }
      else{
        cout << "Plot *" << plotsToDraw[i] << "* is empty for sample " << samplesToDraw[k] 
        << endl; 
      }
    }
    
    hs->Draw("HIST");
    if (i<3){
      hs->SetMinimum(10E-7);
      hs->SetMaximum(10E5);
    }
    else{
      hs->SetMinimum(10E-1);
      hs->SetMaximum(10E4);
    }
    
    hs->GetXaxis()->SetTitle(testHist->GetXaxis()->GetTitle());
    hs->GetYaxis()->SetTitle(testHist->GetYaxis()->GetTitle());
    
    gPad->Update();
    if (i>=3){
      gPad->SetLogx(0);
      gPad->SetLogy(0);
    }
    else{
      gPad->SetLogy();
      gPad->SetLogx();
    }
    
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

// std::vector<std::string> split(const std::string &s, char delim) {
//     std::vector<std::string> elems;
//     split(s, delim, elems);
//     return elems;
// }