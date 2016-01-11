/// EventLoop/xAOD
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

/// ROOT
#include <TSystem.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TColor.h>
#include <Rtypes.h>

/// std C/C++
#include <fstream>

/// private
#include "MyAnalysis/MyxAODAnalysis.h"
#include "HelpFunctions.h"

const char *tmpArr[] = {"h_pt","h_met","h_mt"};

struct histStruct{
  unsigned int nHist = 3;
  const char **histName = tmpArr;
  map<string,TH1D*> histMap;
};

map<string,string> sampleMap;
Color_t colorArr[] = {kBlack, kRed, kGreen, kBlue, kViolet};

void setHistStyle(TH1D* inHist, Color_t kColor);
void setHistStyle(histStruct inHistStruct, Color_t kColor);
TH1D* getSummedHistOverSample(SH::SampleHandler sh, string histName);
void sumUpFirstSampleToOther(vector<histStruct> inVecHistStruct);


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
  
  vector<histStruct> myHists;
  
  histStruct tmpHistStruct;
  
  for (int i=0; i<samples.size(); i++){
    SH::SampleHandler sh;
    cout << "[INFO]\tRead samples from dir: " << samples[i] << endl;
    sh.load (samples[i] + "/hist");
    
    for (int i=0; i<tmpHistStruct.nHist; i++){
      string histName = tmpHistStruct.histName[i];
      tmpHistStruct.histMap[histName] = getSummedHistOverSample(sh,histName);
    }
    
    setHistStyle(tmpHistStruct,colorArr[i]);
    myHists.push_back(tmpHistStruct);
  }

  unsigned int nHistInStruct = myHists[0].nHist;
    
  sumUpFirstSampleToOther(myHists);
  
  TCanvas* tmpCan = new TCanvas("c","c",750,nHistInStruct*500);
  tmpCan->Divide(1,nHistInStruct);
  for (int iSample=0; iSample<myHists.size(); iSample++){
    histStruct tmpHistStruct = myHists[iSample];
    for (int iHistType=0; iHistType<tmpHistStruct.nHist; iHistType++){
      string histName = tmpHistStruct.histName[iHistType];
      TH1D* tmpHist = tmpHistStruct.histMap[histName];
      tmpCan->cd(iHistType+1);
      if (iSample==0){
        tmpHist->Draw();
        tmpHist->GetXaxis()->SetRangeUser(30,10000);
        gPad->SetLogx();
        gPad->SetLogy();
      }
      else
        tmpHist->Draw("same");
    }
  }
  tmpCan->SaveAs("bigTest.png");
  
  
  return 0;
}

TH1D* getSummedHistOverSample(SH::SampleHandler sh, string histName){
  TH1D* outHist = NULL;
  
  for (SH::SampleHandler::iterator iter = sh.begin(); iter != sh.end(); 
      ++ iter){
    TH1D* tmpHist = (TH1D*)(*iter)->readHist (histName.c_str());

    if (tmpHist==NULL){
      cout << "[ERROR]\tCannot read hist " << histName << " of typy TH1D "
              "from sample " << (*iter)->name() << ". Skip sample!" << endl;
      continue;
    }

    if (outHist==NULL)
      outHist = (TH1D*)tmpHist->Clone((histName+"_summed").c_str());
    else
      outHist->Add(tmpHist);
  }
  
  return outHist;
}

void setHistStyle(histStruct inHistStruct, Color_t kColor){
  for (int i=0; i<inHistStruct.nHist; i++){
    string histName = inHistStruct.histName[i];
    setHistStyle(inHistStruct.histMap[histName],kColor);
  }
}

void setHistStyle(TH1D* inHist, Color_t kColor){
  if (inHist==NULL)
    return;
  inHist->SetLineWidth(2);
  inHist->SetLineColor(kColor);
}

void sumUpFirstSampleToOther(vector<histStruct> inVecHistStruct){
  histStruct firstSample = inVecHistStruct[0];
  for (int iSample=1; iSample<inVecHistStruct.size(); iSample++){
    histStruct tmpHistStruct = inVecHistStruct[iSample];
    for (int iHistType=0; iHistType<tmpHistStruct.nHist; iHistType++){
      string histName = tmpHistStruct.histName[iHistType];
      TH1D* firstSampleHist = firstSample.histMap[histName];
      TH1D* tmpHist = tmpHistStruct.histMap[histName];
      tmpHist->Add(firstSampleHist);
    }
  }
}

