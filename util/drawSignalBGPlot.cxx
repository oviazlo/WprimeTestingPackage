/* Copyright 2016 Oleksandr Viazlo */

///*****************************************************************************
///
/// Draw and compare shapes of MC background
/// Script was used mainly for mono-W DM study
/// and for testing...
///
///*****************************************************************************

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

/// std C/C++
#include <fstream>
#include <sstream>

/// private
#include "MyAnalysis/RecoAnalysis.h"
#include "HelpFunctions.h"

// const char *tmpArr[] = {"h_mt","h_mgen","h_pt","h_met"};
const char *tmpArr[] = {"h_mt","h_pt","h_met"};

struct histStruct{
//   unsigned int nHist = 4;
  unsigned int nHist = 3;
  const char **histName = tmpArr;
  map<string,TH1D*> histMap;
};

map<string,string> sampleMap;
Color_t colorArr[] = {kBlack, kBlue, kGreen, kViolet, kOrange, kRed, kYellow, 
  kGreen, kRed, kGreen, kRed, kGreen, kRed, kGreen, kRed, kGreen, kRed, 
  kGreen, kRed, kGreen, kRed, kGreen, kRed, kGreen, kRed, kGreen, kRed,
  kGreen, kRed, kGreen, kRed, kGreen, kRed, kGreen, kRed, kGreen, kRed,
  kGreen, kRed, kGreen, kRed, kGreen, kRed, kGreen, kRed, kGreen, kRed,
};

void setHistStyle(TH1D* inHist, Color_t kColor);
void setHistStyle(histStruct inHistStruct, Color_t kColor);
TH1D* getSummedHistOverSample(SH::SampleHandler sh, string histName);
TH1D* getSummedHistOverIterators(SH::SampleHandler::iterator begin, 
                            SH::SampleHandler::iterator end, string histName);
void sumUpFirstSampleToOther(vector<histStruct> inVecHistStruct);
map<string,TH1D*> getSummedHistMapPerSample(SH::SampleHandler sh, 
                                         string histName);
string getKeyWord(string sampleName);

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
    ("normalize,n", "normalize all histogram to unity 1")
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
  map<string, map<string,TH1D*> > supplementaryMap;
  
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
    
    if (i==0){
      for (int i=0; i<tmpHistStruct.nHist; i++){
        string histName = tmpHistStruct.histName[i];
        supplementaryMap[histName] = getSummedHistMapPerSample(sh,histName);
      }
    }
  }
  
  unsigned int nHistInStruct = myHists[0].nHist;
  
  SetAtlasStyle();
  if (!vm.count("noW"))
    sumUpFirstSampleToOther(myHists);

  TCanvas* tmpCan = new TCanvas("c","c",3508*225/300.,2480*75/300.);
  //TLegend* leg = new TLegend(0.18,0.75,.58,0.99);
  TLegend* leg = new TLegend(0.5,0.75,.95,0.99);
  leg->SetHeader("Samples:");
  tmpCan->Divide(3,1);
  for (int iSample=0; iSample<myHists.size(); iSample++){
    histStruct tmpHistStruct = myHists[iSample];
    for (int iHistType=0; iHistType<tmpHistStruct.nHist; iHistType++){
      string histName = tmpHistStruct.histName[iHistType];
      TH1D* tmpHist = tmpHistStruct.histMap[histName];
      if (vm.count("normalize"))
        tmpHist->Scale(1.0/tmpHist->Integral());
      tmpCan->cd(iHistType+1);
      if (iSample==0){
        gPad->SetLogx();
        gPad->SetLogy();
        tmpHist->GetYaxis()->SetTitle("Entries / pb^{-1}");
        tmpHist->GetXaxis()->SetRangeUser(0,10000);
        if (vm.count("zoom"))
          if (histName=="h_mgen")
            tmpHist->GetXaxis()->SetRangeUser(80,10000);
          else if (histName=="h_mt")
            tmpHist->GetXaxis()->SetRangeUser(40,10000);
          else
            tmpHist->GetXaxis()->SetRangeUser(100,10000);
        if (vm.count("z2"))
          tmpHist->GetXaxis()->SetRangeUser(0,200);
        if (vm.count("yAxis")){
          string tmpStr = vm["yAxis"].as<std::string>();
          stringstream strStream(tmpStr);
          string y1, y2;
          getline (strStream, y1, '-');
          getline (strStream, y2, '-');
          if (y1.size()==0 || y2.size()==0)
            cout << "[WARNING]\tBad input for yAxis: " << tmpStr << endl;
          else
            tmpHist->GetYaxis()->SetRangeUser(atof(y1.c_str()),
                                                   atof(y2.c_str()));
        }
        tmpHist->Draw("H");
        ///********************************************************************
        /// draw supplementary hist
        ///********************************************************************
        if (vm.count("drawSeparately")){
          unsigned int supplementaryHistCounter = 5;
          for (
            map<string,TH1D*>::iterator it = supplementaryMap[histName].begin();
              it!=supplementaryMap[histName].end(); ++it){
            setHistStyle((*it).second,colorArr[supplementaryHistCounter]);
            (*it).second->Draw("H same");
            supplementaryHistCounter++;
          }
        }
        ///********************************************************************
      }
      else
        tmpHist->Draw("H same");
      if (iHistType==0){
        leg->AddEntry(tmpHist,tags[iSample].c_str());
        if (!vm.count("drawSeparately"))
          leg->Draw();
      }
      tmpHist->SaveAs(("pictures/"+tags[iSample]+"_"
      +string(tmpHist->GetName())+".root").c_str());
      
    }
  }
  
  string tmpStr = "";
  stringstream strStream(sampleList);
  getline (strStream, tmpStr, '.');
  
  tmpCan->SaveAs(("pictures/"+tmpStr+".png").c_str());
  tmpCan->SaveAs(("pictures/"+tmpStr+".eps").c_str());
    
  return 0;
}





TH1D* getSummedHistOverSample(SH::SampleHandler sh, string histName){
  return getSummedHistOverIterators(sh.begin(), sh.end(), histName);
}

TH1D* getSummedHistOverIterators(SH::SampleHandler::iterator begin, 
                            SH::SampleHandler::iterator end, string histName){
  TH1D* outHist = NULL;
  bool printDebug = false;
  for (SH::SampleHandler::iterator iter = begin; iter != end; ++ iter){
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

map<string,TH1D*> getSummedHistMapPerSample(SH::SampleHandler sh, 
                                         string histName){
  map<string,TH1D*> outMap;
  SH::SampleHandler::iterator lastKeyWordPosition = sh.begin();
  string lastKeyWord = getKeyWord((*lastKeyWordPosition)->name());
  string currentKeyWord = "";
  
  for (SH::SampleHandler::iterator iter = sh.begin(); iter != sh.end(); 
       ++ iter){
    currentKeyWord = getKeyWord((*iter)->name());
    if (currentKeyWord!=lastKeyWord){
      outMap[lastKeyWord] = getSummedHistOverIterators(lastKeyWordPosition,
                                                       iter, histName);
      lastKeyWord = currentKeyWord;
            lastKeyWordPosition = iter;
    }
    if (iter==--sh.end())
      outMap[lastKeyWord] = getSummedHistOverIterators(lastKeyWordPosition,
                                                       sh.end(), histName);
  }
  return outMap;
}

/// extract keyword from sample name
string getKeyWord(string sampleName){
  string tmpSampleName = "";
  stringstream strStream(sampleName);

  getline (strStream, tmpSampleName, '.'); /// return mc15_13TeV
  getline (strStream, tmpSampleName, '.'); /// return DSID
  
  if ((sampleName.find("MadGraphPythia8EvtGen_A14NNPDF23LO_WWxx_Wlv"))!=
    std::string::npos){
    return tmpSampleName;
  }
  else{
//     ((sampleName.find("PowhegPythia8EvtGen_AZNLOCTEQ6L1_W"))!=
//     std::string::npos){

    getline (strStream, tmpSampleName, '.'); /// return MC generator + mass cut
    
    stringstream strStream2(tmpSampleName);
    string tmpStr = "";
    while(tmpSampleName!=tmpStr){
      tmpStr = tmpSampleName;
      getline (strStream2, tmpSampleName, '_');
    }
    
    /// WARNING FIXME hardcoding to get Wmin and Wplus sample together
    if (tmpStr.find("Wminus")!=string::npos || 
      tmpStr.find("Wplus")!=string::npos){
      tmpStr = "Winclusive"; 
    }
    return tmpStr;
  }
  
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

