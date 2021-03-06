/* Copyright 2016 Oleksandr Viazlo */

///*****************************************************************************
///
/// Plotting script. Plot MCs as THStack. Data - TH1 histogram.
/// Still in developing and testing state
///
///*****************************************************************************

/// WARNING HOWTO run the script
///
/// testHistPlotter -l finalAll_v4.txt -s w,data
/// testHistPlotter -l finalAll_v4.txt
///
/// -l <fileName>
/// list of all directories where histrograms produced by runMainLoop script are
/// can contain two columns: directory and tag. But second column is optional.
/// example:
///  submitDirs/finalData_0027 finalData_0027
///  submitDirs/finalData_0028_m150 finalData_0028_m150
///  submitDirs/finalData_0028_m151 finalData_0028_m151
///  submitDirs/finalMC_diboson finalMC_diboson
///  submitDirs/finalMC_inclusiveZ finalMC_inclusiveZ
///  submitDirs/finalMC_inclusiveW finalMC_inclusiveW
///  submitDirs/finalMC_massBinnedDY finalMC_massBinnedDY
///  submitDirs/finalMC_top finalMC_top
///  submitDirs/finalMC_missingDiboson finalMC_missingDiboson
///
/// -s <sample1,sample2,sample3,...>
/// list of inclusive samples to draw separated by comma
/// allowed sample names to use are specified here in 
/// WprimeMergedSample::WprimeMergedSample() function.
/// if this flag is not used script will use plot default samples: 
/// "diboson","z","top","w","data"

/// TODO
/// 1) add possibility to make plots w/o data
/// 2) draw errors on plots and on ratio plots
/// 3) use tags from input list

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
po::variables_map vm; 
Color_t colorArr[] = {kOrange, kAzure-9, kRed+1, kWhite, kYellow, kBlue, kViolet,
  kGreen, kRed, kGreen, kRed, kGreen, kRed, kGreen, kRed, kGreen, kRed, 
  kGreen, kRed, kGreen, kRed, kGreen, kRed, kGreen, kRed, kGreen, kRed,
  kGreen, kRed, kGreen, kRed, kGreen, kRed, kGreen, kRed, kGreen, kRed,
  kGreen, kRed, kGreen, kRed, kGreen, kRed, kGreen, kRed, kGreen, kRed,
};

void setHistStyle(TH1D* inHist, Color_t kColor);
/// parse name of file from name of input list
TFile* createOutFile(string listName);
vector<string> getSamplesToDraw(WprimeMergedSample *mergedSample);

int main( int argc, char* argv[] ) {

  po::options_description desc("Options"); 
  desc.add_options()
    ("help,h", "Print help messages") 
    ("sampleList,l", po::value<string>(), "file with list of samples and tags")
    ("samplesToDraw,s",po::value<string>()," list of samples to draw, separated by coma")
    ("histFolder,f",po::value<string>()," specify hist folder to use")
    ;
  
  /// get global input arguments:
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
    if (sample == "")
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
  
  TFile* outFile = createOutFile(sampleList);
  
  WprimeMergedSample *mergedSample = new WprimeMergedSample();
  
  vector<string> samplesToDraw = getSamplesToDraw(mergedSample);
  if (samplesToDraw.size()==0)
    return 0;
  
  for (int i=0; i<samples.size(); i++){
    SH::SampleHandler sh;
    cout << "[INFO]\tRead samples from dir: " << samples[i] << endl;
    sh.load (samples[i] + "/hist");
    mergedSample->AddSampleHandler(sh,samples[i]);
  }

  //string prefix = "muon/stage_final_noWeight/hObjDump_";
  string histFolderName = "final";
  if (vm.count("histFolder"))
    histFolderName = vm["histFolder"].as<std::string>();
  
  string prefix = "muon/stage_" + histFolderName + "/hObjDump_";
  vector<string> plotsToDrawWithPrefix = {"pt","met","mt","eta","phi"};
  vector<string> plotsToDrawWoPrefix = {"cutflow_hist"};
  
  map<string,string> prefixMap;
  
  vector<string> plotsToDraw;
  plotsToDraw.reserve(plotsToDrawWithPrefix.size()+plotsToDrawWoPrefix.size());
  plotsToDraw.insert(plotsToDraw.end(),plotsToDrawWithPrefix.begin(),plotsToDrawWithPrefix.end());
  plotsToDraw.insert(plotsToDraw.end(),plotsToDrawWoPrefix.begin(),plotsToDrawWoPrefix.end());
  
  for (int i=0; i<plotsToDrawWithPrefix.size(); i++){
    prefixMap[plotsToDrawWithPrefix[i]] = prefix;
  }
  for (int i=0; i<plotsToDrawWoPrefix.size(); i++){
    prefixMap[plotsToDrawWoPrefix[i]] = "";
  }
  
  SetAtlasStyle();
//   gStyle->SetHistTopMargin(0.);
  TCanvas* tmpCan = new TCanvas("c","c",3508*150/300.,2480*150/300.);
  
  
  for (int i=0; i<plotsToDraw.size(); i++){
      
    THStack *hs = new THStack("hs","Stacked 1D histograms");
    
    TH1D *h2 = NULL;
    TH1D* testHist = NULL;
    TH1D* dataHist = NULL;
    TH1D* multijetHist = NULL;
    
    for (int k=0; k<samplesToDraw.size(); k++){
      if (samplesToDraw[k]=="data"){
        dataHist = mergedSample->GetMergedDataHist(prefixMap[plotsToDraw[i]]+plotsToDraw[i]);
        dataHist->SetName(("mc_"+samplesToDraw[k]+"_" + histFolderName).c_str());
        if(!outFile->cd(plotsToDraw[i].c_str())){
          TDirectory *rdir = outFile->mkdir(plotsToDraw[i].c_str());
          rdir->cd();
        }
        dataHist->Write();
        continue;
      }
      if (samplesToDraw[k]=="multijet"){
        std::string myHistName = prefixMap[plotsToDraw[i]]+plotsToDraw[i];
        std::size_t found = myHistName.find("/hObjDump_");
        myHistName.insert(found,"_QCD");
        testHist = 
//         mergedSample->GetMergedMultijetHist(prefixMap[plotsToDraw[i]]+plotsToDraw[i]);
        mergedSample->GetMergedDataHist(myHistName);
      }
      else
        testHist = mergedSample->GetMergedHist(samplesToDraw[k],prefixMap[plotsToDraw[i]]+plotsToDraw[i]);
      testHist->SetName(("mc_"+samplesToDraw[k]+"_" + histFolderName).c_str());
      
      /// WARNING debug
//       cout << "Print out binflow for hist: " << testHist->GetName() << endl;
//       for (unsigned int iBin=1; iBin<=testHist->GetNbinsX(); iBin++){
//         cout << testHist->GetBinContent(iBin) << endl;
//       }
      
      if (testHist!=NULL){
        setHistStyle(testHist,colorArr[k]);
        hs->Add(testHist);
      }
      else{
        cout << "Plot *" << plotsToDraw[i] << "* is empty for sample " << samplesToDraw[k] 
        << endl; 
      }
      
      if (h2==NULL){
        h2 = (TH1D*)testHist->Clone("h2"); 
      }
      else{
        h2->Add(testHist);
      }
      
      if(!outFile->cd(plotsToDraw[i].c_str())){
        TDirectory *rdir = outFile->mkdir(plotsToDraw[i].c_str());
        rdir->cd();
      }
      testHist->Write();
//       cout << endl;
//       cout << "[DEBUG]\t" << plotsToDraw[i] << " " <<  samplesToDraw[k] << endl;
//       testHist->Print("all");
      
    }
    
    outFile->Write();
    
    TPad *pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1.0);
    pad1->SetBottomMargin(0); /// Upper and lower plot are joined
//     pad1->SetGridx();         /// Vertical grid
    pad1->Draw();             /// Draw the upper pad: pad1
    pad1->cd();               /// pad1 becomes the current pad
    
    std::size_t found = histFolderName.find("noMET_mT_cuts");
    
    hs->Draw("HIST");
    if (dataHist!=NULL)
      dataHist->Draw("Esame");
    if (i<3){
      hs->SetMinimum(10E-7);
      hs->SetMaximum(10E5);
      if (found==std::string::npos){
        hs->GetXaxis()->SetRangeUser(50,2000);
        dataHist->GetXaxis()->SetRangeUser(50,2000);
      }
    }
    else{
      hs->SetMinimum(10E-1);
      hs->SetMaximum(8*10E2);
      
      if (found!=std::string::npos)
        hs->SetMaximum(5*10E3);
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
   
    /// lower plot will be in pad
    tmpCan->cd();          /// Go back to the main canvas before defining pad2
    TPad *pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.3);
    pad2->SetTopMargin(0);
    pad2->SetBottomMargin(0.2);
//     pad2->SetGridx(); /// vertical grid
    pad2->SetGridy(); /// vertical grid
    pad2->Draw();
    pad2->cd();       /// pad2 becomes the current pad

    gPad->Update();
    if (i>=3){
      gPad->SetLogx(0);
    }
    else{
      gPad->SetLogx();
    }
   
    /// Define the ratio plot
    TH1D *h3 = (TH1D*)dataHist->Clone("h3");
    h3->SetLineColor(kBlack);
    h3->SetMinimum(0.8);  /// Define Y ..
    h3->SetMaximum(1.35); /// .. range
    h3->Sumw2();
    h3->SetStats(0);      /// No statistics on lower plot
    h3->Divide(h2);
    h3->SetMarkerStyle(21);
    h3->Draw("ep");       /// Draw the ratio plot

    /// Ratio plot (h3) settings
    h3->SetTitle(""); /// Remove the ratio title

    /// Y axis ratio plot settings
    h3->GetYaxis()->SetTitle("Data/Bkg");
    h3->GetYaxis()->SetNdivisions(505);
    h3->GetYaxis()->SetTitleSize(30);
    h3->GetYaxis()->SetTitleFont(43);
    h3->GetYaxis()->SetTitleOffset(1.55);
    h3->GetYaxis()->SetLabelFont(43); /// Absolute font size in pixel (precision 3)
    h3->GetYaxis()->SetLabelSize(25);
    h3->GetYaxis()->SetRangeUser(0.61,1.49);

    /// X axis ratio plot settings
    h3->GetXaxis()->SetTitleSize(30);
    h3->GetXaxis()->SetTitleFont(43);
    h3->GetXaxis()->SetTitleOffset(4.);
    h3->GetXaxis()->SetLabelFont(43); /// Absolute font size in pixel (precision 3)
    h3->GetXaxis()->SetLabelSize(25);
    
    tmpCan->cd();          /// Go back to the main canvas
    string outFileName = "pictures/" + histFolderName + ".ps";
    if (i==0&&plotsToDraw.size()>1)
        outFileName += "(";
    if (i==(plotsToDraw.size()-1)&&plotsToDraw.size()>1)
        outFileName += ")";
    tmpCan->SaveAs(outFileName.c_str());
    outFileName = "pictures/" + histFolderName + "_" + plotsToDraw[i] + ".ps"; 
    tmpCan->SaveAs(outFileName.c_str());
    outFileName = "pictures/" + histFolderName + "_" + plotsToDraw[i] + ".png"; 
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

TFile* createOutFile(string listName){
  vector<string> tmpStrVec = GetWords(listName,'/');
  string tmpStr = tmpStrVec[tmpStrVec.size()-1];
  tmpStrVec = GetWords(tmpStr,'.');
  (tmpStrVec.size()>1) ? tmpStr = tmpStrVec[tmpStrVec.size()-2] : tmpStrVec[0];
  return new TFile((tmpStr+".root").c_str(),"RECREATE");
}

vector<string> getSamplesToDraw(WprimeMergedSample *mergedSample){
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
      }
    }
    cout << endl;
  }
  else
    samplesToDraw = {"diboson","z","top","w","data"};
  return samplesToDraw;
}



// std::vector<std::string> split(const std::string &s, char delim) {
//     std::vector<std::string> elems;
//     split(s, delim, elems);
//     return elems;
// }
