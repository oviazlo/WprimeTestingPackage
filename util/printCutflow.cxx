/* Copyright 2016 Oleksandr Viazlo */

///*****************************************************************************
///
/// Print cutflow. Used for cutflow comparison.
///
///*****************************************************************************

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
#include <fstream>

/// boost
#include "boost/program_options.hpp"
#include <boost/algorithm/string.hpp>

/// private headers
#include "MyAnalysis/RecoAnalysis.h"
#include "HelpFunctions.h"

map<string,string> sampleMap;
bool MatchPattern(string sampleName, vector<string> patternVec);

int main( int argc, char* argv[] ) {

  po::options_description desc("Options"); 
  desc.add_options()
    ("help,h", "Print help messages") 
    ("folder,f", po::value<string>(), "name of folder to read")
    ("samplePattern,p", po::value<string>(),"specify Sample Pattern")
    ("useOnlyOneFile", "use only one file from folder")
    ("writeCutflowToFile,w", "save cutflow to file")
    ;
  
  /// get global input arguments:
  po::variables_map vm; 
  const size_t returnedMessage = parseOptionsWithBoost(vm,desc,argc,argv);
  if (returnedMessage!=SUCCESS) std::exit(returnedMessage);

  /// Take the submit directory from the input if provided:
  std::string folderPrefix = "submitDirs";
  
  vector <string> strSamplePattertVec;
  if ( vm.count("samplePattern") ){
    string strSamplePattert = vm["samplePattern"].as<std::string>();
    strSamplePattertVec = GetWords(strSamplePattert,'*');
  }
  
  std::vector<std::string> folders;
  if (vm.count("folder")){
    folders = GetWords(vm["folder"].as<std::string>(),',');
  }
  else{
    cout << "[ERROR]\tplease specify input folder with '-f' option" << endl;
    return 0;    
  }

  /// Construct the samples to run on:
  SH::SampleHandler sh;

  for (int i=0; i<folders.size(); i++){
    string folder = folderPrefix + "/" + folders[i];
    cout << "[INFO]\tRead samples from dir: " << folder << endl;
    sh.load (folder + "/hist");
  }
 
  if (sh.size()==0){
    cout << "sh.size() = " << sh.size() << endl;
    cout << "Terminate!!!"<< endl;
    return 0;
  }
  
  SH::Sample* mySample = *(sh.begin());
  TH1I* cutflowHist = NULL;
  cout << "[INFO]\tsh.size() = " << sh.size() << endl;

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
    cout << "[INFO]\tRead following samples:" << endl;
    for (SH::SampleHandler::iterator iter = sh.begin(); iter != sh.end(); 
      ++ iter){
      if (!MatchPattern((*iter)->name(),strSamplePattertVec)) 
        continue;
      cout << (*iter)->name() << endl;
      TH1I* tmpHist = (TH1I*)(*iter)->readHist ("cutflow_hist");
      if (cutflowHist==NULL)
        cutflowHist = tmpHist;
      else
        cutflowHist->Add(tmpHist);
    }
    cout << endl;
  }
    
  if (cutflowHist==NULL){
    cout << "[ERROR]\tno matched pattern found! Terminate!!!" << endl;
    return 0;
  }
    
  
  std::ofstream logFile;
  if (vm.count("writeCutflowToFile")){
    string fileName = "cutflow_" + vm["samplePattern"].as<std::string>() + "_" + vm["folder"].as<std::string>() + ".txt";
    std::ofstream ofs(fileName.c_str(), std::ofstream::out);
    ofs.close();
    logFile.open (fileName.c_str(), std::ofstream::out | std::ofstream::app);
  }
    
  for (int i=1; i<=cutflowHist->GetNbinsX(); i++){
    int binContent = cutflowHist->GetBinContent(i);
    if (binContent<=0) break;
    string binLabel = cutflowHist->GetXaxis()->GetBinLabel(i);
    cout << binLabel << ":\t" << binContent << endl;
    if (vm.count("writeCutflowToFile"))
      logFile << binLabel << ":\t" << binContent << endl;
  }
  
  if (vm.count("writeCutflowToFile"))
    logFile.close();
  
  return 0;
}

bool MatchPattern(string sampleName, vector<string> patternVec){
  for (int i=0; i<patternVec.size(); i++){
    std::size_t found = sampleName.find(patternVec[i]);
    if (found==std::string::npos)
      return false;
  }
  return true;
}
