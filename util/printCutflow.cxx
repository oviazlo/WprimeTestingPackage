#include "xAODRootAccess/Init.h"
#include "SampleHandler/SampleHandler.h"
#include "SampleHandler/Sample.h"
#include "SampleHandler/ToolsDiscovery.h"
#include "EventLoop/Job.h"
#include "EventLoop/ProofDriver.h"
#include "EventLoop/DirectDriver.h"
#include "SampleHandler/DiskListLocal.h"
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

map<string,string> sampleMap;

namespace 
{ 
  const size_t ERROR_IN_COMMAND_LINE = 1; 
  const size_t SUCCESS = 0; 
  const size_t ERROR_UNHANDLED_EXCEPTION = 2; 
  const size_t HELP_CALLED = 3;
  
  namespace po = boost::program_options;
 
} /// namespace 

int parseOptionsWithBoost(po::variables_map &vm, int argc, char* argv[]);
TStyle* AtlasStyle();
void SetAtlasStyle();

int main( int argc, char* argv[] ) {

  /// get global input arguments:
  po::variables_map vm; 
  const size_t returnedMessage = parseOptionsWithBoost(vm,argc,argv);
  if (returnedMessage!=SUCCESS) std::exit(returnedMessage);

  /// Take the submit directory from the input if provided:
  std::string folder = vm["folder"].as<std::string>();  

  // Construct the samples to run on:
  SH::SampleHandler sh;

  sh.load ((folder + "/hist").c_str());

  /// Print what we found:
  sh.print();

  if (sh.size()!=1){
    cout << "Found more(less) than one sample in folder. Abort!" << endl;
    return -1;
  }
  
  SH::Sample* mySample = sh.at(0);
  
  /// get histogram
//   SH::Sample* mySample = sh.get (sampleMap[ vm["sample"].as<std::string>() ]);
  
  //SH::Sample* mySample = sh.get ("mc14_13TeV.203671.Pythia8_AU2MSTW2008LO_Wprime_emutau_2000.merge.DAOD_EXOT9.e3148_s1982_s2008_r5787_r5853_p1846_tid04963913_00");
  //SH::Sample* mySample = sh.get ("mc14_13TeV.158762.Pythia8_AU2MSTW2008LO_Wprime_emutau_3000.merge.DAOD_EXOT9.e3148_s1982_s2008_r5787_r5853_p1816_tid04655450_00");
  //SH::Sample* mySample = sh.findBySource ("mc14_13TeV*Wprime_emutau_*");
  TH1I* cutflowHist = (TH1I*)mySample->readHist ("cutflow_hist");

  //~ cout << "cutflowHist pointer: " << cutflowHist << endl;
  
  for (int i=1; i<=cutflowHist->GetNbinsX(); i++){
    int binContent = cutflowHist->GetBinContent(i);
    if (binContent<=0) break;
    string binLabel = cutflowHist->GetXaxis()->GetBinLabel(i);
    cout << binLabel << ":\t" << binContent << endl;
  }
  

  SetAtlasStyle();
  
  // read pt, MET, Mt histograms, make them pretty and save them
  TH1F* h_pt = (TH1F*)mySample->readHist ("muon/stage2_allCuts/pt");
  TH1F* h_MET = (TH1F*)mySample->readHist ("h_MET_RefFinalFix");
  TH1F* h_Mt = (TH1F*)mySample->readHist ("h_Mt_muonPtCut");
  
  h_pt->GetXaxis()->SetTitle("p_{T}^{#mu} [GeV]");
  h_MET->GetXaxis()->SetTitle("E_{T}^{miss} [GeV]");
  h_Mt->GetXaxis()->SetTitle("m_{T} [GeV]");
  
  h_pt->GetYaxis()->SetTitle("Counts");
  h_MET->GetYaxis()->SetTitle("Counts");
  h_Mt->GetYaxis()->SetTitle("Counts");
  
  h_pt->SetTitle("");
  h_MET->SetTitle("");
  h_Mt->SetTitle("");
  
  TCanvas *can = new TCanvas();
  gPad->SetLogy();
  
  h_pt->Draw();
  h_pt->GetXaxis()->SetRangeUser(0,5000.0);
  can->SaveAs("pT.eps");
  
  h_MET->GetXaxis()->SetRangeUser(0,5000.0);
  h_MET->Draw();
  can->SaveAs("MET.eps");
  
  h_Mt->GetXaxis()->SetRangeUser(0,5000.0);
  h_Mt->Draw();
  can->SaveAs("mT.eps");
  
  return 0;
}



/// parse input arguments
int parseOptionsWithBoost(po::variables_map &vm, int argc, char* argv[]){
  sampleMap["W2"] = "mc14_13TeV.203671.Pythia8_AU2MSTW2008LO_Wprime_emutau_2000.merge.DAOD_EXOT9.e3148_s1982_s2008_r6092_p1872_tid05283006_00";
//   sampleMap["W2"] = "mc14_13TeV.203671.Pythia8_AU2MSTW2008LO_Wprime_emutau_2000.merge.DAOD_EXOT9.e3148_s1982_s2008_r5787_r5853_p1846_tid04963913_00";
  sampleMap["W3"] = "mc14_13TeV.158762.Pythia8_AU2MSTW2008LO_Wprime_emutau_3000.merge.DAOD_EXOT9.e3148_s1982_s2008_r5787_r5853_p1816_tid04655450_00";
  sampleMap["Wplus"] = "mc15_13TeV.301112.PowhegPythia8EvtGen_AZNLOCTEQ6L1_"
  "Wplusmunu_2500M2750.merge.DAOD_EXOT9.e3663_s2608_s2183_r6630_r6264_p2353_tid0"
  "5541743_00";
  
  string availableSamples = "";
  
  typedef map<string,string>::iterator it_type;
  for(it_type iterator = sampleMap.begin(); iterator != sampleMap.end(); iterator++) {
    availableSamples += iterator->first + ", ";
  }
  availableSamples.resize (availableSamples.size()-2);
  
  try 
  { 
    /** Define and parse the program options 
     */ 
    po::options_description desc("Options"); 
    desc.add_options()
      ("help,h", "Print help messages") 
      ("folder,f", po::value<string>()->required(), "name of folder to read")
      ("sample,s", po::value<string>()->required(), ("sample type. Available samples:\n" + availableSamples).c_str())
      ;
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



void SetAtlasStyle()
{
  std::cout << "\nApplying ATLAS style settings...\n" << std::endl ;
  TStyle* atlasStyle = AtlasStyle();
  gROOT->SetStyle("ATLAS");
  gROOT->ForceStyle();
}

TStyle* AtlasStyle() 
{
  TStyle *atlasStyle = new TStyle("ATLAS","Atlas style");

  // use plain black on white colors
  Int_t icol=0; // WHITE
  atlasStyle->SetFrameBorderMode(icol);
  atlasStyle->SetFrameFillColor(icol);
  atlasStyle->SetCanvasBorderMode(icol);
  atlasStyle->SetCanvasColor(icol);
  atlasStyle->SetPadBorderMode(icol);
  atlasStyle->SetPadColor(icol);
  atlasStyle->SetStatColor(icol);
  //atlasStyle->SetFillColor(icol); // don't use: white fill color for *all* objects

  // set the paper & margin sizes
  atlasStyle->SetPaperSize(20,26);

  // set margin sizes
  atlasStyle->SetPadTopMargin(0.08);
  atlasStyle->SetPadRightMargin(0.05);
  atlasStyle->SetPadBottomMargin(0.16);
  atlasStyle->SetPadLeftMargin(0.16);

  // set title offsets (for axis label)
  atlasStyle->SetTitleXOffset(1.4);
  atlasStyle->SetTitleYOffset(1.4);

  // use large fonts
  //Int_t font=72; // Helvetica italics
  Int_t font=42; // Helvetica
  Double_t tsize=0.05;
  atlasStyle->SetTextFont(font);

  atlasStyle->SetTextSize(tsize);
  atlasStyle->SetLabelFont(font,"x");
  atlasStyle->SetTitleFont(font,"x");
  atlasStyle->SetLabelFont(font,"y");
  atlasStyle->SetTitleFont(font,"y");
  atlasStyle->SetLabelFont(font,"z");
  atlasStyle->SetTitleFont(font,"z");
  
  atlasStyle->SetLabelSize(tsize,"x");
  atlasStyle->SetTitleSize(tsize,"x");
  atlasStyle->SetLabelSize(tsize,"y");
  atlasStyle->SetTitleSize(tsize,"y");
  atlasStyle->SetLabelSize(tsize,"z");
  atlasStyle->SetTitleSize(tsize,"z");

  // use bold lines and markers
  atlasStyle->SetMarkerStyle(20);
  atlasStyle->SetMarkerSize(1.2);
  atlasStyle->SetHistLineWidth(2.);
  atlasStyle->SetLineStyleString(2,"[12 12]"); // postscript dashes

  // get rid of X error bars 
  //atlasStyle->SetErrorX(0.001);
  // get rid of error bar caps
  atlasStyle->SetEndErrorSize(0.);

  // do not display any of the standard histogram decorations
  atlasStyle->SetOptTitle(0);
  //atlasStyle->SetOptStat(1111);
  atlasStyle->SetOptStat(0);
  atlasStyle->SetOptFit(1111);
  //atlasStyle->SetOptFit(0);

  // put tick marks on top and RHS of plots
  atlasStyle->SetPadTickX(1);
  atlasStyle->SetPadTickY(1);

  return atlasStyle;

}


