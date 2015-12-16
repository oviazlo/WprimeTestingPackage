#include "xAODRootAccess/Init.h"
#include "SampleHandler/SampleHandler.h"
#include "SampleHandler/Sample.h"
#include "SampleHandler/ToolsDiscovery.h"
#include "EventLoop/Job.h"
#include "EventLoop/ProofDriver.h"
#include "EventLoop/DirectDriver.h"
#include "SampleHandler/DiskListLocal.h"
#include <SampleHandler/ScanDir.h>
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
  std::string folder = "submitDir";
  if (vm.count("folder"))
    folder = vm["folder"].as<std::string>();
  else
    cout << "[INFO]\tread cutflow from default directory: " << folder << endl;

  // Construct the samples to run on:
  SH::SampleHandler sh;

//   sh.load ((folder + "/hist").c_str());

  SH::ScanDir()
  .samplePattern ("hist-*root")
  .scan (sh, folder);
  
  /// Print what we found:
  sh.print();

  if (sh.size()!=1){
    cout << "Found more(less) than one sample in folder. Abort!" << endl;
    return -1;
  }
  
  SH::Sample* mySample = sh.at(0);
  
  /// get histogram
  TH1I* cutflowHist = (TH1I*)mySample->readHist ("cutflow_hist");
  
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



/// parse input arguments
int parseOptionsWithBoost(po::variables_map &vm, int argc, char* argv[]){

  try 
  { 
    /** Define and parse the program options 
     */ 
    po::options_description desc("Options"); 
    desc.add_options()
      ("help,h", "Print help messages") 
      ("folder,f", po::value<string>(), "name of folder to read")
      ("drawHists,d", "draw histograms")
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


