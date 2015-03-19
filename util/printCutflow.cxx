#include "xAODRootAccess/Init.h"
#include "SampleHandler/SampleHandler.h"
#include "SampleHandler/Sample.h"
#include "SampleHandler/ToolsDiscovery.h"
#include "EventLoop/Job.h"
#include "EventLoop/ProofDriver.h"
#include "EventLoop/DirectDriver.h"
#include "SampleHandler/DiskListLocal.h"
#include <TSystem.h>

#include "MyAnalysis/MyxAODAnalysis.h"
#include <EventLoopAlgs/NTupleSvc.h>
#include <EventLoop/OutputStream.h>

/// boost
#include "boost/program_options.hpp"
#include <boost/algorithm/string.hpp>

namespace 
{ 
  const size_t ERROR_IN_COMMAND_LINE = 1; 
  const size_t SUCCESS = 0; 
  const size_t ERROR_UNHANDLED_EXCEPTION = 2; 
  const size_t HELP_CALLED = 3;
  
  namespace po = boost::program_options;
 
} /// namespace 

int parseOptionsWithBoost(po::variables_map &vm, int argc, char* argv[]);

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

	/// get histogram
	SH::Sample* mySample = sh.get ("mc14_13TeV.203671.Pythia8_AU2MSTW2008LO_Wprime_emutau_2000.merge.DAOD_EXOT9.e3148_s1982_s2008_r5787_r5853_p1846_tid04963913_00");
	//SH::Sample* mySample = sh.findBySource ("mc14_13TeV*Wprime_emutau_*");
	TH1I* cutflowHist = (TH1I*)mySample->readHist ("cutflow_hist");

	cout << "cutflowHist pointer: " << cutflowHist << endl;
	
	for (int i=1; i<=cutflowHist->GetNbinsX(); i++){
		int binContent = cutflowHist->GetBinContent(i);
		if (binContent<=0) break;
		string binLabel = cutflowHist->GetXaxis()->GetBinLabel(i);
		cout << binLabel << ":\t" << binContent << endl;
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
      ("folder,f", po::value<string>()->required(), "name of folder to read")
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

