#include "xAODRootAccess/Init.h"
#include "SampleHandler/SampleHandler.h"
#include "SampleHandler/ToolsDiscovery.h"
#include "EventLoop/Job.h"
#include "EventLoop/ProofDriver.h"
#include "SampleHandler/DiskListLocal.h"
#include <TSystem.h>


#include "MyAnalysis/MyxAODAnalysis.h"
#include <EventLoopAlgs/NTupleSvc.h>
#include <EventLoop/OutputStream.h>

int main( int argc, char* argv[] ) {

   // Take the submit directory from the input if provided:
   std::string submitDir = "submitDir";
   if( argc > 1 ) submitDir = argv[ 1 ];

   // Set up the job for xAOD access:
   xAOD::Init().ignore();

   // Construct the samples to run on:
   SH::SampleHandler sh;

   const char* inputFilePath = gSystem->ExpandPathName ("/afs/cern.ch/work/o/oviazlo/Wprime/datasets"); // using $ALRB_TutorialData previously defined
   SH::DiskListLocal list (inputFilePath);
   //SH::scanDir (sh, list, "DAOD_EXOT9.04655450._000002.pool.root.1"); // specifying one particular file for testing
   SH::scanDir (sh, list, "DAOD_EXOT9.*root*");

   // Set the name of the input TTree. It's always "CollectionTree"
   // for xAOD files.
   sh.setMetaString( "nc_tree", "CollectionTree" );

   // Print what we found:
   sh.print();

   // Create an EventLoop job:
   EL::Job job;
   job.sampleHandler( sh );

   // Specify that we only want to run on 1k events
   //job.options()->setDouble(EL::Job::optMaxEvents, 1000);

   // define an output and an ntuple associated to that output 
   EL::OutputStream output  ("myOutput");
   job.outputAdd (output);
   EL::NTupleSvc *ntuple = new EL::NTupleSvc ("myOutput");
   job.algsAdd (ntuple);

   // Add our analysis to the job:
   MyxAODAnalysis* alg = new MyxAODAnalysis();
   job.algsAdd( alg );

   alg->outputName = "myOutput"; // give the name of the output to our algorithm

   // Run the job using the local/direct driver:
   EL::ProofDriver driver;
   driver.numWorkers = 4;
   //EL::DirectDriver driver;
   driver.submit( job, submitDir );

   return 0;
}
