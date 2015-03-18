#include "xAODRootAccess/Init.h"
#include "SampleHandler/SampleHandler.h"
#include "SampleHandler/ToolsDiscovery.h"
#include "EventLoop/Job.h"
#include "EventLoop/DirectDriver.h"
#include "SampleHandler/DiskListLocal.h"
#include <TSystem.h>

#include "MyAnalysis/MyxAODAnalysis.h"
#include <EventLoopAlgs/NTupleSvc.h>
#include <EventLoop/OutputStream.h>

int main( int argc, char* argv[] ) {

   // Take the submit directory from the input if provided:
   std::string submitDir = "submitDir";
   if( argc > 1 ) submitDir = argv[ 1 ];

   int nEvent = -1;
   if( argc > 2 ) nEvent = atoi(argv[ 2 ]);
   
   // Set up the job for xAOD access:
   xAOD::Init().ignore();

   // Construct the samples to run on:
   SH::SampleHandler sh;

   const char* inputFilePath = gSystem->ExpandPathName ("/afs/cern.ch/work/o/oviazlo/Wprime/datasets");
   //const char* inputFilePath = gSystem->ExpandPathName ("$ALRB_TutorialData/r5591"); // using $ALRB_TutorialData previously defined
   SH::DiskListLocal list (inputFilePath);
   SH::scanDir (sh, list, "DAOD_EXOT9*root*");
   //SH::scanDir (sh, list, "DAOD_EXOT9.04655450._000002.pool.root.1");
   //SH::scanDir (sh, list, "AOD.01494882._113691.pool.root.1"); // specifying one particular file for testing

   // Set the name of the input TTree. It's always "CollectionTree"
   // for xAOD files.
   sh.setMetaString( "nc_tree", "CollectionTree" );

   // Print what we found:
   sh.print();

   // Create an EventLoop job:
   EL::Job job;
   job.sampleHandler( sh );

   // Specify that we only want to run on 1k events
   if (nEvent!=-1)
	   job.options()->setDouble(EL::Job::optMaxEvents, nEvent);

   // define an output and an ntuple associated to that output 
   EL::OutputStream output  ("myOutput");
   job.outputAdd (output);
   EL::NTupleSvc *ntuple = new EL::NTupleSvc ("myOutput");
   job.algsAdd (ntuple);

   /// Add our analysis to the job:
   MyxAODAnalysis* alg = new MyxAODAnalysis();
   job.algsAdd( alg );
   /// all configuration to the algorithm should be passed after attachment to the job (practically after this comment)

   alg->outputName = "myOutput"; // give the name of the output to our algorithm
   //alg->m_useHistObjectDumper = false;

   // Run the job using the local/direct driver:
   EL::DirectDriver driver;
   driver.submit( job, submitDir );

   return 0;
}
