#ifndef MyAnalysis_MyxAODAnalysis_H
#define MyAnalysis_MyxAODAnalysis_H

// EventLoop
#include <EventLoop/Algorithm.h>

// Infrastructure includes:
#include "xAODRootAccess/Init.h"
#include "xAODRootAccess/TEvent.h"

#include <TH1.h>
#include <TTree.h>

// GRL
class GoodRunsListSelectionTool;
class xAODJet;
class xAODMuon;
class JetCleaningTool;
class JERTool;

class MyxAODAnalysis : public EL::Algorithm
{
  // put your configuration variables here as public variables.
  // that way they can be set directly from CINT and python.
public:
  // float cutValue;


  // variables that don't get filled at submission time should be
  // protected from being send from the submission node to the worker
  // node (done by the //!)
public:
  // Tree *myTree; //!
  // TH1 *myHist; //!
  xAOD::TEvent *m_event;  //!

  int m_eventCounter; //!

  int m_numCleanEvents; //!

  TH1 *h_jetPt; //!

  // defining the output file name and tree that we will put in the output ntuple, also the one branch that will be in that tree 
  std::string outputName;
  TTree *tree; //!
  int EventNumber; //!

  #ifndef __CINT__
    GoodRunsListSelectionTool *m_grl; //!
    JetCleaningTool *m_jetCleaning; //!  
    JERTool *m_JERTool; //!
  #endif // not __CINT__

  // this is a standard constructor
  MyxAODAnalysis ();

  // these are the functions inherited from Algorithm
  virtual EL::StatusCode setupJob (EL::Job& job);
  virtual EL::StatusCode fileExecute ();
  virtual EL::StatusCode histInitialize ();
  virtual EL::StatusCode changeInput (bool firstFile);
  virtual EL::StatusCode initialize ();
  virtual EL::StatusCode execute ();
  virtual EL::StatusCode postExecute ();
  virtual EL::StatusCode finalize ();
  virtual EL::StatusCode histFinalize ();

  // this is needed to distribute the algorithm to the workers
  ClassDef(MyxAODAnalysis, 1);
};

#endif
