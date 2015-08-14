#ifndef MyAnalysis_MyxAODAnalysis_H
#define MyAnalysis_MyxAODAnalysis_H

/// EventLoop
#include <EventLoop/Algorithm.h>

/// Infrastructure includes:
#include "xAODRootAccess/Init.h"
#include "xAODRootAccess/TEvent.h"

#include <TH1.h>
#include <TH2.h>
#include <TTree.h>

/// shallow copy
#include "xAODCore/ShallowCopy.h"

#include <MyAnalysis/HistObjectDumper.h>
#include <MyAnalysis/BitsetCutflow.h>
#include <METUtilities/METMaker.h>

/// include files for using the trigger tools
#include "TrigConfxAOD/xAODConfigTool.h"
#include "TrigDecisionTool/TrigDecisionTool.h"
#include "JetMomentTools/JetVertexTaggerTool.h"

#include "ElectronPhotonSelectorTools/AsgElectronLikelihoodTool.h"
#include "xAODEgamma/ElectronContainer.h"
#include <ElectronPhotonFourMomentumCorrection/EgammaCalibrationAndSmearingTool.h>

/// IsolationTool
#include "ElectronIsolationSelection/IsolationSelectionTool.h"

#include "ElectronPhotonSelectorTools/AsgElectronLikelihoodTool.h"

/// GRL
class GoodRunsListSelectionTool;
class xAODJet;
class xAODMuon;
class xAODElectron;
class JetCleaningTool;
class JetCalibrationTool;
class JERTool;

class METUtility;

/// muon calibration and smearing tool
namespace CP{
  class MuonSelectionTool;
  class MuonCalibrationAndSmearingTool; /// this tool lives in the namespace CP
  class IsolationSelectionTool;
  //~ class MuonEfficiencyCorrections;
}

// namespace TrigConf{
//   class TrigDecisionTool;
// }
// 
// namespace Trig{
//   class xAODConfigTool;
// }

class MyxAODAnalysis : public EL::Algorithm
{
  /// put your configuration variables here as public variables.
  /// that way they can be set directly from CINT and python.
public:
  bool m_useHistObjectDumper;
  bool m_useBitsetCutflow;
  bool m_useCalibrationAndSmearingTool;
  bool m_runElectronChannel;
  bool m_doWprimeTruthMatching;
  int m_truthoption;
  int truthOption;
  bool m_runElectronChannel;
  bool isMC;
  int count[20];
  
  std::string outputName;
  TTree *tree; //!
  int EventNumber; //!
  int RunNumber; //!
  int LumiBlock; //!
  float MuonPt; //!
  float MuonEta; //!
  float MuonPhi; //!
  float MissingEt; //!
  float MissingEx; //!
  float MissingEy; //!
  float OldMissingEt; //!
  float TransverseMass; //!
  
  /// variables that don't get filled at submission time should be
  /// protected from being send from the submission node to the worker
  /// node (done by the //!)
public:
  // Tree *myTree; //!
  // TH1 *myHist; //!
  xAOD::TEvent *m_event;  //!
  xAOD::TStore *m_store;  //!

  AsgElectronLikelihoodTool* m_LHToolTight2015; //!
  AsgElectronLikelihoodTool* m_LHToolMedium2015; //!
  
  int m_eventCounter; //!

  int m_numCleanEvents; //!
  int m_counter1;
  int m_counter2;

  TH1 *h_ptmuon; //!
  
  TH1 *h_jetPt_LooseBadJets; //!
  TH1 *h_jetPt_VeryLooseBadJets; //!

  TH1 *h_nSelectedMuons; //!

  TH1 *h_MET_RefFinalFix; //!
  TH1 *h_MET_RefJet; //!
  TH1 *h_MET_RefCore; //!
  TH1 *h_MET_MyJet; //!
  TH1 *h_MET_RefTau; //!
  TH1 *h_MET_RefMuon; //!
  TH1 *h_MET_RefGamma; //!
  TH1 *h_MET_RefElectron; //!
  TH1 *h_MET_RefSoft; //!
  TH1 *h_MET_RefTrack; //!
  TH1 *h_MET_RefFinalNew; //!  
  TH1 *h_MET_RefFinalFix_test; //!

  TH1 *h_MET_OldJet; //!
  TH1 *h_MET_OldTau; //!
  TH1 *h_MET_OldMuon; //!
  TH1 *h_MET_OldGamma; //!
  TH1 *h_MET_OldElectron; //!
  TH1 *h_MET_OldSoft; //!
  TH1 *h_MET_OldTrack; //!
  
  TH1 *h_Mt; //!
  TH1 *h_Mt_muonPtCut; //!
  
  TH1 *h_zPrimVtx; //!
  
  TH1 *h_truthMET; //!

  /// defining the output file name and tree that we will put in the output ntuple, also the one branch that will be in that tree 
//   std::string outputName;
//   TTree *tree; //!
  int EventNumber; //!

  #ifndef __CINT__
    GoodRunsListSelectionTool *m_grl; //!
    JetCleaningTool *m_jetCleaning; //!  
    JetCalibrationTool *m_jetCalibrationTool; //!
    JERTool *m_JERTool; //!
    JetVertexTaggerTool* m_jvtTool; //!
    CP::EgammaCalibrationAndSmearingTool *m_eleCalibrationTool; //! 

    met::METMaker * m_metMaker; //!
    Trig::TrigDecisionTool *m_trigDecisionTool; //!
    TrigConf::xAODConfigTool *m_trigConfigTool; //!
    /// Muson Selector Tool
    CP::MuonSelectionTool* m_muonSelection; //!
    CP::MuonSelectionTool* m_loosemuonSelection; //!
    /// MuonCalibrationAndSmearing
    CP::MuonCalibrationAndSmearingTool *m_muonCalibrationAndSmearingTool; //!
//     CP::IsolationSelectionTool *m_isolationSelectionTool; //!
    //~ CP::MuonEfficiencyScaleFactors *m_effi_corr; //!
    //METUtility *m_METUtil;
    
    CP::IsolationSelectionTool* iso_1; //!
    CP::IsolationSelectionTool* iso_2; //!
    CP::IsolationSelectionTool* iso_3; //!
    CP::IsolationSelectionTool* iso_4; //! 
    CP::IsolationSelectionTool* iso_5; //!
    CP::IsolationSelectionTool* iso_6; //!
    CP::IsolationSelectionTool* iso_7; //!
    CP::IsolationSelectionTool* iso_8; //!
    CP::IsolationSelectionTool* iso_9; //!
    
    Trig::TrigDecisionTool *m_trigDecisionTool; //!
    TrigConf::xAODConfigTool *m_trigConfigTool; //!
  #endif /// not __CINT__

  HistObjectDumper *m_HistObjectDumper; //!
  BitsetCutflow* m_BitsetCutflow; //!

  /// this is a standard constructor
  MyxAODAnalysis ();

  /// these are the functions inherited from Algorithm
  virtual EL::StatusCode setupJob (EL::Job& job);
  virtual EL::StatusCode fileExecute ();
  virtual EL::StatusCode histInitialize ();
  virtual EL::StatusCode changeInput (bool firstFile);
  virtual EL::StatusCode initialize ();
  virtual EL::StatusCode execute ();
  virtual EL::StatusCode postExecute ();
  virtual EL::StatusCode finalize ();
  virtual EL::StatusCode histFinalize ();

  /// Custom made functions
//   xAOD::Muon* SelectMuon(const xAOD::MuonContainer* muons, 
//                          xAOD::Vertex* primVertex, 
//                          bool lookForVetoMuon = false);
  
  ConstDataVector<xAOD::MuonContainer>  SelectMuon(const xAOD::MuonContainer* muons,
                         xAOD::Vertex* primVertex,
                         bool lookForVetoMuon = false);
  
  ConstDataVector<xAOD::MuonContainer>  VetoMuon(const xAOD::MuonContainer* muons,
                         xAOD::Vertex* primVertex,
                         bool lookForVetoMuon = false);
  
  ConstDataVector<xAOD::ElectronContainer>  SelectElectron(const xAOD::ElectronContainer* electrons, xAOD::Vertex* primVertex, bool lookForVetoElectron);
  
  bool passMuonSelection(const xAOD::Muon* mu,
               const xAOD::EventInfo* eventInfo,
               xAOD::Vertex* primVertex,
               bool lookForVetoMuon);
  
//   xAOD::Electron* SelectElectron(const xAOD::ElectronContainer* electrons, 
//                                          xAOD::Vertex* primVertex, 
//                                          bool lookForVetoElectron);
  
  /// this is needed to distribute the algorithm to the workers
  ClassDef(MyxAODAnalysis, 1);
};

#endif
