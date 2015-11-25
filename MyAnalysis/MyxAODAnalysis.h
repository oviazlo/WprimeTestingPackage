#ifndef MyAnalysis_MyxAODAnalysis_H
#define MyAnalysis_MyxAODAnalysis_H

#include <EventLoop/Job.h>
#include <EventLoop/StatusCode.h>
#include <EventLoop/Worker.h>
//#include <MyAnalysis/HistObjectDumper.h>
#include "CPAnalysisExamples/errorcheck.h"

/// EDM includes:
#include "xAODEventInfo/EventInfo.h"
#include "xAODTruth/TruthParticle.h"
#include "xAODTruth/TruthVertex.h"
#include "xAODTruth/TruthEventContainer.h"
#include "xAODCore/ShallowCopy.h"
#include "xAODBase/IParticleHelpers.h"

/// Isolation
#include "IsolationSelection/IsolationSelectionTool.h"
#include "JetMomentTools/JetVertexTaggerTool.h"

/// GRL
#include "GoodRunsLists/GoodRunsListSelectionTool.h"

#include "xAODJet/JetContainer.h"
#include "xAODTau/TauJetContainer.h"
#include "xAODTau/TauJetAuxContainer.h"
#include "xAODMuon/MuonContainer.h"
#include "xAODEgamma/ElectronContainer.h"
#include "xAODEgamma/PhotonContainer.h"
#include "xAODEgamma/EgammaContainer.h"
#include "JetSelectorTools/JetCleaningTool.h"
#include "JetResolution/JERTool.h"
#include "JetCalibTools/JetCalibrationTool.h"
#include <TSystem.h> /// used to define JERTool calibration path 
#include "xAODTracking/VertexContainer.h"
#include "xAODParticleEvent/IParticleLink.h"

#include "xAODTracking/TrackParticlexAODHelpers.h"

/// Muons
#include "MuonSelectorTools/MuonSelectionTool.h"
#include "MuonMomentumCorrections/MuonCalibrationAndSmearingTool.h"

/// to check the return correction code status of tools
#include "PATInterfaces/CorrectionCode.h" 

/// MET
#include "METUtilities/METMaker.h" /// METUtils
#include "xAODMissingET/MissingETContainer.h"
#include "xAODMissingET/MissingETAuxContainer.h"
#include "xAODMissingET/MissingETComposition.h"
#include "xAODMissingET/MissingETAssociationMap.h"
#include "METInterface/IMETMaker.h"
#include "METUtilities/CutsMETMaker.h"

#include <TFile.h>
#include <TMath.h>
#include <TLorentzVector.h>

/// OutputStream
#include "EventLoop/OutputStream.h"

#include "xAODRootAccess/tools/Message.h"

/// Helper macro for checking xAOD::TReturnCode return values
#define EL_RETURN_CHECK( CONTEXT, EXP )                     \
   do {                                                     \
      if( ! EXP.isSuccess() ) {                             \
         Error( CONTEXT,                                    \
                XAOD_MESSAGE( "Failed to execute: %s" ),    \
                #EXP );                                     \
         return EL::StatusCode::FAILURE;                    \
      }                                                     \
   } while( false )



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
#include <MyAnalysis/WprimeHist.h>
#include <MyAnalysis/WprimeSample.h>

#include <METUtilities/METMaker.h>

/// include files for using the trigger tools
#include "TrigConfxAOD/xAODConfigTool.h"
#include "TrigDecisionTool/TrigDecisionTool.h"
#include "JetMomentTools/JetVertexTaggerTool.h"

#include "ElectronPhotonSelectorTools/AsgElectronLikelihoodTool.h"
#include "xAODEgamma/ElectronContainer.h"
#include <ElectronPhotonFourMomentumCorrection/EgammaCalibrationAndSmearingTool.h>

/// IsolationTool
// #include "ElectronIsolationSelection/IsolationSelectionTool.h"

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
  bool m_doNotApplyTriggerCuts;
  int m_truthoption;
  int truthOption;
   
  int count[20]; //!
  bool isMC; //!
   
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

  xAOD::TEvent *m_event;  //!
  xAOD::TStore *m_store;  //!

  AsgElectronLikelihoodTool* m_LHToolTight2015; //!
  AsgElectronLikelihoodTool* m_LHToolMedium2015; //!
  
  int m_eventCounter; //!

  int m_numCleanEvents; //!
  int m_counter1;
  int m_counter2;

  ///
  /// h_<object>_<variable>_<specialCutSet>
  ///
  /// possible objects:
  /// - event
  /// - lepton
  /// - vetoElectron
  /// - vetoMuon
  /// - MET
  /// - leptonMET
  /// - 
  
  /// event hists
  TH1F *h_event_nJets; //!
  TH1F *h_event_nJets_Pt40; //!

  /// lepton hists
  TH1F *h_lepton_d0Significance; //!
  TH1F *h_lepton_z0SinTheta; //!
  TH1F *h_lepton_Pt; //!
  TH1F *h_lepton_eta; //!
  TH1F *h_lepton_phi; //!

  /// leptonMET hists
  TH1F *h_leptonMET_deltaPhi; //!
  TH1F *h_leptonMET_Mt; //!  
  
  /// MET hists
  TH1F *h_MET_Et; //!
  TH1F *h_MET_phi; //!

  /// Magnar histogramming way
  TH1F* hMu_pt_off; //!
  TH1F* hMu_mt_off; //!
  TH1F* hMu_MET_Muons_off; //!
  
  /// defining the output file name and tree that we will put in the output 
  /// ntuple, also the one branch that will be in that tree 
  std::string outputName; //!
//   TTree *tree; //!

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
    CP::MuonSelectionTool* m_looseMuonSelection; //!
    /// MuonCalibrationAndSmearing
    CP::MuonCalibrationAndSmearingTool *m_muonCalibrationAndSmearingTool; //!

    CP::IsolationSelectionTool *m_muonisolationSelectionTool; //!
    CP::IsolationSelectionTool *m_eleisolationSelectionTool; //!
    
    const xAOD::EventInfo* eventInfo = 0; //!
    
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
std::pair<unsigned int, unsigned int> SelectMuons(
                                  xAOD::MuonContainer* muons,
                                  xAOD::Vertex* primVertex,
                                  bool fillInCutflow);

std::pair<unsigned int, unsigned int> SelectElectrons(
                                  xAOD::ElectronContainer* electrons,
                                  bool fillInCutflow);
  
  /// this is needed to distribute the algorithm to the workers
  ClassDef(MyxAODAnalysis, 1);
};

#endif
