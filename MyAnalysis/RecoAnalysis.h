/* Copyright 2016 Oleksandr Viazlo */

#ifndef MyAnalysis_RecoAnalysis_H
#define MyAnalysis_RecoAnalysis_H

/// total data lumi in fb-1
#define TOTAL_LUMI 3.20905

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
#include "xAODMuon/MuonAuxContainer.h"
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
#include "MuonEfficiencyCorrections/MuonTriggerScaleFactors.h"
#include "MuonEfficiencyCorrections/MuonEfficiencyScaleFactors.h"

 
/// to check the return correction code status of tools
#include "PATInterfaces/CorrectionCode.h" 
#include "PATInterfaces/SystematicRegistry.h"
#include "PATInterfaces/SystematicVariation.h" 
#include "PATInterfaces/SystematicsUtil.h"

/// MET
#include "METUtilities/METMaker.h" /// METUtils
#include "xAODMissingET/MissingETContainer.h"
#include "xAODMissingET/MissingETAuxContainer.h"
#include "xAODMissingET/MissingETComposition.h"
#include "xAODMissingET/MissingETAssociationMap.h"
#include "METInterface/IMETMaker.h"
#include "METUtilities/CutsMETMaker.h"


/// LPXKfactorTool tool
#include "LPXKfactorTool/LPXKfactorTool.h"

#include <MyAnalysis/DataSetInfo.h>

#include <TFile.h>
#include <TMath.h>
#include <TLorentzVector.h>

/// OutputStream
#include "EventLoop/OutputStream.h"
#include "xAODRootAccess/tools/Message.h"

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

#include "PileupReweighting/PileupReweightingTool.h"
   
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
namespace CP{ /// these tools live in the namespace CP
  class MuonSelectionTool;
  class MuonCalibrationAndSmearingTool; 
  class IsolationSelectionTool;
//   class PileupReweightingTool;
  //~ class MuonEfficiencyCorrections;
}

// namespace TrigConf{
//   class TrigDecisionTool;
// }
// 
// namespace Trig{
//   class xAODConfigTool;
// }

class RecoAnalysis : public EL::Algorithm
{
  /// put your configuration variables here as public variables.
  /// that way they can be set directly from CINT and python.
  
  ///********************************** WARNING ***************************
  /// variables that don't get filled at submission time should be
  /// protected from being send from the submission node to the worker
  /// node (done by the //!)
  ///********************************** WARNING ***************************
public:
  
  /// Analysis Tools
  met::METMaker * m_metMaker; //!
  AsgElectronLikelihoodTool* m_LHToolTight2015; //!
  AsgElectronLikelihoodTool* m_LHToolMedium2015; //!
  GoodRunsListSelectionTool *m_grl; //!
  JetCleaningTool *m_jetCleaning; //!  
  JetCalibrationTool *m_jetCalibrationTool; //!
  JERTool *m_JERTool; //!
  JetVertexTaggerTool* m_jvtTool; //!
  CP::EgammaCalibrationAndSmearingTool *m_eleCalibrationTool; //! 
  Trig::TrigDecisionTool *m_trigDecisionTool; //!
  TrigConf::xAODConfigTool *m_trigConfigTool; //!
  CP::MuonSelectionTool* m_muonSelection; //!
  CP::MuonSelectionTool* m_looseMuonSelection; //!
  CP::MuonCalibrationAndSmearingTool *m_muonCalibrationAndSmearingTool; //!
  CP::IsolationSelectionTool *m_muonisolationSelectionTool; //!
  CP::IsolationSelectionTool *m_eleisolationSelectionTool; //!
  CP::PileupReweightingTool *m_pileupReweightingTool; //!
  CP::MuonTriggerScaleFactors *m_trig_sf; //!
  CP::MuonEfficiencyScaleFactors *m_effi_corr; //!
  CP::MuonEfficiencyScaleFactors *m_effi_corr_iso; //!
  // list of systematics
  std::vector<CP::SystematicSet> m_sysList; //!
  
  LPXKfactorTool* m_LPXKfactorTool; //!
  
  /// Custom classes
  HistObjectDumper *m_HistObjectDumper; //!
  BitsetCutflow* m_BitsetCutflow; //!
  DataSetInfo* m_DataSetInfo; //!
  
  /// 
  xAOD::TEvent *m_event;  //!
  xAOD::TStore *m_store;  //!
  const xAOD::EventInfo* m_eventInfo = 0; //!
  
  /// Bool Flags
  bool m_isMC; //!
  bool m_isMC15b; //!
  bool m_inclusiveWorZ; //!
  /// Bool Flags, initialized at setupJob - on submit stage and are the same for workers
  bool m_useHistObjectDumper;
  bool m_useBitsetCutflow;
  bool m_useCalibrationAndSmearingTool;
  bool m_runElectronChannel;
  bool m_doWprimeTruthMatching;
  bool m_doNotApplyTriggerCuts;
  
  /// Integer Variables
  int EventNumber; //!
  int m_eventCounter; //!
  unsigned int m_DSID; //!
  unsigned int m_pdgIdOfMother; //!

  float MuonPt; //!
  float MuonEta; //!
  float MuonPhi; //!
  float MissingEt; //!
  float MissingEx; //!
  float MissingEy; //!
  float OldMissingEt; //!
  float TransverseMass; //!
  float weighfilterEfficiency; //!
  float weightkFactor; //!
  float sampleLumi; //!

  /// defining the output file name and tree that we will put in the output 
  /// ntuple, also the one branch that will be in that tree 
  string outputName; //!
  string m_sampleName; //!

  ///********************************************
  /// [Histograms]
  /// h_<object>_<variable>_<specialCutSet>
  ///
  /// possible objects:
  /// - event
  /// - lepton
  /// - vetoElectron
  /// - vetoMuon
  /// - MET
  /// - leptonMET
  ///********************************************
  
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
  TH1D* hMu_pt_off; //!
  TH1D* hMu_mt_off; //!
  TH1D* hMu_MET_Muons_off; //!
  
  /// EventInfo
  TH1D* h_event_sampleLumi; //!
  TH1D* h_event_kFactor; //!
  TH1D* h_event_filterEfficiency; //!
  TH1D* h_event_SFWeight; //!
  TH1D* h_event_totalWeight; //!

  TH1D *h_mgen_all; //!
  TH1D *h_mgen; //!
  
  /// this is a standard constructor
  RecoAnalysis ();

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
  ClassDef(RecoAnalysis, 1);
};

#endif
