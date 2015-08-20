#include <EventLoop/Job.h>
#include <EventLoop/StatusCode.h>
#include <EventLoop/Worker.h>
#include <MyAnalysis/MyxAODAnalysis.h>
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

// #include "Histograms.cxx"

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

/// this is needed to distribute the algorithm to the workers
ClassImp(MyxAODAnalysis)

MyxAODAnalysis :: MyxAODAnalysis ()
{
  /// Here you put any code for the base initialization of variables,
  /// e.g. initialize all pointers to 0.  Note that you should only put
  /// the most basic initialization here, since this method will be
  /// called on both the submission and the worker node.  Most of your
  /// initialization code will go into histInitialize() and
  /// initialize().
  
  //m_useHistObjectDumper = true;

}



EL::StatusCode MyxAODAnalysis :: setupJob (EL::Job& job)
{
  /// Here you put code that sets up the job on the submission object
  /// so that it is ready to work with your algorithm, e.g. you can
  /// request the D3PDReader service or add output files.  Any code you
  /// put here could instead also go into the submission script.  The
  /// sole advantage of putting it here is that it gets automatically
  /// activated/deactivated when you add/remove the algorithm from your
  /// job, which may or may not be of value to you.
  
  /// let's initialize the algorithm to use the xAODRootAccess package
  job.useXAOD ();
  xAOD::Init( "MyxAODAnalysis" ).ignore(); /// call before opening first file
  
  m_useHistObjectDumper = true;
  m_useBitsetCutflow = true;
  m_useCalibrationAndSmearingTool = true;
  m_doWprimeTruthMatching = false;
  m_runElectronChannel = false;
  
  m_truthoption = 1;
  cout << "after init"<<endl;
  
  outputName = "outFile";
  
  EL::OutputStream out (outputName);
  job.outputAdd (out);
  
  return EL::StatusCode::SUCCESS;
}



EL::StatusCode MyxAODAnalysis :: histInitialize ()
{
  /// Here you do everything that needs to be done at the very
  /// beginning on each worker node, e.g. create histograms and output
  /// trees.  This method gets called before any input files are
  /// connected.

  /// [Jet histos]
  h_zPrimVtx = new TH1F("h_zPrimVtx", "h_zPrimVtx", 6000, -300.0, 300); 
  wk()->addOutput (h_zPrimVtx);
  h_jetPt_VeryLooseBadJets = new TH1F("h_jetPt_VeryLooseBadJets", 
                                      "VeryLooseBadJets and better with pt>20"
                                      , 300, 0, 3000); /// jet pt [GeV]
  wk()->addOutput (h_jetPt_VeryLooseBadJets);
  h_jetPt_LooseBadJets = new TH1F("h_jetPt_LooseBadJets", 
                                  "LooseBadJets and better with pt>20", 300, 0, 
                                  3000); /// jet pt [GeV]
  wk()->addOutput (h_jetPt_LooseBadJets);
  
  /// [Muon histos]
  h_MET_RefFinalFix = new TH1F("h_MET_RefFinalFix", "h_MET_RefFinalFix", 1000, 
                               0, 10000); 
  wk()->addOutput (h_MET_RefFinalFix);
  h_nSelectedMuons = new TH1I("h_nSelectedMuons", "h_nSelectedMuons", 5, 0, 5); 
  wk()->addOutput (h_nSelectedMuons);
  h_ptmuon = new TH1F("h_ptmuon", "h_ptmuon", 500, 0, 500);
  wk()->addOutput (h_ptmuon);
  
  /// WARNING test
  initHists();
  
//   /// [MET histos]
//   h_MET_RefJet = new TH1F("h_MET_RefJet", "h_MET_RefJet", 2000, 0, 2000);
//   wk()->addOutput (h_MET_RefJet);
//   h_MET_RefMuon = new TH1F("h_MET_RefMuon", "h_MET_RefMuon", 2000, 0, 2000);
//   wk()->addOutput (h_MET_RefMuon);
//   h_MET_RefTau = new TH1F("h_MET_RefTau", "h_MET_RefTau", 2000, 0, 2000);
//   wk()->addOutput (h_MET_RefTau);
//   h_MET_RefSoft = new TH1F("h_MET_RefSoft", "h_MET_RefSoft", 2000, 0, 2000);
//   wk()->addOutput (h_MET_RefSoft);
//   h_MET_RefTrack = new TH1F("h_MET_RefTrack", "h_MET_RefTrack", 2000, 0, 2000);
//   wk()->addOutput (h_MET_RefTrack);
//   h_MET_RefGamma = new TH1F("h_MET_RefGamma", "h_MET_RefGamma", 2000, 0, 2000);
//   wk()->addOutput (h_MET_RefGamma);
//   h_MET_RefElectron = new TH1F("h_MET_RefElectron", "h_MET_RefElectron", 2000, 0, 2000);
//   wk()->addOutput (h_MET_RefElectron);
//   h_MET_RefCore = new TH1F("h_MET_RefCore", "h_MET_RefCore", 2000, 0, 2000);
//   wk()->addOutput (h_MET_RefCore);
//   
//   h_MET_OldJet = new TH1F("h_MET_OldJet", "h_MET_OldJet", 2000, 0, 2000);
//   wk()->addOutput (h_MET_OldJet);
//   h_MET_OldMuon = new TH1F("h_MET_OldMuon", "h_MET_OldMuon", 2000, 0, 2000);
//   wk()->addOutput (h_MET_OldMuon);
//   h_MET_OldTau = new TH1F("h_MET_OldTau", "h_MET_OldTau", 2000, 0, 2000);
//   wk()->addOutput (h_MET_OldTau);
//   h_MET_OldSoft = new TH1F("h_MET_OldSoft", "h_MET_OldSoft", 2000, 0, 2000);
//   wk()->addOutput (h_MET_OldSoft);
//   h_MET_OldTrack = new TH1F("h_MET_OldTrack", "h_MET_OldTrack", 2000, 0, 2000);
//   wk()->addOutput (h_MET_OldTrack);
//   h_MET_OldGamma = new TH1F("h_MET_OldGamma", "h_MET_OldGamma", 2000, 0, 2000);
//   wk()->addOutput (h_MET_OldGamma);
//   h_MET_OldElectron = new TH1F("h_MET_OldElectron", "h_MET_OldElectron", 2000, 0, 2000);
//   wk()->addOutput (h_MET_OldElectron);
    
  h_MET_RefFinalFix_test = new TH1F("h_MET_RefFinalFix_test", 
                                    "h_MET_RefFinalFix_test", 1000, 0, 10000); 
  wk()->addOutput (h_MET_RefFinalFix_test);
  h_truthMET = new TH1F("h_truthMET", "Truth MET", 1000, 0, 10000); 
  wk()->addOutput (h_truthMET);
  
  /// [Mt histos]
  h_Mt = new TH1F("h_Mt", "h_Mt", 1000, 0, 10000); 
  wk()->addOutput (h_Mt);

  h_Mt_muonPtCut = new TH1F("h_Mt_muonPtCut", "h_Mt_muonPtCut", 1000, 0, 10000);
  wk()->addOutput (h_Mt_muonPtCut);

  /// get the output file, create a new TTree and connect it to that output
  /// define what branches will go in that tree

  TFile *outputFile = wk()->getOutputFile (outputName);
  tree = new TTree ("tree", "tree");
  tree->SetDirectory (outputFile);
  tree->Branch("EventNumber", &EventNumber);
  tree->Branch("RunNumber", &RunNumber);
  tree->Branch("LumiBlock", &LumiBlock);
  tree->Branch("MuonPt", &MuonPt);
  tree->Branch("MuonEta", &MuonEta);
  tree->Branch("MuonPhi", &MuonPhi);
  tree->Branch("MissingEt", &MissingEt);
  tree->Branch("MissingEx", &MissingEx);
  tree->Branch("MissingEy", &MissingEy);
  tree->Branch("OldMissingEt", &OldMissingEt);
  tree->Branch("TransverseMass", &TransverseMass);

  return EL::StatusCode::SUCCESS;
}



EL::StatusCode MyxAODAnalysis :: fileExecute ()
{
  /// Here you do everything that needs to be done exactly once for every
  /// single file, e.g. collect a list of all lumi-blocks processed
  return EL::StatusCode::SUCCESS;
}



EL::StatusCode MyxAODAnalysis :: changeInput (bool firstFile)
{
  /// Here you do everything you need to do when we change input files,
  /// e.g. resetting branch addresses on trees.  If you are using
  /// D3PDReader or a similarin that jet loop let's fill this histogram, so the 
  /// loop looks something like: service this method is not needed.
  return EL::StatusCode::SUCCESS;
}



EL::StatusCode MyxAODAnalysis :: initialize ()
{
  /// Here you do everything that you need to do after the first input
  /// file has been connected and before the first event is processed,
  /// e.g. create additional histograms based on which variables are
  /// available in the input files.  You can also create all of your
  /// histograms and trees in here, but be aware that this method
  /// doesn't get called if no events are processed.  So any objects
  /// you create here won't be available in the output if you have no
  /// input events.

  const char* APP_NAME = "MyxAODAnalysis";

  m_event = wk()->xaodEvent();
  m_store = new xAOD::TStore();
  
  /// as a check, let's see the number of events in our xAOD
  /// print long long int
  Info("initialize()", "Number of events = %lli", m_event->getEntries() ); 

  /// Event Info
  const xAOD::EventInfo* eventInfo = 0;
  if( ! m_event->retrieve( eventInfo, "EventInfo").isSuccess() ){
    Error("execute()", "Failed to retrieve event info collection in initialise. Exiting." );
    return EL::StatusCode::FAILURE;
  }  
  
  /// fill the branches of our trees
  bool isData = true;
  if(eventInfo->eventType( xAOD::EventInfo::IS_SIMULATION ) ){
    isData = false;
  }
  cout << "isData flag = " << isData<<endl;
  
  /// count number of events
  m_eventCounter = 0;

  /// number of clean events
  m_numCleanEvents = 0;

  /// GRL
  m_grl = new GoodRunsListSelectionTool("GoodRunsListSelectionTool");
  std::vector<std::string> vecStringGRL;
  vecStringGRL.push_back("/afs/cern.ch/user/a/atlasdqm/grlgen/All_Good/data15_13TeV.periodAllYear_DetStatus-v63-pro18-01_DQDefects-00-01-02_PHYS_StandardGRL_All_Good.xml");
//   vecStringGRL.push_back("/afs/cern.ch/user/a/atlasdqm/grlgen/All_Good/"
//   "data12_8TeV.periodAllYear_DetStatus-v61-pro14-02_DQDefects-00-01-00_PHYS_"
//   "StandardGRL_All_Good.xml");
  CHECK(m_grl->setProperty( "GoodRunsListVec", vecStringGRL));
  /// if true (default) will ignore result of GRL and will just pass all events
  CHECK(m_grl->setProperty("PassThrough", false)); 
  if (!m_grl->initialize().isSuccess()) { /// check this isSuccess
    Error("initialize()", "Failed to properly initialize the GRL. Exiting." );
    return EL::StatusCode::FAILURE;
  }

  /// initialize and configure the jet cleaning tool
//   m_jetCleaning = new JetCleaningTool("JetCleaning");
//   m_jetCleaning->msg().setLevel( MSG::INFO ); 
//   CHECK(m_jetCleaning->setProperty( "CutLevel", "VeryLooseBad"));
//   m_jetCleaning->initialize();
  
  /// Jet Calibration Tool
  cout << "initialise jet calib"<<endl;
  const std::string name = "MyxAODAnalysis"; //string describing the current thread, for logging
  TString jetAlgo = "AntiKt4EMTopo"; //String describing your jet collection, for example AntiKt4EMTopo or AntiKt4EMTopo (see above)
  TString config = "JES_MC15Prerecommendation_April2015.config"; //Path to global config used to initialize the tool (see above)
  TString calibSeq ="JetArea_Residual_Origin_EtaJES_GSC" ; //String describing the calibration sequence to apply (see above)
//   m_jetCalibrationTool = new JetCalibrationTool();
//   CHECK (m_jetCalibrationTool->setProperty("JetCollection",jetAlgo));
//   CHECK (m_jetCalibrationTool->setProperty("ConfigFile",config));
//   CHECK (m_jetCalibrationTool->setProperty("CalibSequence",calibSeq));
//   CHECK (m_jetCalibrationTool->setProperty("IsData",isData));
//   
//   if (! m_jetCalibrationTool->initialize().isSuccess() ){
//     Error("initialize()", "Failed to properly initialize the JES Tool. Exiting." );
//   }
  cout << "initialise jet calib done"<<endl;
  m_metMaker = new met::METMaker("METMakerTool");
  //m_metMaker->msg().setLevel( MSG::DEBUG ); // or DEBUG or VERBOSE
  if ( ! m_metMaker->initialize().isSuccess() ){
    Error("initialize()", "Failed to properly initialize the METMakerTool. Exiting." );
    return EL::StatusCode::FAILURE;
  }
  
  /// Configure the JVT tool.
  m_jvtTool = 0;
  ToolHandle<IJetUpdateJvt> hjvtagup;
  m_jvtTool = new JetVertexTaggerTool("jvtag");
  hjvtagup = ToolHandle<IJetUpdateJvt>("jvtag");
  m_jvtTool->msg().setLevel( MSG::DEBUG ); // or DEBUG or VERBOSE
  CHECK (m_jvtTool->setProperty("JVTFileName","JetMomentTools/JVTlikelihood_20140805.root"));
  if ( ! m_jvtTool->initialize().isSuccess() ){
    Error("initialize()", "Failed to properly initialize the JVTTool. Exiting." );
    return EL::StatusCode::FAILURE;
  }
  
  /// initialize JER 
  /*
  const char* jerFilePath = "$ROOTCOREBIN/data/JetResolution/JERProviderPlots"
  "_2012.root";
  const char* fullJERFilePath = gSystem->ExpandPathName (jerFilePath);
  m_JERTool = new JERTool("JERTool",fullJERFilePath,"AntiKt4LCTopoJES");

  if (! m_JERTool->initialize().isSuccess() ){
    Error("initialize()", "Failed to properly initialize the JER Tool. "
    "Exiting.");
    return EL::StatusCode::FAILURE;
  }
 */
  /// Get Muon Selector Tool
  m_muonSelection = new CP::MuonSelectionTool("MuonSelection");
  // m_muonSelection->msg().setLevel( MSG::INFO );
  m_muonSelection->setProperty( "MaxEta", 2.5 );
  m_muonSelection->setProperty( "MuQuality", 4);
  m_muonSelection->msg().setLevel( MSG::ERROR );
  CHECK (m_muonSelection->initialize().isSuccess());
  m_loosemuonSelection = new CP::MuonSelectionTool("MuonSelection");
  // m_muonSelection->msg().setLevel( MSG::INFO );
  m_loosemuonSelection->setProperty( "MaxEta", 2.5 );
  m_loosemuonSelection->setProperty( "MuQuality", 2);
  m_loosemuonSelection->msg().setLevel( MSG::ERROR );
  CHECK (m_loosemuonSelection->initialize().isSuccess());


  /// initialize the muon calibration and smearing tool
  m_muonCalibrationAndSmearingTool = new CP::MuonCalibrationAndSmearingTool(
    "MuonCorrectionTool");
  //m_muonCalibrationAndSmearingTool->msg().setLevel( MSG::DEBUG );
  if (! m_muonCalibrationAndSmearingTool->initialize().isSuccess() ){
    Error("initialize()", "Failed to properly initialize the "
    "MonCalibrationAndSmearingTool Tool. Exiting." );
    return EL::StatusCode::FAILURE;
  }
  
  //~ m_effi_corr = new CP::MuonEfficiencyScaleFactors m_effi_corr("SomeName");
  //~ m_effi_corr->setProperty("WorkingPoint","CBandST");
  //~ m_effi_corr->setProperty("DataPeriod","2012");
  //~ CHECK (m_effi_corr->initialize().isSuccess());
  
  //m_METUtil = new METUtility;
  //m_METUtil->setVerbosity(true);
  //m_util->setSoftJetCut(20);

  /// Initialize and configure trigger tools
  m_trigConfigTool = new TrigConf::xAODConfigTool("xAODConfigTool"); // gives us access to the meta-data
  EL_RETURN_CHECK( "initialize", m_trigConfigTool->initialize() );
  ToolHandle< TrigConf::ITrigConfigTool > trigConfigHandle( m_trigConfigTool );
  m_trigDecisionTool = new Trig::TrigDecisionTool("TrigDecisionTool");
  EL_RETURN_CHECK( "initialize", m_trigDecisionTool->setProperty( "ConfigTool", trigConfigHandle ) ); // connect the TrigDecisionTool to the ConfigTool
  EL_RETURN_CHECK( "initialize", m_trigDecisionTool->setProperty( "TrigDecisionKey", "xTrigDecision" ) );
  EL_RETURN_CHECK( "initialize", m_trigDecisionTool->initialize() );
    
  /// electrons
  cout << "do e/g calibration" << endl;
  m_eleCalibrationTool = new CP::EgammaCalibrationAndSmearingTool("eletool");
  m_eleCalibrationTool->setProperty("ESModel", "es2015PRE");  // see below for options
  m_eleCalibrationTool->setProperty("decorrelationModel", "1NP_v1");  // see below for options
  m_eleCalibrationTool->initialize();

  cout << "e/g calibration done" << endl;
  
  /// LH electron identification ...
  m_LHToolTight2015    = new AsgElectronLikelihoodTool ("m_LHToolTight2015");
  m_LHToolMedium2015   = new AsgElectronLikelihoodTool ("m_LHToolMedium2015"); 
  
  EL_RETURN_CHECK( "initialize", m_LHToolTight2015->setProperty(
    "primaryVertexContainer","PrimaryVertices"));
  EL_RETURN_CHECK( "initialize", m_LHToolMedium2015->setProperty(
    "primaryVertexContainer","PrimaryVertices"));
  
  std::string confDir = "ElectronPhotonSelectorTools/offline/mc15_20150408/";
  EL_RETURN_CHECK( "initialize", m_LHToolTight2015->setProperty(
    "ConfigFile",confDir+"ElectronLikelihoodTightOfflineConfig2015.conf"));
  EL_RETURN_CHECK( "initialize", m_LHToolMedium2015->setProperty(
    "ConfigFile",confDir+"ElectronLikelihoodMediumOfflineConfig2015.conf"));
  
  EL_RETURN_CHECK( "initialize", m_LHToolTight2015->initialize() );
  EL_RETURN_CHECK( "initialize", m_LHToolMedium2015->initialize() );
  /// ... LH electron identification
  cout << "e/g LH done"<<endl;

  /// Muon Isolation Tool
  m_muonisolationSelectionTool = new CP::IsolationSelectionTool("muoniso");
  EL_RETURN_CHECK( "initialize",m_muonisolationSelectionTool->setProperty( "MuonWP","LooseTrackOnly") );
  EL_RETURN_CHECK( "initialize",m_muonisolationSelectionTool->initialize());
  m_eleisolationSelectionTool = new CP::IsolationSelectionTool("eleiso");
  EL_RETURN_CHECK( "initialize",m_eleisolationSelectionTool->setProperty( "ElectronWP","Loose") );
  EL_RETURN_CHECK( "initialize",m_eleisolationSelectionTool->initialize());
  
  if (m_useHistObjectDumper)
    m_HistObjectDumper = new HistObjectDumper(wk());
  
  if (m_BitsetCutflow)
    m_BitsetCutflow = new BitsetCutflow(wk());
  
  /// TODO do we need this?
  for (int i=0; i<20; i++) {
    count[i]=0;
  }  

  return EL::StatusCode::SUCCESS;
}



EL::StatusCode MyxAODAnalysis :: execute ()
{
  /// Here you do everything that needs to be done on every single
  /// events, e.g. read input variables, apply cuts, and fill
  /// histograms and trees.  This is where most of your actual analysis
  /// code will go.
  
  const char* APP_NAME = "MyxAODAnalysis";

  /// push cutflow bitset to cutflow hist
  m_BitsetCutflow->PushBitSet();
  
  if( (m_eventCounter % 1000) ==0 ) Info("execute()", "Event number = %i", 
    m_eventCounter );
  m_eventCounter++;
  
  count[0]+=1;
  
  ///----------------------------
  /// Event information
  ///--------------------------- 

  const xAOD::EventInfo* eventInfo = 0;
  if( ! m_event->retrieve( eventInfo, "EventInfo").isSuccess() ){
    Error("execute()", "Failed to retrieve event info collection. Exiting." );
    return EL::StatusCode::FAILURE;
  }

  /// fill the branches of our trees
  EventNumber = eventInfo->eventNumber();  

  /// check if the event is data or MC
  /// (many tools are applied either to data or MC)
  bool isMC = false;
  /// check if the event is MC
  if(eventInfo->eventType( xAOD::EventInfo::IS_SIMULATION ) ){
    isMC = true; /// can do something with this later
  }
  
  /// Muon Truth matching. Check do we have muon from W' decay
  if(isMC && m_doWprimeTruthMatching){
    bool foundMuonFromWprimeDecay = false;
    
    /// Create truth vertice container
    const xAOD::TruthVertexContainer* truthVertices = 0; 
    /// retrieve arguments: container type, container key
    if ( !m_event->retrieve( truthVertices, "TruthVertices" ).isSuccess() ){ 
      Error("execute()","Failed to retrieve TruthVertices container. Exiting.");
      return EL::StatusCode::FAILURE;
    }
    
    /// Start iterating over truth container
    xAOD::TruthVertexContainer::const_iterator truthV_itr; 
    for (truthV_itr = truthVertices->begin(); truthV_itr != truthVertices->end();
         ++truthV_itr ) {
      if (foundMuonFromWprimeDecay)
        break;
      for (unsigned int iIn=0; iIn < (*truthV_itr)->nIncomingParticles(); iIn++)
      {
        if (TMath::Abs((*truthV_itr)->incomingParticle(iIn)->pdgId()) == 34) {
          for (unsigned int iOut=0; iOut < (*truthV_itr)->nOutgoingParticles(); 
               iOut++) {
            if (TMath::Abs((*truthV_itr)->outgoingParticle(iOut)->pdgId()) 
              == 13)
              foundMuonFromWprimeDecay = true;
            break;
          }
        }
        if (foundMuonFromWprimeDecay)
          break;
      }
    }
    
    if (!foundMuonFromWprimeDecay) return EL::StatusCode::SUCCESS;
    m_BitsetCutflow->FillCutflow("Truth muon decay");
  }
  
  ///------------------------------------------------------------
  /// Apply event cleaning to remove events due to 
  /// problematic regions of the detector, and incomplete events.
  /// Apply to data.
  ///------------------------------------------------------------
  /// reject event if:
  if(!isMC){
    if(   (eventInfo->errorState(xAOD::EventInfo::LAr)==xAOD::EventInfo::Error )
      || (eventInfo->errorState(xAOD::EventInfo::Tile)==xAOD::EventInfo::Error )
      || (eventInfo->isEventFlagBitSet(xAOD::EventInfo::Core, 18) )  )
    {
      return EL::StatusCode::SUCCESS; /// go to the next event
    } /// end if event flags check
  } /// end if the event is data
  m_numCleanEvents++;
  m_BitsetCutflow->FillCutflow("LAr_Tile_Core");
  count[1]+=1;
  
  
  /// examine the HLT_mu50* chains, see if they passed/failed and their total 
  /// prescale
  bool passTrigger = true;
  auto chainGroup = m_trigDecisionTool->getChainGroup("HLT_mu10.*");
//   auto chainGroup = m_trigDecisionTool->getChainGroup("HLT_mu50.*");
  std::map<std::string,int> triggerCounts;
  for(auto &trig : chainGroup->getListOfTriggers()) {
    auto cg = m_trigDecisionTool->getChainGroup(trig);
    std::string thisTrig = trig;
//     Info( "execute()", "%30s chain passed(1)/failed(0): %d total chain prescale "
//     "(L1*HLT): %.1f", thisTrig.c_str(), cg->isPassed(), cg->getPrescale() );
    if (cg->isPassed()==false)
      passTrigger = false;
  } /// end for loop (c++11 style) over chain group matching "HLT_mu50*" 
  
  /// TODO do we need such chain of triggers?
  chainGroup = m_trigDecisionTool->getChainGroup("HLT_noalg_L1MU10.*");
  for(auto &trig : chainGroup->getListOfTriggers()) {
    auto cg = m_trigDecisionTool->getChainGroup(trig);
    std::string thisTrig = trig;
//     Info( "execute()", "%30s chain passed(1)/failed(0): %d total chain prescale "
//     "(L1*HLT): %.1f", thisTrig.c_str(), cg->isPassed(), cg->getPrescale() );
    if (cg->isPassed())
      passTrigger = true;
  }
  chainGroup = m_trigDecisionTool->getChainGroup("HLT_mu18.*");
  for(auto &trig : chainGroup->getListOfTriggers()) {
    auto cg = m_trigDecisionTool->getChainGroup(trig);
    std::string thisTrig = trig;
//     Info( "execute()", "%30s chain passed(1)/failed(0): %d total chain prescale "
//     "(L1*HLT): %.1f", thisTrig.c_str(), cg->isPassed(), cg->getPrescale() );
    if (cg->isPassed())
      passTrigger = true;
  }
  
  /// TODO do we need !isMC requirement?
  if(!isMC){
    if (passTrigger==false)
      return EL::StatusCode::SUCCESS;
    m_BitsetCutflow->FillCutflow("Trigger");
    count[2]+=1;
  }
  
  /// if data check if event passes GRL
  if(!isMC){ /// it's data!
    if(!m_grl->passRunLB(*eventInfo)){
      return EL::StatusCode::SUCCESS; // go to next event
    }
  } /// end if not MC
  m_BitsetCutflow->FillCutflow("GRL");

  std::vector<const xAOD::IParticle*> uniques;
  /// loop over electrons (needed for MET)
  const xAOD::ElectronContainer* electrons(0);
  m_event->retrieve( electrons, "Electrons");
  if ( !m_event->retrieve( electrons, "Electrons" ).isSuccess() ){ /// retrieve arguments: container$
    Error("execute()", "Failed to retrieve Electrons container. Exiting." );
    return EL::StatusCode::FAILURE;
  }
  
  const xAOD::PhotonContainer* photons(0);
  m_event->retrieve( photons, "Photons");
  if ( !m_event->retrieve( photons, "Photons" ).isSuccess() ){ /// retrieve arguments: container$
    Error("execute()", "Failed to retrieve Photons container. Exiting." );
    return EL::StatusCode::FAILURE;
  }
  
  uniques.clear();
  ConstDataVector<xAOD::PhotonContainer> metPhotons(SG::VIEW_ELEMENTS);
  //for(const auto& ph : *photons) {
  for(const auto& ph : *photons) {
    if(CutsMETMaker::accept(ph)) {
      metPhotons.push_back(ph);
    }
  }
  
  /// full copy 
  /// get muon container of interest
  const xAOD::MuonContainer* muons = 0;
  if ( !m_event->retrieve( muons, "Muons" ).isSuccess() ){ /// retrieve arguments: container$
    Error("execute()", "Failed to retrieve Muons container. Exiting." );
    return EL::StatusCode::FAILURE;
  }

  std::pair<xAOD::MuonContainer*,xAOD::ShallowAuxContainer*> muoncopy = xAOD::shallowCopyContainer(*muons);
  xAOD::setOriginalObjectLink(*muons, *muoncopy.first); // This line is very important ! Overlap removal will fail if line is missing
  for(const auto& muon : *muoncopy.first) {
    if (fabs(muon->eta())<2.5) {
      CP::CorrectionCode result = m_muonCalibrationAndSmearingTool->applyCorrection(*muon);
      if(result != CP::CorrectionCode::Ok){
    cout << "muon pt " << muon->pt()<< " eta " << muon->eta()<<endl;
    throw std::runtime_error("Error when calibrating muons. Exiting." );
      }
    }
  }
  m_store->record(muoncopy.first,  "MuonCopy");
  m_store->record(muoncopy.second, "MuonCopyAux");
  
  uniques.clear();
  ConstDataVector<xAOD::MuonContainer> metMuons(SG::VIEW_ELEMENTS);
  for(const auto& muon : *muoncopy.first) {

    if(CutsMETMaker::accept(muon)) {
      metMuons.push_back(muon);
    }
  }
  
  ConstDataVector<xAOD::MuonContainer> noMuons(SG::VIEW_ELEMENTS);
  
  const xAOD::TauJetContainer* taus(0);
  m_event->retrieve( taus, "TauJets");
  if ( !m_event->retrieve( taus, "TauJets" ).isSuccess() ){ /// retrieve arguments: container$
    Error("execute()", "Failed to retrieve Taus container. Exiting." );
    return EL::StatusCode::FAILURE;
  }
  
  uniques.clear();
  ConstDataVector<xAOD::TauJetContainer> metTaus(SG::VIEW_ELEMENTS);
  for(const auto& tau : *taus) {
    if(CutsMETMaker::accept(tau)) {
      metTaus.push_back(tau);
    }
  }
  
  /// LOOP OVER JETS
  const xAOD::JetContainer* jets(0);
  m_event->retrieve(jets, "AntiKt4EMTopoJets");
  uniques.clear();

  std::pair<xAOD::JetContainer*,xAOD::ShallowAuxContainer*> jetcopy = xAOD::shallowCopyContainer(*jets);
  xAOD::setOriginalObjectLink(*jets, *jetcopy.first); // This line is very important ! Overlap removal will fail if line is missing
  for(const auto& jet : *jetcopy.first) {

    //cout << "original jet pt " << (*jet).pt()/1000. <<endl;
    /*
    CP::CorrectionCode result = m_jetCalibrationTool->applyCorrection(*jet);
    if(result != CP::CorrectionCode::Ok){
      throw std::runtime_error("Error when calibrating jets. Exiting." );
    }
    */
    bool hasjvt = jet->isAvailable<float>("Jvt");
    float oldjvt = -999;
    if (hasjvt) oldjvt = jet->auxdata<float>("Jvt");
    //cout << "calib jet pt " << (*jet).pt()/1000. << " calib jet pt " <<endl;
    float newjvt = m_jvtTool->updateJvt(*jet);
    jet->auxdecor<float>("Jvt") = newjvt;
    //cout << "oldjvt " << oldjvt << " newjvt " << newjvt<<endl;
  }

  m_store->record(jetcopy.first,  "JetCopy");
  m_store->record(jetcopy.second, "JetCopyAux");

  std::pair<xAOD::ElectronContainer*,xAOD::ShallowAuxContainer*> elecopy = xAOD::shallowCopyContainer(*electrons);
  xAOD::setOriginalObjectLink(*electrons, *elecopy.first); // This line is very important ! Overlap removal will fail if line is missing
  for(const auto& ele : *elecopy.first) {

    //cout << "original jet pt " << (*jet).pt()/1000. <<endl;
    CP::CorrectionCode result = m_eleCalibrationTool->applyCorrection(*ele);
    if(result != CP::CorrectionCode::Ok){
      throw std::runtime_error("Error when calibrating jets. Exiting." );
    }
  }
  m_store->record(elecopy.first,  "EleCopy");
  m_store->record(elecopy.second, "EleCopyAux");
  
  const xAOD::VertexContainer* vertices = 0;
  /// retrieve arguments: container type, container key
  if ( !m_event->retrieve( vertices, "PrimaryVertices" ).isSuccess() ){ 
    Error("execute()","Failed to retrieve PrimaryVertices container. Exiting.");
    return EL::StatusCode::FAILURE;
  }
  
  xAOD::VertexContainer::const_iterator vtx_itr = vertices->begin();
  xAOD::VertexContainer::const_iterator vtx_end = vertices->end();
  xAOD::Vertex* primVertex = 0;
  int nGoodVtx = 0;
  for( ; vtx_itr != vtx_end; ++vtx_itr ) {
    if ((*vtx_itr)->vertexType()==xAOD::VxType::PriVtx){
      primVertex = (*vtx_itr);
      nGoodVtx++;
    }
  }

  /// LOOP OVER JETS

  /// Loop over all jets in the event
  /// get jet container of interest
//   const xAOD::JetContainer* jets = 0;

  /// retrieve arguments: container type, container key
  if ( !m_event->retrieve( jets, "AntiKt4EMTopoJets" ).isSuccess() ){ 
    Error("execute()","Failed to retrieve AntiKt4EMTopoJets container. "
    "Exiting.");
    return EL::StatusCode::FAILURE;
  }

  /// loop over the jets in the container
  /// not used in current cutflow
  /*
  bool foundMassyJet = false;
  xAOD::JetContainer::const_iterator jet_itr = jets->begin();
  xAOD::JetContainer::const_iterator jet_end = jets->end();
  for( ; jet_itr != jet_end; ++jet_itr ) {
    
    /// check for Jet pt (not interested in Jets with pt<=20)
    if (( (*jet_itr)->pt()) * 0.001 <= 20.0) continue;
    
    /// check if jet satisfy at least VeryLooseBad jet requirement
    /// if no - discar event
    if( !m_jetCleaning->accept( **jet_itr )){
      return EL::StatusCode::SUCCESS;
    }
    
  } // end for loop over jets
  
  m_BitsetCutflow->FillCutflow("JetCleaning");
  */
  
  if (nGoodVtx>1)
    cout << "WARNING!!!! Found more then one prim.vertex: nGoodVtx = " 
    << nGoodVtx << endl;
  if (nGoodVtx==0)
    return EL::StatusCode::SUCCESS;
  m_BitsetCutflow->FillCutflow("Primary vertex");
  
  
  // select electrons and muon candidate using W' cut flow
  xAOD::Electron* signalEl = 0;

  static SG::AuxElement::Accessor< xAOD::IParticleLink > accSetOriginLink ("originalObjectLink");
  ConstDataVector<xAOD::MuonContainer> signalMuon(SG::VIEW_ELEMENTS); // Create a new muon container
  // not really used right now
  if (!m_runElectronChannel){
    /// loop over the muons in the container
    /// signal selection
    //signalMuon = SelectMuon(muons,primVertex);
    for(const auto& mu : *muoncopy.first) {
      if (passMuonSelection(mu, eventInfo, primVertex, false)) {
    //const xAOD::IParticleLink originLink( *muons, (*mu).index() );
    //accSetOriginLink(*mu) = originLink;
    signalMuon.push_back(mu);
    //cout << "pass muon with " << mu->pt()<<endl;
      }
    }

    if (signalMuon.size()==0)
      return EL::StatusCode::SUCCESS;
  }
  else{
    // to be filled
  }
  
  ConstDataVector<xAOD::MuonContainer> newMuons(SG::VIEW_ELEMENTS);
  
  for(const auto& mu : *muoncopy.first) {
    // put your cuts here
      if ((*mu).muonType()!=xAOD::Muon_v1::Combined) continue;
      double muPt = ((*mu).pt()) * 0.001;
      if (muPt < 55 || muPt>=9999999) continue; /// veto muon
      if(!m_muonSelection->accept(mu)) continue;
      /// do significance 
      double d0_sig = TMath::Abs((*mu).primaryTrackParticle()->d0()) /
    TMath::Sqrt((*mu).primaryTrackParticle()->definingParametersCovMatrix()(0,0)
            + eventInfo->beamPosSigmaX()*eventInfo->beamPosSigmaX() );
      if (d0_sig>3.0) continue;
      /// zo cut
      double z0_vrtPVx = (*mu).primaryTrackParticle()->z0() +
    (*mu).primaryTrackParticle()->vz() - primVertex->z();
      double sintheta = 1.0/TMath::CosH((*mu).eta());
      if (abs( z0_vrtPVx*sintheta )>10.0) continue;
      if (!m_muonisolationSelectionTool->accept(*mu)) continue;
      const xAOD::IParticleLink originLink( *muons, (*mu).index() );
      accSetOriginLink(*mu) = originLink;
      newMuons.push_back(mu);
  }
  
  //cout << "no of signal muons " <<newMuons.size()<<endl; 
  /// look for veto muon, add this for the moment
  //if (signalMuon.size()>1) return EL::StatusCode::SUCCESS;
  if (signalMuon.size()==0 || signalMuon.size()>1) return EL::StatusCode::SUCCESS;

  ConstDataVector<xAOD::MuonContainer> vetoMu(SG::VIEW_ELEMENTS); // Create a new muon container
  vetoMu = VetoMuon(muoncopy.first,primVertex,true);
  //cout << "Run/Event " << EventNumber << " signal " << signalMuon.size() << ", veto " << vetoMu.size()<<endl;
  if (vetoMu.size()>1) return EL::StatusCode::SUCCESS;
  if (signalMuon.size()>1) return EL::StatusCode::SUCCESS;
  m_BitsetCutflow->FillCutflow("Veto muon");
  count[12]+=1;
  ConstDataVector<xAOD::ElectronContainer> vetoEl(SG::VIEW_ELEMENTS);
  vetoEl = SelectElectron(elecopy.first,primVertex,true);
  if (vetoEl.size()!=0) return EL::StatusCode::SUCCESS;
  m_BitsetCutflow->FillCutflow("Veto electron");
  count[13]+=1;
  if (nGoodVtx>1)
    cout << "WARNING!!!! Found more then one prim.vertex: nGoodVtx = " << nGoodVtx << endl;
  if (nGoodVtx==0)  return EL::StatusCode::SUCCESS;
  m_BitsetCutflow->FillCutflow("Primary vertex");
  count[14]+=1;
  // move met to last step as it needs our selected muons as input
  // calculate MET
  bool doJVFCut = false;
  std::string softTerm = "PVSoftTrk";
  std::string finalTerm = "FinalTrk";

  // comment for the moment as this does not work on EXOT9
  xAOD::MissingETContainer* met = new xAOD::MissingETContainer;
  xAOD::MissingETAuxContainer* met_aux = new xAOD::MissingETAuxContainer;
  met->setStore(met_aux);

  const xAOD::MissingETAssociationMap* metMap(0);
  m_event->retrieve( metMap, "METAssoc_AntiKt4EMTopo" );
  metMap->resetObjSelectionFlags() ;

  const xAOD::MissingETContainer* metcore(0);
  m_event->retrieve( metcore, "MET_Core_AntiKt4EMTopo" );

  const xAOD::MissingETContainer* metold(0);
  m_event->retrieve( metold, "MET_Reference_AntiKt4EMTopo" );

  if ( !m_event->retrieve( metold, "MET_Reference_AntiKt4EMTopo" ).isSuccess() ){ /// retrieve arguments: container$
    
    Error("execute()", "Failed to retrieve old met container. Exiting." );
    return EL::StatusCode::FAILURE;
  }

  bool doTracks = true;
  // not needed any more as now I have my selected muons
  //cout << "input electrons"<<endl;
  //m_metMaker->rebuildMET("NewRefEle", xAOD::Type::Electron, met, metElectrons.asDataVector(), metMap);
  //m_metMaker->rebuildMET("NewRefEle", xAOD::Type::Electron, met, vetoEl.asDataVector(), metMap);

  //m_metMaker->rebuildMET("NewRefGamma", xAOD::Type::Photon, met, metPhotons.asDataVector(), metMap);
  //m_metMaker->rebuildMET("NewRefTau", xAOD::Type::Tau, met, metTaus.asDataVector(), metMap);
  //cout << "input jets"<<endl;
  m_metMaker->rebuildMET("NewMuons", xAOD::Type::Muon, met, signalMuon.asDataVector(), metMap);
  //m_metMaker->rebuildMET("NewMuons", xAOD::Type::Muon, met, newMuons.asDataVector(), metMap);
  //m_metMaker->rebuildMET("NewMuons", xAOD::Type::Muon, met, metMuons.asDataVector(), metMap);
  //m_metMaker->rebuildMET("NewMuons", xAOD::Type::Muon, met, noMuons.asDataVector(), metMap);

  m_metMaker->rebuildJetMET("NewRefJet", "PVSoftTrk", met, jetcopy.first, metcore, metMap,false);
  // m_metMaker->rebuildJetMET("NewRefJet",  "SoftClus", "PVSoftTrk", met, jetcopy.first, metcore, metMap,false); 
 //m_metMaker->rebuildJetMET("RefJet",  "SoftClus", softTerm, met, metJets.asDataVector(), metcore, metMap, doJVFCut); 

  // m_metMaker->rebuildJetMET("RefJet",  "SoftClus", softTerm, met, jets, metcore, metMap, doJVFCut); 
 //cout << "input final term"<<endl;
  m_metMaker->buildMETSum("NewFinalTrk", met, (*met)[softTerm]->source());

  double mpx = (*met)["NewFinalTrk"]->mpx();
  double mpy = (*met)["NewFinalTrk"]->mpy();

  TLorentzVector *metVec = new TLorentzVector();
  metVec->SetPx(mpx);
  metVec->SetPy(mpy);
  metVec->SetPz(0.0);
  metVec->SetE(sqrt(mpx*mpx + mpy*mpy));
  double phi_met = metVec->Phi();
  //cout << "MET new " << mpx << ", " << mpy << ", " << sqrt(mpx*mpx + mpy*mpy)<< " phi " << phi_met<<endl;
  
  // uncalibrated met
  /// get MET_RefFinalFix container of interest
  //cout << "get old met values"<<endl;
  const xAOD::MissingETContainer* metcontainer = 0;

  /// we need to use MET_RefFinalFix, according to:
  /// https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/Run2xAODMissingET
  if ( !m_event->retrieve(metcontainer, "MET_Reference_AntiKt4EMTopo" ).isSuccess() ){ /// retrieve arguments: container$
    Error("execute()", "Failed to retrieve MET_Reference_AntiKt4EMTopo container. Exiting." );
    return EL::StatusCode::FAILURE;
  }

  /// We want to check default Final MET, wo any recalibration first
  /*
  xAOD::MissingETContainer::const_iterator met_itold = metcontainer->find("FinalTrk"); 
  
  if (met_itold == metcontainer->end()) {
    Error("execute()", "No RefFinal inside MET container" );
  }
  
  double mpx2 = (*met_itold)->mpx();
  double mpy2 = (*met_itold)->mpy();
  cout << "MET old " << mpx2 << ", " << mpy2 << ", " << sqrt(mpx2*mpx2 + mpy2*mpy2)<<endl;
  */
  // easier
  double mpx2 = (*met)["NewRefJet"]->mpx();
  double mpy2 = (*met)["NewRefJet"]->mpy();
  double refjet = sqrt(pow(mpx2,2)+pow(mpy2,2));
  //cout << "new met refjet " << mpx2/1000.<< ", " << refjet<<endl;
  mpx2 = (*met)["NewMuons"]->mpx();
  mpy2 = (*met)["NewMuons"]->mpy();
  //double refmuon = 0.;
  double refmuon = sqrt(pow(mpx2,2)+pow(mpy2,2));
  
  double refsoft = sqrt(pow(mpx2,2)+pow(mpy2,2));
  mpx2 = (*met)["PVSoftTrk"]->mpx();
  mpy2 = (*met)["PVSoftTrk"]->mpy();
  double reftrk = sqrt(pow(mpx2,2)+pow(mpy2,2));
  
  mpx2 = (*metcontainer)["RefJet"]->mpx();
  mpy2 = (*metcontainer)["RefJet"]->mpy();
  double refjetold = sqrt(pow(mpx2,2)+pow(mpy2,2));
  //cout << "old met refjet " << mpx2/1000.<<endl;

  mpx2 = (*metcontainer)["Muons"]->mpx();
  mpy2 = (*metcontainer)["Muons"]->mpy();
  //double refmuon = 0.;
  double refmuonold = sqrt(pow(mpx2,2)+pow(mpy2,2));
  //cout << "old met refmuon " << mpx2/1000.<<endl;

  mpx2 = (*metcontainer)["RefGamma"]->mpx();
  mpy2 = (*metcontainer)["RefGamma"]->mpy();
  double refgammaold = sqrt(pow(mpx2,2)+pow(mpy2,2));
  //cout << "old met refgam " << mpx2/1000.<<endl;

  mpx2 = (*metcontainer)["RefEle"]->mpx();
  mpy2 = (*metcontainer)["RefEle"]->mpy();
  double refeleold = sqrt(pow(mpx2,2)+pow(mpy2,2));
  //cout << "old met refele " << mpx2/1000.<<endl;

  mpx2 = (*metcontainer)["RefTau"]->mpx();
  mpy2 = (*metcontainer)["RefTau"]->mpy();
  double reftauold = sqrt(pow(mpx2,2)+pow(mpy2,2));

  //cout << "soft term"<<endl;
  mpx2 = (*metcontainer)["SoftClus"]->mpx();
  mpy2 = (*metcontainer)["SoftClus"]->mpy();
  double refsoftold = sqrt(pow(mpx2,2)+pow(mpy2,2));
  //cout << "soft: "<<refsoft<<endl;
  //cout << "old met refsoft " << mpx2/1000.<<endl;

  mpx2 = (*metcontainer)["PVSoftTrk"]->mpx();
  mpy2 = (*metcontainer)["PVSoftTrk"]->mpy();
  double reftrkold = sqrt(pow(mpx2,2)+pow(mpy2,2));
  //cout << "rtrk: "<<reftrk<<endl;
  //cout << "old met reftrk " << mpx2/1000.<< " " << reftrkold<<endl;
 
  //cout << "METcore result " << refcore<<endl;
  mpx2 = (*metcontainer)["FinalTrk"]->mpx();
  mpy2 = (*metcontainer)["FinalTrk"]->mpy();
  double reftotal = sqrt(pow(mpx,2)+pow(mpy,2));
  double oldtotal = sqrt(pow(mpx2,2)+pow(mpy2,2));
  //cout << "new MET: " << reftotal/1000. << " old MET: " << oldtotal/1000.<<endl;
  //cout << "MET components: jet " << refjet/1000 << " muon " << refmuon/1000 << " gamma " << refgamma/1000 << " electron " << refele/1000 << " tau " << reftau/1000.<< " total " <<reftotal/1000.<< endl;
  //cout << "MET old " << mpx2 << ", " << mpy2 << ", " << sqrt(mpx2*mpx2 + mpy2*mpy2)<<endl;

  //cout << "MET new " << mpx << ", " << mpy << ", " << sqrt(mpx*mpx + mpy*mpy)<<endl;
  TLorentzVector *metVecold = new TLorentzVector();
  metVecold->SetPx(mpx2);
  metVecold->SetPy(mpy2);
  metVecold->SetPz(0.0);
  metVecold->SetE(sqrt(mpx2*mpx2 + mpy2*mpy2));
  double phi_metold = metVecold->Phi();
  //cout << "new MET: "<<metVec->Pt()*0.001<< " old MET: " << metVecold->Pt()* 0.001<<endl; 
  //if (metVec->Pt()<55000)
  //if (metVecold->Pt()<20000)
  //return EL::StatusCode::SUCCESS;
  //m_BitsetCutflow->FillCutflow("MET cut");
  count[15]+=1;
  // calculate MT
  xAOD::MuonContainer::const_iterator muon_itr = signalMuon.begin();
  xAOD::MuonContainer::const_iterator muon_end = signalMuon.end();
  for (muon_itr=signalMuon.begin(); muon_itr != muon_end; ++muon_itr ){
    xAOD::Muon* mu = 0;
    mu = const_cast<xAOD::Muon*> (*muon_itr);
    double phi_mu = mu->phi();
    double pt_mu = mu->pt();
    //double phi_mu = signalMuon->phi();
    double Mt = sqrt( 2*mu->pt()*sqrt(mpx*mpx + mpy*mpy) * (1.0 - TMath::Cos( phi_mu - phi_met )) );
    //double Mt = sqrt( 2*pt_mu*sqrt(mpx2*mpx2 + mpy2*mpy2) * (1.0 - TMath::Cos( phi_mu - phi_metold )) );
    //cout << "MT = " << Mt<<endl;
    h_Mt_muonPtCut->Fill(Mt * 0.001);
    //if (Mt<60000.)
    //return EL::StatusCode::SUCCESS;
    m_BitsetCutflow->FillCutflow("MT cut");
    count[16]+=1;

    EventNumber = eventInfo->eventNumber();
    RunNumber = eventInfo->runNumber();
    LumiBlock = eventInfo->lumiBlock();
    MuonPt = mu->pt();
  MuonEta = mu->eta();
    MuonPhi = mu->phi();
    MissingEt = metVec->Pt();
    MissingEx = metVec->Px();
    MissingEy = metVec->Py();
    OldMissingEt = metVecold->Pt();
    TransverseMass = Mt;
    tree->Fill();
  }

  /// ************************************************
  /// Fill other event-based distributions (e.g. MET)
  /// ************************************************
  h_MET_RefJet->Fill(refjet* 0.001);
  //h_MET_RefCore->Fill(refcore* 0.001);
  //h_MET_RefTau->Fill(reftau* 0.001);
  //h_MET_RefGamma->Fill(refgamma* 0.001);
  //h_MET_RefElectron->Fill(refele* 0.001);
  h_MET_RefMuon->Fill(refmuon* 0.001);
  h_MET_RefSoft->Fill(refsoft* 0.001);
  h_MET_RefTrack->Fill(reftrk* 0.001);

  h_MET_OldJet->Fill(refjetold* 0.001);
  h_MET_OldTau->Fill(reftauold* 0.001);
  h_MET_OldGamma->Fill(refgammaold* 0.001);
  h_MET_OldElectron->Fill(refeleold* 0.001);
  h_MET_OldMuon->Fill(refmuonold* 0.001);
  h_MET_OldSoft->Fill(refsoftold* 0.001);
  h_MET_OldTrack->Fill(reftrkold* 0.001);
  
//   h_MET_RefFinalNew->Fill(sqrt(mpx*mpx + mpy*mpy)* 0.001);
//   h_MET_RefFinalFix->Fill(sqrt(mpx2*mpx2 + mpy2*mpy2)* 0.001);
  h_MET_RefFinalFix_test->Fill(reftotal * 0.001);

  m_BitsetCutflow->FillCutflow("End");
  /*
  /// get MET_RefFinalFix container of interest
  const xAOD::MissingETContainer* metcontainer = 0;
  
  /// we need to use MET_Core_AntiKt4EMTopo, according to:
  /// https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/Run2xAODMissingET
  /// FIXME use MET_Reference_AntiKt4EMTopo instead of MET_Core_AntiKt4EMTopo
  if ( !m_event->retrieve( metcontainer, "MET_Reference_AntiKt4EMTopo" ).isSuccess())
  { /// retrieve arguments: container$
    Error("execute()","Failed to retrieve MET_RefFinalFix container. Exiting.");
    return EL::StatusCode::FAILURE;
  }

  
  
  
  /// We want to check default Final MET, wo any recalibration first
  xAOD::MissingETContainer::const_iterator met_it = metcontainer->find("FinalTrk"); 

  if (met_it == metcontainer->end()) {
    Error("execute()", "No RefFinal inside MET container" );
  }

  double mpx = (*met_it)->mpx();  
  double mpy = (*met_it)->mpy();
  TLorentzVector *metVec = new TLorentzVector();
  metVec->SetPx(mpx);
  metVec->SetPy(mpy);
  metVec->SetPz(0.0);
  metVec->SetE(sqrt(mpx*mpx + mpy*mpy));
  
  double phi_met = metVec->Phi();
  
  //~ // create a shallow copy of the muons container
  //~ std::pair< xAOD::MuonContainer*, xAOD::ShallowAuxContainer* > muons_shallowCopy = xAOD::shallowCopyContainer( *muons );
//~ 
  //~ // iterate over our shallow copy
  //~ xAOD::MuonContainer::iterator muonSC_itr = (muons_shallowCopy.first)->begin();
  //~ xAOD::MuonContainer::iterator muonSC_end = (muons_shallowCopy.first)->end();
//~ 
  //~ for( ; muonSC_itr != muonSC_end; ++muonSC_itr ) {
   //~ if(m_muonCalibrationAndSmearingTool->applyCorrection(**muonSC_itr) == CP::CorrectionCode::Error){ // apply correction and check return code
       //~ // Can have CorrectionCode values of Ok, OutOfValidityRange, or Error. Here only checking for Error.
       //~ // If OutOfValidityRange is returned no modification is made and the original muon values are taken.
       //~ Error("execute()", "MuonCalibrationAndSmearingTool returns Error CorrectionCode");
   //~ }
   //~ //if(!m_muonSelection->accept(**muonSC_itr)) continue;
   //~ Info("execute()", "  corrected muon pt = %.2f GeV", ((*muonSC_itr)->pt() * 0.001));  
  //~ } // end for loop over shallow copied muons
  //~ delete muons_shallowCopy.first;
  //~ delete muons_shallowCopy.second;
  
  /// get muon container of interest
  const xAOD::MuonContainer* muons = 0;
  if ( !m_event->retrieve( muons, "Muons" ).isSuccess() ){ /// retrieve arguments: container$
    Error("execute()", "Failed to retrieve Muons container. Exiting." );
    return EL::StatusCode::FAILURE;
  }  
  
  const xAOD::ElectronContainer* electrons = 0;
  if ( !m_event->retrieve( electrons, "Electrons" ).isSuccess() ){ /// retrieve arguments: container$
    Error("execute()", "Failed to retrieve Electrons container. Exiting." );
    return EL::StatusCode::FAILURE;
  }  
  
  xAOD::Muon* signalMuon = 0;
  xAOD::Muon* signalEl = 0;
  
  if (!m_runElectronChannel){
    /// loop over the muons in the container
    /// signal selection
    signalMuon = SelectMuon(muons,primVertex);
    if (signalMuon==0)
      return EL::StatusCode::SUCCESS;
  }
  else{
    
  }



  /// look for veto muon
  xAOD::Muon* vetoMu = SelectMuon(muons,primVertex,true);
  if (vetoMu!=0)
    return EL::StatusCode::SUCCESS;
  m_BitsetCutflow->FillCutflow("Veto muon");

  xAOD::Electron* vetoEl = SelectElectron(electrons,primVertex,true);
  if (vetoEl!=0)
    return EL::StatusCode::SUCCESS;
  m_BitsetCutflow->FillCutflow("Veto electron");
  
  
  double phi_mu;
  double Mt;

  if (m_runElectronChannel){
    phi_mu = signalEl->phi();
    Mt = sqrt( 2*signalEl->pt()*sqrt(mpx*mpx + mpy*mpy) * 
    (1.0 - TMath::Cos( phi_mu - phi_met )) );
  }
  else{
    phi_mu = signalMuon->phi();
    Mt = sqrt( 2*signalMuon->pt()*sqrt(mpx*mpx + mpy*mpy) * 
    (1.0 - TMath::Cos( phi_mu - phi_met )) );
  }
  
  
  
  h_Mt_muonPtCut->Fill(Mt * 0.001);
      
//   if (m_useHistObjectDumper) m_HistObjectDumper->plotMuon(mu,"allCuts");
  
  
  /// ************************************************
  /// Fill other event-based distributions (e.g. MET)
  /// ************************************************
  h_MET_RefFinalFix->Fill(sqrt(mpx*mpx + mpy*mpy)* 0.001);
  h_MET_RefFinalFix_test->Fill(metVec->Pt() * 0.001);

  m_BitsetCutflow->FillCutflow("End");
 */
  return EL::StatusCode::SUCCESS;
}



EL::StatusCode MyxAODAnalysis :: postExecute ()
{
  /// Here you do everything that needs to be done after the main event
  /// processing.  This is typically very rare, particularly in user
  /// code.  It is mainly used in implementing the NTupleSvc.
  return EL::StatusCode::SUCCESS;
}



EL::StatusCode MyxAODAnalysis :: finalize ()
{
  /// This method is the mirror image of initialize(), meaning it gets
  /// called after the last event has been processed on the worker node
  /// and allows you to finish up any objects you created in
  /// initialize() before they are written to disk.  This is actually
  /// fairly rare, since this happens separately for each worker node.
  /// Most of the time you want to do your post-processing on the
  /// submission node after all your histogram outputs have been
  /// merged.  This is different from histFinalize() in that it only
  /// gets called on worker nodes that processed input events.
  
  const char* APP_NAME = "MyxAODAnalysis";  

  /// push cutflow for last event
  /// FIXME probably this breaks possibility to use PROOF
  m_BitsetCutflow->PushBitSet();
  
  ///*************************
  ///
  /// deleting of all tools
  ///
  /// ************************
  
  for (int i=0; i<17; i++) {
      cout << "cut " << i << " value: " << count[i]<<endl;
    }
  
  /// GRL
  if( m_grl ) {
    delete m_grl;
    m_grl = 0;
  }
  /// Jet cleaning tool
//   if( m_jetCleaning ) {
//     delete m_jetCleaning;
//     m_jetCleaning = 0;
//   }
//   
//   if( m_jetCalibrationTool ) {
//     delete m_jetCalibrationTool;
//     m_jetCalibrationTool = 0;
//   }
  
  
  /*
  if(m_JERTool){
    delete m_JERTool;
    m_JERTool = 0;
  }
  */
  Info("finalize()", "Number of clean events = %i", m_numCleanEvents);
  /// Muon calibration and smearing tool
  if(m_muonCalibrationAndSmearingTool){
    delete m_muonCalibrationAndSmearingTool;
    m_muonCalibrationAndSmearingTool = 0;
  }
  /// Muon selector tool
  if(m_muonSelection){
    delete m_muonSelection;
    m_muonSelection = 0;
  }
  //~ /// Muon eff.corr. tool
  //~ if(m_effi_corr){
    //~ delete m_effi_corr;
    //~ m_effi_corr = 0;
  //~ }
  /*
  if(*m_METUtil){
      delete m_METUtil;
      m_METUtil = 0;
    }
  */
  
  if(m_useHistObjectDumper && m_HistObjectDumper){
    delete m_HistObjectDumper;
    m_HistObjectDumper = 0;
  }
  
  if(m_useBitsetCutflow && m_BitsetCutflow){
    delete m_BitsetCutflow;
    m_BitsetCutflow = 0;
  }
  
  /// cleaning up trigger tools
  if( m_trigConfigTool ) {
      delete m_trigConfigTool;
      m_trigConfigTool = 0;
   }
   if( m_trigDecisionTool ) {
      delete m_trigDecisionTool;
      m_trigDecisionTool = 0;
   }

  return EL::StatusCode::SUCCESS;
}



EL::StatusCode MyxAODAnalysis :: histFinalize ()
{
  /// This method is the mirror image of histInitialize(), meaning it
  /// gets called after the last event has been processed on the worker
  /// node and allows you to finish up any objects you created in
  /// histInitialize() before they are written to disk.  This is
  /// actually fairly rare, since this happens separately for each
  /// worker node.  Most of the time you want to do your
  /// post-processing on the submission node after all your histogram
  /// outputs have been merged.  This is different from finalize() in
  /// that it gets called on all worker nodes regardless of whether
  /// they processed input events.
  return EL::StatusCode::SUCCESS;
}


ConstDataVector<xAOD::MuonContainer> MyxAODAnalysis :: VetoMuon(const xAOD::MuonContainer* muons,
                                         xAOD::Vertex* primVertex,
                                  bool lookForVetoMuon){

  // option 1: signal selection; 2: muon veto; 3: 
  // would need to fill 
  // ConstDataVector<xAOD::MuonContainer> muonfinal(SG::VIEW_ELEMENTS); // Create a new muon container
  const char* APP_NAME = "MyxAODAnalysis";
  //nmuons = 0;

  const xAOD::EventInfo* eventInfo = 0;
  if( ! m_event->retrieve( eventInfo, "EventInfo").isSuccess() ){
    Error("execute()", "Failed to retrieve event info collection. Exiting." );
  }

  //  xAOD::MuonContainer* outMuonContainer = 0;

  ConstDataVector<xAOD::MuonContainer> vetoMuonContainer(SG::VIEW_ELEMENTS);
  xAOD::Muon* vetoMuon = 0;
  int muonCounter = 0;
  int nmuon=0;
  //static SG::AuxElement::Accessor< xAOD::IParticleLink > accSetOriginLink ("originalObjectLink");
  xAOD::MuonContainer::const_iterator muon_itr = muons->begin();
  xAOD::MuonContainer::const_iterator muon_end = muons->end();
  if (!lookForVetoMuon) {
    for( ; muon_itr != muon_end; ++muon_itr ) {
      xAOD::Muon* mu = 0;
      mu = const_cast<xAOD::Muon*> (*muon_itr);
      //cout << "muon type " << mu->muonType()<<endl;
      if (mu->muonType()==xAOD::Muon_v1::Combined) nmuon+=1;
    }
    //cout << "nmuon " << nmuon<< " size " << muons->size()<<endl;
    if (nmuon>0) {
      count[3]+=1;
    } else {
      return vetoMuonContainer;
    }
  }
muon_itr = muons->begin();
  for( ; muon_itr != muon_end; ++muon_itr ) {

    if (m_useHistObjectDumper) m_HistObjectDumper->plotMuon((*muon_itr),"noCuts");

    m_BitsetCutflow->FillCutflow("oneMuon",!lookForVetoMuon);

    xAOD::Muon* mu = 0;
    mu = const_cast<xAOD::Muon*> (*muon_itr);

    // keep lins

    if (mu->muonType()!=xAOD::Muon_v1::Combined) continue;
    m_BitsetCutflow->FillCutflow("Combined",!lookForVetoMuon);
    if (!lookForVetoMuon) count[4]+=1;
    double muPt = (mu->pt()) * 0.001;
    double lowPtCut = 55.0; /// GeV
    //double lowPtCut = 20.0; /// GeV
    double highPtCut = 9999999.9; /// GeV
    if (lookForVetoMuon){
      lowPtCut = 20.0;
      highPtCut = 999999.0;
    }

    if (muPt < lowPtCut || muPt>=highPtCut) {
      continue; /// veto muon
      //cout << "muon rejected in pt cut " << muPt<<endl;
    }
    m_BitsetCutflow->FillCutflow("mu_pt",!lookForVetoMuon);
    if (!lookForVetoMuon) count[5]+=1;

    if (!lookForVetoMuon) {
      if(!m_muonSelection->accept(mu)) continue;
      m_BitsetCutflow->FillCutflow("MCP selector",!lookForVetoMuon);
      count[6]+=1;
    } else {
      if(!m_loosemuonSelection->accept(mu)) continue;
    }
    /// do significance 
    double d0_sig = TMath::Abs(mu->primaryTrackParticle()->d0()) /
    TMath::Sqrt(mu->primaryTrackParticle()->definingParametersCovMatrix()(0,0)
    + eventInfo->beamPosSigmaX()*eventInfo->beamPosSigmaX() );
    if (d0_sig>3.0) continue;
    m_BitsetCutflow->FillCutflow("d0",!lookForVetoMuon);
    if (!lookForVetoMuon) count[8]+=1;


/// zo cut
    double z0_vrtPVx = mu->primaryTrackParticle()->z0() +
    mu->primaryTrackParticle()->vz() - primVertex->z();
    double sintheta = 1.0/TMath::CosH(mu->eta());
    if (abs( z0_vrtPVx*sintheta )>10.0) continue;
    m_BitsetCutflow->FillCutflow("z0",!lookForVetoMuon);
    if (!lookForVetoMuon) count[9]+=1;

    /// Isolation stuff

    if (!lookForVetoMuon) {
      if (!m_muonisolationSelectionTool->accept(*mu)) continue;
      m_BitsetCutflow->FillCutflow("Isolation",!lookForVetoMuon);
      count[10]+=1;
    }

    /// check truth
    if(isMC){

    }
    // save this muon
    //const xAOD::IParticleLink originLink( *muons, mu->index() );
    //accSetOriginLink(*mu) = originLink;
    vetoMuonContainer.push_back(mu);
    //cout << "veto muon found with pt " << mu->pt()/1000.<<endl;
    vetoMuon = mu;
    muonCounter++;

  } /// end for loop over muons

  /// if there are two or more signal muons - return NULL
  /// Veto muon is defined to have pT within 20..55 GeV range
  /// that's why veto muon selection is not overlapping with signal region
  /// if more than two veto muons found - return any of them
  //  if (muonCounter==1 || lookForVetoMuon) 
  return vetoMuonContainer;
//else
//  return 0;
}



ConstDataVector<xAOD::MuonContainer> MyxAODAnalysis :: SelectMuon(const xAOD::MuonContainer* muons,
                                         xAOD::Vertex* primVertex,
                                  bool lookForVetoMuon){

  // option 1: signal selection; 2: muon veto; 3: 
  // would need to fill 
  //ConstDataVector<xAOD::MuonContainer> muonfinal(SG::VIEW_ELEMENTS); // Create a new muon container
  const char* APP_NAME = "MyxAODAnalysis";
  //nmuons = 0;

  const xAOD::EventInfo* eventInfo = 0;
  if( ! m_event->retrieve( eventInfo, "EventInfo").isSuccess() ){
    Error("execute()", "Failed to retrieve event info collection. Exiting." );
  }

  //  xAOD::MuonContainer* outMuonContainer = 0;

  ConstDataVector<xAOD::MuonContainer> outMuonContainer(SG::VIEW_ELEMENTS);
  xAOD::Muon* outMuon = 0;
  int muonCounter = 0;
  int nmuon=0;
  static SG::AuxElement::Accessor< xAOD::IParticleLink > accSetOriginLink ("originalObjectLink");
  xAOD::MuonContainer::const_iterator muon_itr = muons->begin();
  xAOD::MuonContainer::const_iterator muon_end = muons->end();
  if (!lookForVetoMuon) {
    for( ; muon_itr != muon_end; ++muon_itr ) {
      xAOD::Muon* mu = 0;
      mu = const_cast<xAOD::Muon*> (*muon_itr);
      //cout << "muon type " << mu->muonType()<<endl;
      if (mu->muonType()==xAOD::Muon_v1::Combined) nmuon+=1;
    }
    //cout << "nmuon " << nmuon<< " size " << muons->size()<<endl;
    if (nmuon>0) {
      count[3]+=1;
    } else {
      return outMuonContainer;
    }
  }


  muon_itr = muons->begin();
  for( ; muon_itr != muon_end; ++muon_itr ) {

    if (m_useHistObjectDumper) m_HistObjectDumper->plotMuon((*muon_itr),"noCuts");

    m_BitsetCutflow->FillCutflow("oneMuon",!lookForVetoMuon);

    xAOD::Muon* mu = 0;
    /// TODO return m_useMuonCalibrationAndSmearingTool flag
//     if (m_useMuonCalibrationAndSmearingTool){
//       if( !m_muonCalibrationAndSmearingTool->correctedCopy( **muon_itr, mu ) ) {
//         Error(APP_NAME, "Cannot really apply calibration nor smearing");
//         continue;
//       }
//     }
//     else{
//       mu = const_cast<xAOD::Muon*> (*muon_itr);
//     }
    mu = const_cast<xAOD::Muon*> (*muon_itr);

    // keep lins

    if (mu->muonType()!=xAOD::Muon_v1::Combined) continue;
    m_BitsetCutflow->FillCutflow("Combined",!lookForVetoMuon);
    if (!lookForVetoMuon) count[4]+=1;
    double muPt = (mu->pt()) * 0.001;
    //double lowPtCut = 55.0; /// GeV
    double lowPtCut = 55.0; /// GeV
    double highPtCut = 99999.9; /// GeV
    if (lookForVetoMuon){
      lowPtCut = 20.0;
      highPtCut = 999999.0;
    }

    if (muPt < lowPtCut || muPt>=highPtCut) continue; /// veto muon
    m_BitsetCutflow->FillCutflow("mu_pt",!lookForVetoMuon);
    if (!lookForVetoMuon) count[5]+=1;

    if (!lookForVetoMuon) {
      if(!m_muonSelection->accept(mu)) continue;
      m_BitsetCutflow->FillCutflow("MCP selector",!lookForVetoMuon);
      count[6]+=1;
    } else {
      if(!m_loosemuonSelection->accept(mu)) continue;
    }
    /// do significance 
    double d0_sig = TMath::Abs(mu->primaryTrackParticle()->d0()) /
    TMath::Sqrt(mu->primaryTrackParticle()->definingParametersCovMatrix()(0,0)
    + eventInfo->beamPosSigmaX()*eventInfo->beamPosSigmaX() );
    if (d0_sig>3.0) continue;
    m_BitsetCutflow->FillCutflow("d0",!lookForVetoMuon);
    if (!lookForVetoMuon) count[8]+=1;

    /// zo cut
    double z0_vrtPVx = mu->primaryTrackParticle()->z0() +
    mu->primaryTrackParticle()->vz() - primVertex->z();
    double sintheta = 1.0/TMath::CosH(mu->eta());
    if (abs( z0_vrtPVx*sintheta )>10.0) continue;
    m_BitsetCutflow->FillCutflow("z0",!lookForVetoMuon);
    if (!lookForVetoMuon) count[9]+=1;

    /// Isolation stuff

    if (!lookForVetoMuon) {
      if (!m_muonisolationSelectionTool->accept(*mu)) continue;
      m_BitsetCutflow->FillCutflow("Isolation",!lookForVetoMuon);
      count[10]+=1;
    }

    // check truth
    if(isMC){
    }
    // save this muon
    const xAOD::IParticleLink originLink( *muons, mu->index() );
    accSetOriginLink(*mu) = originLink;
    outMuonContainer.push_back(mu);
    //cout << "signal muon found with pt " << mu->pt()/1000.<<endl;
    outMuon = mu;
    muonCounter++;

  } /// end for loop over muons

  /// if there are two or more signal muons - return NULL
  /// Veto muon is defined to have pT within 20..55 GeV range
  /// that's why veto muon selection is not overlapping with signal region
  /// if more than two veto muons found - return any of them
  //  if (muonCounter==1 || lookForVetoMuon) 
  return outMuonContainer;
//else
//  return 0;
}

bool MyxAODAnalysis :: passMuonSelection(const xAOD::Muon* mu,
                     const xAOD::EventInfo* eventInfo,
                                         xAOD::Vertex* primVertex,
                     bool lookForVetoMuon){

  const char* APP_NAME = "MyxAODAnalysis";
  bool pass = false;

  if (mu->muonType()!=xAOD::Muon_v1::Combined) return pass;
  m_BitsetCutflow->FillCutflow("Combined",!lookForVetoMuon);
  if (!lookForVetoMuon) count[4]+=1;
  double muPt = (mu->pt()) * 0.001;
  //double lowPtCut = 55.0; /// GeV
  double lowPtCut = 20.0; /// GeV
  double highPtCut = 99999.9; /// GeV
  if (lookForVetoMuon){
    lowPtCut = 20.0;
    highPtCut = 999999.0;
  }

  if (muPt < lowPtCut || muPt>=highPtCut) return pass; /// veto muon
  m_BitsetCutflow->FillCutflow("mu_pt",!lookForVetoMuon);
  if (!lookForVetoMuon) count[5]+=1;

  if (!lookForVetoMuon) {
    if(!m_muonSelection->accept(mu)) return pass;
    m_BitsetCutflow->FillCutflow("MCP selector",!lookForVetoMuon);
    count[6]+=1;
  } else {
    if(!m_loosemuonSelection->accept(mu)) return pass;
  }
  /// do significance 
  double d0_sig = TMath::Abs(mu->primaryTrackParticle()->d0()) /
    TMath::Sqrt(mu->primaryTrackParticle()->definingParametersCovMatrix()(0,0)
        + eventInfo->beamPosSigmaX()*eventInfo->beamPosSigmaX() );
  if (d0_sig>3.0) return pass;
  m_BitsetCutflow->FillCutflow("d0",!lookForVetoMuon);
  if (!lookForVetoMuon) count[8]+=1;

  /// zo cut
  double z0_vrtPVx = mu->primaryTrackParticle()->z0() +
    mu->primaryTrackParticle()->vz() - primVertex->z();
  double sintheta = 1.0/TMath::CosH(mu->eta());
  if (abs( z0_vrtPVx*sintheta )>10.0) return pass;
  m_BitsetCutflow->FillCutflow("z0",!lookForVetoMuon);
  if (!lookForVetoMuon) count[9]+=1;
  /// Isolation stuff

  if (!lookForVetoMuon) {
    if (!m_muonisolationSelectionTool->accept(*mu)) return pass;
    m_BitsetCutflow->FillCutflow("Isolation",!lookForVetoMuon);
    count[10]+=1;
  }

  // check truth
  if(isMC){
  }
  
  pass = true;
  //cout << "signal muon found with pt " << mu->pt()/1000.<<endl;

  return pass;
}

ConstDataVector<xAOD::ElectronContainer> MyxAODAnalysis :: SelectElectron(const xAOD::ElectronContainer* electrons, xAOD::Vertex* primVertex, bool lookForVetoElectron){

  const char* APP_NAME = "MyxAODAnalysis";

    const xAOD::EventInfo* eventInfo = 0;
  if( ! m_event->retrieve( eventInfo, "EventInfo").isSuccess() ){
    Error("execute()", "Failed to retrieve event info collection. Exiting." );
  }

  //cout << "in selectElectron"<<endl;
  xAOD::Electron* outElectron = 0;
  ConstDataVector<xAOD::ElectronContainer> outElectronContainer(SG::VIEW_ELEMENTS);
  int electronCounter = 0;
  xAOD::ElectronContainer::const_iterator electron_itr = electrons->begin();
  xAOD::ElectronContainer::const_iterator electron_end = electrons->end();
  int i = 0;
  for( ; electron_itr != electron_end; ++electron_itr ) {

    m_BitsetCutflow->FillCutflow("oneElectron",!lookForVetoElectron);

    /// Calibration and corrections ...
    xAOD::Electron* el = 0;
    xAOD::Electron* el2 = 0;
    el2 = const_cast<xAOD::Electron*> (*electron_itr);
    float oldpt = el2->pt()/1000.;

    el = const_cast<xAOD::Electron*> (*electron_itr);
      //}

    i++;
    /// ... Calibration and corrections
    /// Eta
    double Eta = el->caloCluster()->eta();
    if ( abs(Eta) > 2.47 || (abs(Eta) > 1.37 && abs(Eta) < 1.52)) continue;
    m_BitsetCutflow->FillCutflow("Eta",!lookForVetoElectron);

    /// OQ
    if (el->isGoodOQ(xAOD::EgammaParameters::BADCLUSELECTRON)==false) continue;
    m_BitsetCutflow->FillCutflow("OQ",!lookForVetoElectron);

    /// pT cut ...
    double elPt = (el->pt()) * 0.001;
    double lowPtCut = 65.0; /// GeV
    if (lookForVetoElectron){
      lowPtCut = 20.0;
    }
//cout << "pt " << elPt<< " old " << oldpt<<endl;
    if (elPt < lowPtCut) continue; /// veto electron
    m_BitsetCutflow->FillCutflow("el_pt",!lookForVetoElectron);
    /// ... pT cut

    /// d0 significance ...

    const xAOD::TrackParticle *trk = el->trackParticle();
    float d0_sig =  TMath::Abs(trk->d0())/TMath::Sqrt(trk->definingParametersCovMatrix()(0,0)
    + eventInfo->beamPosSigmaX()*eventInfo->beamPosSigmaX() );

   if (d0_sig>5.0) continue;
    m_BitsetCutflow->FillCutflow("d0",!lookForVetoElectron);
    /// ... d0 significance

    /// ID ...
    bool goodID = false;
    if (lookForVetoElectron)
      goodID = m_LHToolMedium2015->accept(el) || (m_LHToolTight2015->accept(el));
    //goodID = m_LHToolMedium2015->accept(el);
    else
      goodID = m_LHToolMedium2015->accept(el) &&
      m_LHToolTight2015->accept(el);

    if (!goodID) continue;
    m_BitsetCutflow->FillCutflow("ID",!lookForVetoElectron);
    /// ... ID

    /// Isolation
    if (!m_eleisolationSelectionTool->accept(*el)) continue;
    m_BitsetCutflow->FillCutflow("Isolation",!lookForVetoElectron);

    outElectron = el;
    outElectronContainer.push_back(el);
    electronCounter++;

  } /// end for loop over electrons

  /// if there are two or more signal electrons - return NULL
  /// Veto electron is defined to have pT within 20..55 GeV range
  /// that's why veto electron selection is not overlapping with signal region
  /// if more than two veto electrons found - return any of them
  return outElectronContainer;

}



/*

xAOD::Muon* MyxAODAnalysis :: SelectMuon(const xAOD::MuonContainer* muons, 
                                         xAOD::Vertex* primVertex, 
                                         bool lookForVetoMuon){

  const char* APP_NAME = "MyxAODAnalysis";
  
  const xAOD::EventInfo* eventInfo = 0;
  if( ! m_event->retrieve( eventInfo, "EventInfo").isSuccess() ){
    Error("execute()", "Failed to retrieve event info collection. Exiting." );
  }
  
  xAOD::Muon* outMuon = 0;
  int muonCounter = 0;
  xAOD::MuonContainer::const_iterator muon_itr = muons->begin();
  xAOD::MuonContainer::const_iterator muon_end = muons->end();
  for( ; muon_itr != muon_end; ++muon_itr ) {
  
    if (m_useHistObjectDumper) m_HistObjectDumper->plotMuon((*muon_itr),"noCuts");
  
    m_BitsetCutflow->FillCutflow("oneMuon",!lookForVetoMuon);
  
    xAOD::Muon* mu = 0;
    if (m_useCalibrationAndSmearingTool){
      if( !m_muonCalibrationAndSmearingTool->correctedCopy( **muon_itr, mu ) ) {
        Error(APP_NAME, "Cannot really apply calibration nor smearing");
        continue;
      }
    }
    else{
      mu = const_cast<xAOD::Muon*> (*muon_itr);
    }
    
    if (mu->muonType()!=xAOD::Muon_v1::Combined) continue;
    m_BitsetCutflow->FillCutflow("Combined",!lookForVetoMuon);
    
    double muPt = (mu->pt()) * 0.001;
    double lowPtCut = 55.0; /// GeV
    double highPtCut = 99999.9; /// GeV
    if (lookForVetoMuon){
      lowPtCut = 20.0;
      highPtCut = 55.0;
    }
      
    if (muPt < lowPtCut || muPt>=highPtCut) continue; /// veto muon
    m_BitsetCutflow->FillCutflow("mu_pt",!lookForVetoMuon);
    
    if(!m_muonSelection->accept(mu)) continue;
    m_BitsetCutflow->FillCutflow("MCP selector",!lookForVetoMuon);
        
    if (!m_muonSelection->passedHighPtCuts(*mu)) continue;
    m_BitsetCutflow->FillCutflow("MS Hits",!lookForVetoMuon);
    
    /// do significance 
    double d0_sig = TMath::Abs(mu->primaryTrackParticle()->d0()) / 
    TMath::Sqrt(mu->primaryTrackParticle()->definingParametersCovMatrix()(0,0)
    + eventInfo->beamPosSigmaX()*eventInfo->beamPosSigmaX() );
    if (d0_sig>3.0) continue;
    m_BitsetCutflow->FillCutflow("d0",!lookForVetoMuon);
    
    /// zo cut
    double z0_vrtPVx = mu->primaryTrackParticle()->z0() + 
    mu->primaryTrackParticle()->vz() - primVertex->z(); 
    double sintheta = 1.0/TMath::CosH(mu->eta());
    if (abs( z0_vrtPVx*sintheta )>10.0) continue;
    m_BitsetCutflow->FillCutflow("z0",!lookForVetoMuon);
    
    /// Isolation stuff
    if (!m_isolationSelectionTool->accept(*mu)) continue;
    m_BitsetCutflow->FillCutflow("Isolation",!lookForVetoMuon);
    
    outMuon = mu;
    muonCounter++;
      
    
  
  } /// end for loop over muons

  /// if there are two or more signal muons - return NULL
  /// Veto muon is defined to have pT within 20..55 GeV range
  /// that's why veto muon selection is not overlapping with signal region
  /// if more than two veto muons found - return any of them
  if (muonCounter==1 || lookForVetoMuon) return outMuon;
  else
    return 0;

}



xAOD::Electron* MyxAODAnalysis :: SelectElectron(const xAOD::ElectronContainer* electrons, 
                                         xAOD::Vertex* primVertex, 
                                         bool lookForVetoElectron){

  const char* APP_NAME = "MyxAODAnalysis";
  
  const xAOD::EventInfo* eventInfo = 0;
  if( ! m_event->retrieve( eventInfo, "EventInfo").isSuccess() ){
    Error("execute()", "Failed to retrieve event info collection. Exiting." );
  }
  
  xAOD::Electron* outElectron = 0;
  int electronCounter = 0;
  xAOD::ElectronContainer::const_iterator electron_itr = electrons->begin();
  xAOD::ElectronContainer::const_iterator electron_end = electrons->end();
  for( ; electron_itr != electron_end; ++electron_itr ) {
  
    m_BitsetCutflow->FillCutflow("oneElectron",!lookForVetoElectron);
  
    /// Calibration and corrections ...
    xAOD::Electron* el = 0;
    if (m_useCalibrationAndSmearingTool){
//       if( !m_electronCalibrationAndSmearingTool->correctedCopy( **electron_itr, el ) ) {
//         Error(APP_NAME, "Cannot really apply calibration nor smearing");
//         continue;
//       }
    }
    else{
      el = const_cast<xAOD::Electron*> (*electron_itr);
    }
    /// ... Calibration and corrections
    
    /// Eta
    double Eta = el->caloCluster()->eta();
    if ( abs(Eta) > 2.47 || (1.37 > abs(Eta) < 1.52)) continue;
    m_BitsetCutflow->FillCutflow("Eta",!lookForVetoElectron);
    
    /// OQ
    if (el->isGoodOQ(xAOD::EgammaParameters::BADCLUSELECTRON)==false) continue;
    m_BitsetCutflow->FillCutflow("OQ",!lookForVetoElectron);
    
    /// pT cut ...
    double elPt = (el->pt()) * 0.001;
    double lowPtCut = 65.0; /// GeV
    if (lookForVetoElectron){
      lowPtCut = 20.0;
    }
    if (elPt < lowPtCut) continue; /// veto electron
    m_BitsetCutflow->FillCutflow("el_pt",!lookForVetoElectron);
    /// ... pT cut
    
    /// d0 significance ...
    
    const xAOD::TrackParticle *trk = el->trackParticle();
    float d0_sig = fabs(trk->d0())/sqrt(trk->definingParametersCovMatrix()(0,0)
    + eventInfo->beamPosSigmaX()*eventInfo->beamPosSigmaX() );
    if (d0_sig>5.0) continue;
    m_BitsetCutflow->FillCutflow("d0",!lookForVetoElectron);
    /// ... d0 significance
    
    /// ID ...
    bool goodID = false;
    if (lookForVetoElectron)
      goodID = m_LHToolMedium2015->accept(el) && (!m_LHToolTight2015->accept(el));
    else
      goodID = m_LHToolMedium2015->accept(el) && 
      m_LHToolTight2015->accept(el);
    
    if (!goodID) continue;
    m_BitsetCutflow->FillCutflow("ID",!lookForVetoElectron);
    /// ... ID
    
    /// Isolation
    if (!m_isolationSelectionTool->accept(*el)) continue;
    m_BitsetCutflow->FillCutflow("Isolation",!lookForVetoElectron);
    
    outElectron = el;
    electronCounter++;
    
  } /// end for loop over electrons

  /// if there are two or more signal electrons - return NULL
  /// Veto electron is defined to have pT within 20..55 GeV range
  /// that's why veto electron selection is not overlapping with signal region
  /// if more than two veto electrons found - return any of them
  if (electronCounter==1 || lookForVetoElectron) return outElectron;
  else
    return 0;

}
*/

