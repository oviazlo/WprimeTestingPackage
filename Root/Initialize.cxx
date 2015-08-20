#include <MyAnalysis/MyxAODAnalysis.h>


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
