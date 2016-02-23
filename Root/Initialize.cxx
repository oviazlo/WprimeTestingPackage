///*****************************************************************************
///
/// Implementation of all functions which are called BEFORE execute:
///     - constructor RecoAnalysis
///     - initialize:
///             * retrieve and configure all tools
///     - setupJob:
///             * initialization of flags
///     - fileExecute
///     - changeInput
///
///*****************************************************************************
#include <MyAnalysis/RecoAnalysis.h>

/// this is needed to distribute the algorithm to the workers
ClassImp(RecoAnalysis)

EL::StatusCode RecoAnalysis :: initialize ()
{
  /// Here you do everything that you need to do after the first input
  /// file has been connected and before the first event is processed,
  /// e.g. create additional histograms based on which variables are
  /// available in the input files.  You can also create all of your
  /// histograms and trees in here, but be aware that this method
  /// doesn't get called if no events are processed.  So any objects
  /// you create here won't be available in the output if you have no
  /// input events.

  const char* APP_NAME = "RecoAnalysis";

  m_event = wk()->xaodEvent();
  m_store = new xAOD::TStore();
  
  /// as a check, let's see the number of events in our xAOD
  /// print long long int
  Info("initialize()", "Number of events = %lli", m_event->getEntries() ); 

  /// Event Info
  const xAOD::EventInfo* eventInfo = 0;
  if( ! m_event->retrieve( eventInfo, "EventInfo").isSuccess() ){
    Error("execute()", 
          "Failed to retrieve event info collection in initialise. Exiting." );
    return EL::StatusCode::FAILURE;
  }  
  
  /// TODO Initialize all flag
  m_isMC = false;
  if(eventInfo->eventType( xAOD::EventInfo::IS_SIMULATION ) ){
    m_isMC = true;
  }
  
  /// count number of events
  m_eventCounter = 0;

  /// GRL
  m_grl = new GoodRunsListSelectionTool("GoodRunsListSelectionTool");
  std::vector<std::string> vecStringGRL; 
  vecStringGRL.push_back(
    "$ROOTCOREBIN/data/MyAnalysis/data15_13TeV.periodAllYear_DetStatus"
    "-v73-pro19-08_DQDefects-00-01-02_PHYS_StandardGRL_All_Good_25ns.xml");

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
  /// Recommendation taken from:
  /// https://twiki.cern.ch/twiki/bin/view/AtlasProtected/ApplyJetCalibration2015
  cout << "initialise jet calib"<<endl;
  /// String describing the current thread, for logging
  const std::string name = "RecoAnalysis"; 
  /// String describing your jet collection, 
  /// for example AntiKt4EMTopo or AntiKt4EMTopo (see above)
  TString jetAlgo = "AntiKt4EMTopo"; 
  /// Path to global config used to initialize the tool (see above)
//   TString config = "JES_MC15Prerecommendation_April2015.config"; 
  TString config = "JES_2015dataset_recommendation_Feb2016.config";
  /// String describing the calibration sequence to apply (see above)
  TString calibSeq ="JetArea_Residual_Origin_EtaJES_GSC" ; 
  
  m_jetCalibrationTool = 
  new JetCalibrationTool(name, jetAlgo, config, calibSeq, !m_isMC);
//   CHECK (m_jetCalibrationTool->setProperty("JetCollection",jetAlgo));
//   CHECK (m_jetCalibrationTool->setProperty("ConfigFile",config));
//   CHECK (m_jetCalibrationTool->setProperty("CalibSequence",calibSeq));
//   CHECK (m_jetCalibrationTool->setProperty("IsData",!m_isMC));
  
  EL_RETURN_CHECK("initialize jetCalibrationTool",
                  m_jetCalibrationTool->initializeTool(name));
  
  m_metMaker = new met::METMaker("METMakerTool");
  //m_metMaker->msg().setLevel( MSG::DEBUG ); // or DEBUG or VERBOSE
  if ( ! m_metMaker->initialize().isSuccess() ){
    Error("initialize()", "Failed to properly initialize the METMakerTool. "
    "Exiting." );
    return EL::StatusCode::FAILURE;
  }
  
  /// Configure the JVT tool.
  m_jvtTool = 0;
  ToolHandle<IJetUpdateJvt> hjvtagup;
  m_jvtTool = new JetVertexTaggerTool("jvtag");
  hjvtagup = ToolHandle<IJetUpdateJvt>("jvtag");
  m_jvtTool->msg().setLevel( MSG::DEBUG ); // or DEBUG or VERBOSE
  CHECK (m_jvtTool->setProperty("JVTFileName","JetMomentTools/"
  "JVTlikelihood_20140805.root")); /// TODO is it the latest file?
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
  /// Source: https://twiki.cern.ch/twiki/bin/view/Atlas/MuonSelectionTool
  m_muonSelection = new CP::MuonSelectionTool("MuonSelection");
  // m_muonSelection->msg().setLevel( MSG::INFO );
  m_muonSelection->setProperty( "MaxEta", 2.5 );
  m_muonSelection->setProperty( "MuQuality", 4); /// according to source: Tight -> MuQuality=0. But then cutflow doesn't match. While MuQuality=4 leads to matching of cutflows... FIXME TODO Cross check with Magnard
  m_muonSelection->msg().setLevel( MSG::ERROR );
  CHECK (m_muonSelection->initialize().isSuccess());

  m_looseMuonSelection = new CP::MuonSelectionTool("MuonLooseSelection");
  // m_looseMuonSelection->msg().setLevel( MSG::INFO );
  m_looseMuonSelection->setProperty( "MaxEta", 2.5 );
  m_looseMuonSelection->setProperty( "MuQuality", 1); /// WARNING before was 2
  m_looseMuonSelection->msg().setLevel( MSG::ERROR );
  CHECK (m_looseMuonSelection->initialize().isSuccess());
  
  
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
  /// gives us access to the meta-data
  m_trigConfigTool = new TrigConf::xAODConfigTool("xAODConfigTool"); 
  EL_RETURN_CHECK( "initialize", m_trigConfigTool->initialize() );
  ToolHandle< TrigConf::ITrigConfigTool > trigConfigHandle( m_trigConfigTool );
  m_trigDecisionTool = new Trig::TrigDecisionTool("TrigDecisionTool");
  /// connect the TrigDecisionTool to the ConfigTool
  EL_RETURN_CHECK( "initialize", m_trigDecisionTool->setProperty( 
  "ConfigTool", trigConfigHandle ) ); 
  EL_RETURN_CHECK( "initialize", m_trigDecisionTool->setProperty( 
  "TrigDecisionKey", "xTrigDecision" ) );
  EL_RETURN_CHECK( "initialize", m_trigDecisionTool->initialize() );
  
  /// electrons
  /// Recommendation are taken from:
  /// https://twiki.cern.ch/twiki/bin/view/AtlasProtected/ElectronPhotonFourMomentumCorrection
  cout << "do e/g calibration" << endl;
  m_eleCalibrationTool = new CP::EgammaCalibrationAndSmearingTool("eletool");
  /// see below for options
  m_eleCalibrationTool->setProperty("ESModel", "es2015PRE"); /// seems OK
  /// see below for options
  /// FIXME what to use here??? Ask Magnar/Markus!!!
  m_eleCalibrationTool->setProperty("decorrelationModel", "1NP_v1");  
  m_eleCalibrationTool->initialize();
  
  cout << "e/g calibration done" << endl;
  
  /// LH electron identification ...
  m_LHToolTight2015    = new AsgElectronLikelihoodTool ("m_LHToolTight2015");
  m_LHToolMedium2015   = new AsgElectronLikelihoodTool ("m_LHToolMedium2015"); 
  
  EL_RETURN_CHECK( "initialize", m_LHToolTight2015->setProperty(
    "primaryVertexContainer","PrimaryVertices"));
  EL_RETURN_CHECK( "initialize", m_LHToolMedium2015->setProperty(
    "primaryVertexContainer","PrimaryVertices"));

  /// FIXME do we need to use mc15_20160113 instead?
  /// Source: https://twiki.cern.ch/twiki/bin/view/AtlasProtected/EGammaIdentificationRun2#Electron_d0_and_z0_cut_definitio  
  std::string confDir = "ElectronPhotonSelectorTools/offline/mc15_20151012/";
  EL_RETURN_CHECK( "m_LHToolTight2015 setProperty ConfigFile", 
                   m_LHToolTight2015->setProperty(
    "ConfigFile",confDir+"ElectronLikelihoodTightOfflineConfig2015.conf"));
  EL_RETURN_CHECK( "m_LHToolMedium2015  setProperty ConfigFile", 
                   m_LHToolMedium2015->setProperty(
    "ConfigFile",confDir+"ElectronLikelihoodMediumOfflineConfig2015.conf"));
  
  EL_RETURN_CHECK( "initialize", m_LHToolTight2015->initialize() );
  EL_RETURN_CHECK( "initialize", m_LHToolMedium2015->initialize() );
  /// ... LH electron identification
  cout << "e/g LH done"<<endl;
  
  /// Muon Isolation Tool
  m_muonisolationSelectionTool = new CP::IsolationSelectionTool("muoniso");
  EL_RETURN_CHECK( "initialize",m_muonisolationSelectionTool
  ->setProperty( "MuonWP","LooseTrackOnly") );
  EL_RETURN_CHECK( "initialize",m_muonisolationSelectionTool->initialize());
  m_eleisolationSelectionTool = new CP::IsolationSelectionTool("eleiso");
  EL_RETURN_CHECK( "initialize",m_eleisolationSelectionTool
  ->setProperty( "ElectronWP","Loose") );
  EL_RETURN_CHECK( "initialize",m_eleisolationSelectionTool->initialize());
  
  if (m_useHistObjectDumper)
    m_HistObjectDumper = new HistObjectDumper(wk());
  
  if (m_BitsetCutflow)
    m_BitsetCutflow = new BitsetCutflow(wk());
 
  weightkFactor = 1.0;
  weighfilterEfficiency = 1.0;
  weightCrossSection = 1.0;
  
  m_LPXKfactorTool = NULL;
  if (m_isMC){
    m_LPXKfactorTool = new LPXKfactorTool("LPXKfactorTool");
    CHECK(m_LPXKfactorTool->setProperty("isMC15", true)); 
    CHECK(m_LPXKfactorTool->setProperty("applyEWCorr", true)); 
    CHECK(m_LPXKfactorTool->setProperty("applyPICorr", true)); 
    
    EL_RETURN_CHECK( "m_LPXKfactorTool initialize",m_LPXKfactorTool->initialize());
  }
  
  /// Source: https://twiki.cern.ch/twiki/bin/view/AtlasProtected/ExtendedPileupReweighting
  /// Example: https://twiki.cern.ch/twiki/bin/view/AtlasProtected/HiggsZZllllPreparationRunII2014#Data_Period_Fraction_and_PileUpR
  /// FIXME do I have to extract weight as below?
  /// weight *= eventInfo->auxdata< double >("PileupWeight");
  m_pileupReweightingTool = new CP::PileupReweightingTool("PileupReweightingTool");
  std::vector<std::string> confFiles;
  std::vector<std::string> lcalcFiles;
  confFiles.push_back("$ROOTCOREBIN/data/MyAnalysis/prwFileWplusenuMC15b_25ns.root");
  lcalcFiles.push_back("$ROOTCOREBIN/data/MyAnalysis/ilumicalc_histograms_None_276262-284484.root");
  m_pileupReweightingTool->setProperty("ConfigFiles",confFiles);
  m_pileupReweightingTool->setProperty("LumiCalcFiles",lcalcFiles); 
  /// FIXME cross check with Markus should one setup DefaultChannel or not
  
//   m_pileupReweightingTool->setProperty("DefaultChannel", 361106).ignore(); /// Monika line
  m_pileupReweightingTool->setProperty("DefaultChannel", 361100).ignore(); /// from ROOT file
  m_pileupReweightingTool->initialize();
  cout << "pileup reweighting initialised" << endl;
  
  m_sampleName = wk()->metaData()->getString ("sample_name");
  string tmpSampleName = "";
  stringstream strStream(m_sampleName);
  getline (strStream, tmpSampleName, '.'); /// return mc15_13TeV
  getline (strStream, tmpSampleName, '.'); /// return DSID
  m_DSID = atoi(tmpSampleName.c_str());
  
  m_inclusiveWorZ = false;
  if (m_DSID>=361100 && m_DSID<=361108){
    m_inclusiveWorZ = true;
    cout << "[WARNING]\tm_cut120GeVForInclusiveW is activated!!!" << endl;
    if (m_DSID>=361106) /// Z sample
      m_pdgIdOfMother = 23; /// FIXME but gamma is 22... 
                            /// do we have gamma in Z inc. samples?
                            /// UPDATE: looks like it's OK to use only pdgId=23
                            /// for all Z samples...
    else /// W sample
      m_pdgIdOfMother = 24;
  }
  
  return EL::StatusCode::SUCCESS;
}

RecoAnalysis :: RecoAnalysis ()
{
  /// Here you put any code for the base initialization of variables,
  /// e.g. initialize all pointers to 0.  Note that you should only put
  /// the most basic initialization here, since this method will be
  /// called on both the submission and the worker node.  Most of your
  /// initialization code will go into histInitialize() and
  /// initialize().
  
  //m_useHistObjectDumper = true;

}



EL::StatusCode RecoAnalysis :: setupJob (EL::Job& job)
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
  xAOD::Init( "RecoAnalysis" ).ignore(); /// call before opening first file
  
  m_useHistObjectDumper = true;
  m_useBitsetCutflow = true;
  m_useCalibrationAndSmearingTool = true;
  m_doWprimeTruthMatching = false;
  m_runElectronChannel = false;
  
  cout << "after init"<<endl;
  
  outputName = "outFile";
  
  EL::OutputStream out (outputName);
  job.outputAdd (out);
  
  return EL::StatusCode::SUCCESS;
}

EL::StatusCode RecoAnalysis :: fileExecute ()
{
  /// Here you do everything that needs to be done exactly once for every
  /// single file, e.g. collect a list of all lumi-blocks processed
  return EL::StatusCode::SUCCESS;
}



EL::StatusCode RecoAnalysis :: changeInput (bool firstFile)
{
  /// Here you do everything you need to do when we change input files,
  /// e.g. resetting branch addresses on trees.  If you are using
  /// D3PDReader or a similarin that jet loop let's fill this histogram, so the 
  /// loop looks something like: service this method is not needed.
  return EL::StatusCode::SUCCESS;
}
