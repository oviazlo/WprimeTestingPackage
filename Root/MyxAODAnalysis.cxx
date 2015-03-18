#include <EventLoop/Job.h>
#include <EventLoop/StatusCode.h>
#include <EventLoop/Worker.h>
#include <MyAnalysis/MyxAODAnalysis.h>
//#include <MyAnalysis/HistObjectDumper.h>
#include "CPAnalysisExamples/errorcheck.h"

/// EDM includes:
#include "xAODEventInfo/EventInfo.h"

/// GRL
#include "GoodRunsLists/GoodRunsListSelectionTool.h"

#include "xAODJet/JetContainer.h"
#include "xAODMuon/MuonContainer.h"
#include "JetSelectorTools/JetCleaningTool.h"
#include "JetResolution/JERTool.h"
#include <TSystem.h> /// used to define JERTool calibration path 

/// Muons
#include "MuonSelectorTools/MuonSelectionTool.h"
#include "MuonMomentumCorrections/MuonCalibrationAndSmearingTool.h"

#include "PATInterfaces/CorrectionCode.h" /// to check the return correction code status of tools

//#include "MissingETUtility/METUtility.h" /// METUtils
#include "xAODMissingET/MissingETContainer.h"

#include <TFile.h>
#include <TMath.h>
#include <TLorentzVector.h>

/// this is needed to distribute the algorithm to the workers
ClassImp(MyxAODAnalysis)



MyxAODAnalysis :: MyxAODAnalysis ()
{
	// Here you put any code for the base initialization of variables,
	// e.g. initialize all pointers to 0.  Note that you should only put
	// the most basic initialization here, since this method will be
	// called on both the submission and the worker node.  Most of your
	// initialization code will go into histInitialize() and
	// initialize().
	
	//m_useHistObjectDumper = true;

}



EL::StatusCode MyxAODAnalysis :: setupJob (EL::Job& job)
{
	// Here you put code that sets up the job on the submission object
	// so that it is ready to work with your algorithm, e.g. you can
	// request the D3PDReader service or add output files.  Any code you
	// put here could instead also go into the submission script.  The
	// sole advantage of putting it here is that it gets automatically
	// activated/deactivated when you add/remove the algorithm from your
	// job, which may or may not be of value to you.
	
	// let's initialize the algorithm to use the xAODRootAccess package
	job.useXAOD ();
	xAOD::Init( "MyxAODAnalysis" ).ignore(); // call before opening first file

	m_useHistObjectDumper = true;

	return EL::StatusCode::SUCCESS;
}



EL::StatusCode MyxAODAnalysis :: histInitialize ()
{
	// Here you do everything that needs to be done at the very
	// beginning on each worker node, e.g. create histograms and output
	// trees.  This method gets called before any input files are
	// connected.



	h_jetPt = new TH1F("h_jetPt", "h_jetPt", 100, 0, 500); // jet pt [GeV]
	wk()->addOutput (h_jetPt);

	/// Muon Pt histograms:
	h_muPt_uncorr_wSelector = new TH1F("h_muPt_uncorr_wSelector", "h_muPt_uncorr_wSelector", 300, 0, 3000); // jet pt [GeV]
	wk()->addOutput (h_muPt_uncorr_wSelector);

	h_muPt_uncorr_woSelector = new TH1F("h_muPt_uncorr_woSelector", "h_muPt_uncorr_woSelector", 300, 0, 3000); // jet pt [GeV]
	wk()->addOutput (h_muPt_uncorr_woSelector);
	
	h_muPt_corr_wSelector = new TH1F("h_muPt_corr_wSelector", "h_muPt_corr_wSelector", 300, 0, 3000); // jet pt [GeV]
	wk()->addOutput (h_muPt_corr_wSelector);
	
	h_muPt_corr_woSelector = new TH1F("h_muPt_corr_woSelector", "h_muPt_corr_woSelector", 300, 0, 3000); // jet pt [GeV]
	wk()->addOutput (h_muPt_corr_woSelector);

	h_MET_RefFinalFix = new TH1F("h_MET_RefFinalFix", "h_MET_RefFinalFix", 500, 0, 5000); 
	wk()->addOutput (h_MET_RefFinalFix);
	
	h_MET_RefFinalFix_test = new TH1F("h_MET_RefFinalFix_test", "h_MET_RefFinalFix_test", 500, 0, 5000); 
	wk()->addOutput (h_MET_RefFinalFix_test);
	
	h_Mt = new TH1F("h_Mt", "h_Mt", 500, 0, 5000); 
	wk()->addOutput (h_Mt);

	h_Mt_muonPtCut = new TH1F("h_Mt_muonPtCut", "h_Mt_muonPtCut", 500, 0, 5000); 
	wk()->addOutput (h_Mt_muonPtCut);

	// get the output file, create a new TTree and connect it to that output
	// define what branches will go in that tree
	TFile *outputFile = wk()->getOutputFile (outputName);
	tree = new TTree ("tree", "tree");
	tree->SetDirectory (outputFile);
	tree->Branch("EventNumber", &EventNumber);

	return EL::StatusCode::SUCCESS;
}



EL::StatusCode MyxAODAnalysis :: fileExecute ()
{
	// Here you do everything that needs to be done exactly once for every
	// single file, e.g. collect a list of all lumi-blocks processed
	return EL::StatusCode::SUCCESS;
}



EL::StatusCode MyxAODAnalysis :: changeInput (bool firstFile)
{
	// Here you do everything you need to do when we change input files,
	// e.g. resetting branch addresses on trees.  If you are using
	// D3PDReader or a similarin that jet loop let's fill this histogram, so the loop looks something like: service this method is not needed.
	return EL::StatusCode::SUCCESS;
}



EL::StatusCode MyxAODAnalysis :: initialize ()
{
	// Here you do everything that you need to do after the first input
	// file has been connected and before the first event is processed,
	// e.g. create additional histograms based on which variables are
	// available in the input files.  You can also create all of your
	// histograms and trees in here, but be aware that this method
	// doesn't get called if no events are processed.  So any objects
	// you create here won't be available in the output if you have no
	// input events.

	const char* APP_NAME = "MyxAODAnalysis";

	m_event = wk()->xaodEvent();

	// as a check, let's see the number of events in our xAOD
	Info("initialize()", "Number of events = %lli", m_event->getEntries() ); // print long long int

	// count number of events
	m_eventCounter = 0;

	// number of clean events
	m_numCleanEvents = 0;

	// GRL
	m_grl = new GoodRunsListSelectionTool("GoodRunsListSelectionTool");
	std::vector<std::string> vecStringGRL;
	vecStringGRL.push_back("/afs/cern.ch/user/a/atlasdqm/grlgen/All_Good/data12_8TeV.periodAllYear_DetStatus-v61-pro14-02_DQDefects-00-01-00_PHYS_StandardGRL_All_Good.xml");
	CHECK(m_grl->setProperty( "GoodRunsListVec", vecStringGRL));
	CHECK(m_grl->setProperty("PassThrough", false)); // if true (default) will ignore result of GRL and will just pass all events
	if (!m_grl->initialize().isSuccess()) { // check this isSuccess
		Error("initialize()", "Failed to properly initialize the GRL. Exiting." );
		return EL::StatusCode::FAILURE;
	}

	// initialize and configure the jet cleaning tool
	m_jetCleaning = new JetCleaningTool("JetCleaning");
	m_jetCleaning->msg().setLevel( MSG::DEBUG ); 
	CHECK(m_jetCleaning->setProperty( "CutLevel", "MediumBad"));
	m_jetCleaning->initialize();

	// initialize JER 
	const char* jerFilePath = "$ROOTCOREBIN/data/JetResolution/JERProviderPlots_2012.root";
	const char* fullJERFilePath = gSystem->ExpandPathName (jerFilePath);
	m_JERTool = new JERTool("JERTool",fullJERFilePath,"AntiKt4LCTopoJES");

	if (! m_JERTool->initialize().isSuccess() ){
		Error("initialize()", "Failed to properly initialize the JER Tool. Exiting." );
		return EL::StatusCode::FAILURE;
	}

	/// Get Muon Selector Tool
	m_muonSelection = new CP::MuonSelectionTool("MuonSelection");
	m_muonSelection->msg().setLevel( MSG::ERROR );
	m_muonSelection->setProperty( "MaxEta", 2.4 );
	m_muonSelection->setProperty( "MuQuality", 1);
	CHECK (m_muonSelection->initialize().isSuccess());


	/// initialize the muon calibration and smearing tool
	m_muonCalibrationAndSmearingTool = new CP::MuonCalibrationAndSmearingTool( "MuonCorrectionTool" );
	//m_muonCalibrationAndSmearingTool->msg().setLevel( MSG::DEBUG );
	if (! m_muonCalibrationAndSmearingTool->initialize().isSuccess() ){
		Error("initialize()", "Failed to properly initialize the MonCalibrationAndSmearingTool Tool. Exiting." );
		return EL::StatusCode::FAILURE;
	}
	
	//~ m_effi_corr = new CP::MuonEfficiencyScaleFactors m_effi_corr("SomeName");
	//~ m_effi_corr->setProperty("WorkingPoint","CBandST");
	//~ m_effi_corr->setProperty("DataPeriod","2012");
	//~ CHECK (m_effi_corr->initialize().isSuccess());
	
	//m_METUtil = new METUtility;
	//m_METUtil->setVerbosity(true);
	//m_util->setSoftJetCut(20);

	if (m_useHistObjectDumper)
		m_HistObjectDumper = new HistObjectDumper(wk());

	return EL::StatusCode::SUCCESS;
}



EL::StatusCode MyxAODAnalysis :: execute ()
{
	/// Here you do everything that needs to be done on every single
	/// events, e.g. read input variables, apply cuts, and fill
	/// histograms and trees.  This is where most of your actual analysis
	/// code will go.
	
	const char* APP_NAME = "MyxAODAnalysis";

	if( (m_eventCounter % 100) ==0 ) Info("execute()", "Event number = %i", m_eventCounter );
	m_eventCounter++;

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

	/// if data check if event passes GRL
	if(!isMC){ /// it's data!
		if(!m_grl->passRunLB(*eventInfo)){
			return EL::StatusCode::SUCCESS; // go to next event
		}
	} /// end if not MC

	///------------------------------------------------------------
	/// Apply event cleaning to remove events due to 
	/// problematic regions of the detector, and incomplete events.
	/// Apply to data.
	///------------------------------------------------------------
	/// reject event if:
	if(!isMC){
		if(   (eventInfo->errorState(xAOD::EventInfo::LAr)==xAOD::EventInfo::Error ) || (eventInfo->errorState(xAOD::EventInfo::Tile)==xAOD::EventInfo::Error ) || (eventInfo->isEventFlagBitSet(xAOD::EventInfo::Core, 18) )  )
		{
			return EL::StatusCode::SUCCESS; /// go to the next event
		} /// end if event flags check
	} /// end if the event is data
	m_numCleanEvents++;

  int numGoodJets = 0;

	/// LOOP OVER JETS
	
	// Loop over all jets in the event
	// get jet container of interest
	const xAOD::JetContainer* jets = 0;
	if ( !m_event->retrieve( jets, "AntiKt4LCTopoJets" ).isSuccess() ){ // retrieve arguments: container type, container key
		Error("execute()", "Failed to retrieve AntiKt4LCTopoJets container. Exiting." );
		return EL::StatusCode::FAILURE;
	}

	// loop over the jets in the container
	xAOD::JetContainer::const_iterator jet_itr = jets->begin();
	xAOD::JetContainer::const_iterator jet_end = jets->end();
	for( ; jet_itr != jet_end; ++jet_itr ) {
		if( !m_jetCleaning->accept( **jet_itr )) continue; //only keep good clean jets
		numGoodJets++;
		//Info("execute()", " clean jet pt = %.2f GeV", ((*jet_itr)->pt() * 0.001)); // just to print out something
		h_jetPt->Fill( ( (*jet_itr)->pt()) * 0.001); // GeV
		
		//~ // JER and uncert
		//~ if(isMC){ // assuming isMC flag has been set based on eventInfo->eventType( xAOD::EventInfo::IS_SIMULATION ) 
			//~ // Get the MC resolution
			//~ double mcRes = m_JERTool->getRelResolutionMC((*jet_itr));
			//~ // Get the resolution uncertainty
			//~ double uncert = m_JERTool->getUncertainty((*jet_itr), false, false); // getUncertainty(const xAOD::Jet* jet, bool alt2, bool isAFII)
			//~ Info("execute()", "jet mcRes = %f , uncert = %f", mcRes, uncert);
		//~ } // end if MC
	} // end for loop over jets

	//Info("execute()", "  number of jets = %lu; number of clean jets = %lu", jets->size(), numGoodJets);
  
  	/// get MET_RefFinalFix container of interest
	const xAOD::MissingETContainer* metcontainer = 0;
	
	/// we need to use MET_RefFinalFix, according to:
	/// https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/Run2xAODMissingET
	if ( !m_event->retrieve( metcontainer, "MET_RefFinalFix" ).isSuccess() ){ /// retrieve arguments: container$
		Error("execute()", "Failed to retrieve MET_RefFinalFix container. Exiting." );
		return EL::StatusCode::FAILURE;
	}

	/// We want to check default Final MET, wo any recalibration first
	xAOD::MissingETContainer::const_iterator met_it = metcontainer->find("Final"); 

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
	
	h_MET_RefFinalFix->Fill(sqrt(mpx*mpx + mpy*mpy)* 0.001);
	h_MET_RefFinalFix_test->Fill(metVec->Pt() * 0.001);

	double phi_met = metVec->Phi();
  
	/// get muon container of interest
	const xAOD::MuonContainer* muons = 0;
	if ( !m_event->retrieve( muons, "Muons" ).isSuccess() ){ /// retrieve arguments: container$
		Error("execute()", "Failed to retrieve Muons container. Exiting." );
		return EL::StatusCode::FAILURE;
	}
  
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
  
	/// loop over the muons in the container
	xAOD::MuonContainer::const_iterator muon_itr = muons->begin();
	xAOD::MuonContainer::const_iterator muon_end = muons->end();
	for( ; muon_itr != muon_end; ++muon_itr ) {
	
		if (m_useHistObjectDumper) m_HistObjectDumper->plotMuon((*muon_itr),"noCuts");
	
		xAOD::Muon* mu = 0;
		if( !m_muonCalibrationAndSmearingTool->correctedCopy( **muon_itr, mu ) ) {
			Error(APP_NAME, "Cannot really apply calibration nor smearing");
			continue;
		}

		if(m_muonSelection->accept(mu)){
			double phi_mu = mu->phi();
			double Mt = sqrt( 2*mu->pt()*sqrt(mpx*mpx + mpy*mpy) * (1.0 - TMath::Cos( phi_mu - phi_met )) );
			h_Mt->Fill(Mt * 0.001);

			if (( mu->pt()) * 0.001 >= 50.0){
				h_Mt_muonPtCut->Fill(Mt * 0.001);
				if (m_useHistObjectDumper) m_HistObjectDumper->plotMuon(mu,"allCuts");
			}
		}
		else{
			if (m_useHistObjectDumper) m_HistObjectDumper->plotMuon(mu,"rejectedByMuonSelectorToot");
		}
	
	} /// end for loop over muons


	


	tree->Fill();

	return EL::StatusCode::SUCCESS;
}



EL::StatusCode MyxAODAnalysis :: postExecute ()
{
	// Here you do everything that needs to be done after the main event
	// processing.  This is typically very rare, particularly in user
	// code.  It is mainly used in implementing the NTupleSvc.
	return EL::StatusCode::SUCCESS;
}



EL::StatusCode MyxAODAnalysis :: finalize ()
{
	// This method is the mirror image of initialize(), meaning it gets
	// called after the last event has been processed on the worker node
	// and allows you to finish up any objects you created in
	// initialize() before they are written to disk.  This is actually
	// fairly rare, since this happens separately for each worker node.
	// Most of the time you want to do your post-processing on the
	// submission node after all your histogram outputs have been
	// merged.  This is different from histFinalize() in that it only
	// gets called on worker nodes that processed input events.
	
	const char* APP_NAME = "MyxAODAnalysis";  
	
	///*************************
	///
	/// deleting of all tools
	///
	/// ************************
	
	/// GRL
	if( m_grl ) {
		delete m_grl;
		m_grl = 0;
	}
	/// Jet cleaning tool
	if( m_jetCleaning ) {
		delete m_jetCleaning;
		m_jetCleaning = 0;
	}
	if(m_JERTool){
		delete m_JERTool;
		m_JERTool = 0;
	}
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

	return EL::StatusCode::SUCCESS;
}



EL::StatusCode MyxAODAnalysis :: histFinalize ()
{
	// This method is the mirror image of histInitialize(), meaning it
	// gets called after the last event has been processed on the worker
	// node and allows you to finish up any objects you created in
	// histInitialize() before they are written to disk.  This is
	// actually fairly rare, since this happens separately for each
	// worker node.  Most of the time you want to do your
	// post-processing on the submission node after all your histogram
	// outputs have been merged.  This is different from finalize() in
	// that it gets called on all worker nodes regardless of whether
	// they processed input events.
	return EL::StatusCode::SUCCESS;
}
