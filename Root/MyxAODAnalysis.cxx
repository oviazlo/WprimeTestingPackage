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

/// GRL
#include "GoodRunsLists/GoodRunsListSelectionTool.h"

#include "xAODJet/JetContainer.h"
#include "xAODMuon/MuonContainer.h"
#include "JetSelectorTools/JetCleaningTool.h"
#include "JetResolution/JERTool.h"
#include <TSystem.h> /// used to define JERTool calibration path 
#include "xAODTracking/VertexContainer.h"

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
	m_useBitsetCutflow = true;
	m_useMuonCalibrationAndSmearingTool = true;

	return EL::StatusCode::SUCCESS;
}



EL::StatusCode MyxAODAnalysis :: histInitialize ()
{
	// Here you do everything that needs to be done at the very
	// beginning on each worker node, e.g. create histograms and output
	// trees.  This method gets called before any input files are
	// connected.

	h_zPrimVtx = new TH1F("h_zPrimVtx", "h_zPrimVtx", 6000, -300.0, 300); // jet pt [GeV]
	wk()->addOutput (h_zPrimVtx);

	h_jetPt_VeryLooseBadJets = new TH1F("h_jetPt_VeryLooseBadJets", "VeryLooseBadJets and better with pt>20", 300, 0, 3000); // jet pt [GeV]
	wk()->addOutput (h_jetPt_VeryLooseBadJets);

	h_jetPt_LooseBadJets = new TH1F("h_jetPt_LooseBadJets", "LooseBadJets and better with pt>20", 300, 0, 3000); // jet pt [GeV]
	wk()->addOutput (h_jetPt_LooseBadJets);
	
	/// Muon Pt histograms:
	h_MET_RefFinalFix = new TH1F("h_MET_RefFinalFix", "h_MET_RefFinalFix", 1000, 0, 10000); 
	wk()->addOutput (h_MET_RefFinalFix);
	
	h_MET_RefFinalFix_test = new TH1F("h_MET_RefFinalFix_test", "h_MET_RefFinalFix_test", 1000, 0, 10000); 
	wk()->addOutput (h_MET_RefFinalFix_test);
	
	h_Mt = new TH1F("h_Mt", "h_Mt", 1000, 0, 10000); 
	wk()->addOutput (h_Mt);

	h_Mt_muonPtCut = new TH1F("h_Mt_muonPtCut", "h_Mt_muonPtCut", 1000, 0, 10000); 
	wk()->addOutput (h_Mt_muonPtCut);

	h_nSelectedMuons = new TH1I("h_nSelectedMuons", "h_nSelectedMuons", 5, 0, 5); 
	wk()->addOutput (h_nSelectedMuons);
	
	h_truthMET = new TH1F("h_truthMET", "Truth MET", 1000, 0, 10000); 
	wk()->addOutput (h_truthMET);
	
	// get the output file, create a new TTree and connect it to that output
	// define what branches will go in that tree
// 	TFile *outputFile = wk()->getOutputFile (outputName);
// 	tree = new TTree ("tree", "tree");
// 	tree->SetDirectory (outputFile);
// 	tree->Branch("EventNumber", &EventNumber);

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
	m_jetCleaning->msg().setLevel( MSG::INFO ); 
	CHECK(m_jetCleaning->setProperty( "CutLevel", "VeryLooseBad"));
	m_jetCleaning->initialize();
	
	m_jetCleaning_LooseBad = new JetCleaningTool("JetCleaning_LooseBad");
	m_jetCleaning_LooseBad->msg().setLevel( MSG::INFO ); 
	CHECK(m_jetCleaning_LooseBad->setProperty( "CutLevel", "LooseBad"));
	m_jetCleaning_LooseBad->initialize();
	
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
	m_muonSelection->msg().setLevel( MSG::INFO );
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
	
	if (m_BitsetCutflow)
		m_BitsetCutflow = new BitsetCutflow(wk());
	

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
	
	if( (m_eventCounter % 1000) ==0 ) Info("execute()", "Event number = %i", m_eventCounter );
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
	m_BitsetCutflow->FillCutflow("GRL");
	

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
	m_BitsetCutflow->FillCutflow("LAr_Tile_Core");

	const xAOD::VertexContainer* vertices = 0;
	if ( !m_event->retrieve( vertices, "PrimaryVertices" ).isSuccess() ){ // retrieve arguments: container type, container key
		Error("execute()", "Failed to retrieve PrimaryVertices container. Exiting." );
		return EL::StatusCode::FAILURE;
	}
	
	xAOD::VertexContainer::const_iterator vtx_itr = vertices->begin();
	xAOD::VertexContainer::const_iterator vtx_end = vertices->end();
	int nGoodVtx = 0;
	for( ; vtx_itr != vtx_end; ++vtx_itr ) {
		//~ h_zPrimVtx->Fill((*vtx_itr)->z());
		//~ if (((*vtx_itr)->vertexType()==xAOD::VxType::PriVtx)&&(abs((*vtx_itr)->z())<200.0))
		if ((abs((*vtx_itr)->z())<200.0)&&((*vtx_itr)->nTrackParticles()>=3))
			nGoodVtx++;
	}

	/// LOOP OVER JETS
	int nVeryLooseBadJets = 0;
	int nLooseBadJets = 0;
	
	/// Loop over all jets in the event
	/// get jet container of interest
	const xAOD::JetContainer* jets = 0;
	if ( !m_event->retrieve( jets, "AntiKt4LCTopoJets" ).isSuccess() ){ // retrieve arguments: container type, container key
		Error("execute()", "Failed to retrieve AntiKt4LCTopoJets container. Exiting." );
		return EL::StatusCode::FAILURE;
	}

	/// loop over the jets in the container
	xAOD::JetContainer::const_iterator jet_itr = jets->begin();
	xAOD::JetContainer::const_iterator jet_end = jets->end();
	for( ; jet_itr != jet_end; ++jet_itr ) {
		
		/// check for Jet pt (not interested in Jets with pt<=20)
		if (( (*jet_itr)->pt()) * 0.001 <= 20.0) continue;
		
		/// check for VeryLooseBad jets
		if( !m_jetCleaning->accept( **jet_itr )){
			m_BitsetCutflow->FillCutflow("VeryBadJet");
			continue;
		}
		nVeryLooseBadJets++;
		h_jetPt_VeryLooseBadJets->Fill( ( (*jet_itr)->pt()) * 0.001); // GeV
		
		/// check for LooseBad jets
		if( !m_jetCleaning_LooseBad->accept( **jet_itr )) continue; //only keep good clean jets
		nLooseBadJets++;
		h_jetPt_LooseBadJets->Fill( ( (*jet_itr)->pt()) * 0.001); // GeV

	} // end for loop over jets

	if (nVeryLooseBadJets!=nLooseBadJets){
		return EL::StatusCode::SUCCESS;
	}
	
	m_BitsetCutflow->FillCutflow("JetCleaning");
	
	if (nGoodVtx==0)
		return EL::StatusCode::SUCCESS;
	m_BitsetCutflow->FillCutflow("Primary vertex");
	
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
	int muonCounter = 0;
	xAOD::MuonContainer::const_iterator muon_itr = muons->begin();
	xAOD::MuonContainer::const_iterator muon_end = muons->end();
	for( ; muon_itr != muon_end; ++muon_itr ) {
	
		if (m_useHistObjectDumper) m_HistObjectDumper->plotMuon((*muon_itr),"noCuts");
	
		m_BitsetCutflow->FillCutflow("oneMuon");
	
		xAOD::Muon* mu = 0;
		if (m_useMuonCalibrationAndSmearingTool){
			if( !m_muonCalibrationAndSmearingTool->correctedCopy( **muon_itr, mu ) ) {
				Error(APP_NAME, "Cannot really apply calibration nor smearing");
				continue;
			}
		}
		else{
			mu = const_cast<xAOD::Muon*> (*muon_itr);
		}
		
		if (mu->muonType()!=xAOD::Muon_v1::Combined) continue;
			m_BitsetCutflow->FillCutflow("Combined");
		
		if (( mu->pt()) * 0.001 < 55.0) continue;

		m_BitsetCutflow->FillCutflow("mu_pt");
		
		if(m_muonSelection->accept(mu)){

			m_BitsetCutflow->FillCutflow("MCP selector");
			
			//~ uint8_t nMSPrecLayers = -1;
			//~ mu->primaryTrackParticle()->summaryValue(nMSPrecLayers, xAOD::numberOfPrecisionLayers);	/// < layers with at least 3 hits [unit8_t].
			//~ if (nMSPrecLayers<3) continue;
			//~ m_BitsetCutflow->FillCutflow("3 hits in 3 MS layers");
			
			//~ uint8_t nLayersWithPhiHit = -1;
			//~ mu->primaryTrackParticle()->summaryValue(nLayersWithPhiHit, xAOD::numberOfPhiLayers);		/// < layers with a trigger phi hit [unit8_t].
			//~ if (nLayersWithPhiHit<2) continue;
			//~ m_BitsetCutflow->FillCutflow("2 phi layers");
			
			uint8_t n_innerSmallHits = -1;
			mu->summaryValue(n_innerSmallHits, xAOD::innerSmallHits);
			
			uint8_t n_innerLargeHits = -1;
			mu->summaryValue(n_innerLargeHits, xAOD::innerLargeHits);
			
			uint8_t n_middleSmallHits = -1;
			mu->summaryValue(n_middleSmallHits, xAOD::middleSmallHits);
			
			uint8_t n_middleLargeHits = -1;
			mu->summaryValue(n_middleLargeHits, xAOD::middleLargeHits);
			
			uint8_t n_outerSmallHits = -1;
			mu->summaryValue(n_outerSmallHits, xAOD::outerSmallHits);
			
			uint8_t n_outerLargeHits = -1;
			mu->summaryValue(n_outerLargeHits, xAOD::outerLargeHits);
			
			uint8_t n_extendedSmallHits = -1;
			mu->summaryValue(n_extendedSmallHits, xAOD::extendedSmallHits);
			
			uint8_t n_extendedLargeHits = -1;
			mu->summaryValue(n_extendedLargeHits, xAOD::extendedLargeHits);
			
			uint8_t n_innerMDTLayerHits = n_innerSmallHits + n_innerLargeHits;
			uint8_t n_middleMDTLayerHits = n_middleSmallHits + n_middleLargeHits;
			uint8_t n_outerMDTLayerHits = n_outerSmallHits + n_outerLargeHits;
			uint8_t n_extendedMDTLayerHits = n_extendedSmallHits + n_extendedLargeHits;
			
			if ((n_innerMDTLayerHits<3)||(n_middleMDTLayerHits<3)||(n_outerMDTLayerHits<3)) continue;
			m_BitsetCutflow->FillCutflow("3 hits in all 3 MS layers");
						
			uint8_t n_phiLayer1Hits = -1;
			mu->summaryValue(n_phiLayer1Hits, xAOD::phiLayer1Hits);
			
			uint8_t n_phiLayer2Hits = -1;
			mu->summaryValue(n_phiLayer2Hits, xAOD::phiLayer2Hits);
			
			uint8_t n_phiLayer3Hits = -1;
			mu->summaryValue(n_phiLayer3Hits, xAOD::phiLayer3Hits);
			
			uint8_t n_phiLayer4Hits = -1;
			mu->summaryValue(n_phiLayer4Hits, xAOD::phiLayer4Hits);
						
			int nPhiLayers = 0;
			if (n_phiLayer1Hits>0) nPhiLayers++;
			if (n_phiLayer2Hits>0) nPhiLayers++;
			if (n_phiLayer3Hits>0) nPhiLayers++;
			if (n_phiLayer4Hits>0) nPhiLayers++;
			
			if (nPhiLayers<2) continue;
			m_BitsetCutflow->FillCutflow("2 phi layers");
			
			/// Isolation stuff
			float muPtCone30 = 0.; // your variable that will be filled after calling the isolation function
			mu->isolation(muPtCone30, xAOD::Iso::ptcone30);  // second arg is an enum defined in xAODPrimitives/IsolationType.h

			if (muPtCone30/mu->pt() >= 0.05) continue;
			m_BitsetCutflow->FillCutflow("Isolation");
			
			double phi_mu = mu->phi();
			double Mt = sqrt( 2*mu->pt()*sqrt(mpx*mpx + mpy*mpy) * (1.0 - TMath::Cos( phi_mu - phi_met )) );
			
			//~ if (Mt * 0.001<700.0) continue;
			//~ m_BitsetCutflow->FillCutflow("Dummy Mt cut");
			
			h_Mt_muonPtCut->Fill(Mt * 0.001);
			
			if (m_useHistObjectDumper) m_HistObjectDumper->plotMuon(mu,"allCuts");
			muonCounter++;
		}
		else{
			if (m_useHistObjectDumper) m_HistObjectDumper->plotMuon(mu,"rejectedByMuonSelectorToot");
		}
	
	} /// end for loop over muons

	h_nSelectedMuons->Fill(muonCounter);
	
	if (muonCounter==0) return EL::StatusCode::SUCCESS;
	/// ************************************************
	/// Fill other event-based distributions (e.g. MET)
	/// ************************************************
	h_MET_RefFinalFix->Fill(sqrt(mpx*mpx + mpy*mpy)* 0.001);
	h_MET_RefFinalFix_test->Fill(metVec->Pt() * 0.001);

	m_BitsetCutflow->FillCutflow("End");
	/*
	const xAOD::TruthEventContainer* xTruthEventContainer = NULL;
	CHECK(m_event->retrieve( xTruthEventContainer, "TruthEvent"));
		
		
	cout << "Event number: " << m_eventCounter << endl;
	if (m_eventCounter==4960) return EL::StatusCode::SUCCESS;
	if (m_eventCounter==9484) return EL::StatusCode::SUCCESS;
	
	xAOD::TruthEventContainer::const_iterator itr;
	for (itr = xTruthEventContainer->begin(); itr!=xTruthEventContainer->end(); ++itr) {
		int nVert = (*itr)->nTruthVertices();
		int nPart = (*itr)->nTruthParticles();
		
		//~ cout << "nVert = " << nVert << "\tnPart = " << nPart << endl;
		
		for (int iVtx=0; iVtx<nVert; iVtx++){
			const xAOD::TruthVertex* vertex = (*itr)->truthVertex(iVtx);
			int nIn = vertex->nIncomingParticles();
			bool particleIsFound = false;
			for (int j=0;j<nIn;j++){
					const xAOD::TruthParticle* inPart = vertex->incomingParticle(j);
					int pdgId = inPart->pdgId();
					if (abs(pdgId)==34){
						int nDecPart = vertex->nOutgoingParticles();
						if (nDecPart<2) continue;
						particleIsFound = true;
						//~ if (pdgId==34)
							//~ cout << "W'+ decays to " << nDecPart << " particles: ";
						//~ else
							//~ cout << "W'- decays to " << nDecPart << " particles: ";
							
						int lowestEnergeticTauIndex = -1;
						int tauNeutrinoIndex = -1;
						double tauPt = 9999.9;
						for (int j=0; j<nDecPart; j++){
							const xAOD::TruthParticle* decayPart = vertex->outgoingParticle(j);
							if (abs(decayPart->pdgId())==14){ // muon neutrino
								double truthMET = TMath::Sqrt(decayPart->px()*decayPart->px() + decayPart->py()*decayPart->py())*0.001;
								h_truthMET->Fill(truthMET);
								break;
							}
							if (abs(decayPart->pdgId())==16){ // tau neutrino
								tauNeutrinoIndex = j;
							}
							if (abs(decayPart->pdgId())==15){ // tau
								double currentTauPt = TMath::Sqrt(decayPart->px()*decayPart->px() + decayPart->py()*decayPart->py())*0.001;
								if (tauPt>currentTauPt){
									tauPt = currentTauPt;
									lowestEnergeticTauIndex = j;
								}
							}
							//~ cout << endl << " " << j << ": " << decayPart->pdgId() << " (pT = " << decayPartPt << " GeV)";
						}
						if (lowestEnergeticTauIndex!=-1){
							double met_x = vertex->outgoingParticle(tauNeutrinoIndex)->px();
							double met_y = vertex->outgoingParticle(tauNeutrinoIndex)->py();
							const xAOD::TruthParticle* lowestEnergeticTau = vertex->outgoingParticle(lowestEnergeticTauIndex);
							const xAOD::TruthVertex* tauDecayVtx = lowestEnergeticTau->decayVtx();
							int n_tauDaughters = tauDecayVtx->nOutgoingParticles();
							cout << endl;
							for (int k=0; k<n_tauDaughters; k++){
								const xAOD::TruthParticle* tauDaughter = tauDecayVtx->outgoingParticle(k);
								int tauDaughterAbsPdgId = abs(tauDaughter->pdgId());
								cout << "tau daughter number " << k << ":\t" << tauDaughterAbsPdgId << endl;
								if (tauDaughterAbsPdgId==12 || tauDaughterAbsPdgId==14 || tauDaughterAbsPdgId==16){
									met_x += tauDaughter->px();
									met_y += tauDaughter->py();
								}
							}
							h_truthMET->Fill(TMath::Sqrt(met_x*met_x + met_y*met_y)*0.001);
							cout << "truth MET = " << TMath::Sqrt(met_x*met_x + met_y*met_y)*0.001 << "\t Reco MET = " << metVec->Pt() * 0.001 << endl;
							break;
						}
						//~ if (particleIsFound == false){
							//~ if (pdgId==34)
								//~ cout << "W'+ decays to " << nDecPart << " particles: ";
							//~ else
								//~ cout << "W'- decays to " << nDecPart << " particles: ";
							//~ for (int j=0; j<nDecPart; j++){
								//~ const xAOD::TruthParticle* decayPart = vertex->outgoingParticle(j);
								//~ double decayPartPt = TMath::Sqrt(decayPart->px()*decayPart->px() + decayPart->py()*decayPart->py())*0.001;
								//~ cout << endl << " " << j << ": " << decayPart->pdgId() << " (pT = " << decayPartPt << " GeV)";
							//~ }
							//~ cout << endl;
						//~ }
						//~ particleIsFound = true;
						//~ cout << endl;
					}
					if (particleIsFound) break;
			}
			if (particleIsFound) break;
		}
		
		//~ for (int iPart=0; iPart<nPart; iPart++){
			//~ const xAOD::TruthParticle* particle = (*itr)->truthParticle(iPart);
			//~ int pdgId = particle->pdgId();
			//~ if (abs(pdgId)==34){
				//~ const xAOD::TruthVertex* decayVtx = particle->decayVtx();
				//~ int nDecPart = decayVtx->nOutgoingParticles();
				//~ if (nDecPart<2) continue;
				//~ cout << endl;
				//~ if (pdgId==34)
					//~ cout << "W'+ decays to " << nDecPart << " particles: ";
				//~ else
					//~ cout << "W'- decays to " << nDecPart << " particles: ";
				//~ for (int j=0; j<nDecPart; j++){
					//~ const xAOD::TruthParticle* decayPart = decayVtx->outgoingParticle(j);
					//~ double decayPartPt = TMath::Sqrt(decayPart->px()*decayPart->px() + decayPart->py()*decayPart->py())*0.001;
					//~ cout << endl << " " << j << ": " << decayPart->pdgId() << " (pT = " << decayPartPt << " GeV)";
				//~ }
				//~ cout << endl;
				//~ cout << "Reco MET: " << sqrt(mpx*mpx + mpy*mpy)*0.001 << endl;
				//~ break; // break loop cause we found our particle
			//~ }
		//~ }
	}
	
// 	tree->Fill();
*/
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

	/// push cutflow for last event
	m_BitsetCutflow->PushBitSet();
	
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
	if( m_jetCleaning_LooseBad ) {
		delete m_jetCleaning_LooseBad;
		m_jetCleaning_LooseBad = 0;
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
	
	if(m_useBitsetCutflow && m_BitsetCutflow){
		delete m_BitsetCutflow;
		m_BitsetCutflow = 0;
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
