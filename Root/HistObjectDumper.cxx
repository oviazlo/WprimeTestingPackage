#include <MyAnalysis/HistObjectDumper.h>
#include <xAODTracking/TrackingPrimitives.h>
#include <xAODPrimitives/IsolationType.h>

/// this is needed to distribute the algorithm to the workers
ClassImp(HistObjectDumper)

HistObjectDumper::HistObjectDumper(EL::Worker *wk){
/// **************************************************************
///
/// Init all histograms here...
///
/// **************************************************************
	/// copy pointer to EventLoop worker
	m_wk = wk;

	/// Prepare reference histograms
	/// (they will not be saved to the output)
	TH1* muon_pt_original = new TH1F("pt","muon_pt", 300, 0, 3000);
	map<string,TH1*> map_muon_pt;
	map_muon_pt["reference"] = muon_pt_original;
	m_muonHistMap["pt"] = map_muon_pt;
	
	TH1* muon_eta_original = new TH1F("eta","muon_eta", 48, -3.0,3.0);
	map<string,TH1*> map_muon_eta;
	map_muon_eta["reference"] = muon_eta_original;
	m_muonHistMap["eta"] = map_muon_eta;
	
	TH1* muon_phi_original = new TH1F("phi","muon_phi", 60, -TMath::Pi(),TMath::Pi());
	map<string,TH1*> map_muon_phi;
	map_muon_phi["reference"] = muon_phi_original;
	m_muonHistMap["phi"] = map_muon_phi;
	
	TH1* muon_quality_original = new TH1F("quality","muon_quality", 4, -0.5, 3.5);
	map<string,TH1*> map_muon_quality;
	map_muon_quality["reference"] = muon_quality_original;
	m_muonHistMap["quality"] = map_muon_quality;
	
/// 	Muon_v1::Combined - 0
/// 	Muon_v1::MuonStandAlone - 1
/// 	Muon_v1::SegmentTagged - 2
/// 	Muon_v1::CaloTagged - 3
/// 	Muon_v1::SiliconAssociatedForwardMuon - 4

	TH1* muon_type_original = new TH1F("type","muon_type", 5, -0.5, 4.5);
	muon_type_original->GetXaxis()->SetBinLabel(muon_type_original->FindBin(xAOD::Muon_v1::Combined),"Combined");
	muon_type_original->GetXaxis()->SetBinLabel(muon_type_original->FindBin(xAOD::Muon_v1::SegmentTagged),"SegmentTagged");
	muon_type_original->GetXaxis()->SetBinLabel(muon_type_original->FindBin(xAOD::Muon_v1::CaloTagged),"CaloTagged");
	muon_type_original->GetXaxis()->SetBinLabel(muon_type_original->FindBin(xAOD::Muon_v1::SiliconAssociatedForwardMuon),"SiliconAssociatedForwardMuon");
	muon_type_original->GetXaxis()->SetBinLabel(muon_type_original->FindBin(xAOD::Muon_v1::MuonStandAlone),"MuonStandAlone");
	map<string,TH1*> map_muon_type;
	map_muon_type["reference"] = muon_type_original;
	m_muonHistMap["type"] = map_muon_type;

/// 	xAOD::Muon_v1::MuidCo = 1
/// 	xAOD::Muon_v1::STACO = 2
/// 	xAOD::Muon_v1::MuTag = 3
	
	TH1* muon_author_original = new TH1F("author","muon_author", 13, -1.5, 11.5);
	muon_author_original->GetXaxis()->SetBinLabel(muon_author_original->FindBin(-1),"non-defined");
	muon_author_original->GetXaxis()->SetBinLabel(muon_author_original->FindBin(xAOD::Muon_v1::unknown),"unknown");
	muon_author_original->GetXaxis()->SetBinLabel(muon_author_original->FindBin(xAOD::Muon_v1::MuidCo),"MuidCo");
	muon_author_original->GetXaxis()->SetBinLabel(muon_author_original->FindBin(xAOD::Muon_v1::STACO),"STACO");
	muon_author_original->GetXaxis()->SetBinLabel(muon_author_original->FindBin(xAOD::Muon_v1::MuTag),"MuTag");
	muon_author_original->GetXaxis()->SetBinLabel(muon_author_original->FindBin(xAOD::Muon_v1::MuTagIMO),"MuTagIMO");
	muon_author_original->GetXaxis()->SetBinLabel(muon_author_original->FindBin(xAOD::Muon_v1::MuidSA),"MuidSA");
	muon_author_original->GetXaxis()->SetBinLabel(muon_author_original->FindBin(xAOD::Muon_v1::MuGirl),"MuGirl");
	muon_author_original->GetXaxis()->SetBinLabel(muon_author_original->FindBin(xAOD::Muon_v1::MuGirlLowBeta),"MuGirlLowBeta");
	muon_author_original->GetXaxis()->SetBinLabel(muon_author_original->FindBin(xAOD::Muon_v1::CaloTag),"CaloTag");
	muon_author_original->GetXaxis()->SetBinLabel(muon_author_original->FindBin(xAOD::Muon_v1::CaloLikelihood),"CaloLikelihood");
	muon_author_original->GetXaxis()->SetBinLabel(muon_author_original->FindBin(xAOD::Muon_v1::ExtrapolateMuonToIP),"ExtrapolateMuonToIP");
	map<string,TH1*> map_muon_author;
	map_muon_author["reference"] = muon_author_original;
	m_muonHistMap["author"] = map_muon_author;
	
	TH1* muon_precMSLayers_original = new TH1F("precMSLayers","number of MS layers with at least 3 hits", 11, -0.5, 10.5);
	map<string,TH1*> map_muon_precMSLayers;
	map_muon_precMSLayers["reference"] = muon_precMSLayers_original;
	m_muonHistMap["precMSLayers"] = map_muon_precMSLayers;
	
	TH1* muon_phiMSLayers_original = new TH1F("phiMSLayers","number of MS layers with a trigger phi hit", 11, -0.5, 10.5);
	map<string,TH1*> map_muon_phiMSLayers;
	map_muon_phiMSLayers["reference"] = muon_phiMSLayers_original;
	m_muonHistMap["phiMSLayers"] = map_muon_phiMSLayers;
	
	TH1* muon_ptCone30OverPt_original = new TH1F("ptCone30OverPt","ptCone30/pt", 100, 0.0, 1.0);
	map<string,TH1*> map_muon_ptCone30OverPt;
	map_muon_ptCone30OverPt["reference"] = muon_ptCone30OverPt_original;
	m_muonHistMap["ptCone30OverPt"] = map_muon_ptCone30OverPt;
	
	TH1* muon_ptCone30_original = new TH1F("ptCone30","ptCone30", 30000, 0, 3000);
	map<string,TH1*> map_muon_ptCone30;
	map_muon_ptCone30["reference"] = muon_ptCone30_original;
	m_muonHistMap["ptCone30"] = map_muon_ptCone30;
}

HistObjectDumper::~HistObjectDumper(){
	
}

int HistObjectDumper::InitNewStageHists(map<string,map<string,TH1*> >& inMap, string folderName, string stage_tag){

	int iStage;

	typedef map<string,map<string,TH1*> >::iterator it_type;
	for(it_type iterator = inMap.begin(); iterator != inMap.end(); iterator++) {
		
		/// convert int to string
		iStage = iterator->second.size();
		std::ostringstream tmpStream;
		tmpStream << iStage;
		
		/// construct histName
		string histName = folderName + "/stage" + tmpStream.str();
		if (stage_tag!="")
			histName += "_" + stage_tag;
		histName += "/" + string(iterator->second["reference"]->GetName());
		
		/// Clone referance hist with adding to suffix to the name of new hist and push it back to vector
		iterator->second[stage_tag] = (TH1*)iterator->second["reference"]->Clone(histName.c_str());
		
		/// register new hist in EventLoop worker
		m_wk->addOutput(iterator->second[stage_tag]);
		
	}

	m_StageNameDict[stage_tag] = true;
		
	return iStage;
	
}

void HistObjectDumper::plotMuon(const xAOD::Muon* mu, string stage_tag){
/// **************************************************************
///
/// Filling all histograms...
///	
/// **************************************************************
	
	if (m_StageNameDict[stage_tag]==false)
		InitNewStageHists(m_muonHistMap,"muon",stage_tag);
	
	m_muonHistMap["pt"][stage_tag]->Fill(mu->pt()*0.001);
	m_muonHistMap["eta"][stage_tag]->Fill(mu->eta());
	m_muonHistMap["phi"][stage_tag]->Fill(mu->phi());
	m_muonHistMap["quality"][stage_tag]->Fill(mu->quality());
	m_muonHistMap["type"][stage_tag]->Fill(mu->muonType());
	m_muonHistMap["author"][stage_tag]->Fill(mu->author());
	
	uint8_t nMSPrecLayers = -1;
	uint8_t nLayersWithPhiHit = -1;
	
	mu->primaryTrackParticle()->summaryValue(nMSPrecLayers, xAOD::numberOfPrecisionLayers);	/// < layers with at least 3 hits [unit8_t].
	mu->primaryTrackParticle()->summaryValue(nLayersWithPhiHit, xAOD::numberOfPhiLayers);		/// < layers with a trigger phi hit [unit8_t].
	
	m_muonHistMap["precMSLayers"][stage_tag]->Fill(nMSPrecLayers);
	m_muonHistMap["phiMSLayers"][stage_tag]->Fill(nLayersWithPhiHit);
	
	float muPtCone30 = 0.; /// your variable that will be filled after calling the isolation function
	mu->isolation(muPtCone30, xAOD::Iso::ptcone30);  /// second arg is an enum defined in xAODPrimitives/IsolationType.h
	
	m_muonHistMap["ptCone30"][stage_tag]->Fill(muPtCone30*0.001);
	m_muonHistMap["ptCone30OverPt"][stage_tag]->Fill(muPtCone30/mu->pt());
}


void HistObjectDumper::printMuon(const xAOD::Muon* mu){
/// **************************************************************
///
/// Print info about selected muon
///     
/// **************************************************************
  
  float muPtCone30 = 0.; /// your variable that will be filled after calling the isolation function
  mu->isolation(muPtCone30, xAOD::Iso::ptcone30);  /// second arg is an enum defined in xAODPrimitives/IsolationType.h
  
  cout << "pt:\t" << mu->pt()*0.001 << endl;
  cout << "eta:\t" << mu->eta() << endl;
  cout << "phi:\t" << mu->phi() << endl;
  cout << "qual:\t" << mu->quality() << endl;
  cout << "ptCone/pt:\t" << muPtCone30/mu->pt() << endl;
   
}