#include <MyAnalysis/HistObjectDumper.h>
#include <xAODTracking/TrackingPrimitives.h>

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
	
	TH1* muon_type_original = new TH1F("type","muon_type", 4, -0.5, 4.5);
	muon_type_original->GetXaxis()->SetBinLabel(xAOD::Muon_v1::Combined,"Combined");
	muon_type_original->GetXaxis()->SetBinLabel(xAOD::Muon_v1::SegmentTagged,"SegmentTagged");
	muon_type_original->GetXaxis()->SetBinLabel(xAOD::Muon_v1::CaloTagged,"CaloTagged");
	muon_type_original->GetXaxis()->SetBinLabel(xAOD::Muon_v1::SiliconAssociatedForwardMuon,"SiliconAssociatedForwardMuon");
	muon_type_original->GetXaxis()->SetBinLabel(xAOD::Muon_v1::MuonStandAlone,"MuonStandAlone");
	map<string,TH1*> map_muon_type;
	map_muon_type["reference"] = muon_type_original;
	m_muonHistMap["type"] = map_muon_type;
	
	TH1* muon_author_original = new TH1F("author","muon_author", 12, -1.5, 10.5);
	muon_author_original->GetXaxis()->SetBinLabel(-1,"non-defined");
	muon_author_original->GetXaxis()->SetBinLabel(xAOD::Muon_v1::unknown,"unknown");
	muon_author_original->GetXaxis()->SetBinLabel(xAOD::Muon_v1::MuidCo,"MuidCo");
	muon_author_original->GetXaxis()->SetBinLabel(xAOD::Muon_v1::STACO,"STACO");
	muon_author_original->GetXaxis()->SetBinLabel(xAOD::Muon_v1::MuTag,"MuTag");
	muon_author_original->GetXaxis()->SetBinLabel(xAOD::Muon_v1::MuTagIMO,"MuTagIMO");
	muon_author_original->GetXaxis()->SetBinLabel(xAOD::Muon_v1::MuidSA,"MuidSA");
	muon_author_original->GetXaxis()->SetBinLabel(xAOD::Muon_v1::MuGirl,"MuGirl");
	muon_author_original->GetXaxis()->SetBinLabel(xAOD::Muon_v1::MuGirlLowBeta,"MuGirlLowBeta");
	muon_author_original->GetXaxis()->SetBinLabel(xAOD::Muon_v1::CaloTag,"CaloTag");
	muon_author_original->GetXaxis()->SetBinLabel(xAOD::Muon_v1::CaloLikelihood,"CaloLikelihood");
	muon_author_original->GetXaxis()->SetBinLabel(xAOD::Muon_v1::ExtrapolateMuonToIP,"ExtrapolateMuonToIP");
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
	
}
