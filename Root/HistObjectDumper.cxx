#include <MyAnalysis/HistObjectDumper.h>

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
}

HistObjectDumper::~HistObjectDumper(){
	
}

int HistObjectDumper::InitNewStageHists(map<string,map<string,TH1*> > inMap, string folderName, string stage_tag){

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
}
