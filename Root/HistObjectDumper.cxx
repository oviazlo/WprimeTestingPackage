#include <MyAnalysis/HistObjectDumper.h>

/// this is needed to distribute the algorithm to the workers
ClassImp(HistObjectDumper)

HistObjectDumper::HistObjectDumper(EL::Worker *wk){
	
	m_wk = wk;
	
	TH1* muon_pt_original = new TH1F("pt","muon_pt", 300, 0, 3000);
	TH1* muon_eta_original = new TH1F("eta","muon_eta", 48, -3.0,3.0);
	
	map<string,TH1*> map_muon_pt;
	map<string,TH1*> map_muon_eta;
	
	map_muon_pt["reference"] = muon_pt_original;
	map_muon_eta["reference"] = muon_eta_original;
	
	m_muonHistMap["pt"] = map_muon_pt;
	m_muonHistMap["eta"] = map_muon_eta;
}

HistObjectDumper::~HistObjectDumper(){
	
}

int HistObjectDumper::InitNewStageHists(string stage_tag){

    int iStage;

	typedef map<string,map<string,TH1*> >::iterator it_type;
	for(it_type iterator = m_muonHistMap.begin(); iterator != m_muonHistMap.end(); iterator++) {
		
		/// convert int to string
		iStage = iterator->second.size();
		std::ostringstream tmpStream;
		tmpStream << iStage;
		
		/// construct histName
		string histName = "muon/stage" + tmpStream.str();
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

void HistObjectDumper::plotMuon(xAOD::Muon* mu, string stage_tag){
	
	if (m_StageNameDict[stage_tag]==false)
		InitNewStageHists(stage_tag);
	
	m_muonHistMap["pt"][stage_tag]->Fill(mu->pt()*0.001);
	m_muonHistMap["eta"][stage_tag]->Fill(mu->eta());
	
}
