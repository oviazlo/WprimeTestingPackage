#include <MyAnalysis/HistObjectDumper.h>

/// this is needed to distribute the algorithm to the workers
ClassImp(HistObjectDumper)

HistObjectDumper::HistObjectDumper(){
	
	TH1* muon_pt_original = new TH1F("muon_pt","muon_pt", 300, 0, 3000);
	TH1* muon_eta_original = new TH1F("muon_eta","muon_eta", 300, 0, 3000);
	
	muon_pt_.push_back(muon_pt_original);
	muon_eta_.push_back(muon_eta_original);
	
	muonHistMap_["pt"] = muon_pt_;
	muonHistMap_["eta"] = muon_eta_;
}

HistObjectDumper::~HistObjectDumper(){
	
}

int HistObjectDumper::InitNewStageHists(string stage_tag){

    int iStage;

	typedef map<string,vector<TH1*> >::iterator it_type;
	for(it_type iterator = muonHistMap_.begin(); iterator != muonHistMap_.end(); iterator++) {
		
		/// convert int to string
		iStage = iterator->second.size();
		std::ostringstream tmpStream;
		tmpStream << iStage;
		
		/// construct suffix
		string suffix = "_stage" + tmpStream.str();
		if (stage_tag!="")
			suffix + "_" + stage_tag;
		
		/// Clone referance hist with adding to suffix to the name of new hist and push it back to vector
		iterator->second.push_back((TH1*)iterator->second[0]->Clone((iterator->second[0]->GetName() + suffix).c_str()));
		
	}

    return iStage;
	
}

void plotMuon(xAOD::Muon* mu, string stage_tag){
	
	InitNewStageHists(stage_tag);
	muonHistMap_["pt"]->back()->Fill(mu->pt());
	muonHistMap_["eta"]->back()->Fill(mu->eta());
	
}