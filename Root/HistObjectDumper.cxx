#include <MyAnalysis/HistObjectDumper.h>

/// this is needed to distribute the algorithm to the workers
ClassImp(HistObjectDumper)

HistObjectDumper::HistObjectDumper(){
	
	TH1* muon_pt_original = new TH1F("muon_pt","muon_pt", 300, 0, 3000);
	TH1* muon_eta_original = new TH1F("muon_eta","muon_eta", 300, 0, 3000);
	
	muon_pt.push_back(muon_pt_original);
	muon_eta.push_back(muon_eta_original);
	
	muonHistMap["pt"] = muon_pt;
	muonHistMap["eta"] = muon_eta;
}

HistObjectDumper::~HistObjectDumper(){
	
}

int HistObjectDumper::InitNewStageHists(string stage_tag){

    int iStage;

	typedef map<string,vector<TH1*> >::iterator it_type;
	for(it_type iterator = muonHistMap.begin(); iterator != muonHistMap.end(); iterator++) {
		
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

