#ifndef HistObjectDumper_H
#define HistObjectDumper_H

#include "xAODMuon/MuonContainer.h"
#include <TH1.h>
#include <iostream>
#include <map>

using namespace std;

class HistObjectDumper
{

public:
	HistObjectDumper();
	~HistObjectDumper();

// 	plotMuon(xAOD::Muon* mu, string stage_tag = "");
	
private:
	
	int InitNewStageHists(string stage_tag = "");
	
	vector<TH1*> muon_pt; //!
	vector<TH1*> muon_eta; //!

	///
	/// string - variable name
	///
	map<string,vector<TH1*> > muonHistMap; //!
	
	/// this is needed to distribute the algorithm to the workers
	ClassDef(HistObjectDumper, 1);

};

#endif
