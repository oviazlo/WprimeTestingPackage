#ifndef HistObjectDumper_H
#define HistObjectDumper_H

#include "EventLoop/Worker.h"
#include "xAODMuon/MuonContainer.h"

/// EDM include(s):
#include "xAODMuon/Muon.h"

#include <TH1.h>
#include <iostream>
#include <map>
#include <TMath.h>

using namespace std;

class HistObjectDumper
{

public:
	HistObjectDumper(EL::Worker *wk);
	~HistObjectDumper();

	void plotMuon(const xAOD::Muon* mu, string stage_tag);
        void printMuon(const xAOD::Muon* mu);
	
private:
	
	int InitNewStageHists(map<string,map<string,TH1*> >& inMap, string folderName, string stage_tag);
	
	/// link to EventLoop worker;
	EL::Worker *m_wk; //!
	
	/// string - variable name
	map<string,map<string,TH1*> > m_muonHistMap; //!
	
	map<string,bool> m_StageNameDict; //!
	
	/// this is needed to distribute the algorithm to the workers
	ClassDef(HistObjectDumper, 1);

};

#endif
