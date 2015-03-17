#ifndef HistObjectDumper_H
#define HistObjectDumper_H

#include "EventLoop/Algorithm.h"
#include "xAODMuon/MuonContainer.h"

#include <TH1.h>
#include <iostream>
#include <map>

using namespace std;

class HistObjectDumper
{

public:
	HistObjectDumper(Worker *wk);
	~HistObjectDumper();

	void plotMuon(xAOD::Muon* mu, string stage_tag);
	
private:
	
	int InitNewStageHists(string stage_tag);
	
	/// link to EventLoop worker;
	Worker *m_wk; //!
	
	/// string - variable name
	map<string,map<string,TH1*> > m_muonHistMap; //!
	
	map<string,bool> m_StageNameDict; //!
	
	/// this is needed to distribute the algorithm to the workers
	ClassDef(HistObjectDumper, 1);

};

#endif