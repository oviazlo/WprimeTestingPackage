#ifndef HistObjectDumper_H
#define HistObjectDumper_H

/// EventLoop/xAOD
#include "EventLoop/Worker.h"
#include "xAODMuon/MuonContainer.h"
#include "xAODMuon/Muon.h"
#include "xAODMissingET/MissingETContainer.h"
#include "xAODMissingET/MissingETAuxContainer.h"

/// custom
#include <MyAnalysis/WprimeHist.h>

/// ROOT
#include <TMath.h>
#include <TH1.h>

/// std c++
#include <iostream>
#include <map>


using namespace std;

class HistObjectDumper
{

public:
	HistObjectDumper(EL::Worker *wk);
	~HistObjectDumper();

	void plotMuon(const xAOD::Muon* mu, string stage_tag, double weight = 1.0);
    void plotMtAndMet(const xAOD::Muon* mu, xAOD::MissingET* finalTrkMet,
                      string stage_tag, double weight = 1.0);
      
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
