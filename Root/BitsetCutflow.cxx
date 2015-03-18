#include <MyAnalysis/BitsetCutflow.h>

/// this is needed to distribute the algorithm to the workers
ClassImp(BitsetCutflow)

BitsetCutflow::BitsetCutflow(EL::Worker *wk){
	m_wk = wk;
	m_cutflowHist = new TH1I("cutflow_hist","Cutflow",99,0,99);
	m_wk->addOutput(m_cutflowHist);
	m_cutflowStepCounter = 0;
}

BitsetCutflow::~BitsetCutflow(){
	
}

void BitsetCutflow::PushBitSet(){
	
	/// All events
	m_cutflowHist->Fill(0);
	
	for (unsigned i=1; i<=m_StageNameDict.size(); i++){
		/// first bit is always zero
		if (m_bitset[i]!=m_bitset[0]) m_cutflowHist->Fill(i); 
	}
	
	/// clean bit for next event
	m_bitset.reset()
}

BitsetCutflow::void FillCutflow(string cutflowStepName){
	/// check if there is a slot already created for this cutflow step
	int slotPosition = m_StageNameDict[cutflowStepName];
	if (slotPosition==0){
		m_cutflowStepCounter++;
		m_StageNameDict[cutflowStepName] = m_cutflowStepCounter;
		slotPosition = m_cutflowStepCounter;
		m_cutflowHist->GetXaxis()->SetBinLabel(m_cutflowStepCounter,cutflowStepName.c_str());
	}
	m_bitset[slotPosition]=1;
}

