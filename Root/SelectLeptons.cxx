/* Copyright 2016 Oleksandr Viazlo */

///*****************************************************************************
///
/// Custom user-defined functions, which are not standart in EventLoop:
///     - SelectMuons
///     - SelectElectrons
///
/// Both functions return pair(nSignalLeptons,nVetoLeptons)
///
///*****************************************************************************
#include <MyAnalysis/RecoAnalysis.h>

/// this is needed to distribute the algorithm to the workers
ClassImp(RecoAnalysis)

std::map<string, unsigned int> RecoAnalysis :: SelectMuons(
                                  xAOD::MuonContainer* muons,
                                  xAOD::Vertex* primVertex,
                                  bool fillInCutflow){

  xAOD::MuonContainer::iterator muon_itr = muons->begin();
  xAOD::MuonContainer::iterator muon_end = muons->end();
  
  unsigned int nVetoLeptons = 0;
  unsigned int nSignalLeptons = 0;
  unsigned int nSignalLeptons_failIso = 0;
  unsigned int nVetoLeptons_failIso = 0;
  
  for( ; muon_itr != muon_end; ++muon_itr ) {
    
    bool fillCutflow = true;
    
    /// Apply loosened vet selection first, fill in veto flags
    /// Left cut for signal selection apply in the end and fill signal flags
    ( *muon_itr )->auxdata< bool >( "veto" ) = false;
    ( *muon_itr )->auxdata< bool >( "signal" ) = false;
    ( *muon_itr )->auxdata< bool >( "veto_failIso" ) = false;
    ( *muon_itr )->auxdata< bool >( "signal_failIso" ) = false;
    
    if ( ( *muon_itr )->muonType()!=xAOD::Muon_v1::Combined ) continue;
    if (fillCutflow) 
      m_BitsetCutflow->FillCutflow( "Combined",fillInCutflow );
    
    double muPt = ( ( *muon_itr )->pt() ) * 0.001;
    double lowPtCutVeto = 20.0; /// GeV
    double lowPtCutSignal = 55.0; /// GeV
    double highPtCut = 999999.0; /// GeV
    
    /// fill cutflow only if pass highPtCut
    bool passHighPtCut = false;
    if ( muPt < lowPtCutVeto ) continue;
    if ( muPt > lowPtCutSignal )
      passHighPtCut = true;
    fillCutflow = passHighPtCut;
    
    if (fillCutflow){
      m_BitsetCutflow->FillCutflow( "mu_pt",fillInCutflow );
//       cout << "[CUTFLOW_DEBUG_PT]\t" << m_eventInfo->runNumber() << "\t" << m_eventInfo->eventNumber() <<  endl;
    }
    
    /// medium OR high-pT
    if( !(( m_looseMuonSelection->accept( ( *muon_itr ) ) ) ||
        ( m_muonSelection->accept( ( *muon_itr ) ) ))
    ) continue;

    if (passHighPtCut)
      m_BitsetCutflow->FillCutflow( "MCP (medium OR high-pT) selector",fillInCutflow );
    
    /// high-pT selector
    bool passHighPtSelection = false;
    if( m_muonSelection->accept( ( *muon_itr ) ) ){
      passHighPtSelection = true;
    }
    fillCutflow = (passHighPtCut && passHighPtSelection);
    
    /// do significance 
    const xAOD::TrackParticle *tp = ( *muon_itr )->primaryTrackParticle();
    double d0_sig = xAOD::TrackingHelpers::d0significance
      ( tp, m_eventInfo->beamPosSigmaX(), m_eventInfo->beamPosSigmaY(), 
        m_eventInfo->beamPosSigmaXY() );
    if (TMath::Abs(d0_sig)>3.0) continue; 
    if (passHighPtCut)
      m_BitsetCutflow->FillCutflow("d0",fillInCutflow);
    
    /// zo cut
    double z0_vrtPVx = ( *muon_itr )->primaryTrackParticle()->z0() +
    ( *muon_itr )->primaryTrackParticle()->vz() - primVertex->z();
    double sintheta = 1.0/TMath::CosH(( *muon_itr )->eta());
    if (abs( z0_vrtPVx*sintheta )>10.0) continue;
    if (passHighPtCut)
      m_BitsetCutflow->FillCutflow( "z0",fillInCutflow );
    
    if (fillCutflow){ 
      m_BitsetCutflow->FillCutflow( "MCP selector",fillInCutflow );
//       cout << "[CUTFLOW_DEBUG_MCP_SELECTOR]\t" << m_eventInfo->runNumber() << "\t" << m_eventInfo->eventNumber() <<  endl;
    }
    
    if (m_muonisolationSelectionTool->accept(( **muon_itr ))){
      if (fillCutflow)
        m_BitsetCutflow->FillCutflow("Isolation",fillInCutflow);
      
      /// check veto lepton if it satisfy signal requirements
      if (passHighPtCut && passHighPtSelection){
        ( *muon_itr )->auxdata< bool >( "signal" ) = true;
        nSignalLeptons++;
      }
      else{
        ( *muon_itr )->auxdata< bool >( "veto" ) = true;
        nVetoLeptons++;
      }
    }
    else{ /// lepton fail isolation
      /// check veto lepton if it satisfy signal requirements
      if (passHighPtCut && passHighPtSelection){
        ( *muon_itr )->auxdata< bool >( "signal_failIso" ) = true;
        nSignalLeptons_failIso++;
      }
      else{
        ( *muon_itr )->auxdata< bool >( "veto_failIso" ) = true;
        nVetoLeptons_failIso++;
      }
    }
  }
  
  std::map<string, unsigned int> outMap;
  outMap["nSignalLeptons"] = nSignalLeptons;
  outMap["nVetoLeptons"] = nVetoLeptons;
  outMap["nSignalLeptons_failIso"] = nSignalLeptons_failIso;
  outMap["nVetoLeptons_failIso"] = nVetoLeptons_failIso;
  
  return outMap;
 
}

std::map<string, unsigned int> RecoAnalysis :: SelectElectrons(
                                  xAOD::ElectronContainer* electrons,
                                  bool fillInCutflow){
  
  xAOD::ElectronContainer::iterator el_itr = electrons->begin();
  xAOD::ElectronContainer::iterator el_end = electrons->end();
  
  unsigned int nVetoLeptons = 0;
  unsigned int nSignalLeptons = 0;
  unsigned int nSignalLeptons_failIso = 0;
  unsigned int nVetoLeptons_failIso = 0;
  
  for( ; el_itr != el_end; ++el_itr ) {

    bool fillCutflow = true;

    /// Apply loosened vet selection first, fill in veto flags
    /// Left cut for signal selection apply in the end and fill signal flags
    ( *el_itr )->auxdata< bool >( "veto" ) = false;
    ( *el_itr )->auxdata< bool >( "signal" ) = false;
    ( *el_itr )->auxdata< bool >( "veto_failIso" ) = false;
    ( *el_itr )->auxdata< bool >( "signal_failIso" ) = false;
    
    if (fillCutflow) 
      m_BitsetCutflow->FillCutflow("oneElectron",fillInCutflow);
//     cout << "pt=" << (*el_itr)->pt() * 0.001 << "\tpass oneElectron" << endl;
    
    double elEta = (*el_itr)->caloCluster()->etaBE(2);
    if ( abs(elEta) > 2.47 || (abs(elEta) > 1.37 && abs(elEta) < 1.52)) continue;
    if (fillCutflow)
      m_BitsetCutflow->FillCutflow("Eta",fillInCutflow);
//     cout << "pt=" << (*el_itr)->pt() * 0.001 << "\tpass Eta" << endl;
    
    bool objQual = (*el_itr)->isGoodOQ(xAOD::EgammaParameters::BADCLUSELECTRON);
    if (!objQual) continue;
    if (fillCutflow)
      m_BitsetCutflow->FillCutflow("OQ",fillInCutflow);
//     cout << "pt=" << (*el_itr)->pt() * 0.001 << "\tpass OQ" << endl;
    
    double elPt = (*el_itr)->pt() * 0.001;

    double lowPtCutVeto = 20.0; /// GeV
    double lowPtCutSignal = 55.0; /// GeV
    double highPtCut = 999999.0; /// GeV
    
    /// fill cutflow only if pass highPtCut
    bool passHighPtCut = false;
    if ( elPt < lowPtCutVeto || elPt>=highPtCut ) continue;
    if ( ( elPt > lowPtCutSignal ) && ( elPt<highPtCut ) ){
      passHighPtCut = true;
    }
    fillCutflow = passHighPtCut;
    
    if (fillCutflow) 
      m_BitsetCutflow->FillCutflow( "el_pt",fillInCutflow );
//     cout << "pt=" << (*el_itr)->pt() * 0.001 << "\tpass el_pt" << endl;
    
    /// do significance 
    const xAOD::TrackParticle *tp = (*el_itr)->trackParticle();
    double d0_sig = xAOD::TrackingHelpers::d0significance
      ( tp, m_eventInfo->beamPosSigmaX(), m_eventInfo->beamPosSigmaY(), 
        m_eventInfo->beamPosSigmaXY() );
    if (TMath::Abs(d0_sig)>5.0) continue;
    if (fillCutflow)
      m_BitsetCutflow->FillCutflow("d0",fillInCutflow);
//     cout << "pt=" << (*el_itr)->pt() * 0.001 << "\tpass d0" << endl;
    
    bool passHighPtSelection = false;
    if ( !m_LHToolMedium2015->accept((*el_itr)) ) continue;
    if (m_LHToolTight2015->accept((*el_itr)))
      passHighPtSelection = true;
      
    fillCutflow = (passHighPtCut && passHighPtSelection);
    
    if (fillCutflow)
      m_BitsetCutflow->FillCutflow("ID",fillInCutflow);
    
    if (m_eleisolationSelectionTool->accept(**el_itr)){
      if (fillCutflow)
        m_BitsetCutflow->FillCutflow("Isolation",fillInCutflow);
      
      /// check veto lepton if it satisfy signal requirements
      if (passHighPtCut && passHighPtSelection){
        ( *el_itr )->auxdata< bool >( "signal" ) = true;
        nSignalLeptons++;
      }
      else{
        ( *el_itr )->auxdata< bool >( "veto" ) = true;
        nVetoLeptons++;
      }
    }
    else{ /// lepton fail isolation
      /// check veto lepton if it satisfy signal requirements
      if (passHighPtCut && passHighPtSelection){
        ( *el_itr )->auxdata< bool >( "signal_failIso" ) = true;
        nSignalLeptons_failIso++;
      }
      else{
        ( *el_itr )->auxdata< bool >( "veto_failIso" ) = true;
        nVetoLeptons_failIso++;
      }
    }
  }
  
  std::map<string, unsigned int> outMap;
  outMap["nSignalLeptons"] = nSignalLeptons;
  outMap["nVetoLeptons"] = nVetoLeptons;
  outMap["nSignalLeptons_failIso"] = nSignalLeptons_failIso;
  outMap["nVetoLeptons_failIso"] = nVetoLeptons_failIso;
  
  return outMap;
 
}

void RecoAnalysis :: plotPtBinnedHists(const xAOD::Muon* mu, xAOD::MissingET* finalTrkMet, string mainRootFileDirName, double weight){
 
  const int nPtBins = 11;
  double ptBinEdge[nPtBins] = {55,60,65,70,75,80,90,100,120,200,2000};
    
  ///********************************************************************
  /// save histos for each pT-bin
  double muonPt = mu->pt()*0.001;
  string ptPostfix = "";
  stringstream ss;
  
  for (int iBin=0; iBin<(nPtBins-1); iBin++){
    if (muonPt>=ptBinEdge[iBin] && muonPt<ptBinEdge[iBin+1]){
      
      
      ss.str(std::string());
      ss << ptBinEdge[iBin] << "-" << ptBinEdge[iBin+1];
            
      ptPostfix = "_ptBin_" + ss.str();
      m_HistObjectDumper->plotMuon(mu,mainRootFileDirName + ptPostfix,weight);
      m_HistObjectDumper->plotMtAndMet(mu,finalTrkMet,
                                       mainRootFileDirName + ptPostfix,weight);
    }
  }
  if (muonPt>=ptBinEdge[nPtBins-1]){
    ss.str(std::string());
    ss << (int)ptBinEdge[nPtBins-1];
    ptPostfix = "_ptBin_" + ss.str();
    m_HistObjectDumper->plotMuon(mu,mainRootFileDirName + ptPostfix,weight);
    m_HistObjectDumper->plotMtAndMet(mu,finalTrkMet,
                                     mainRootFileDirName + ptPostfix,weight);
  }
  ///******************************************************************** 
  
}
