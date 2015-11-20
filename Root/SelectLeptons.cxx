///*****************************************************************************
///
/// Custom user-defined functions, which are not standart in EventLoop:
///     - SelectMuons
///     - SelectElectrons
///
/// Both functions return pair(nSignalLeptons,nVetoLeptons)
///
///*****************************************************************************
#include <MyAnalysis/MyxAODAnalysis.h>

/// this is needed to distribute the algorithm to the workers
ClassImp(MyxAODAnalysis)

std::pair<unsigned int, unsigned int> MyxAODAnalysis :: SelectMuons(
                                  xAOD::MuonContainer* muons,
                                  xAOD::Vertex* primVertex,
                                  bool fillInCutflow){

  xAOD::MuonContainer::iterator muon_itr = muons->begin();
  xAOD::MuonContainer::iterator muon_end = muons->end();
  
  int nVetoLeptons = 0;
  int nSignalLeptons = 0;
  
  for( ; muon_itr != muon_end; ++muon_itr ) {
    
    bool fillCutflow = true;
    
    /// Apply loosened vet selection first, fill in veto flags
    /// Left cut for signal selection apply in the end and fill signal flags
    ( *muon_itr )->auxdata< bool >( "veto" ) = false;
    ( *muon_itr )->auxdata< bool >( "signal" ) = false;
    
    if ( ( *muon_itr )->muonType()!=xAOD::Muon_v1::Combined ) continue;
    if (fillCutflow) 
      m_BitsetCutflow->FillCutflow( "Combined",fillInCutflow );
    
    double muPt = ( ( *muon_itr )->pt() ) * 0.001;
    double lowPtCutVeto = 20.0; /// GeV
    double lowPtCutSignal = 55.0; /// GeV
    double highPtCut = 999999.0; /// GeV
    
    /// fill cutflow only if pass highPtCut
    bool passHighPtCut = false;
    if ( muPt < lowPtCutVeto || muPt>=highPtCut ) continue;
    if ( ( muPt > lowPtCutSignal ) && ( muPt<highPtCut ) ){
      passHighPtCut = true;
    }
    fillCutflow = passHighPtCut;
    
    if (fillCutflow) 
      m_BitsetCutflow->FillCutflow( "mu_pt",fillInCutflow );
    
    /// medium OR high-pT
    bool passHighPtSelection = false;
    if( ( !m_looseMuonSelection->accept( ( *muon_itr ) ) ) ||
      ( !m_muonSelection->accept( ( *muon_itr ) ) )
    ) continue;
    
    if( m_muonSelection->accept( ( *muon_itr ) ) ){
      passHighPtSelection = true;
    }
    bool fillCutflow = (passHighPtCut && passHighPtSelection);
    
    if (fillCutflow) 
      m_BitsetCutflow->FillCutflow( "MCP selector",fillInCutflow );
    
    /// do significance 
    const xAOD::TrackParticle *tp = ( *muon_itr )->primaryTrackParticle();
    double d0_sig = xAOD::TrackingHelpers::d0significance
      ( tp, eventInfo->beamPosSigmaX(), eventInfo->beamPosSigmaY(), 
        eventInfo->beamPosSigmaXY() );
    if (d0_sig>3.0) continue;
    if (fillCutflow)
      m_BitsetCutflow->FillCutflow("d0",fillInCutflow);
    
    /// zo cut
    double z0_vrtPVx = ( *muon_itr )->primaryTrackParticle()->z0() +
    ( *muon_itr )->primaryTrackParticle()->vz() - primVertex->z();
    double sintheta = 1.0/TMath::CosH(( *muon_itr )->eta());
    if (abs( z0_vrtPVx*sintheta )>10.0) continue;
    if (fillCutflow)
      m_BitsetCutflow->FillCutflow( "z0",fillInCutflow );
    
    if (!m_muonisolationSelectionTool->accept(( **muon_itr ))) continue;
    if (fillCutflow)
      m_BitsetCutflow->FillCutflow( "Isolation",fillInCutflow );
    
    /// all requirement above satisfy veto lepton
    ( *muon_itr )->auxdata< bool >( "veto" ) = true;
    nVetoLeptons++;
    
    /// check veto lepton if it satisfy signal requirements
    if (passHighPtCut && passHighPtSelection){
      ( *muon_itr )->auxdata< bool >( "signal" ) = true;
      nSignalLeptons++;
    }
  }
  
  /// all signal leptons are counted as veto...
  /// subtract them to get exclusive veto leptons
  std::pair<unsigned int, unsigned int> outPair (nSignalLeptons,
                                            nVetoLeptons-nSignalLeptons);
  
  return outPair;
  
}



std::pair<unsigned int, unsigned int> MyxAODAnalysis :: SelectElectrons(
                                  xAOD::ElectronContainer* electrons,
                                  bool fillInCutflow){
  
  xAOD::ElectronContainer::iterator el_itr = electrons->begin();
  xAOD::ElectronContainer::iterator el_end = electrons->end();
  
  int nVetoLeptons = 0;
  int nSignalLeptons = 0;
  
  for( ; el_itr != el_end; ++el_itr ) {

    bool fillCutflow = true;

    /// Apply loosened vet selection first, fill in veto flags
    /// Left cut for signal selection apply in the end and fill signal flags
    ( *el_itr )->auxdata< bool >( "veto" ) = false;
    ( *el_itr )->auxdata< bool >( "signal" ) = false;
    
    if (fillCutflow) 
      m_BitsetCutflow->FillCutflow("oneElectron",fillInCutflow);
    
    double elEta = (*el_itr)->caloCluster()->etaBE(2);
    if ( abs(elEta) > 2.47 || (abs(elEta) > 1.37 && abs(elEta) < 1.52)) continue;
    if (fillCutflow)
      m_BitsetCutflow->FillCutflow("Eta",fillInCutflow);
    
    bool objQual = (*el_itr)->isGoodOQ(xAOD::EgammaParameters::BADCLUSELECTRON);
    if (!objQual) continue;
    if (fillCutflow)
      m_BitsetCutflow->FillCutflow("OQ",fillInCutflow);
    
    double elPt = (*el_itr)->pt() * 0.001;

    double lowPtCutVeto = 20.0; /// GeV
    double lowPtCutSignal = 65.0; /// GeV
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
    
    /// do significance 
    const xAOD::TrackParticle *tp = (*el_itr)->trackParticle();
    double d0_sig = xAOD::TrackingHelpers::d0significance
      ( tp, eventInfo->beamPosSigmaX(), eventInfo->beamPosSigmaY(), 
        eventInfo->beamPosSigmaXY() );
    if (d0_sig>5.0) continue;
    if (fillCutflow)
      m_BitsetCutflow->FillCutflow("d0",fillInCutflow);
    
    bool passHighPtSelection = false;
    if ( !m_LHToolMedium2015->accept((*el_itr)) ) continue;
    if (m_LHToolTight2015->accept((*el_itr)))
      passHighPtSelection = true;
      
    bool fillCutflow = (passHighPtCut && passHighPtSelection);
    
    if (fillCutflow)
      m_BitsetCutflow->FillCutflow("ID",fillInCutflow);
    
    if (!m_eleisolationSelectionTool->accept(**el_itr)) continue;
    m_BitsetCutflow->FillCutflow("Isolation",fillInCutflow);
    
    /// all requirement above satisfy veto lepton
    ( *el_itr )->auxdata< bool >( "veto" ) = true;
    nVetoLeptons++;
    
    /// check veto lepton if it satisfy signal requirements
    if (passHighPtCut && passHighPtSelection){
      ( *el_itr )->auxdata< bool >( "signal" ) = true;
      nSignalLeptons++;
    }
  }
  
  /// all signal leptons are counted as veto...
  /// subtract them to get exclusive veto leptons
  std::pair<unsigned int, unsigned int> outPair (nSignalLeptons,
                                            nVetoLeptons-nSignalLeptons);
  
  return outPair; 
  
}


