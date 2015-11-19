///*****************************************************************************
///
/// Custom user-defined functions, which are not standart in EventLoop:
///     - VetoMuon
///     - SelectMuon
///     - passMuonSelection
///     - SelectElectron
///
///*****************************************************************************
#include <MyAnalysis/MyxAODAnalysis.h>

/// this is needed to distribute the algorithm to the workers
ClassImp(MyxAODAnalysis)

ConstDataVector<xAOD::MuonContainer> MyxAODAnalysis :: VetoMuon(const xAOD::MuonContainer* muons,
                                         xAOD::Vertex* primVertex,
                                  bool lookForVetoMuon){

  /// option 1: signal selection; 2: muon veto; 3: 
  /// would need to fill 

  const char* APP_NAME = "MyxAODAnalysis";

  const xAOD::EventInfo* eventInfo = 0;
  if( ! m_event->retrieve( eventInfo, "EventInfo").isSuccess() ){
    Error("execute()", "Failed to retrieve event info collection. Exiting." );
  }

  ConstDataVector<xAOD::MuonContainer> vetoMuonContainer(SG::VIEW_ELEMENTS);
  xAOD::Muon* vetoMuon = 0;
  int muonCounter = 0;
  int nmuon=0;
  //static SG::AuxElement::Accessor< xAOD::IParticleLink > accSetOriginLink ("originalObjectLink");
  xAOD::MuonContainer::const_iterator muon_itr = muons->begin();
  xAOD::MuonContainer::const_iterator muon_end = muons->end();
  if (!lookForVetoMuon) {
    for( ; muon_itr != muon_end; ++muon_itr ) {
      xAOD::Muon* mu = 0;
      mu = const_cast<xAOD::Muon*> (*muon_itr);
      //cout << "muon type " << mu->muonType()<<endl;
      if (mu->muonType()==xAOD::Muon_v1::Combined) nmuon+=1;
    }
    //cout << "nmuon " << nmuon<< " size " << muons->size()<<endl;
    if (nmuon>0) {
      count[3]+=1;
    } else {
      return vetoMuonContainer;
    }
  }
  muon_itr = muons->begin();
  for( ; muon_itr != muon_end; ++muon_itr ) {

    if (m_useHistObjectDumper) m_HistObjectDumper->plotMuon((*muon_itr),"noCuts");

    m_BitsetCutflow->FillCutflow("oneMuon",!lookForVetoMuon);

    xAOD::Muon* mu = 0;
    mu = const_cast<xAOD::Muon*> (*muon_itr);

    /// keep lins

    if (mu->muonType()!=xAOD::Muon_v1::Combined) continue;
    m_BitsetCutflow->FillCutflow("Combined",!lookForVetoMuon);
    if (!lookForVetoMuon) count[4]+=1;
    double muPt = (mu->pt()) * 0.001;
    double lowPtCut = 55.0; /// GeV
    //double lowPtCut = 20.0; /// GeV
    double highPtCut = 9999999.9; /// GeV
    if (lookForVetoMuon){
      lowPtCut = 20.0;
      highPtCut = 999999.0;
    }

    if (muPt < lowPtCut || muPt>=highPtCut) {
      continue; /// veto muon
      //cout << "muon rejected in pt cut " << muPt<<endl;
    }
    m_BitsetCutflow->FillCutflow("mu_pt",!lookForVetoMuon);
    if (!lookForVetoMuon) count[5]+=1;

    if (!lookForVetoMuon) {
      if(!m_muonSelection->accept(mu)) continue;
      m_BitsetCutflow->FillCutflow("MCP selector",!lookForVetoMuon);
      count[6]+=1;
    } else {
      if(!m_loosemuonSelection->accept(mu)) continue;
    }
    /// do significance 
//     double d0_sig = TMath::Abs(mu->primaryTrackParticle()->d0()) /
//     TMath::Sqrt(mu->primaryTrackParticle()->definingParametersCovMatrix()(0,0)
//     + eventInfo->beamPosSigmaX()*eventInfo->beamPosSigmaX() );
//     if (d0_sig>3.0) continue;

    xAOD::TrackParticle *tp = mu->primaryTrackParticle();
    double d0_sig = xAOD::TrackingHelpers::d0significance
    ( tp, eventInfo->beamPosSigmaX(), eventInfo->beamPosSigmaY(), 
      eventInfo->beamPosSigmaXY() );
    if (d0_sig>3.0) continue;
    m_BitsetCutflow->FillCutflow("d0",!lookForVetoMuon);
    if (!lookForVetoMuon) count[8]+=1;


    /// zo cut
    double z0_vrtPVx = mu->primaryTrackParticle()->z0() +
    mu->primaryTrackParticle()->vz() - primVertex->z();
    double sintheta = 1.0/TMath::CosH(mu->eta());
    if (abs( z0_vrtPVx*sintheta )>10.0) continue;
    m_BitsetCutflow->FillCutflow("z0",!lookForVetoMuon);
    if (!lookForVetoMuon) count[9]+=1;

    /// Isolation stuff
    if (!lookForVetoMuon) {
      if (!m_muonisolationSelectionTool->accept(*mu)) continue;
      m_BitsetCutflow->FillCutflow("Isolation",!lookForVetoMuon);
      count[10]+=1;
    }

    /// check truth
    if(isMC){

    }
    /// save this muon
    //const xAOD::IParticleLink originLink( *muons, mu->index() );
    //accSetOriginLink(*mu) = originLink;
    vetoMuonContainer.push_back(mu);
    //cout << "veto muon found with pt " << mu->pt()/1000.<<endl;
    vetoMuon = mu;
    muonCounter++;

  } /// end for loop over muons

  /// if there are two or more signal muons - return NULL
  /// Veto muon is defined to have pT within 20..55 GeV range
  /// that's why veto muon selection is not overlapping with signal region
  /// if more than two veto muons found - return any of them
  //  if (muonCounter==1 || lookForVetoMuon) 
  return vetoMuonContainer;
//else
//  return 0;
}



ConstDataVector<xAOD::MuonContainer> MyxAODAnalysis :: SelectMuon(const xAOD::MuonContainer* muons,
                                         xAOD::Vertex* primVertex,
                                  bool lookForVetoMuon){

  // option 1: signal selection; 2: muon veto; 3: 
  // would need to fill 
  //ConstDataVector<xAOD::MuonContainer> muonfinal(SG::VIEW_ELEMENTS); // Create a new muon container
  const char* APP_NAME = "MyxAODAnalysis";
  //nmuons = 0;

  const xAOD::EventInfo* eventInfo = 0;
  if( ! m_event->retrieve( eventInfo, "EventInfo").isSuccess() ){
    Error("execute()", "Failed to retrieve event info collection. Exiting." );
  }

  //  xAOD::MuonContainer* outMuonContainer = 0;

  ConstDataVector<xAOD::MuonContainer> outMuonContainer(SG::VIEW_ELEMENTS);
  xAOD::Muon* outMuon = 0;
  int muonCounter = 0;
  int nmuon=0;
  static SG::AuxElement::Accessor< xAOD::IParticleLink > accSetOriginLink ("originalObjectLink");
  xAOD::MuonContainer::const_iterator muon_itr = muons->begin();
  xAOD::MuonContainer::const_iterator muon_end = muons->end();
  if (!lookForVetoMuon) {
    for( ; muon_itr != muon_end; ++muon_itr ) {
      xAOD::Muon* mu = 0;
      mu = const_cast<xAOD::Muon*> (*muon_itr);
      //cout << "muon type " << mu->muonType()<<endl;
      if (mu->muonType()==xAOD::Muon_v1::Combined) nmuon+=1;
    }
    //cout << "nmuon " << nmuon<< " size " << muons->size()<<endl;
    if (nmuon>0) {
      count[3]+=1;
    } else {
      return outMuonContainer;
    }
  }


  muon_itr = muons->begin();
  for( ; muon_itr != muon_end; ++muon_itr ) {

    if (m_useHistObjectDumper) m_HistObjectDumper->plotMuon((*muon_itr),"noCuts");

    m_BitsetCutflow->FillCutflow("oneMuon",!lookForVetoMuon);

    xAOD::Muon* mu = 0;
    /// TODO return m_useMuonCalibrationAndSmearingTool flag
//     if (m_useMuonCalibrationAndSmearingTool){
//       if( !m_muonCalibrationAndSmearingTool->correctedCopy( **muon_itr, mu ) ) {
//         Error(APP_NAME, "Cannot really apply calibration nor smearing");
//         continue;
//       }
//     }
//     else{
//       mu = const_cast<xAOD::Muon*> (*muon_itr);
//     }
    mu = const_cast<xAOD::Muon*> (*muon_itr);

    // keep lins

    if (mu->muonType()!=xAOD::Muon_v1::Combined) continue;
    m_BitsetCutflow->FillCutflow("Combined",!lookForVetoMuon);
    if (!lookForVetoMuon) count[4]+=1;
    double muPt = (mu->pt()) * 0.001;
    //double lowPtCut = 55.0; /// GeV
    double lowPtCut = 55.0; /// GeV
    double highPtCut = 99999.9; /// GeV
    if (lookForVetoMuon){
      lowPtCut = 20.0;
      highPtCut = 999999.0;
    }

    if (muPt < lowPtCut || muPt>=highPtCut) continue; /// veto muon
    m_BitsetCutflow->FillCutflow("mu_pt",!lookForVetoMuon);
    if (!lookForVetoMuon) count[5]+=1;

    if (!lookForVetoMuon) {
      if(!m_muonSelection->accept(mu)) continue;
      m_BitsetCutflow->FillCutflow("MCP selector",!lookForVetoMuon);
      count[6]+=1;
    } else {
      if(!m_loosemuonSelection->accept(mu)) continue;
    }
    /// do significance 
    xAOD::TrackParticle *tp = mu->primaryTrackParticle();
    double d0_sig = xAOD::TrackingHelpers::d0significance
    ( tp, eventInfo->beamPosSigmaX(), eventInfo->beamPosSigmaY(), 
      eventInfo->beamPosSigmaXY() );
    if (d0_sig>3.0) continue;
    m_BitsetCutflow->FillCutflow("d0",!lookForVetoMuon);
    if (!lookForVetoMuon) count[8]+=1;

    /// zo cut
    double z0_vrtPVx = mu->primaryTrackParticle()->z0() +
    mu->primaryTrackParticle()->vz() - primVertex->z();
    double sintheta = 1.0/TMath::CosH(mu->eta());
    if (abs( z0_vrtPVx*sintheta )>10.0) continue;
    m_BitsetCutflow->FillCutflow("z0",!lookForVetoMuon);
    if (!lookForVetoMuon) count[9]+=1;

    /// Isolation stuff

    if (!lookForVetoMuon) {
      if (!m_muonisolationSelectionTool->accept(*mu)) continue;
      m_BitsetCutflow->FillCutflow("Isolation",!lookForVetoMuon);
      count[10]+=1;
    }

    // check truth
    if(isMC){
    }
    // save this muon
    const xAOD::IParticleLink originLink( *muons, mu->index() );
    accSetOriginLink(*mu) = originLink;
    outMuonContainer.push_back(mu);
    //cout << "signal muon found with pt " << mu->pt()/1000.<<endl;
    outMuon = mu;
    muonCounter++;

  } /// end for loop over muons

  /// if there are two or more signal muons - return NULL
  /// Veto muon is defined to have pT within 20..55 GeV range
  /// that's why veto muon selection is not overlapping with signal region
  /// if more than two veto muons found - return any of them
  //  if (muonCounter==1 || lookForVetoMuon) 
  return outMuonContainer;
//else
//  return 0;
}

bool MyxAODAnalysis :: passMuonSelection(const xAOD::Muon* mu,
                     const xAOD::EventInfo* eventInfo,
                                         xAOD::Vertex* primVertex,
                     bool lookForVetoMuon){

  const char* APP_NAME = "MyxAODAnalysis";
  bool pass = false;

  if (mu->muonType()!=xAOD::Muon_v1::Combined) return pass;
  m_BitsetCutflow->FillCutflow("Combined",!lookForVetoMuon);
  if (!lookForVetoMuon) count[4]+=1;
  double muPt = (mu->pt()) * 0.001;
  double lowPtCut = 55.0; /// GeV
  double highPtCut = 99999.9; /// GeV
  if (lookForVetoMuon){
    lowPtCut = 20.0;
    highPtCut = 999999.0;
  }

  if (muPt < lowPtCut || muPt>=highPtCut) return pass; /// veto muon
  m_BitsetCutflow->FillCutflow("mu_pt",!lookForVetoMuon);
  if (!lookForVetoMuon) count[5]+=1;

  if (!lookForVetoMuon) {
    if(!m_muonSelection->accept(mu)) return pass;
    m_BitsetCutflow->FillCutflow("MCP selector",!lookForVetoMuon);
    count[6]+=1;
  } else {
    if(!m_loosemuonSelection->accept(mu)) return pass;
  }
  /// do significance 
  xAOD::TrackParticle *tp = mu->primaryTrackParticle();
  double d0_sig = xAOD::TrackingHelpers::d0significance
    ( tp, eventInfo->beamPosSigmaX(), eventInfo->beamPosSigmaY(), 
      eventInfo->beamPosSigmaXY() );
  if (d0_sig>3.0) return pass;
  m_BitsetCutflow->FillCutflow("d0",!lookForVetoMuon);
  if (!lookForVetoMuon) count[8]+=1;

  /// zo cut
  double z0_vrtPVx = mu->primaryTrackParticle()->z0() +
    mu->primaryTrackParticle()->vz() - primVertex->z();
  double sintheta = 1.0/TMath::CosH(mu->eta());
  if (abs( z0_vrtPVx*sintheta )>10.0) return pass;
  m_BitsetCutflow->FillCutflow("z0",!lookForVetoMuon);
  if (!lookForVetoMuon) count[9]+=1;
  /// Isolation stuff

  if (!lookForVetoMuon) {
    if (!m_muonisolationSelectionTool->accept(*mu)) return pass;
    m_BitsetCutflow->FillCutflow("Isolation",!lookForVetoMuon);
    count[10]+=1;
  }

  // check truth
  if(isMC){
  }
  
  pass = true;
  //cout << "signal muon found with pt " << mu->pt()/1000.<<endl;

  return pass;
}

ConstDataVector<xAOD::ElectronContainer> MyxAODAnalysis :: SelectElectron(const xAOD::ElectronContainer* electrons, xAOD::Vertex* primVertex, bool lookForVetoElectron){

  const char* APP_NAME = "MyxAODAnalysis";

    const xAOD::EventInfo* eventInfo = 0;
  if( ! m_event->retrieve( eventInfo, "EventInfo").isSuccess() ){
    Error("execute()", "Failed to retrieve event info collection. Exiting." );
  }

  //cout << "in selectElectron"<<endl;
  xAOD::Electron* outElectron = 0;
  ConstDataVector<xAOD::ElectronContainer> outElectronContainer(SG::VIEW_ELEMENTS);
  int electronCounter = 0;
  xAOD::ElectronContainer::const_iterator electron_itr = electrons->begin();
  xAOD::ElectronContainer::const_iterator electron_end = electrons->end();
  int i = 0;
  for( ; electron_itr != electron_end; ++electron_itr ) {

    m_BitsetCutflow->FillCutflow("oneElectron",!lookForVetoElectron);

    /// Calibration and corrections ...
    xAOD::Electron* el = 0;
    xAOD::Electron* el2 = 0;
    el2 = const_cast<xAOD::Electron*> (*electron_itr);
    float oldpt = el2->pt()/1000.;

    el = const_cast<xAOD::Electron*> (*electron_itr);
      //}

    i++;
    /// ... Calibration and corrections
    /// Eta
    double Eta = el->caloCluster()->eta();
    if ( abs(Eta) > 2.47 || (abs(Eta) > 1.37 && abs(Eta) < 1.52)) continue;
    m_BitsetCutflow->FillCutflow("Eta",!lookForVetoElectron);

    /// OQ
    if (el->isGoodOQ(xAOD::EgammaParameters::BADCLUSELECTRON)==false) continue;
    m_BitsetCutflow->FillCutflow("OQ",!lookForVetoElectron);

    /// pT cut ...
    double elPt = (el->pt()) * 0.001;
    double lowPtCut = 65.0; /// GeV
    if (lookForVetoElectron){
      lowPtCut = 20.0;
    }
//cout << "pt " << elPt<< " old " << oldpt<<endl;
    if (elPt < lowPtCut) continue; /// veto electron
    m_BitsetCutflow->FillCutflow("el_pt",!lookForVetoElectron);
    /// ... pT cut

    /// d0 significance ...

    xAOD::TrackParticle *tp = el->primaryTrackParticle();
    double d0_sig = xAOD::TrackingHelpers::d0significance
    ( tp, eventInfo->beamPosSigmaX(), eventInfo->beamPosSigmaY(), 
      eventInfo->beamPosSigmaXY() );

   if (d0_sig>5.0) continue;
    m_BitsetCutflow->FillCutflow("d0",!lookForVetoElectron);
    /// ... d0 significance

    /// ID ...
    bool goodID = false;
    if (lookForVetoElectron)
      goodID = m_LHToolMedium2015->accept(el) || (m_LHToolTight2015->accept(el));
    //goodID = m_LHToolMedium2015->accept(el);
    else
      goodID = m_LHToolMedium2015->accept(el) &&
      m_LHToolTight2015->accept(el);

    if (!goodID) continue;
    m_BitsetCutflow->FillCutflow("ID",!lookForVetoElectron);
    /// ... ID

    /// Isolation
    if (!m_eleisolationSelectionTool->accept(*el)) continue;
    m_BitsetCutflow->FillCutflow("Isolation",!lookForVetoElectron);

    outElectron = el;
    outElectronContainer.push_back(el);
    electronCounter++;

  } /// end for loop over electrons

  /// if there are two or more signal electrons - return NULL
  /// Veto electron is defined to have pT within 20..55 GeV range
  /// that's why veto electron selection is not overlapping with signal region
  /// if more than two veto electrons found - return any of them
  return outElectronContainer;

}


/*

xAOD::Muon* MyxAODAnalysis :: SelectMuon(const xAOD::MuonContainer* muons, 
                                         xAOD::Vertex* primVertex, 
                                         bool lookForVetoMuon){

  const char* APP_NAME = "MyxAODAnalysis";
  
  const xAOD::EventInfo* eventInfo = 0;
  if( ! m_event->retrieve( eventInfo, "EventInfo").isSuccess() ){
    Error("execute()", "Failed to retrieve event info collection. Exiting." );
  }
  
  xAOD::Muon* outMuon = 0;
  int muonCounter = 0;
  xAOD::MuonContainer::const_iterator muon_itr = muons->begin();
  xAOD::MuonContainer::const_iterator muon_end = muons->end();
  for( ; muon_itr != muon_end; ++muon_itr ) {
  
    if (m_useHistObjectDumper) m_HistObjectDumper->plotMuon((*muon_itr),"noCuts");
  
    m_BitsetCutflow->FillCutflow("oneMuon",!lookForVetoMuon);
  
    xAOD::Muon* mu = 0;
    if (m_useCalibrationAndSmearingTool){
      if( !m_muonCalibrationAndSmearingTool->correctedCopy( **muon_itr, mu ) ) {
        Error(APP_NAME, "Cannot really apply calibration nor smearing");
        continue;
      }
    }
    else{
      mu = const_cast<xAOD::Muon*> (*muon_itr);
    }
    
    if (mu->muonType()!=xAOD::Muon_v1::Combined) continue;
    m_BitsetCutflow->FillCutflow("Combined",!lookForVetoMuon);
    
    double muPt = (mu->pt()) * 0.001;
    double lowPtCut = 55.0; /// GeV
    double highPtCut = 99999.9; /// GeV
    if (lookForVetoMuon){
      lowPtCut = 20.0;
      highPtCut = 55.0;
    }
      
    if (muPt < lowPtCut || muPt>=highPtCut) continue; /// veto muon
    m_BitsetCutflow->FillCutflow("mu_pt",!lookForVetoMuon);
    
    if(!m_muonSelection->accept(mu)) continue;
    m_BitsetCutflow->FillCutflow("MCP selector",!lookForVetoMuon);
        
    if (!m_muonSelection->passedHighPtCuts(*mu)) continue;
    m_BitsetCutflow->FillCutflow("MS Hits",!lookForVetoMuon);
    
    /// do significance 
    double d0_sig = TMath::Abs(mu->primaryTrackParticle()->d0()) / 
    TMath::Sqrt(mu->primaryTrackParticle()->definingParametersCovMatrix()(0,0)
    + eventInfo->beamPosSigmaX()*eventInfo->beamPosSigmaX() );
    if (d0_sig>3.0) continue;
    m_BitsetCutflow->FillCutflow("d0",!lookForVetoMuon);
    
    /// zo cut
    double z0_vrtPVx = mu->primaryTrackParticle()->z0() + 
    mu->primaryTrackParticle()->vz() - primVertex->z(); 
    double sintheta = 1.0/TMath::CosH(mu->eta());
    if (abs( z0_vrtPVx*sintheta )>10.0) continue;
    m_BitsetCutflow->FillCutflow("z0",!lookForVetoMuon);
    
    /// Isolation stuff
    if (!m_isolationSelectionTool->accept(*mu)) continue;
    m_BitsetCutflow->FillCutflow("Isolation",!lookForVetoMuon);
    
    outMuon = mu;
    muonCounter++;
      
    
  
  } /// end for loop over muons

  /// if there are two or more signal muons - return NULL
  /// Veto muon is defined to have pT within 20..55 GeV range
  /// that's why veto muon selection is not overlapping with signal region
  /// if more than two veto muons found - return any of them
  if (muonCounter==1 || lookForVetoMuon) return outMuon;
  else
    return 0;

}



xAOD::Electron* MyxAODAnalysis :: SelectElectron(const xAOD::ElectronContainer* electrons, 
                                         xAOD::Vertex* primVertex, 
                                         bool lookForVetoElectron){

  const char* APP_NAME = "MyxAODAnalysis";
  
  const xAOD::EventInfo* eventInfo = 0;
  if( ! m_event->retrieve( eventInfo, "EventInfo").isSuccess() ){
    Error("execute()", "Failed to retrieve event info collection. Exiting." );
  }
  
  xAOD::Electron* outElectron = 0;
  int electronCounter = 0;
  xAOD::ElectronContainer::const_iterator electron_itr = electrons->begin();
  xAOD::ElectronContainer::const_iterator electron_end = electrons->end();
  for( ; electron_itr != electron_end; ++electron_itr ) {
  
    m_BitsetCutflow->FillCutflow("oneElectron",!lookForVetoElectron);
  
    /// Calibration and corrections ...
    xAOD::Electron* el = 0;
    if (m_useCalibrationAndSmearingTool){
//       if( !m_electronCalibrationAndSmearingTool->correctedCopy( **electron_itr, el ) ) {
//         Error(APP_NAME, "Cannot really apply calibration nor smearing");
//         continue;
//       }
    }
    else{
      el = const_cast<xAOD::Electron*> (*electron_itr);
    }
    /// ... Calibration and corrections
    
    /// Eta
    double Eta = el->caloCluster()->eta();
    if ( abs(Eta) > 2.47 || (1.37 > abs(Eta) < 1.52)) continue;
    m_BitsetCutflow->FillCutflow("Eta",!lookForVetoElectron);
    
    /// OQ
    if (el->isGoodOQ(xAOD::EgammaParameters::BADCLUSELECTRON)==false) continue;
    m_BitsetCutflow->FillCutflow("OQ",!lookForVetoElectron);
    
    /// pT cut ...
    double elPt = (el->pt()) * 0.001;
    double lowPtCut = 65.0; /// GeV
    if (lookForVetoElectron){
      lowPtCut = 20.0;
    }
    if (elPt < lowPtCut) continue; /// veto electron
    m_BitsetCutflow->FillCutflow("el_pt",!lookForVetoElectron);
    /// ... pT cut
    
    /// d0 significance ...
    
    const xAOD::TrackParticle *trk = el->trackParticle();
    float d0_sig = fabs(trk->d0())/sqrt(trk->definingParametersCovMatrix()(0,0)
    + eventInfo->beamPosSigmaX()*eventInfo->beamPosSigmaX() );
    if (d0_sig>5.0) continue;
    m_BitsetCutflow->FillCutflow("d0",!lookForVetoElectron);
    /// ... d0 significance
    
    /// ID ...
    bool goodID = false;
    if (lookForVetoElectron)
      goodID = m_LHToolMedium2015->accept(el) && (!m_LHToolTight2015->accept(el));
    else
      goodID = m_LHToolMedium2015->accept(el) && 
      m_LHToolTight2015->accept(el);
    
    if (!goodID) continue;
    m_BitsetCutflow->FillCutflow("ID",!lookForVetoElectron);
    /// ... ID
    
    /// Isolation
    if (!m_isolationSelectionTool->accept(*el)) continue;
    m_BitsetCutflow->FillCutflow("Isolation",!lookForVetoElectron);
    
    outElectron = el;
    electronCounter++;
    
  } /// end for loop over electrons

  /// if there are two or more signal electrons - return NULL
  /// Veto electron is defined to have pT within 20..55 GeV range
  /// that's why veto electron selection is not overlapping with signal region
  /// if more than two veto electrons found - return any of them
  if (electronCounter==1 || lookForVetoElectron) return outElectron;
  else
    return 0;

}
*/

