///*****************************************************************************
///
/// Implementation of execute() function.
///
///*****************************************************************************

#include <MyAnalysis/MyxAODAnalysis.h>

/// this is needed to distribute the algorithm to the workers
ClassImp(MyxAODAnalysis)

EL::StatusCode MyxAODAnalysis :: execute ()
{
  
  const char* APP_NAME = "MyxAODAnalysis";

  /// push cutflow bitset to cutflow hist
  m_BitsetCutflow->PushBitSet();
  
  if( (m_eventCounter % 1000) ==0 ) Info("execute()", "Event number = %i", 
    m_eventCounter );
  m_eventCounter++;
  
  count[0]+=1;
  
  ///----------------------------
  /// Event information
  ///--------------------------- 

//   const xAOD::EventInfo* eventInfo = 0;
  if( ! m_event->retrieve( eventInfo, "EventInfo").isSuccess() ){
    Error("execute()", "Failed to retrieve event info collection. Exiting." );
    return EL::StatusCode::FAILURE;
  }

  EventNumber = eventInfo->eventNumber();  

  bool isMC = false;
  if(eventInfo->eventType( xAOD::EventInfo::IS_SIMULATION ) ){
    isMC = true;
  }
  
  /// Muon Truth matching. Check do we have muon from W' decay
  /// WARNING do we need it for all MCs?
  if(isMC && m_doWprimeTruthMatching){
    bool foundMuonFromWprimeDecay = false;
    
    /// Create truth vertice container
    const xAOD::TruthVertexContainer* truthVertices = 0;
    if ( !m_event->retrieve( truthVertices, "TruthVertices" ).isSuccess() ){ 
      Error("execute()","Failed to retrieve TruthVertices container. Exiting.");
      return EL::StatusCode::FAILURE;
    }
    
    /// Start iterating over truth container
    xAOD::TruthVertexContainer::const_iterator truthV_itr; 
    for (truthV_itr = truthVertices->begin(); truthV_itr != truthVertices->end();
         ++truthV_itr ) {
      if (foundMuonFromWprimeDecay)
        break;
      for (unsigned int iIn=0; iIn < (*truthV_itr)->nIncomingParticles(); iIn++)
      {
        /// 34 - Wprime
        if (TMath::Abs((*truthV_itr)->incomingParticle(iIn)->pdgId()) == 34) { 
          
          for (unsigned int iOut=0; iOut < (*truthV_itr)->nOutgoingParticles(); 
               iOut++) {
            if (TMath::Abs((*truthV_itr)->outgoingParticle(iOut)->pdgId()) 
              == 13) /// 13 - muon
              foundMuonFromWprimeDecay = true;
            break;
          }
        }
        if (foundMuonFromWprimeDecay)
          break;
      }
    }
    
    if (!foundMuonFromWprimeDecay) return EL::StatusCode::SUCCESS;
    m_BitsetCutflow->FillCutflow("Truth muon decay");
  }
  
  /// pass GRL?
  if(!isMC){
    if(!m_grl->passRunLB(*eventInfo)){
      return EL::StatusCode::SUCCESS;
    }
  }
  m_BitsetCutflow->FillCutflow("GRL");
  
  ///------------------------------------------------------------
  /// Apply event cleaning to remove events due to 
  /// problematic regions of the detector, and incomplete events.
  /// Apply to data.
  ///------------------------------------------------------------
  /// reject event if:
  if(!isMC){
    if(  (eventInfo->errorState(xAOD::EventInfo::LAr)==xAOD::EventInfo::Error )
      || (eventInfo->errorState(xAOD::EventInfo::Tile)==xAOD::EventInfo::Error )
      || (eventInfo->errorState(xAOD::EventInfo::SCT)==xAOD::EventInfo::Error )
      || (eventInfo->isEventFlagBitSet(xAOD::EventInfo::Core, 18) )  )
    {
      return EL::StatusCode::SUCCESS; 
    }
  }
  m_numCleanEvents++;
  m_BitsetCutflow->FillCutflow("EventCleaning");
  
  /// Primary vertex
  const xAOD::VertexContainer* vertices = 0;
  if ( !m_event->retrieve( vertices, "PrimaryVertices" ).isSuccess() ){ 
    Error("execute()","Failed to retrieve PrimaryVertices container. Exiting.");
    return EL::StatusCode::FAILURE;
  }
  
  xAOD::VertexContainer::const_iterator vtx_itr = vertices->begin();
  xAOD::VertexContainer::const_iterator vtx_end = vertices->end();
  xAOD::Vertex* primVertex = 0;
  int nGoodVtx = 0;
  for( ; vtx_itr != vtx_end; ++vtx_itr ) {
    if ((*vtx_itr)->vertexType()==xAOD::VxType::PriVtx){
      primVertex = (*vtx_itr);
      nGoodVtx++;
    }
  }
  
  if (nGoodVtx>1)
    cout << "WARNING!!!! Found more then one prim.vertex: nGoodVtx = " 
    << nGoodVtx << endl;
  if (nGoodVtx==0)
    return EL::StatusCode::SUCCESS;
  m_BitsetCutflow->FillCutflow("Primary vertex");
  
  /// triggers  
  /// list of triggers to use
  std::vector<std::string> triggerChains = {"HLT_mu50*"};
  if (m_runElectronChannel){
    triggerChains.erase (triggerChains.begin()+1);
    triggerChains.push_back("HLT_e24_lhmedium_L1EM20VH");
    triggerChains.push_back("HLT_e60_lhmedium*");
    triggerChains.push_back("HLT_e120_lhloose*");
  }
  
  bool passTriggerOR = false;
  
  for(std::vector<std::string>::iterator it = triggerChains.begin(); it != 
    triggerChains.end(); ++it) {
    
    auto chainGroup = m_trigDecisionTool->getChainGroup(*it);
    for(auto &trig : chainGroup->getListOfTriggers()) {
      auto cg = m_trigDecisionTool->getChainGroup(trig);
      if (cg->isPassed()==true){
        passTriggerOR = true;
        m_BitsetCutflow->FillCutflow(*it);
      }
    }
  }
  
  if(m_doNotApplyTriggerCuts){
    if (passTriggerOR==false)
      return EL::StatusCode::SUCCESS;
    m_BitsetCutflow->FillCutflow("Trigger");
  }
  

  const xAOD::MuonContainer* muons = 0;
  if ( !m_event->retrieve( muons, "Muons" ).isSuccess() ){ 
    Error("execute()", "Failed to retrieve Muons container. Exiting." );
    return EL::StatusCode::FAILURE;
  }
  
  std::pair<xAOD::MuonContainer*,xAOD::ShallowAuxContainer*> classifiedMuons = 
  xAOD::shallowCopyContainer(*muons);
  xAOD::setOriginalObjectLink(*muons, *classifiedMuons.first); 
  for(const auto& muon : *classifiedMuons.first) { /// loop over muoncopy vector
    CP::CorrectionCode result = 
    m_muonCalibrationAndSmearingTool->applyCorrection(*muon); 
    if(result != CP::CorrectionCode::Ok){
      cout << "muon pt " << muon->pt()<< " eta " << muon->eta()<<endl;
      throw std::runtime_error("Error when calibrating muons. Exiting." );
    }
  }

//   m_store->record(classifiedMuons.first,  "classifiedMuons");
//   m_store->record(classifiedMuons.second, "classifiedMuonsAux");

  std::pair<unsigned int, unsigned int> muPair;
  
  if (!m_runElectronChannel){
    muPair = SelectMuons(classifiedMuons.first, primVertex, true);
    
    xAOD::MuonContainer::iterator muon_itr = classifiedMuons.first->begin();
    xAOD::MuonContainer::iterator muon_end = classifiedMuons.first->end();
    
    for( ; muon_itr != muon_end; ++muon_itr ) {
      if ((*muon_itr)->auxdata< bool >( "signal" )){
        m_HistObjectDumper->plotMuon((*muon_itr),"signal muons");
      }
      if ((*muon_itr)->auxdata< bool >( "veto" ))
        m_HistObjectDumper->plotMuon((*muon_itr),"veto muons");
    }
    
    if (muPair.first!=1 || muPair.second!=0)
      return EL::StatusCode::SUCCESS;
    m_BitsetCutflow->FillCutflow("Muon Veto");
  }

  const xAOD::ElectronContainer* electrons(0);
  m_event->retrieve( electrons, "Electrons");
  if ( !m_event->retrieve( electrons, "Electrons" ).isSuccess() ){ 
    Error("execute()", "Failed to retrieve Electrons container. Exiting." );
    return EL::StatusCode::FAILURE;
  }
  
  std::pair<xAOD::ElectronContainer*,xAOD::ShallowAuxContainer*> 
  classifiedElectrons = xAOD::shallowCopyContainer(*electrons);
  xAOD::setOriginalObjectLink(*electrons, *classifiedElectrons.first); 
  for(const auto& electron : *classifiedElectrons.first) { 
    CP::CorrectionCode result = 
    m_eleCalibrationTool->applyCorrection(*electron); 
    if(result != CP::CorrectionCode::Ok){
      cout << "electron pt " << electron->pt()<< " eta " << electron->eta()<<endl;
      throw std::runtime_error("Error when calibrating muons. Exiting." );
    }
  }

//   m_store->record(classifiedElectrons.first,  "classifiedElectrons");
//   m_store->record(classifiedElectrons.second, "classifiedElectronsAux");
  
  std::pair<unsigned int, unsigned int> elPair;
  elPair = SelectElectrons( classifiedElectrons.first, m_runElectronChannel );

  if (!m_runElectronChannel){
    
    /// check dR for all vetoed electrons with selected muon
    xAOD::MuonContainer::iterator muon_itr = classifiedMuons.first->begin();
    xAOD::MuonContainer::iterator muon_end = classifiedMuons.first->end();
    
    xAOD::ElectronContainer::iterator elec_itr = classifiedElectrons.first->begin();
    xAOD::ElectronContainer::iterator elec_end = classifiedElectrons.first->end();
    
    int nOverlapElec = 0;
    
    for( ; muon_itr != muon_end; ++muon_itr ) {
      if ((*muon_itr)->auxdata< bool >( "signal" )){
        TLorentzVector part_muon = (*muon_itr)->p4();
        for( ; elec_itr != elec_end; ++elec_itr ) {
          if ((*elec_itr)->auxdata< bool >( "veto" ) || 
              (*elec_itr)->auxdata< bool >( "signal" )){
            TLorentzVector part_elec = (*elec_itr)->p4();
            double dR = part_muon.DeltaR(part_elec);
            if (dR<0.1) 
              nOverlapElec++;
          }
        }
        break;
      }
    }
        
    /// TODO verify this cut somehow...
    if ((elPair.second+elPair.first)!=nOverlapElec)
      return EL::StatusCode::SUCCESS;
    m_BitsetCutflow->FillCutflow("Electron Veto");
  }
  else{
    if (elPair.first!=1 || elPair.second!=0)
      return EL::StatusCode::SUCCESS;
    m_BitsetCutflow->FillCutflow("Electron Veto");
    
    muPair = SelectMuons(classifiedMuons.first, primVertex, 
                         !m_runElectronChannel);
    
    if (muPair.first!=0 || muPair.second!=0)
      return EL::StatusCode::SUCCESS;
    m_BitsetCutflow->FillCutflow("Muon Veto");
  }
  
  
  /// calibrate jets for MET
  const xAOD::JetContainer* jets(0);
  m_event->retrieve(jets, "AntiKt4EMTopoJets");
  std::pair<xAOD::JetContainer*,xAOD::ShallowAuxContainer*> metJets = 
  xAOD::shallowCopyContainer(*jets);
  xAOD::setOriginalObjectLink(*jets, *metJets.first); 
  for(const auto& jet : *metJets.first) {
    CP::CorrectionCode result = m_jetCalibrationTool->applyCorrection(*jet);
    if(result != CP::CorrectionCode::Ok){
          throw std::runtime_error("Error when calibrating jets. Exiting." );
    }
  }
//   m_event->record(metJets.first, "CalibAntiKt4EMTopoJets");
//   m_event->record(metJets.second,"CalibAntiKt4EMTopoJetsAux.");
  
  /// metPhotons
  const xAOD::PhotonContainer* photons(0);
  m_event->retrieve( photons, "Photons");
  if ( !m_event->retrieve( photons, "Photons" ).isSuccess() ){ 
    Error("execute()", "Failed to retrieve Photons container. Exiting." );
    return EL::StatusCode::FAILURE;
  }
  
  /// preselect photons for MET
  /// TODO do I need to make hardcopy here?
  /// or maybe I can just preselect food photons?
  /// WARNING implementation with hardcopy
//   xAOD::PhotonContainer* metPhotons = new xAOD::PhotonContainer();
//   xAOD::AuxContainerBase* metPhotonsAux = new xAOD::AuxContainerBase();
//   metPhotons->setStore( metPhotonsAux ); ///< Connect the two
//   
//   for(const auto& ph : *photons) {
//     if( !CutsMETMaker::accept(ph) ) continue;
// 
//     double photonPt = ph->pt() * 0.001;
//     if ( photonPt < 25.0 ) continue;
// 
//     double photonEta = ph->caloCluster()->etaBE(2);
//     if ( abs(photonEta) >= 2.37 ) continue;
//     if ( (abs(photonEta > 1.37)) && (abs(photonEta) < 1.52) ) continue;
//     
//     /// WARNING it's strange as for me. But it's way it's done in tutorial
//     xAOD::Photon* photon = new xAOD::Photon();
//     metPhotons->push_back( photon );
//     *photon= *ph; /// copies auxdata from one auxstore to the other
//   }
  
  /// WARNING implementation with only preselection
  ConstDataVector<xAOD::PhotonContainer> metPhotons(SG::VIEW_ELEMENTS); 
  for(const auto& ph : *photons) {
    if( !CutsMETMaker::accept(ph) ) continue;

    double photonPt = ph->pt() * 0.001;
    if ( photonPt < 25.0 ) continue;

    double photonEta = ph->caloCluster()->etaBE(2);
    if ( abs(photonEta) >= 2.37 ) continue;
    if ( (abs(photonEta > 1.37)) && (abs(photonEta) < 1.52) ) continue;
    metPhotons.push_back(ph); 
  }
  
  /// metTaus
  const xAOD::TauJetContainer* taus(0);
  m_event->retrieve( taus, "Taus");
  if ( !m_event->retrieve( taus, "Taus" ).isSuccess() ){ 
    Error("execute()", "Failed to retrieve Taus container. Exiting." );
    return EL::StatusCode::FAILURE;
  }
  
  /// WARNING implementation with only preselection
  ConstDataVector<xAOD::TauJetContainer> metTaus(SG::VIEW_ELEMENTS); 
  for(const auto& tau : *taus) {
    if( !CutsMETMaker::accept(tau) ) continue;
    metTaus.push_back(tau);
  }
 
  /// metElectrons
  ConstDataVector<xAOD::ElectronContainer> metElectrons(SG::VIEW_ELEMENTS);
  for (const auto& elec : *classifiedElectrons.first) {
    if (elec->auxdata< bool >( "signal" )) metElectrons.push_back(elec);
  }

  /// metMuons
  ConstDataVector<xAOD::MuonContainer> metMuons(SG::VIEW_ELEMENTS);
  for (const auto& muon : *classifiedMuons.first) {
    if (muon->auxdata< bool >( "signal" )) metMuons.push_back(muon);
  }
  
  /// Recalculate MET
  bool doJVTCut = true;
  std::string softTerm = "PVSoftTrk";
  std::string finalTerm = "FinalTrk";

  xAOD::MissingETContainer* met = new xAOD::MissingETContainer;
  xAOD::MissingETAuxContainer* met_aux = new xAOD::MissingETAuxContainer;
  met->setStore(met_aux);

  const xAOD::MissingETAssociationMap* metMap(0);
  m_event->retrieve( metMap, "METAssoc_AntiKt4EMTopo" );
  
  const xAOD::MissingETContainer* metcore(0);
  m_event->retrieve( metcore, "MET_Core_AntiKt4EMTopo" );
  
  m_metMaker = new met::METMaker("METMakerTool");
  EL_RETURN_CHECK("init m_metMaker", m_metMaker->initialize() );

  m_metMaker->rebuildMET("Muons", xAOD::Type::Muon, met, 
                           metMuons.asDataVector(), metMap);
  m_metMaker->rebuildMET("RefEle", xAOD::Type::Electron, met, 
                           metElectrons.asDataVector(), metMap);
  m_metMaker->rebuildMET("RefGamma", xAOD::Type::Photon, met, 
                           metPhotons.asDataVector(), metMap);
  m_metMaker->rebuildMET("RefTau", xAOD::Type::Tau, met, 
                           metTaus.asDataVector(), metMap);
  m_metMaker->rebuildJetMET("RefJet", softTerm, met,
                            metJets.first, metcore, metMap, 
                            doJVTCut);                                              

  m_metMaker->buildMETSum(finalTerm, met, (*met)[softTerm]->source());
  
  cout << "met = " << met << endl;
  
  /// FIXME exit from execute here for debugging purpose
  return EL::StatusCode::SUCCESS;
}

























