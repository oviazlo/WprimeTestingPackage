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

  const xAOD::EventInfo* eventInfo = 0;
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
  std::vector<std::string> triggerChains = {"HLT_mu50.*"};
  bool passTrigger = true;
  
  for(std::vector<std::string>::iterator it = triggerChains.begin(); it != 
    triggerChains.end(); ++it) {
    auto chainGroup = m_trigDecisionTool->getChainGroup(*it);
    for(auto &trig : chainGroup->getListOfTriggers()) {
      auto cg = m_trigDecisionTool->getChainGroup(trig);
      std::string thisTrig = trig;
      if (cg->isPassed()==false){
        passTrigger = false;
      }
      else{
        m_BitsetCutflow->FillCutflow(*it);
      }
    }
  }
  
  if(m_doNotApplyTriggerCuts){
    if (passTrigger==false)
      return EL::StatusCode::SUCCESS;
    m_BitsetCutflow->FillCutflow("Trigger");
  }
  
  std::vector<const xAOD::IParticle*> uniques;
  /// loop over electrons (needed for MET)
  const xAOD::ElectronContainer* electrons(0);
  m_event->retrieve( electrons, "Electrons");
  if ( !m_event->retrieve( electrons, "Electrons" ).isSuccess() ){ 
    Error("execute()", "Failed to retrieve Electrons container. Exiting." );
    return EL::StatusCode::FAILURE;
  }
  
  /*
  uniques.clear();
  ConstDataVector<xAOD::ElectronContainer> metElectrons(SG::VIEW_ELEMENTS);
  for(const auto& electron : *electrons) {
    if(CutsMETMaker::accept(electron)) {
  metElectrons.push_back(electron);
    }
  }
  */
  
  /// complex loop over muon (full copy)
  const xAOD::MuonContainer* muons = 0;
  if ( !m_event->retrieve( muons, "Muons" ).isSuccess() ){ 
    Error("execute()", "Failed to retrieve Muons container. Exiting." );
    return EL::StatusCode::FAILURE;
  }

  /// TODO cross-check this functionality
  std::pair<xAOD::MuonContainer*,xAOD::ShallowAuxContainer*> muoncopy = 
  xAOD::shallowCopyContainer(*muons);
  /// This line is very important! Overlap removal will fail if line is missing
  xAOD::setOriginalObjectLink(*muons, *muoncopy.first); 
  for(const auto& muon : *muoncopy.first) { /// loop over muoncopy vector
    if (fabs(muon->eta())<2.5) {
      /// change is valid only for muoncopy vector (I assume)
      CP::CorrectionCode result = 
      m_muonCalibrationAndSmearingTool->applyCorrection(*muon); 
      if(result != CP::CorrectionCode::Ok){
        cout << "muon pt " << muon->pt()<< " eta " << muon->eta()<<endl;
        throw std::runtime_error("Error when calibrating muons. Exiting." );
      }
    }
  }
  /// TODO what's that? whe try to create new container in m_store? 
  /// is this the place where we have ERRORs?
  m_store->record(muoncopy.first,  "MuonCopy");
  m_store->record(muoncopy.second, "MuonCopyAux");
  
  uniques.clear();
  ConstDataVector<xAOD::MuonContainer> metMuons(SG::VIEW_ELEMENTS);
  for(const auto& muon : *muoncopy.first) { 
    /// metMuons are muons after calibration and smearing...
    if(CutsMETMaker::accept(muon)) {
      metMuons.push_back(muon);
    }
  }
  
//   ConstDataVector<xAOD::MuonContainer> noMuons(SG::VIEW_ELEMENTS);
  
  /// simple loop over taus
  const xAOD::TauJetContainer* taus(0);
  m_event->retrieve( taus, "TauJets");
  if ( !m_event->retrieve( taus, "TauJets" ).isSuccess() ){ 
    Error("execute()", "Failed to retrieve Taus container. Exiting." );
    return EL::StatusCode::FAILURE;
  }
  
  uniques.clear();
  ConstDataVector<xAOD::TauJetContainer> metTaus(SG::VIEW_ELEMENTS);
  for(const auto& tau : *taus) {
    if(CutsMETMaker::accept(tau)) {
      metTaus.push_back(tau);
    }
  }
  
  /// complex loop over jets
  const xAOD::JetContainer* jets(0);
  m_event->retrieve(jets, "AntiKt4EMTopoJets");
  uniques.clear();

  std::pair<xAOD::JetContainer*,xAOD::ShallowAuxContainer*> jetcopy = 
  xAOD::shallowCopyContainer(*jets);
  xAOD::setOriginalObjectLink(*jets, *jetcopy.first); 
  for(const auto& jet : *jetcopy.first) {
    /// TODO check it, is it up-to-date way to calibrate jets
    //cout << "original jet pt " << (*jet).pt()/1000. <<endl;
    /*
    CP::CorrectionCode result = m_jetCalibrationTool->applyCorrection(*jet);
    if(result != CP::CorrectionCode::Ok){
      throw std::runtime_error("Error when calibrating jets. Exiting." );
    }
    */
    bool hasjvt = jet->isAvailable<float>("Jvt");
    float oldjvt = -999;
    if (hasjvt) oldjvt = jet->auxdata<float>("Jvt");
    //cout << "calib jet pt " << (*jet).pt()/1000. << " calib jet pt " <<endl;
    float newjvt = m_jvtTool->updateJvt(*jet);
    jet->auxdecor<float>("Jvt") = newjvt;
    //cout << "oldjvt " << oldjvt << " newjvt " << newjvt<<endl;
  }

  /// recond Jet copies to storage...
  m_store->record(jetcopy.first,  "JetCopy");
  m_store->record(jetcopy.second, "JetCopyAux");

  /// complex loop over electrons
  std::pair<xAOD::ElectronContainer*,xAOD::ShallowAuxContainer*> elecopy = 
  xAOD::shallowCopyContainer(*electrons);
  xAOD::setOriginalObjectLink(*electrons, *elecopy.first); 
  for(const auto& ele : *elecopy.first) {

    //cout << "original jet pt " << (*jet).pt()/1000. <<endl;
    CP::CorrectionCode result = m_eleCalibrationTool->applyCorrection(*ele);
    if(result != CP::CorrectionCode::Ok){
      throw std::runtime_error("Error when calibrating jets. Exiting." );
    }
  }
  
  /// recond Electron copies to storage...
  m_store->record(elecopy.first,  "EleCopy");
  m_store->record(elecopy.second, "EleCopyAux");

  /// LOOP OVER JETS

  /// Loop over all jets in the event
  /// get jet container of interest
//   const xAOD::JetContainer* jets = 0;

  /// retrieve arguments: container type, container key
  /// WARNING retieve jets again... why not use previous container?
//   if ( !m_event->retrieve( jets, "AntiKt4EMTopoJets" ).isSuccess() ){ 
//     Error("execute()","Failed to retrieve AntiKt4EMTopoJets container. "
//     "Exiting.");
//     return EL::StatusCode::FAILURE;
//   }

  /// loop over the jets in the container
  /// not used in current cutflow
  /*
  bool foundMassyJet = false;
  xAOD::JetContainer::const_iterator jet_itr = jets->begin();
  xAOD::JetContainer::const_iterator jet_end = jets->end();
  for( ; jet_itr != jet_end; ++jet_itr ) {
    
    /// check for Jet pt (not interested in Jets with pt<=20)
    if (( (*jet_itr)->pt()) * 0.001 <= 20.0) continue;
    
    /// check if jet satisfy at least VeryLooseBad jet requirement
    /// if no - discar event
    if( !m_jetCleaning->accept( **jet_itr )){
      return EL::StatusCode::SUCCESS;
    }
    
  } /// end for loop over jets
  
  m_BitsetCutflow->FillCutflow("JetCleaning");
  */
  
  /// select electrons and muon candidate using W' cut flow
  xAOD::Electron* signalEl = 0;

  static SG::AuxElement::Accessor< xAOD::IParticleLink > 
  accSetOriginLink("originalObjectLink"); /// WARNING what's this?
  ConstDataVector<xAOD::MuonContainer> signalMuon(SG::VIEW_ELEMENTS);
  /// not really used right now
  if (!m_runElectronChannel){
    /// loop over the muons in the container
    /// signal selection
    //signalMuon = SelectMuon(muons,primVertex);
    for(const auto& mu : *muoncopy.first) {
      if (passMuonSelection(mu, eventInfo, primVertex, false)) {
        //const xAOD::IParticleLink originLink( *muons, (*mu).index() );
        //accSetOriginLink(*mu) = originLink;
        signalMuon.push_back(mu);
        //cout << "pass muon with " << mu->pt()<<endl;
      }
    }

    if (signalMuon.size()==0)
      return EL::StatusCode::SUCCESS;
  }
  else{
    /// to be filled: electron selection
  }

  /// look for veto muon, add this for the moment
  //if (signalMuon.size()>1) return EL::StatusCode::SUCCESS;
  if (signalMuon.size()==0 || signalMuon.size()>1) 
    return EL::StatusCode::SUCCESS;

  /// Create a new muon container
  ConstDataVector<xAOD::MuonContainer> vetoMu(SG::VIEW_ELEMENTS); 
  vetoMu = VetoMuon(muoncopy.first,primVertex,true);
  //cout << "Run/Event " << EventNumber << " signal " << signalMuon.size() 
  //<< ", veto " << vetoMu.size()<<endl;
  if (vetoMu.size()>1) return EL::StatusCode::SUCCESS;
  if (signalMuon.size()>1) return EL::StatusCode::SUCCESS;
  m_BitsetCutflow->FillCutflow("Veto muon");
  count[12]+=1;
  ConstDataVector<xAOD::ElectronContainer> vetoEl(SG::VIEW_ELEMENTS);
  vetoEl = SelectElectron(elecopy.first,primVertex,true);
  if (vetoEl.size()!=0) return EL::StatusCode::SUCCESS;
  m_BitsetCutflow->FillCutflow("Veto electron");
  count[13]+=1;
  if (nGoodVtx>1)
    cout << "WARNING!!!! Found more then one prim.vertex: nGoodVtx = " 
    << nGoodVtx << endl;
  if (nGoodVtx==0)  return EL::StatusCode::SUCCESS;
  m_BitsetCutflow->FillCutflow("Primary vertex");
  count[14]+=1;
  
  
  
  /// simple loop over photons
  /// need to make preselection -> use deep copy
  const xAOD::PhotonContainer* photons(0);
  m_event->retrieve( photons, "Photons");
  if ( !m_event->retrieve( photons, "Photons" ).isSuccess() ){ 
    Error("execute()", "Failed to retrieve Photons container. Exiting." );
    return EL::StatusCode::FAILURE;
  }
  
  /// Create the new container and its auxiliary store.
  xAOD::PhotonContainer* metPhotons = new xAOD::PhotonContainer();
  xAOD::AuxContainerBase* metPhotonsAux = new xAOD::AuxContainerBase();
  metPhotons->setStore( metPhotonsAux ); ///< Connect the two
  
  xAOD::PhotonContainer::iterator photon_itr = (photons.first)->begin();
  xAOD::PhotonContainer::iterator photon_end = (photons.first)->end();
  
  for( ; photon_itr != photon_end; ++photon_itr ) {
    if(!CutsMETMaker::accept(**photon_itr)) continue;
    xAOD::Photon* photon = new xAOD::Photon();
        
    double photonPt = ((**photon_itr)->pt()) * 0.001;
    if (photonPt<25.0) continue;
    
    metPhotons->push_back( photon ); /// jet acquires the goodJets auxstore
    *photon= **photon_itr; /// copies auxdata from one auxstore to the other
    
  }
  
  
  
  
  /// move met to last step as it needs our selected muons as input
  /// calculate MET
  bool doJVFCut = false;
  std::string softTerm = "PVSoftTrk";
  std::string finalTerm = "FinalTrk";

  /// comment for the moment as this does not work on EXOT9
  xAOD::MissingETContainer* met = new xAOD::MissingETContainer;
  xAOD::MissingETAuxContainer* met_aux = new xAOD::MissingETAuxContainer;
  met->setStore(met_aux);

  const xAOD::MissingETAssociationMap* metMap(0);
  m_event->retrieve( metMap, "METAssoc_AntiKt4EMTopo" );
  metMap->resetObjSelectionFlags() ;

  const xAOD::MissingETContainer* metcore(0);
  m_event->retrieve( metcore, "MET_Core_AntiKt4EMTopo" );

  const xAOD::MissingETContainer* metold(0);
  m_event->retrieve( metold, "MET_Reference_AntiKt4EMTopo" );

  if ( !m_event->retrieve( metold, "MET_Reference_AntiKt4EMTopo" ).isSuccess() ){
    Error("execute()", "Failed to retrieve old met container. Exiting." );
    return EL::StatusCode::FAILURE;
  }

  bool doTracks = true;
  /// not needed any more as now I have my selected muons
  m_metMaker->rebuildMET("NewMuons", xAOD::Type::Muon, met, 
                         signalMuon.asDataVector(), metMap);
  m_metMaker->rebuildJetMET("NewRefJet", "PVSoftTrk", met, jetcopy.first, 
                            metcore, metMap,false);
  m_metMaker->buildMETSum("NewFinalTrk", met, (*met)[softTerm]->source());

  double mpx = (*met)["NewFinalTrk"]->mpx();
  double mpy = (*met)["NewFinalTrk"]->mpy();

  TLorentzVector *metVec = new TLorentzVector();
  metVec->SetPx(mpx);
  metVec->SetPy(mpy);
  metVec->SetPz(0.0);
  metVec->SetE(sqrt(mpx*mpx + mpy*mpy));
  double phi_met = metVec->Phi();
  //cout << "MET new " << mpx << ", " << mpy << ", " << sqrt(mpx*mpx + mpy*mpy) 
  //<< " phi " << phi_met<<endl;
  
  /// uncalibrated met
  /// get MET_RefFinalFix container of interest
  //cout << "get old met values"<<endl;
  const xAOD::MissingETContainer* metcontainer = 0;

  /// we need to use MET_RefFinalFix, according to:
  /// https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/Run2xAODMissingET
  if ( !m_event->retrieve(metcontainer, 
    "MET_Reference_AntiKt4EMTopo" ).isSuccess() ){ 
    Error("execute()", 
          "Failed to retrieve MET_Reference_AntiKt4EMTopo container. Exiting." );
    return EL::StatusCode::FAILURE;
  }

  /// We want to check default Final MET, wo any recalibration first
  /*
  xAOD::MissingETContainer::const_iterator met_itold = 
  metcontainer->find("FinalTrk"); 
  
  if (met_itold == metcontainer->end()) {
    Error("execute()", "No RefFinal inside MET container" );
  }
  
  double mpx2 = (*met_itold)->mpx();
  double mpy2 = (*met_itold)->mpy();
  cout << "MET old " << mpx2 << ", " << mpy2 << ", " 
  << sqrt(mpx2*mpx2 + mpy2*mpy2)<<endl;
  */
  /// easier
  double mpx2 = (*met)["NewRefJet"]->mpx();
  double mpy2 = (*met)["NewRefJet"]->mpy();
  double refjet = sqrt(pow(mpx2,2)+pow(mpy2,2));
  //cout << "new met refjet " << mpx2/1000.<< ", " << refjet<<endl;
  mpx2 = (*met)["NewMuons"]->mpx();
  mpy2 = (*met)["NewMuons"]->mpy();
  //double refmuon = 0.;
  double refmuon = sqrt(pow(mpx2,2)+pow(mpy2,2));
  
  double refsoft = sqrt(pow(mpx2,2)+pow(mpy2,2));
  mpx2 = (*met)["PVSoftTrk"]->mpx();
  mpy2 = (*met)["PVSoftTrk"]->mpy();
  double reftrk = sqrt(pow(mpx2,2)+pow(mpy2,2));
  
  mpx2 = (*metcontainer)["RefJet"]->mpx();
  mpy2 = (*metcontainer)["RefJet"]->mpy();
  double refjetold = sqrt(pow(mpx2,2)+pow(mpy2,2));
  //cout << "old met refjet " << mpx2/1000.<<endl;

  mpx2 = (*metcontainer)["Muons"]->mpx();
  mpy2 = (*metcontainer)["Muons"]->mpy();
  //double refmuon = 0.;
  double refmuonold = sqrt(pow(mpx2,2)+pow(mpy2,2));
  //cout << "old met refmuon " << mpx2/1000.<<endl;

  mpx2 = (*metcontainer)["RefGamma"]->mpx();
  mpy2 = (*metcontainer)["RefGamma"]->mpy();
  double refgammaold = sqrt(pow(mpx2,2)+pow(mpy2,2));
  //cout << "old met refgam " << mpx2/1000.<<endl;

  mpx2 = (*metcontainer)["RefEle"]->mpx();
  mpy2 = (*metcontainer)["RefEle"]->mpy();
  double refeleold = sqrt(pow(mpx2,2)+pow(mpy2,2));
  //cout << "old met refele " << mpx2/1000.<<endl;

  mpx2 = (*metcontainer)["RefTau"]->mpx();
  mpy2 = (*metcontainer)["RefTau"]->mpy();
  double reftauold = sqrt(pow(mpx2,2)+pow(mpy2,2));

  //cout << "soft term"<<endl;
  mpx2 = (*metcontainer)["SoftClus"]->mpx();
  mpy2 = (*metcontainer)["SoftClus"]->mpy();
  double refsoftold = sqrt(pow(mpx2,2)+pow(mpy2,2));
  //cout << "soft: "<<refsoft<<endl;
  //cout << "old met refsoft " << mpx2/1000.<<endl;

  mpx2 = (*metcontainer)["PVSoftTrk"]->mpx();
  mpy2 = (*metcontainer)["PVSoftTrk"]->mpy();
  double reftrkold = sqrt(pow(mpx2,2)+pow(mpy2,2));
  //cout << "rtrk: "<<reftrk<<endl;
  //cout << "old met reftrk " << mpx2/1000.<< " " << reftrkold<<endl;
 
  //cout << "METcore result " << refcore<<endl;
  mpx2 = (*metcontainer)["FinalTrk"]->mpx();
  mpy2 = (*metcontainer)["FinalTrk"]->mpy();
  double reftotal = sqrt(pow(mpx,2)+pow(mpy,2));
  double oldtotal = sqrt(pow(mpx2,2)+pow(mpy2,2));
  //cout << "new MET: " << reftotal/1000. << " old MET: " << oldtotal/1000.<<endl;
  //cout << "MET components: jet " << refjet/1000 << " muon " 
  //<< refmuon/1000 << " gamma " << refgamma/1000 << " electron " 
  //<< refele/1000 << " tau " << reftau/1000.<< " total " <<reftotal/1000.<< endl;
  //cout << "MET old " << mpx2 << ", " << mpy2 << ", " 
  //<< sqrt(mpx2*mpx2 + mpy2*mpy2)<<endl;

  //cout << "MET new " << mpx << ", " << mpy << ", " << sqrt(mpx*mpx + mpy*mpy) 
  //<<endl;
  TLorentzVector *metVecold = new TLorentzVector();
  metVecold->SetPx(mpx2);
  metVecold->SetPy(mpy2);
  metVecold->SetPz(0.0);
  metVecold->SetE(sqrt(mpx2*mpx2 + mpy2*mpy2));
  double phi_metold = metVecold->Phi();
  //cout << "new MET: "<<metVec->Pt()*0.001<< " old MET: " 
  //<< metVecold->Pt()* 0.001<<endl; 
  //if (metVec->Pt()<55000)
  //if (metVecold->Pt()<20000)
  //return EL::StatusCode::SUCCESS;
  //m_BitsetCutflow->FillCutflow("MET cut");
  count[15]+=1;
  /// calculate MT
  xAOD::MuonContainer::const_iterator muon_itr = signalMuon.begin();
  xAOD::MuonContainer::const_iterator muon_end = signalMuon.end();
  for (muon_itr=signalMuon.begin(); muon_itr != muon_end; ++muon_itr ){
    xAOD::Muon* mu = 0;
    mu = const_cast<xAOD::Muon*> (*muon_itr);
    double phi_mu = mu->phi();
    double pt_mu = mu->pt();
    //double phi_mu = signalMuon->phi();
    double Mt = sqrt( 2*mu->pt()*sqrt(mpx*mpx + mpy*mpy) * 
    (1.0 - TMath::Cos( phi_mu - phi_met )) );
    //double Mt = sqrt( 2*pt_mu*sqrt(mpx2*mpx2 + mpy2*mpy2) * 
    //(1.0 - TMath::Cos( phi_mu - phi_metold )) );
    //cout << "MT = " << Mt<<endl;
    h_Mt_muonPtCut->Fill(Mt * 0.001);
    //if (Mt<60000.)
    //return EL::StatusCode::SUCCESS;
    m_BitsetCutflow->FillCutflow("MT cut");
    count[16]+=1;

    EventNumber = eventInfo->eventNumber();
    RunNumber = eventInfo->runNumber();
    LumiBlock = eventInfo->lumiBlock();
    MuonPt = mu->pt();
    MuonEta = mu->eta();
    MuonPhi = mu->phi();
    MissingEt = metVec->Pt();
    MissingEx = metVec->Px();
    MissingEy = metVec->Py();
    OldMissingEt = metVecold->Pt();
    TransverseMass = Mt;
//     tree->Fill();
    
    if (m_useHistObjectDumper) m_HistObjectDumper->plotMuon(mu,"allCuts");
    
  }

  /// ************************************************
  /// Fill other event-based distributions (e.g. MET)
  /// ************************************************
  h_MET_RefJet->Fill(refjet* 0.001);
  //h_MET_RefCore->Fill(refcore* 0.001);
  //h_MET_RefTau->Fill(reftau* 0.001);
  //h_MET_RefGamma->Fill(refgamma* 0.001);
  //h_MET_RefElectron->Fill(refele* 0.001);
  h_MET_RefMuon->Fill(refmuon* 0.001);
  h_MET_RefSoft->Fill(refsoft* 0.001);
  h_MET_RefTrack->Fill(reftrk* 0.001);

  h_MET_OldJet->Fill(refjetold* 0.001);
  h_MET_OldTau->Fill(reftauold* 0.001);
  h_MET_OldGamma->Fill(refgammaold* 0.001);
  h_MET_OldElectron->Fill(refeleold* 0.001);
  h_MET_OldMuon->Fill(refmuonold* 0.001);
  h_MET_OldSoft->Fill(refsoftold* 0.001);
  h_MET_OldTrack->Fill(reftrkold* 0.001);
  
//   h_MET_RefFinalNew->Fill(sqrt(mpx*mpx + mpy*mpy)* 0.001);
//   h_MET_RefFinalFix->Fill(sqrt(mpx2*mpx2 + mpy2*mpy2)* 0.001);
  h_MET_RefFinalFix_test->Fill(reftotal * 0.001);

  m_BitsetCutflow->FillCutflow("End");
  /*
  /// get MET_RefFinalFix container of interest
  const xAOD::MissingETContainer* metcontainer = 0;
  
  /// we need to use MET_Core_AntiKt4EMTopo, according to:
  /// https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/Run2xAODMissingET
  /// FIXME use MET_Reference_AntiKt4EMTopo instead of MET_Core_AntiKt4EMTopo
  if ( !m_event->retrieve( metcontainer, 
    "MET_Reference_AntiKt4EMTopo" ).isSuccess())
  { /// retrieve arguments: container$
    Error("execute()","Failed to retrieve MET_RefFinalFix container. Exiting.");
    return EL::StatusCode::FAILURE;
  }

  
  
  
  /// We want to check default Final MET, wo any recalibration first
  xAOD::MissingETContainer::const_iterator met_it = metcontainer->find("FinalTrk"); 

  if (met_it == metcontainer->end()) {
    Error("execute()", "No RefFinal inside MET container" );
  }

  double mpx = (*met_it)->mpx();  
  double mpy = (*met_it)->mpy();
  TLorentzVector *metVec = new TLorentzVector();
  metVec->SetPx(mpx);
  metVec->SetPy(mpy);
  metVec->SetPz(0.0);
  metVec->SetE(sqrt(mpx*mpx + mpy*mpy));
  
  double phi_met = metVec->Phi();
  
  //~ // create a shallow copy of the muons container
  //~ std::pair< xAOD::MuonContainer*, xAOD::ShallowAuxContainer* > 
  muons_shallowCopy = xAOD::shallowCopyContainer( *muons );
//~ 
  //~ // iterate over our shallow copy
  //~ xAOD::MuonContainer::iterator muonSC_itr = 
  //~ (muons_shallowCopy.first)->begin();
  //~ xAOD::MuonContainer::iterator muonSC_end = 
  //~ (muons_shallowCopy.first)->end();
//~ 
  //~ for( ; muonSC_itr != muonSC_end; ++muonSC_itr ) {
   //~ if(m_muonCalibrationAndSmearingTool->applyCorrection(**muonSC_itr) == 
   CP::CorrectionCode::Error){ // apply correction and check return code
       //~ // Can have CorrectionCode values of Ok, OutOfValidityRange, or Error.
       //~ // Here only checking for Error.
       //~ // If OutOfValidityRange is returned no modification is made and the 
       //~ // original muon values are taken.
       //~ Error("execute()", 
       //~ "MuonCalibrationAndSmearingTool returns Error CorrectionCode");
   //~ }
   //~ //if(!m_muonSelection->accept(**muonSC_itr)) continue;
   //~ Info("execute()", "  corrected muon pt = %.2f GeV", ((*muonSC_itr)->pt() 
   * 0.001));  
  //~ } // end for loop over shallow copied muons
  //~ delete muons_shallowCopy.first;
  //~ delete muons_shallowCopy.second;
  
  /// get muon container of interest
  const xAOD::MuonContainer* muons = 0;
  if ( !m_event->retrieve( muons, "Muons" ).isSuccess() ){ 
    Error("execute()", "Failed to retrieve Muons container. Exiting." );
    return EL::StatusCode::FAILURE;
  }  
  
  const xAOD::ElectronContainer* electrons = 0;
  if ( !m_event->retrieve( electrons, "Electrons" ).isSuccess() ){ 
    Error("execute()", "Failed to retrieve Electrons container. Exiting." );
    return EL::StatusCode::FAILURE;
  }  
  
  xAOD::Muon* signalMuon = 0;
  xAOD::Muon* signalEl = 0;
  
  if (!m_runElectronChannel){
    /// loop over the muons in the container
    /// signal selection
    signalMuon = SelectMuon(muons,primVertex);
    if (signalMuon==0)
      return EL::StatusCode::SUCCESS;
  }
  else{
    
  }



  /// look for veto muon
  xAOD::Muon* vetoMu = SelectMuon(muons,primVertex,true);
  if (vetoMu!=0)
    return EL::StatusCode::SUCCESS;
  m_BitsetCutflow->FillCutflow("Veto muon");

  xAOD::Electron* vetoEl = SelectElectron(electrons,primVertex,true);
  if (vetoEl!=0)
    return EL::StatusCode::SUCCESS;
  m_BitsetCutflow->FillCutflow("Veto electron");
  
  
  double phi_mu;
  double Mt;

  if (m_runElectronChannel){
    phi_mu = signalEl->phi();
    Mt = sqrt( 2*signalEl->pt()*sqrt(mpx*mpx + mpy*mpy) * 
    (1.0 - TMath::Cos( phi_mu - phi_met )) );
  }
  else{
    phi_mu = signalMuon->phi();
    Mt = sqrt( 2*signalMuon->pt()*sqrt(mpx*mpx + mpy*mpy) * 
    (1.0 - TMath::Cos( phi_mu - phi_met )) );
  }
  
  
  
  h_Mt_muonPtCut->Fill(Mt * 0.001);
      
//   if (m_useHistObjectDumper) m_HistObjectDumper->plotMuon(mu,"allCuts");
  
  
  /// ************************************************
  /// Fill other event-based distributions (e.g. MET)
  /// ************************************************
  h_MET_RefFinalFix->Fill(sqrt(mpx*mpx + mpy*mpy)* 0.001);
  h_MET_RefFinalFix_test->Fill(metVec->Pt() * 0.001);

  m_BitsetCutflow->FillCutflow("End");
 */
  return EL::StatusCode::SUCCESS;
}

