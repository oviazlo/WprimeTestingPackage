/* Copyright 2016 Oleksandr Viazlo */

///*****************************************************************************
///
/// Implementation of execute() function.
///
///*****************************************************************************

#include <MyAnalysis/RecoAnalysis.h>

/// this is needed to distribute the algorithm to the workers
ClassImp(RecoAnalysis)

EL::StatusCode RecoAnalysis :: execute ()
{
  
  const char* APP_NAME = "RecoAnalysis";

  /// push cutflow bitset to cutflow hist
  m_BitsetCutflow->PushBitSet();
  
  if( (m_eventCounter % 1000) ==0 ) Info("execute()", "Event number = %i", 
    m_eventCounter );
  m_eventCounter++;

  EL_RETURN_CHECK("execute()", m_event->retrieve( m_eventInfo, "EventInfo") );
  EventNumber = m_eventInfo->eventNumber();
  
  /// WARNING cutflow
//   static const int arr[] = {646, 1951, 158, 1279, 720, 1884, 1019, 241, 1238, 43};
//   vector<int> runNumbersToPlay (arr, arr + sizeof(arr) / sizeof(arr[0]) );
//   vector<int>::iterator it;
//   it = std::find(runNumbersToPlay.begin(), runNumbersToPlay.end(), EventNumber);
//   if (it==runNumbersToPlay.end())
//     return EL::StatusCode::SUCCESS; 
  
  /// use PileUp Reweighting Tool
  m_pileupReweightingTool->apply( *m_eventInfo );
  
  /// Muon Truth matching. Check do we have muon from W' decay
  /// WARNING do we need it for all MCs?
  if(m_isMC && m_doWprimeTruthMatching){
    bool foundMuonFromWprimeDecay = false;
    
    /// Create truth vertice container
    const xAOD::TruthVertexContainer* truthVertices = 0;
    EL_RETURN_CHECK("retrieve TruthVertices", 
                    m_event->retrieve( truthVertices, "TruthVertices" ));
    
    /// Start iterating over truth container
    xAOD::TruthVertexContainer::const_iterator truthV_itr; 
    for (truthV_itr = truthVertices->begin(); 
         truthV_itr != truthVertices->end(); ++truthV_itr )
    {
      if (foundMuonFromWprimeDecay)
        break;
      for (unsigned int iIn=0; iIn < (*truthV_itr)->nIncomingParticles(); 
           iIn++)
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
  
  if(m_inclusiveWorZ){
    bool massAbove120Gev = false;
    /// Create truth vertice container
    const xAOD::TruthVertexContainer* truthVertices = 0;
    EL_RETURN_CHECK("retrieve TruthVertices", 
                    m_event->retrieve( truthVertices, "TruthVertices" ));
    
    /// Start iterating over truth container
    xAOD::TruthVertexContainer::const_iterator truthV_itr; 
    for (truthV_itr = truthVertices->begin(); 
         truthV_itr != truthVertices->end(); ++truthV_itr )
    {
      if (massAbove120Gev)
        break;
      for (unsigned int iIn=0; iIn < (*truthV_itr)->nIncomingParticles(); 
           iIn++)
      {
        /// 34 - Wprime
        if (TMath::Abs((*truthV_itr)->incomingParticle(iIn)->pdgId()) == 
          m_pdgIdOfMother) { 
          
          double m_trueWmass = 
          TMath::Sqrt(TMath::Power((*truthV_itr)->incomingParticle(iIn)->e(),2)
              - TMath::Power((*truthV_itr)->incomingParticle(iIn)->pz(),2)
              - TMath::Power((*truthV_itr)->incomingParticle(iIn)->py(),2)
              - TMath::Power((*truthV_itr)->incomingParticle(iIn)->px(),2)
                   )*0.001; 

          h_mgen_all->Fill(m_trueWmass);
          if (m_trueWmass > 120.0){
            massAbove120Gev = true;
            h_mgen->Fill(m_trueWmass);
          }
        }
        if (massAbove120Gev)
          break;
      }
    }
    
    if (massAbove120Gev) return EL::StatusCode::SUCCESS;
    m_BitsetCutflow->FillCutflow("massBelow120Gev");
    
  }
  
  /// pass GRL?
  if(!m_isMC){
    if(!m_grl->passRunLB(*m_eventInfo)){
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
  if(!m_isMC){
    if(  (m_eventInfo->errorState(xAOD::EventInfo::LAr)==xAOD::EventInfo::Error )
      || (m_eventInfo->errorState(xAOD::EventInfo::Tile)==xAOD::EventInfo::Error )
      || (m_eventInfo->errorState(xAOD::EventInfo::SCT)==xAOD::EventInfo::Error )
      || (m_eventInfo->isEventFlagBitSet(xAOD::EventInfo::Core, 18) )  )
    {
      return EL::StatusCode::SUCCESS; 
    }
  }
  m_BitsetCutflow->FillCutflow("EventCleaning");
  
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
  
    /// Primary vertex
  const xAOD::VertexContainer* vertices = 0;
  EL_RETURN_CHECK("retrieve PrimaryVertices", 
                  m_event->retrieve( vertices, "PrimaryVertices" ));
  
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

  ///***************************************************************************
  /// Apply corrections for muons and electrons
  const xAOD::MuonContainer* muons = 0;
  EL_RETURN_CHECK("retrieve Muons",
                  m_event->retrieve( muons, "Muons" ));
  
  std::pair<xAOD::MuonContainer*,xAOD::ShallowAuxContainer*> classifiedMuons = 
  xAOD::shallowCopyContainer(*muons);
  xAOD::setOriginalObjectLink(*muons, *classifiedMuons.first); 
  
  for(const auto& muon : *classifiedMuons.first) { /// loop over muoncopy vector
    if (abs(muon->eta())>2.5) continue; /// WARNING FIXME temprorary hack
    CP::CorrectionCode result = 
    m_muonCalibrationAndSmearingTool->applyCorrection(*muon); 
    if(result != CP::CorrectionCode::Ok){
      cout << "muon pt " << muon->pt()<< " eta " << muon->eta()<<endl;
      throw std::runtime_error("Error when calibrating muons. Exiting." );
    }
  }

  const xAOD::ElectronContainer* electrons(0);
  EL_RETURN_CHECK("retrieve Electrons",
                  m_event->retrieve( electrons, "Electrons" ));
  
  std::pair<xAOD::ElectronContainer*,xAOD::ShallowAuxContainer*> 
  classifiedElectrons = xAOD::shallowCopyContainer(*electrons);
  xAOD::setOriginalObjectLink(*electrons, *classifiedElectrons.first);

  for(const auto& electron : *classifiedElectrons.first) { 
    CP::CorrectionCode result = 
    m_eleCalibrationTool->applyCorrection(*electron); 
    if(result != CP::CorrectionCode::Ok){
      cout << "electron pt " << electron->pt()<< " eta " << electron->eta()<<endl;
      throw std::runtime_error("Error when calibrating electrons. Exiting." );
    }
  }
  ///***************************************************************************
  
  std::map<string, unsigned int> muCounter;
  muCounter = SelectMuons(classifiedMuons.first, primVertex, !m_runElectronChannel);
  
  std::map<string, unsigned int> elCounter;
  elCounter = SelectElectrons( classifiedElectrons.first, m_runElectronChannel );

  ConstDataVector<xAOD::TauJetContainer> metTaus(SG::VIEW_ELEMENTS);
  ConstDataVector<xAOD::PhotonContainer> metPhotons(SG::VIEW_ELEMENTS);
  std::pair<xAOD::JetContainer*,xAOD::ShallowAuxContainer*> metJets;
  
  if (!m_runElectronChannel){ /// muon channel
    
    xAOD::MuonContainer::iterator muon_itr = classifiedMuons.first->begin();
    xAOD::MuonContainer::iterator muon_end = classifiedMuons.first->end();
    
    for( ; muon_itr != muon_end; ++muon_itr ) {
      if ((*muon_itr)->auxdata< bool >( "signal" )){
        m_HistObjectDumper->plotMuon((*muon_itr),"signal_muons");
      }
      if ((*muon_itr)->auxdata< bool >( "veto" ))
        m_HistObjectDumper->plotMuon((*muon_itr),"veto_muons");
    }
    
    ///*************************************************************************
    /// fake background inplementation for data only
    if (!m_isMC){
     
      /// calibrate jets, photons and taus. Need to be done here. They will be used also
      /// for signal leptons which pass isolation!
      
      ///***************************************************************************
      /// calibrate jets for MET
      const xAOD::JetContainer* jets(0);
      EL_RETURN_CHECK("retrieve AntiKt4EMTopoJets",
                      m_event->retrieve(jets, "AntiKt4EMTopoJets"));
      
      metJets = xAOD::shallowCopyContainer(*jets);
      xAOD::setOriginalObjectLink(*jets, *metJets.first); 
      for(const auto& jet : *metJets.first) {
        CP::CorrectionCode result = m_jetCalibrationTool->applyCorrection(*jet);
        if(result != CP::CorrectionCode::Ok){
              throw std::runtime_error("Error when calibrating jets. Exiting." );
        }
        /// Correcting the xAOD JVT value
        /// source: https://twiki.cern.ch/twiki/bin/view/AtlasProtected/JetVertexTaggerTool
        bool hasjvt = jet->isAvailable<float>("Jvt");
        float oldjvt = -999;
        if (hasjvt) oldjvt = jet->auxdata<float>("Jvt");
        float newjvt = m_jvtTool->updateJvt(*jet);
        jet->auxdecor<float>("Jvt") = newjvt;
      }
      isMetJetsFilled = true;
      ///***************************************************************************
      
      ///***************************************************************************
      /// metPhotons
      const xAOD::PhotonContainer* photons(0);
      EL_RETURN_CHECK("retrieve Photons",
                      m_event->retrieve( photons, "Photons" ));

      for(const auto& ph : *photons) {
        if( !CutsMETMaker::accept(ph) ) continue;

        double photonPt = ph->pt() * 0.001;
        if ( photonPt < 25.0 ) continue;

        double photonEta = ph->caloCluster()->etaBE(2);
        if ( abs(photonEta) >= 2.37 ) continue;
        if ( (abs(photonEta > 1.37)) && (abs(photonEta) < 1.52) ) continue;
        metPhotons.push_back(ph); 
      }
      isMetPhotonsFilled = true;
      ///***************************************************************************
      
      ///***************************************************************************
      /// metTaus
      const xAOD::TauJetContainer* taus(0);
      EL_RETURN_CHECK("retrieve TauRecContainer",
                      m_event->retrieve( taus, "TauJets" ));
      
      /// WARNING implementation with only preselection
      for(const auto& tau : *taus) {
        if( !CutsMETMaker::accept(tau) ) continue;
        metTaus.push_back(tau);
      }
      isMetTausFilled = true;
      ///***************************************************************************
    
      ///***************************************************************************
      /// metElectrons_failIso
      ConstDataVector<xAOD::ElectronContainer> metElectrons_failIso(SG::VIEW_ELEMENTS);
      for (const auto& elec : *classifiedElectrons.first) {
        if (elec->auxdata< bool >( "signal_failIso" )){
          if (elec->auxdata< bool >( "overlap_failIso" )==false)
            metElectrons_failIso.push_back(elec);
        }
      }
      ///***************************************************************************

      ///***************************************************************************
      /// metMuons_failIso
      ConstDataVector<xAOD::MuonContainer> metMuons_failIso(SG::VIEW_ELEMENTS);
      for (const auto& muon : *classifiedMuons.first) {
        if (muon->auxdata< bool >( "signal_failIso" )) metMuons_failIso.push_back(muon);
      }
      ///***************************************************************************                
      
      for (int i=0; i<1; i++){ /// hack make this loop just to enable "break" functionality
        
        /// failIso muon veto
        if (muCounter["nSignalLeptons_failIso"]!=1 || muCounter["nVetoLeptons_failIso"]!=0)
          break;
        
        /// check dR for all vetoed electrons with selected muon
        xAOD::MuonContainer::iterator muon_itr = classifiedMuons.first->begin();
        xAOD::MuonContainer::iterator muon_end = classifiedMuons.first->end();
        
        xAOD::ElectronContainer::iterator elec_itr = classifiedElectrons.first->begin();
        xAOD::ElectronContainer::iterator elec_end = classifiedElectrons.first->end();

        int nOverlapElec_failIso = 0;
        for( ; muon_itr != muon_end; ++muon_itr ) {
          if ((*muon_itr)->auxdata< bool >( "signal_failIso" )){
            TLorentzVector part_muon = (*muon_itr)->p4();
            for( ; elec_itr != elec_end; ++elec_itr ) {
              if ((*elec_itr)->auxdata< bool >( "veto_failIso" ) || 
                  (*elec_itr)->auxdata< bool >( "signal_failIso" )){
                (*elec_itr)->auxdata< bool >( "overlap__failIso" ) = false;
                TLorentzVector part_elec = (*elec_itr)->p4();
                double dR = part_muon.DeltaR(part_elec);
                if (dR<0.1){
                  nOverlapElec_failIso++;
                  (*elec_itr)->auxdata< bool >( "overlap_failIso" ) = true;
                }
              }
            }
            break;
          }
        }
        
        /// failIso electron veto
        if ((elCounter["nVetoLeptons_failIso"] + elCounter["nSignalLeptons_failIso"])
          != nOverlapElec_failIso)
          break;
        
        ///***************************************************************************
        /// Recalculate and get MET
        bool doJVTCut = true;
        std::string softTerm = "PVSoftTrk";
        std::string finalTerm = "FinalTrk";

        xAOD::MissingETContainer* met_failIso = new xAOD::MissingETContainer;
        xAOD::MissingETAuxContainer* met_aux_failIso = new xAOD::MissingETAuxContainer;
        met_failIso->setStore(met_aux_failIso);

        const xAOD::MissingETAssociationMap* metMap_failIso(0);
        EL_RETURN_CHECK("retrieve METAssoc_AntiKt4EMTopo",
                        m_event->retrieve( metMap_failIso, "METAssoc_AntiKt4EMTopo" ));
        
        const xAOD::MissingETContainer* metcore_failIso(0);
        EL_RETURN_CHECK("retrieve MET_Core_AntiKt4EMTopo",
                        m_event->retrieve( metcore_failIso, "MET_Core_AntiKt4EMTopo" ));
        
        m_metMaker->rebuildMET("RefEle", xAOD::Type::Electron, met_failIso, 
                                metElectrons_failIso.asDataVector(), metMap_failIso);
        m_metMaker->rebuildMET("RefGamma", xAOD::Type::Photon, met_failIso, 
                                metPhotons.asDataVector(), metMap_failIso);
        m_metMaker->rebuildMET("RefTau", xAOD::Type::Tau, met_failIso, 
                                metTaus.asDataVector(), metMap_failIso);
        m_metMaker->rebuildMET("Muons", xAOD::Type::Muon, met_failIso, 
                                metMuons_failIso.asDataVector(), metMap_failIso);
        m_metMaker->rebuildJetMET("RefJet", softTerm, met_failIso,
                                  metJets.first, metcore_failIso, metMap_failIso, 
                                  doJVTCut);        

        m_metMaker->buildMETSum(finalTerm, met_failIso, (*met_failIso)[softTerm]->source());
        
        xAOD::MissingET * finalTrkMet_failIso = (*met_failIso)["FinalTrk"];
        if ( finalTrkMet_failIso == 0) {
            throw std::runtime_error("Pointer finalTrkMet_failIso is NULL. Exiting." );
        }
        
        double missingEt    = finalTrkMet_failIso->met()/1000.;
        double m_met_sumet  = finalTrkMet_failIso->sumet()/1000.;
        double missingEtPhi = finalTrkMet_failIso->phi();
        
        /// fill in only for data, so no weights
        m_HistObjectDumper->plotMuon(metMuons_failIso[0],"final_noMET_mT_cuts_failIso",1.0);
        m_HistObjectDumper->plotMtAndMet
        (metMuons_failIso[0],finalTrkMet_failIso,"final_noMET_mT_cuts_failIso",1.0);
       
        plotPtBinnedHists(metMuons_failIso[0],finalTrkMet_failIso,
                          "final_noMET_mT_cuts_failIso",1.0);
        
        ///***************************************************************************
        /// calculate QCD nkg. hist
        if (!m_isMC){
          double QCD_weight = getQCDWeight(metMuons_failIso[0]->pt() * 0.001,false);
          m_HistObjectDumper->plotMuon(metMuons_failIso[0],"final_noMET_mT_cuts_QCD",QCD_weight);
          m_HistObjectDumper->plotMtAndMet
          (metMuons_failIso[0],finalTrkMet_failIso,"final_noMET_mT_cuts_QCD",QCD_weight);
        }
        ///***************************************************************************
        
        ///*********************************************************************
        /// calculate mT
        double leptonEt;
        double leptonPhi;
        
        leptonEt = metMuons_failIso[0]->pt() * 0.001;
        leptonPhi = metMuons_failIso[0]->phi();
        
        float deltaPhi = TMath::Abs(leptonPhi-missingEtPhi);
        if(deltaPhi > TMath::Pi())
              deltaPhi = TMath::TwoPi() - deltaPhi;
        
        double mT = sqrt( 2 * missingEt * leptonEt * 
        (1 - TMath::Cos( deltaPhi ) ) );
        ///*********************************************************************
        
        if (missingEt>metCut && mT>mtCut){
          m_HistObjectDumper->plotMuon(metMuons_failIso[0],"final_failIso",1.0);
          m_HistObjectDumper->plotMtAndMet
          (metMuons_failIso[0],finalTrkMet_failIso,"final_failIso",1.0);
        
          plotPtBinnedHists(metMuons_failIso[0],finalTrkMet_failIso,
                            "final_failIso",1.0);
          
          ///***************************************************************************
          /// calculate QCD nkg. hist
          if (!m_isMC){
            double QCD_weight = getQCDWeight(metMuons_failIso[0]->pt() * 0.001,false);
            m_HistObjectDumper->plotMuon(metMuons_failIso[0],"final_QCD",QCD_weight);
            m_HistObjectDumper->plotMtAndMet
            (metMuons_failIso[0],finalTrkMet_failIso,"final_QCD",QCD_weight);
          }
          ///***************************************************************************

        }
        
        /// see presentation by Saminder Dhaliwal in Lepton+X meeting on March 22.
        /// https://indico.cern.ch/event/512374/
        double missingVetoRegion = false;
        double absMuPhi = abs(metMuons_failIso[0]->phi());
        double absMuEta = abs(metMuons_failIso[0]->eta());
        if (absMuEta>1.05 && absMuEta<1.3){
          if ((absMuEta>0.21 && absMuEta<0.57) || (absMuEta>1.00 && absMuEta<1.33) || 
              (absMuEta>1.78 && absMuEta<2.14) || (absMuEta>2.57 && absMuEta<2.93))
          {
            missingVetoRegion = true;
          }
        }
        
        if (missingVetoRegion==false){
          m_HistObjectDumper->plotMuon(metMuons_failIso[0],"final_noMET_mT_cuts_updatedHighPtVeto_failIso",
                                       1);
          m_HistObjectDumper->plotMtAndMet(metMuons_failIso[0],finalTrkMet_failIso,
                                           "final_noMET_mT_cuts_updatedHighPtVeto_failIso",
                                           1);
          
          plotPtBinnedHists(metMuons_failIso[0],finalTrkMet_failIso,
                            "final_noMET_mT_cuts_updatedHighPtVeto_failIso",1);
        }
        
      }
      
    }
    ///*************************************************************************
    
    ///*************************************************************************
    /// Muon and electron veto
    if (muCounter["nSignalLeptons"]!=1 || muCounter["nVetoLeptons"]!=0)
      return EL::StatusCode::SUCCESS;
    m_BitsetCutflow->FillCutflow("Muon Veto");
    
    /// check dR for all vetoed electrons with selected muon
    muon_itr = classifiedMuons.first->begin();
    muon_end = classifiedMuons.first->end();
    
    xAOD::ElectronContainer::iterator elec_itr = classifiedElectrons.first->begin();
    xAOD::ElectronContainer::iterator elec_end = classifiedElectrons.first->end();
    
    int nOverlapElec = 0;
    int nOverlapElec_failIso = 0;
    
    for( ; muon_itr != muon_end; ++muon_itr ) {
      if ((*muon_itr)->auxdata< bool >( "signal" )){
        TLorentzVector part_muon = (*muon_itr)->p4();
        for( ; elec_itr != elec_end; ++elec_itr ) {
          if ((*elec_itr)->auxdata< bool >( "veto" ) || 
              (*elec_itr)->auxdata< bool >( "signal" )){
            (*elec_itr)->auxdata< bool >( "overlap" ) = false;
            TLorentzVector part_elec = (*elec_itr)->p4();
            double dR = part_muon.DeltaR(part_elec);
            if (dR<0.1){
              nOverlapElec++;
              (*elec_itr)->auxdata< bool >( "overlap" ) = true;
            }
          }
        }
        break;
      }
    }
    
    /// TODO verify this cut somehow...
    if ((elCounter["nVetoLeptons"] + elCounter["nSignalLeptons"])!=nOverlapElec)
      return EL::StatusCode::SUCCESS;
    m_BitsetCutflow->FillCutflow("Electron Veto");

//     cout << "[CUTFLOW_DEBUG]\t" << m_eventInfo->runNumber() << "\t" << EventNumber <<  endl;
  } /// FIXME no implementation of fake background for electron channel!!!
    /// FIXME no implementation of overlap check in electron channel!!!
  else{ /// m_runElectronChannel == true
    if (elCounter["nSignalLeptons"]!=1 || elCounter["nVetoLeptons"]!=0)
      return EL::StatusCode::SUCCESS;
    m_BitsetCutflow->FillCutflow("Electron Veto");
    
    if (muCounter["nSignalLeptons"]!=0 || muCounter["nVetoLeptons"]!=0)
      return EL::StatusCode::SUCCESS;
    m_BitsetCutflow->FillCutflow("Muon Veto");
  }
  ///*************************************************************************
  
  ///***************************************************************************
  /// calibrate jets for MET
  if (!isMetJetsFilled){
    const xAOD::JetContainer* jets(0);
    EL_RETURN_CHECK("retrieve AntiKt4EMTopoJets",
                    m_event->retrieve(jets, "AntiKt4EMTopoJets"));
    
    metJets = xAOD::shallowCopyContainer(*jets);
    xAOD::setOriginalObjectLink(*jets, *metJets.first); 
    for(const auto& jet : *metJets.first) {
      CP::CorrectionCode result = m_jetCalibrationTool->applyCorrection(*jet);
      if(result != CP::CorrectionCode::Ok){
            throw std::runtime_error("Error when calibrating jets. Exiting." );
      }
      /// Correcting the xAOD JVT value
      /// source: https://twiki.cern.ch/twiki/bin/view/AtlasProtected/JetVertexTaggerTool
      bool hasjvt = jet->isAvailable<float>("Jvt");
      float oldjvt = -999;
      if (hasjvt) oldjvt = jet->auxdata<float>("Jvt");
      float newjvt = m_jvtTool->updateJvt(*jet);
      jet->auxdecor<float>("Jvt") = newjvt;

    }
  }
  ///***************************************************************************
  
  ///***************************************************************************
  /// metPhotons
  if (!isMetPhotonsFilled){
    const xAOD::PhotonContainer* photons(0);
    EL_RETURN_CHECK("retrieve Photons",
                    m_event->retrieve( photons, "Photons" ));

    for(const auto& ph : *photons) {
      if( !CutsMETMaker::accept(ph) ) continue;

      double photonPt = ph->pt() * 0.001;
      if ( photonPt < 25.0 ) continue;

      double photonEta = ph->caloCluster()->etaBE(2);
      if ( abs(photonEta) >= 2.37 ) continue;
      if ( (abs(photonEta > 1.37)) && (abs(photonEta) < 1.52) ) continue;
      metPhotons.push_back(ph); 
    }
  }
  ///***************************************************************************
  
  ///***************************************************************************
  /// metTaus
  if (!isMetTausFilled){
    const xAOD::TauJetContainer* taus(0);
    EL_RETURN_CHECK("retrieve TauRecContainer",
                    m_event->retrieve( taus, "TauJets" ));
    
    /// WARNING implementation with only preselection
    for(const auto& tau : *taus) {
      if( !CutsMETMaker::accept(tau) ) continue;
      metTaus.push_back(tau);
    }
  }
  ///***************************************************************************
 
  ///***************************************************************************
  /// metElectrons
  ConstDataVector<xAOD::ElectronContainer> metElectrons(SG::VIEW_ELEMENTS);
  for (const auto& elec : *classifiedElectrons.first) {
    if (elec->auxdata< bool >( "signal" )){
      if (elec->auxdata< bool >( "overlap" )==false)
        metElectrons.push_back(elec);
    }
  }
  ///***************************************************************************

  ///***************************************************************************
  /// metMuons
  ConstDataVector<xAOD::MuonContainer> metMuons(SG::VIEW_ELEMENTS);
  for (const auto& muon : *classifiedMuons.first) {
    /// FIXME if run electron channel we also need to exclude muons which overlaps with electrons
    /// from metMuons (which are usid in MET rebuilding)
    if (muon->auxdata< bool >( "signal" )) metMuons.push_back(muon);
  }
  ///***************************************************************************
  
  ///***************************************************************************
  /// Recalculate and get MET
  bool doJVTCut = true;
  std::string softTerm = "PVSoftTrk";
  std::string finalTerm = "FinalTrk";

  xAOD::MissingETContainer* met = new xAOD::MissingETContainer;
  xAOD::MissingETAuxContainer* met_aux = new xAOD::MissingETAuxContainer;
  met->setStore(met_aux);

  const xAOD::MissingETAssociationMap* metMap(0);
  EL_RETURN_CHECK("retrieve METAssoc_AntiKt4EMTopo",
                  m_event->retrieve( metMap, "METAssoc_AntiKt4EMTopo" ));
  
  const xAOD::MissingETContainer* metcore(0);
  EL_RETURN_CHECK("retrieve MET_Core_AntiKt4EMTopo",
                  m_event->retrieve( metcore, "MET_Core_AntiKt4EMTopo" ));
  
  m_metMaker->rebuildMET("RefEle", xAOD::Type::Electron, met, 
                           metElectrons.asDataVector(), metMap);
  m_metMaker->rebuildMET("RefGamma", xAOD::Type::Photon, met, 
                           metPhotons.asDataVector(), metMap);
  m_metMaker->rebuildMET("RefTau", xAOD::Type::Tau, met, 
                           metTaus.asDataVector(), metMap);
  m_metMaker->rebuildMET("Muons", xAOD::Type::Muon, met, 
                           metMuons.asDataVector(), metMap);
  m_metMaker->rebuildJetMET("RefJet", softTerm, met,
                            metJets.first, metcore, metMap, 
                            doJVTCut);        

  m_metMaker->buildMETSum(finalTerm, met, (*met)[softTerm]->source());
  
  xAOD::MissingET * finalTrkMet = (*met)["FinalTrk"];
  if ( finalTrkMet == 0) {
      throw std::runtime_error("Pointer finalTrkMet is NULL. Exiting." );
  }
  
  double missingEt    = finalTrkMet->met()/1000.;
  double m_met_sumet  = finalTrkMet->sumet()/1000.;
  double missingEtPhi = finalTrkMet->phi();
  ///***************************************************************************
  
  
  /// get MC weights
  double SFWeight = 1.0;
  double trgSF = -999.0;
  float recoEffSF = -999.0;
  float isoSF = -999.0;
  double puWeight = -999.0;
  
  if (m_isMC){
	puWeight = m_pileupReweightingTool->getCombinedWeight( *m_eventInfo );

    m_LPXKfactorTool->applySystematicVariation(CP::SystematicSet()); //nominal
    m_LPXKfactorTool->execute();
    weightkFactor = m_eventInfo->auxdecor<double>("KfactorWeight");
    weighfilterEfficiency = m_LPXKfactorTool->getMCFilterEfficiency();
    sampleLumi = m_DataSetInfo->getnEventsPerSample(m_DSID)/m_LPXKfactorTool->getMCCrossSection()/1000000; /// nb --> fb

    unsigned int randomRunNumber = m_pileupReweightingTool->getRandomRunNumber( *m_eventInfo, false );
    m_trig_sf->setRunNumber(randomRunNumber);

    xAOD::MuonContainer *SelectedMuon = new xAOD::MuonContainer;
    xAOD::MuonAuxContainer *SelectedMuonAux = new xAOD::MuonAuxContainer;
    SelectedMuon->setStore(SelectedMuonAux);
    xAOD::Muon* newMuon = new xAOD::Muon;
    newMuon->makePrivateStore(*(metMuons[0]));
    SelectedMuon->push_back(newMuon);

    m_trig_sf->applySystematicVariation(CP::SystematicSet()); //nominal
    if (m_trig_sf->getTriggerScaleFactor(*SelectedMuon, trgSF, "HLT_mu50") != CP::CorrectionCode::Ok) {
      Error("execute()", "m_trig_sf returns not Ok CorrectionCode");
      return EL::StatusCode::FAILURE;
    }
      
    m_effi_corr->applySystematicVariation(CP::SystematicSet()); //nominal
    if (m_effi_corr->getEfficiencyScaleFactor(*(metMuons[0]),  recoEffSF ) != CP::CorrectionCode::Ok) {
      if (m_effi_corr->getEfficiencyScaleFactor(*(metMuons[0]), recoEffSF) == CP::CorrectionCode::OutOfValidityRange)
        recoEffSF = 1.0;
      else {
        Error("execute()", "m_effi_corr returns not Ok CorrectionCode");
        return EL::StatusCode::FAILURE;
      }
    }
      
    m_effi_corr_iso->applySystematicVariation(CP::SystematicSet()); //nominal
    if (m_effi_corr_iso->getEfficiencyScaleFactor(*(metMuons[0]),  isoSF ) != CP::CorrectionCode::Ok) {
      Error("execute()", "m_effi_corr_iso returns not Ok CorrectionCode");
      return EL::StatusCode::FAILURE;
    }

    delete SelectedMuon;
    delete SelectedMuonAux; 
    
    SFWeight = trgSF * recoEffSF * isoSF * puWeight;
  
  }
  double totalWeight = 1.0;
  double totalWeight_wo_trgSF = 1.0;
  double totalWeight_wo_recoEffSF = 1.0;
  double totalWeight_wo_isoSF = 1.0;
  if (m_isMC){ /// WARNING is that correct? SF and weights only for MCs?
    totalWeight = weighfilterEfficiency*weightkFactor*(TOTAL_LUMI/sampleLumi);
    totalWeight_wo_trgSF =      totalWeight         * recoEffSF * isoSF * puWeight;
    totalWeight_wo_recoEffSF =  totalWeight * trgSF             * isoSF * puWeight;
    totalWeight_wo_isoSF =      totalWeight * trgSF * recoEffSF         * puWeight;
    totalWeight *= SFWeight;
  }
  
  m_HistObjectDumper->plotMuon(metMuons[0],"final_noMET_mT_cuts",totalWeight);
  m_HistObjectDumper->plotMtAndMet(metMuons[0],finalTrkMet,"final_noMET_mT_cuts",totalWeight);
  
  ///***************************************************************************
  /// calculate QCD nkg. hist
  if (!m_isMC){
    double QCD_weight = getQCDWeight(metMuons[0]->pt() * 0.001,true);
    m_HistObjectDumper->plotMuon(metMuons[0],"final_noMET_mT_cuts_QCD",QCD_weight);
    m_HistObjectDumper->plotMtAndMet
    (metMuons[0],finalTrkMet,"final_noMET_mT_cuts_QCD",QCD_weight);
  }
  ///***************************************************************************
  
  if (!m_isMC)
    plotPtBinnedHists(metMuons[0],finalTrkMet,"final_noMET_mT_cuts",totalWeight);
  
  m_HistObjectDumper->plotMuon(metMuons[0],"final_noMET_mT_cuts_noTrgSF",totalWeight_wo_trgSF);
  m_HistObjectDumper->plotMtAndMet(metMuons[0],finalTrkMet,"final_noMET_mT_cuts_noTrgSF",totalWeight_wo_trgSF);
  
  m_HistObjectDumper->plotMuon(metMuons[0],"final_noMET_mT_cuts_noRecoEffSF",totalWeight_wo_recoEffSF);
  m_HistObjectDumper->plotMtAndMet(metMuons[0],finalTrkMet,"final_noMET_mT_cuts_noRecoEffSF",totalWeight_wo_recoEffSF);
  
  m_HistObjectDumper->plotMuon(metMuons[0],"final_noMET_mT_cuts_noIsoSF",totalWeight_wo_isoSF);
  m_HistObjectDumper->plotMtAndMet(metMuons[0],finalTrkMet,"final_noMET_mT_cuts_noIsoSF",totalWeight_wo_isoSF);
  
  /// see presentation by Saminder Dhaliwal in Lepton+X meeting on March 22.
  /// https://indico.cern.ch/event/512374/
  double missingVetoRegion = false;
  double absMuPhi = abs(metMuons[0]->phi());
  double absMuEta = abs(metMuons[0]->eta());
  if (absMuEta>1.05 && absMuEta<1.3){
    if ((absMuEta>0.21 && absMuEta<0.57) || (absMuEta>1.00 && absMuEta<1.33) || 
        (absMuEta>1.78 && absMuEta<2.14) || (absMuEta>2.57 && absMuEta<2.93))
    {
      missingVetoRegion = true;
    }
  }
  
  if (missingVetoRegion==false){
    m_HistObjectDumper->plotMuon(metMuons[0],"final_noMET_mT_cuts_updatedHighPtVeto",totalWeight);
    m_HistObjectDumper->plotMtAndMet(metMuons[0],finalTrkMet,"final_noMET_mT_cuts_updatedHighPtVeto",totalWeight);
  }
  
  double metCut = 55.0;
  double mtCut = 110.0;
  
  if (missingEt<metCut)
    return EL::StatusCode::SUCCESS;
  m_BitsetCutflow->FillCutflow("MET");
  
  double leptonEt;
  double leptonPhi;
  
  if (m_runElectronChannel){
    leptonEt = metElectrons[0]->pt() * 0.001;
    leptonPhi = metElectrons[0]->phi();
  }
  else{
    leptonEt = metMuons[0]->pt() * 0.001;
    leptonPhi = metMuons[0]->phi();
  }
  
  float deltaPhi = TMath::Abs(leptonPhi-missingEtPhi);
  if(deltaPhi > TMath::Pi())
        deltaPhi = TMath::TwoPi() - deltaPhi;
  
  double mT = sqrt( 2 * missingEt * leptonEt * 
  (1 - TMath::Cos( deltaPhi ) ) );
  
  if (mT<mtCut)
    return EL::StatusCode::SUCCESS;
  m_BitsetCutflow->FillCutflow("mT");
  
  
  
  h_event_sampleLumi->Fill(sampleLumi);
  h_event_kFactor->Fill(weightkFactor);
  h_event_filterEfficiency->Fill(weighfilterEfficiency);
  h_event_SFWeight->Fill(SFWeight);
  h_event_totalWeight->Fill(totalWeight);
  
  hMu_pt_off->Fill(leptonEt,totalWeight);
  hMu_mt_off->Fill(mT,totalWeight);
  hMu_MET_Muons_off->Fill(missingEt,totalWeight);

  m_HistObjectDumper->plotMuon(metMuons[0],"final",totalWeight);
  m_HistObjectDumper->plotMtAndMet(metMuons[0],finalTrkMet,"final",totalWeight);
  
  if (!m_isMC)
    plotPtBinnedHists(metMuons[0],finalTrkMet,"final",totalWeight);
  
  ///***************************************************************************
  /// calculate QCD nkg. hist
  if (!m_isMC){
    double QCD_weight = getQCDWeight(metMuons[0]->pt() * 0.001,true);
    m_HistObjectDumper->plotMuon(metMuons[0],"final_QCD",QCD_weight);
    m_HistObjectDumper->plotMtAndMet
    (metMuons[0],finalTrkMet,"final_QCD",QCD_weight);
  }
  ///***************************************************************************
  
  m_HistObjectDumper->plotMuon(metMuons[0],"final_noWeight",1);
  m_HistObjectDumper->plotMtAndMet(metMuons[0],finalTrkMet,"final_noWeight",1);
  
  m_HistObjectDumper->plotMuon(metMuons[0],"final_noTrgSF",totalWeight_wo_trgSF);
  m_HistObjectDumper->plotMtAndMet(metMuons[0],finalTrkMet,"final_noTrgSF",totalWeight_wo_trgSF);
  
  m_HistObjectDumper->plotMuon(metMuons[0],"final_noRecoEffSF",totalWeight_wo_recoEffSF);
  m_HistObjectDumper->plotMtAndMet(metMuons[0],finalTrkMet,"final_noRecoEffSF",totalWeight_wo_recoEffSF);
  
  m_HistObjectDumper->plotMuon(metMuons[0],"final_noIsoSF",totalWeight_wo_isoSF);
  m_HistObjectDumper->plotMtAndMet(metMuons[0],finalTrkMet,"final_noIsoSF",totalWeight_wo_isoSF);
  
  if (missingVetoRegion==false){
    m_HistObjectDumper->plotMuon(metMuons[0],"final_updatedHighPtVeto",totalWeight);
    m_HistObjectDumper->plotMtAndMet(metMuons[0],finalTrkMet,"final_updatedHighPtVeto",totalWeight);
  }
  
  /// WARNING cutflow
//   cout << "[CUTFLOW]\t" << endl;
//   cout << "[CUTFLOW]\t Scale Factor cross checks" << endl;
//   cout << "[CUTFLOW]\t EventNumber = " << EventNumber << endl;
//   cout << "[CUTFLOW]\t mT = " << mT << endl;
//   cout << "[CUTFLOW]\t leptonEt = " << leptonEt << endl;
//   cout << "[CUTFLOW]\t missingEt = " << missingEt << endl;
//   cout << "[CUTFLOW]\t metMuons[0]->eta() = " << metMuons[0]->eta() << endl;
//   cout << "[CUTFLOW]\t trgSF = " << trgSF << endl;
//   cout << "[CUTFLOW]\t recoEffSF = " << recoEffSF << endl;
//   cout << "[CUTFLOW]\t isoSF = " << isoSF << endl;
//   cout << "[CUTFLOW]\t puWeight = " << puWeight << endl;
//   cout << "[CUTFLOW]\t weightkFactor = " << weightkFactor << endl;
//   cout << "[CUTFLOW]\t xSec = " << m_LPXKfactorTool->getMCCrossSection() << endl;
//   cout << "[CUTFLOW]\t nEvents = " << m_DataSetInfo->getnEventsPerSample(m_DSID) << endl;
//   cout << "[CUTFLOW]\t sampleLumi = " << sampleLumi << endl;
//   cout << "[CUTFLOW]\t" << endl;
  
  return EL::StatusCode::SUCCESS;
}



