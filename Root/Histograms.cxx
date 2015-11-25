///*****************************************************************************
///
/// Implementation of histInitialize function:
///     - histInitialize:
///             * create and register in worker all output histograms
///
///*****************************************************************************
#include <MyAnalysis/MyxAODAnalysis.h>
#include <TCollection.h>
#include <TKey.h>

EL::StatusCode MyxAODAnalysis :: histInitialize ()
{
  /// Here you do everything that needs to be done at the very
  /// beginning on each worker node, e.g. create histograms and output
  /// trees.  This method gets called before any input files are
  /// connected.

    /// event hists
//   h_event_nJets = new TH1F("h_event_nJets", "h_event_nJets", 13, 0, 12);
//   h_event_nJets_Pt40 = new TH1F("h_event_nJets_Pt40", "h_event_nJets_Pt40", 
//                                 13, 0, 12);
// 
//   /// lepton hists
//   h_lepton_d0Significance = 
//   new TH1F("h_lepton_d0Significance", "h_lepton_d0Significance", 50, 0, 3);
//   h_lepton_z0SinTheta =
//   new TH1F("h_lepton_z0SinTheta", "h_lepton_z0SinTheta", 50, 0, 1);
//   h_lepton_Pt = new TH1F("h_lepton_Pt", "h_lepton_Pt", 
//   h_lepton_eta = new TH1F("h_lepton_eta", "h_lepton_eta", 
//   h_lepton_phi = new TH1F("h_lepton_phi", "h_lepton_phi", 
// 
//   /// leptonMET hists
//   h_leptonMET_deltaPhi = new TH1F("h_leptonMET_deltaPhi", "h_leptonMET_deltaPhi", 
//                                   50, 0, TMath::Pi());
//   h_leptonMET_Mt = new TH1F("h_leptonMET_Mt", "h_leptonMET_Mt", 
//   
//   /// MET hists
//   h_MET_Et = new TH1F("h_MET_Et", "h_MET_Et", 
//   h_MET_phi = new TH1F("h_MET_phi", "h_MET_phi", 
    
  hMu_pt_off = (TH1D*)WprimeHist::standard("pt","","","");
//   hMu_pt_off->SetName("hMu_pt_off");
  wk()->addOutput(hMu_pt_off); 
  
  hMu_mt_off = (TH1D*)WprimeHist::standard("mt","","","");
//   hMu_mt_off->SetName("hMu_mt_off");
  wk()->addOutput(hMu_mt_off);

  hMu_MET_Muons_off = (TH1D*)WprimeHist::standard("met","","","");
//   hMu_MET_Muons_off->SetName("hMu_MET_Muons_off");
  wk()->addOutput(hMu_MET_Muons_off); 
  
  TIter next(gDirectory->GetListOfKeys());
  TKey *key;
  cout << "List of created hists:" << endl;
  while ((key = (TKey*)next())) {
    cout << key->GetName() << endl;
    TClass *cl = gROOT->GetClass(key->GetClassName());
    if (!cl->InheritsFrom("TH1")) continue;
    TH1 *h = (TH1*)key->ReadObj();
    cout << h->GetName() << endl;
  }
  
  
  
  
  
  /// [Jet histos]
//   h_zPrimVtx = new TH1F("h_zPrimVtx", "h_zPrimVtx", 6000, -300.0, 300); 
//   wk()->addOutput (h_zPrimVtx);
//   h_jetPt_VeryLooseBadJets = new TH1F("h_jetPt_VeryLooseBadJets", 
//                                       "VeryLooseBadJets and better with pt>20"
//                                       , 300, 0, 3000); /// jet pt [GeV]
//   wk()->addOutput (h_jetPt_VeryLooseBadJets);
//   h_jetPt_LooseBadJets = new TH1F("h_jetPt_LooseBadJets", 
//                                   "LooseBadJets and better with pt>20", 300, 0, 
//                                   3000); /// jet pt [GeV]
//   wk()->addOutput (h_jetPt_LooseBadJets);
//   
//   /// [Muon histos]
//   h_MET_RefFinalFix = new TH1F("h_MET_RefFinalFix", "h_MET_RefFinalFix", 1000, 
//                                0, 10000); 
//   wk()->addOutput (h_MET_RefFinalFix);
//   h_nSelectedMuons = new TH1I("h_nSelectedMuons", "h_nSelectedMuons", 5, 0, 5); 
//   wk()->addOutput (h_nSelectedMuons);
//   h_ptmuon = new TH1F("h_ptmuon", "h_ptmuon", 500, 0, 500);
//   wk()->addOutput (h_ptmuon);
//   
//   /// [MET histos]
//   h_MET_RefJet = new TH1F("h_MET_RefJet", "h_MET_RefJet", 2000, 0, 2000);
//   wk()->addOutput (h_MET_RefJet);
//   h_MET_RefMuon = new TH1F("h_MET_RefMuon", "h_MET_RefMuon", 2000, 0, 2000);
//   wk()->addOutput (h_MET_RefMuon);
//   h_MET_RefTau = new TH1F("h_MET_RefTau", "h_MET_RefTau", 2000, 0, 2000);
//   wk()->addOutput (h_MET_RefTau);
//   h_MET_RefSoft = new TH1F("h_MET_RefSoft", "h_MET_RefSoft", 2000, 0, 2000);
//   wk()->addOutput (h_MET_RefSoft);
//   h_MET_RefTrack = new TH1F("h_MET_RefTrack", "h_MET_RefTrack", 2000, 0, 2000);
//   wk()->addOutput (h_MET_RefTrack);
//   h_MET_RefGamma = new TH1F("h_MET_RefGamma", "h_MET_RefGamma", 2000, 0, 2000);
//   wk()->addOutput (h_MET_RefGamma);
//   h_MET_RefElectron = new TH1F("h_MET_RefElectron", "h_MET_RefElectron", 2000, 0, 2000);
//   wk()->addOutput (h_MET_RefElectron);
//   h_MET_RefCore = new TH1F("h_MET_RefCore", "h_MET_RefCore", 2000, 0, 2000);
//   wk()->addOutput (h_MET_RefCore);
//   
//   h_MET_OldJet = new TH1F("h_MET_OldJet", "h_MET_OldJet", 2000, 0, 2000);
//   wk()->addOutput (h_MET_OldJet);
//   h_MET_OldMuon = new TH1F("h_MET_OldMuon", "h_MET_OldMuon", 2000, 0, 2000);
//   wk()->addOutput (h_MET_OldMuon);
//   h_MET_OldTau = new TH1F("h_MET_OldTau", "h_MET_OldTau", 2000, 0, 2000);
//   wk()->addOutput (h_MET_OldTau);
//   h_MET_OldSoft = new TH1F("h_MET_OldSoft", "h_MET_OldSoft", 2000, 0, 2000);
//   wk()->addOutput (h_MET_OldSoft);
//   h_MET_OldTrack = new TH1F("h_MET_OldTrack", "h_MET_OldTrack", 2000, 0, 2000);
//   wk()->addOutput (h_MET_OldTrack);
//   h_MET_OldGamma = new TH1F("h_MET_OldGamma", "h_MET_OldGamma", 2000, 0, 2000);
//   wk()->addOutput (h_MET_OldGamma);
//   h_MET_OldElectron = new TH1F("h_MET_OldElectron", "h_MET_OldElectron", 2000, 0, 2000);
//   wk()->addOutput (h_MET_OldElectron);
//     
//   h_MET_RefFinalFix_test = new TH1F("h_MET_RefFinalFix_test", 
//                                     "h_MET_RefFinalFix_test", 1000, 0, 10000); 
//   wk()->addOutput (h_MET_RefFinalFix_test);
//   h_truthMET = new TH1F("h_truthMET", "Truth MET", 1000, 0, 10000); 
//   wk()->addOutput (h_truthMET);
//   
//   /// [Mt histos]
//   h_Mt = new TH1F("h_Mt", "h_Mt", 1000, 0, 10000); 
//   wk()->addOutput (h_Mt);
// 
//   h_Mt_muonPtCut = new TH1F("h_Mt_muonPtCut", "h_Mt_muonPtCut", 1000, 0, 10000);
//   wk()->addOutput (h_Mt_muonPtCut);

  /// get the output file, create a new TTree and connect it to that output
  /// define what branches will go in that tree

//   TFile *outputFile = wk()->getOutputFile (outputName);
//   tree = new TTree ("tree", "tree");
//   tree->SetDirectory (outputFile);
//   tree->Branch("EventNumber", &EventNumber);
//   tree->Branch("RunNumber", &RunNumber);
//   tree->Branch("LumiBlock", &LumiBlock);
//   tree->Branch("MuonPt", &MuonPt);
//   tree->Branch("MuonEta", &MuonEta);
//   tree->Branch("MuonPhi", &MuonPhi);
//   tree->Branch("MissingEt", &MissingEt);
//   tree->Branch("MissingEx", &MissingEx);
//   tree->Branch("MissingEy", &MissingEy);
//   tree->Branch("OldMissingEt", &OldMissingEt);
//   tree->Branch("TransverseMass", &TransverseMass);

  return EL::StatusCode::SUCCESS;
}