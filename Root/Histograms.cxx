#include <MyAnalysis/MyxAODAnalysis.h>

EL::StatusCode MyxAODAnalysis::initHists (){
  
  /// [MET histos]
  h_MET_RefJet = new TH1F("h_MET_RefJet", "h_MET_RefJet", 2000, 0, 2000);
  wk()->addOutput (h_MET_RefJet);
  h_MET_RefMuon = new TH1F("h_MET_RefMuon", "h_MET_RefMuon", 2000, 0, 2000);
  wk()->addOutput (h_MET_RefMuon);
  h_MET_RefTau = new TH1F("h_MET_RefTau", "h_MET_RefTau", 2000, 0, 2000);
  wk()->addOutput (h_MET_RefTau);
  h_MET_RefSoft = new TH1F("h_MET_RefSoft", "h_MET_RefSoft", 2000, 0, 2000);
  wk()->addOutput (h_MET_RefSoft);
  h_MET_RefTrack = new TH1F("h_MET_RefTrack", "h_MET_RefTrack", 2000, 0, 2000);
  wk()->addOutput (h_MET_RefTrack);
  h_MET_RefGamma = new TH1F("h_MET_RefGamma", "h_MET_RefGamma", 2000, 0, 2000);
  wk()->addOutput (h_MET_RefGamma);
  h_MET_RefElectron = new TH1F("h_MET_RefElectron", "h_MET_RefElectron", 2000, 0, 2000);
  wk()->addOutput (h_MET_RefElectron);
  h_MET_RefCore = new TH1F("h_MET_RefCore", "h_MET_RefCore", 2000, 0, 2000);
  wk()->addOutput (h_MET_RefCore);
  
  h_MET_OldJet = new TH1F("h_MET_OldJet", "h_MET_OldJet", 2000, 0, 2000);
  wk()->addOutput (h_MET_OldJet);
  h_MET_OldMuon = new TH1F("h_MET_OldMuon", "h_MET_OldMuon", 2000, 0, 2000);
  wk()->addOutput (h_MET_OldMuon);
  h_MET_OldTau = new TH1F("h_MET_OldTau", "h_MET_OldTau", 2000, 0, 2000);
  wk()->addOutput (h_MET_OldTau);
  h_MET_OldSoft = new TH1F("h_MET_OldSoft", "h_MET_OldSoft", 2000, 0, 2000);
  wk()->addOutput (h_MET_OldSoft);
  h_MET_OldTrack = new TH1F("h_MET_OldTrack", "h_MET_OldTrack", 2000, 0, 2000);
  wk()->addOutput (h_MET_OldTrack);
  h_MET_OldGamma = new TH1F("h_MET_OldGamma", "h_MET_OldGamma", 2000, 0, 2000);
  wk()->addOutput (h_MET_OldGamma);
  h_MET_OldElectron = new TH1F("h_MET_OldElectron", "h_MET_OldElectron", 2000, 0, 2000);
  wk()->addOutput (h_MET_OldElectron); 
  
}