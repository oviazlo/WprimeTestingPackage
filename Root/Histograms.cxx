///*****************************************************************************
///
/// Implementation of histInitialize function:
///     - histInitialize:
///             * create and register in worker all output histograms
///
///*****************************************************************************
#include <MyAnalysis/RecoAnalysis.h>

EL::StatusCode RecoAnalysis :: histInitialize ()
{
  /// Here you do everything that needs to be done at the very
  /// beginning on each worker node, e.g. create histograms and output
  /// trees.  This method gets called before any input files are
  /// connected.
 
  h_mgen_all = (TH1D*)WprimeHist::standard("mgen","hAll","","");
  wk()->addOutput(h_mgen_all); 
  
  h_mgen = (TH1D*)WprimeHist::standard("mgen","h","","");
  wk()->addOutput(h_mgen); 
  
  hMu_pt_off = (TH1D*)WprimeHist::standard("pt","h","","");
  wk()->addOutput(hMu_pt_off); 
  
  hMu_mt_off = (TH1D*)WprimeHist::standard("mt","h","","");
  wk()->addOutput(hMu_mt_off);

  hMu_MET_Muons_off = (TH1D*)WprimeHist::standard("met","h","","");
  wk()->addOutput(hMu_MET_Muons_off); 
  
  h_event_crossSectionWeight = (TH1D*)WprimeHist::standard("evtwt_xSec","h","",
                                                           "");
  wk()->addOutput(h_event_crossSectionWeight);
  
  h_event_kFactor = (TH1D*)WprimeHist::standard("evtwt_kFactor","h","","");
  wk()->addOutput(h_event_kFactor);
  
  h_event_filterEfficiency = (TH1D*)WprimeHist::standard("evtwt_filterEff","h",
                                                         "","");
  wk()->addOutput(h_event_filterEfficiency);
  
  h_event_totalWeight = (TH1D*)WprimeHist::standard("evtwt","h","","");
  wk()->addOutput(h_event_totalWeight);
  
  return EL::StatusCode::SUCCESS;
}