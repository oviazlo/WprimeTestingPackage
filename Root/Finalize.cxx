///*****************************************************************************
///
/// Implementation of all functions which are called AFTER execute:
///     - postExecute
///     - finalize:
///             * deleting of all tools
///     - histFinalize
///
///*****************************************************************************
#include <MyAnalysis/MyxAODAnalysis.h>

/// this is needed to distribute the algorithm to the workers
ClassImp(MyxAODAnalysis)

EL::StatusCode MyxAODAnalysis :: postExecute ()
{
  /// Here you do everything that needs to be done after the main event
  /// processing.  This is typically very rare, particularly in user
  /// code.  It is mainly used in implementing the NTupleSvc.
  return EL::StatusCode::SUCCESS;
}



EL::StatusCode MyxAODAnalysis :: finalize ()
{
  /// This method is the mirror image of initialize(), meaning it gets
  /// called after the last event has been processed on the worker node
  /// and allows you to finish up any objects you created in
  /// initialize() before they are written to disk.  This is actually
  /// fairly rare, since this happens separately for each worker node.
  /// Most of the time you want to do your post-processing on the
  /// submission node after all your histogram outputs have been
  /// merged.  This is different from histFinalize() in that it only
  /// gets called on worker nodes that processed input events.
  
  const char* APP_NAME = "MyxAODAnalysis";  

  /// push cutflow for last event
  /// FIXME probably this breaks possibility to use PROOF
  m_BitsetCutflow->PushBitSet();
  
  ///*************************
  ///
  /// deleting of all tools
  ///
  /// ************************
  
  for (int i=0; i<17; i++) {
      cout << "cut " << i << " value: " << count[i]<<endl;
    }
  
  /// GRL
  if( m_grl ) {
    delete m_grl;
    m_grl = 0;
  }
  /// Jet cleaning tool
//   if( m_jetCleaning ) {
//     delete m_jetCleaning;
//     m_jetCleaning = 0;
//   }
//   
//   if( m_jetCalibrationTool ) {
//     delete m_jetCalibrationTool;
//     m_jetCalibrationTool = 0;
//   }
  
  
  /*
  if(m_JERTool){
    delete m_JERTool;
    m_JERTool = 0;
  }
  */
  Info("finalize()", "Number of clean events = %i", m_numCleanEvents);
  /// Muon calibration and smearing tool
  if(m_muonCalibrationAndSmearingTool){
    delete m_muonCalibrationAndSmearingTool;
    m_muonCalibrationAndSmearingTool = 0;
  }
  /// Muon selector tool
  if(m_muonSelection){
    delete m_muonSelection;
    m_muonSelection = 0;
  }
  //~ /// Muon eff.corr. tool
  //~ if(m_effi_corr){
    //~ delete m_effi_corr;
    //~ m_effi_corr = 0;
  //~ }
  /*
  if(*m_METUtil){
      delete m_METUtil;
      m_METUtil = 0;
    }
  */
  
  if(m_useHistObjectDumper && m_HistObjectDumper){
    delete m_HistObjectDumper;
    m_HistObjectDumper = 0;
  }
  
  if(m_useBitsetCutflow && m_BitsetCutflow){
    delete m_BitsetCutflow;
    m_BitsetCutflow = 0;
  }
  
  /// cleaning up trigger tools
  if( m_trigConfigTool ) {
      delete m_trigConfigTool;
      m_trigConfigTool = 0;
   }
   if( m_trigDecisionTool ) {
      delete m_trigDecisionTool;
      m_trigDecisionTool = 0;
   }

  return EL::StatusCode::SUCCESS;
}



EL::StatusCode MyxAODAnalysis :: histFinalize ()
{
  /// This method is the mirror image of histInitialize(), meaning it
  /// gets called after the last event has been processed on the worker
  /// node and allows you to finish up any objects you created in
  /// histInitialize() before they are written to disk.  This is
  /// actually fairly rare, since this happens separately for each
  /// worker node.  Most of the time you want to do your
  /// post-processing on the submission node after all your histogram
  /// outputs have been merged.  This is different from finalize() in
  /// that it gets called on all worker nodes regardless of whether
  /// they processed input events.
  return EL::StatusCode::SUCCESS;
}


