/* Copyright 2016 Oleksandr Viazlo */

#ifndef DataSetInfo_H
#define DataSetInfo_H  
  
#include <TROOT.h>
#include <map>
#include <iostream>
  
class DataSetInfo {

public:

  DataSetInfo();
  ~DataSetInfo();
  double getnEventsPerSample(unsigned int datasetID);
  /// return cross-section in nb
  double getCrossSection(unsigned int datasetID);
  
private:

  std::map<unsigned int, double> m_nEventsMap; //!
  std::map<unsigned int, double> m_crossSectionMap; //!
  /// this is needed to distribute the algorithm to the workers
  /// WARNING works ONLY if include section contains some ROOT header!!!
  ClassDef(DataSetInfo, 1);
  
};
  
#endif
