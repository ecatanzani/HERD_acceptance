#ifndef MCGENSPECTRUM_H_
#define MCGENSPECTRUM_H_

#include "algorithm/Algorithm.h"
#include "analysis/AnalysisManager.h"

// HerdSoftware headers
#include "dataobjects/MCTruth.h"

using namespace EA;

class TH1D;

class mcGenSpectrum : public Algorithm {
public:
  mcGenSpectrum(const std::string &name);
  bool Initialize();
  bool Process();
  bool Finalize();

private:

  std::vector<double> axispar;
  std::vector<double> binning;
  bool logaxis;
  std::string title;
  
  void GenerateLogBinning();
  void GenerateBinning();
  

  unsigned int ngen;
  std::vector<double> momrange;
  
  std::shared_ptr<TH1D> histo; // Objects to be pushed on global store must be held by a shared_ptr

  // Utility variables
  observer_ptr<EventDataStore> _evStore; // Pointer to the event data store
};

#endif /* MCGENSPECTRUM */
