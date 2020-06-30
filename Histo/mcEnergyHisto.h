#ifndef MCENERGYHISTO_H_
#define MCENERGYHISTO_H_

#include "algorithm/Algorithm.h"

// HerdSoftware headers
#include "dataobjects/MCTruth.h"

using namespace EA;

class TH1D;

class mcEnergyHisto : public Algorithm
{
public:
  mcEnergyHisto(const std::string &name);
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

  // Created global objects
  std::shared_ptr<TH1D> histo; // Objects to be pushed on global store must be held by a shared_ptr

  // Utility variables
  observer_ptr<EventDataStore> _evStore; // Pointer to the event data store
};

#endif /* ENEHISTOHISTO_H_ */
