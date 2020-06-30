#ifndef CALOTEST_H_
#define CALOTEST_H_

#include "algorithm/Algorithm.h"

// HerdSoftware headers

using namespace EA;

class CaloTestStore;

class CaloTest : public Algorithm {
public:
  CaloTest(const std::string &name);

  bool Initialize();

  bool Process();

  bool Finalize();

private:

  // Utility variables
  observer_ptr<EventDataStore> _evStore; // Pointer to the event data store

};


#endif /* CALOTEST_H_ */
