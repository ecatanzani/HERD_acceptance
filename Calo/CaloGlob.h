#ifndef CALOGLOB_H_
#define CALOGLOB_H_

#include "algorithm/Algorithm.h"

// HerdSoftware headers

using namespace EA;

class CaloGlobStore;

class CaloGlob : public Algorithm {
public:
  CaloGlob(const std::string &name);

  bool Initialize();

  bool Process();

  bool Finalize();

private:
  bool filterenable;
  bool calohitscutmc;


  //Store pointer
  std::shared_ptr<CaloGlobStore> _processstore;
  

  // Utility variables
  observer_ptr<EventDataStore> _evStore; // Pointer to the event data store

};

class CaloGlobStore : public Algorithm {
public:
  
  CaloGlobStore(const std::string &name);
  bool Initialize();

  bool Process();

  bool Finalize();

  bool Reset();

  int calonhits;
  float calototedep;
  int calonclusters;
  
private:
  // Algorithm parameters
  // Created global objects

  // Utility variables
};

#endif /* CALOGLOB_H_ */
