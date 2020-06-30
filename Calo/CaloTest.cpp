// Example headers
#include "CaloTest.h"
#include "dataobjects/CaloHits.h"
#include "dataobjects/CaloClusters.h"
#include "dataobjects/CaloGeoParams.h"
#include "dataobjects/MCTruth.h"

// Root headers

// C/C++ standard headers
#include <numeric>
#include <cmath>

RegisterAlgorithm(CaloTest);


CaloTest::CaloTest(const std::string &name) :
  Algorithm{name}
   {
    
  }

bool CaloTest::Initialize() {
  const std::string routineName("CaloTest::Initialize");
  _evStore = GetDataStoreManager()->GetEventDataStore("evStore"); if (!_evStore) { COUT(ERROR) << "Event data store not found." << ENDL; return false; }
  return true;
}

bool CaloTest::Process() {
  const std::string routineName("CaloTest::Process");

  //Add the ProcessStore object for this event to the event data store
  auto globStore = GetDataStoreManager()->GetGlobalDataStore("globStore");
  if (!globStore) {COUT(ERROR) << "Global data store not found." << ENDL;}
  auto caloHits = _evStore->GetObject<Herd::CaloHits>("caloHitsMC");
  if (!caloHits) { COUT(DEBUG) << "CaloHitsMC not present for event " << GetEventLoopProxy()->GetCurrentEvent() << ENDL; return false; }
  auto caloGeoParams = globStore->GetObject<Herd::CaloGeoParams>("caloGeoParams");
  if (!caloGeoParams) { COUT(DEBUG) << "caloGeoParams not present for event " << GetEventLoopProxy()->GetCurrentEvent() << ENDL; return false; }
  auto caloClusters = _evStore->GetObject<Herd::CaloClusters>("caloClusters");
  if (!caloClusters) { COUT(DEBUG) << "caloClusters not present for event " << GetEventLoopProxy()->GetCurrentEvent() << ENDL; return false; }


  float calototedep = std::accumulate(caloHits->begin(), caloHits->end(), 0.,[](float sum, const Herd::Hit &hit) { return sum + hit.EDep(); });
  int calonhits =     std::accumulate(caloHits->begin(), caloHits->end(), 0.,[](int n, const Herd::Hit &hit) { if( hit.EDep()>0) return n+1; });
  int calonclusters = std::accumulate(caloClusters->begin(), caloClusters->end(), 0.,[](int n, const Herd::CaloHits &calohit) { return n+1; });
  
  COUT(INFO)<<calototedep<<" "<<calonhits<<" "<<caloClusters->size()<<ENDL;

return true;
}

bool CaloTest::Finalize() {
  const std::string routineName("CaloTest::Finalize");
  return true;
}


