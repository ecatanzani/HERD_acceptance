// Example headers
#include "CaloGlob.h"
#include "dataobjects/CaloHits.h"
#include "dataobjects/CaloClusters.h"
#include "dataobjects/CaloGeoParams.h"
#include "dataobjects/MCTruth.h"

// Root headers

// C/C++ standard headers
#include <numeric>
#include <cmath>

RegisterAlgorithm(CaloGlob);
RegisterAlgorithm(CaloGlobStore);


CaloGlob::CaloGlob(const std::string &name) :
  Algorithm{name},
  filterenable{true},
  calohitscutmc{false}
   {
    DefineParameter("filterenable",  filterenable); 
    DefineParameter("calohitscutmc", calohitscutmc);
  }

bool CaloGlob::Initialize() {
  const std::string routineName("CaloGlob::Initialize");

  _evStore = GetDataStoreManager()->GetEventDataStore("evStore"); if (!_evStore) { COUT(ERROR) << "Event data store not found." << ENDL; return false; }

  _processstore = std::make_shared<CaloGlobStore>("caloGlobStore");
  
  // Setup the filter                                                                                                                                                                                                                       
  if (filterenable) SetFilterStatus(FilterStatus::ENABLED); else SetFilterStatus(FilterStatus::DISABLED);

  return true;
}

bool CaloGlob::Process() {
  const std::string routineName("CaloGlob::Process");

  //Add the ProcessStore object for this event to the event data store
  //_processstore->Reset();
  _processstore->Reset();
  _evStore->AddObject("caloGlobStore",_processstore);

  //Set Filter Status
  SetFilterResult(FilterResult::ACCEPT);

  auto globStore = GetDataStoreManager()->GetGlobalDataStore("globStore");
  if (!globStore) {COUT(ERROR) << "Global data store not found." << ENDL;}
  //auto calotrack = _evStore->GetObject<Herd::TrackInfoForCalo>("trackInfoForCaloMC");
  //if (!calotrack) { COUT(DEBUG) << "TrackInfoForCalo  not present for event " << GetEventLoopProxy()->GetCurrentEvent() << ENDL; return false; }
  auto caloHits = _evStore->GetObject<Herd::CaloHits>("caloHitsMC");
  if (!caloHits) { COUT(DEBUG) << "CaloHitsMC not present for event " << GetEventLoopProxy()->GetCurrentEvent() << ENDL; return false; }
  auto caloGeoParams = globStore->GetObject<Herd::CaloGeoParams>("caloGeoParams");
  if (!caloGeoParams) { COUT(DEBUG) << "caloGeoParams not present for event " << GetEventLoopProxy()->GetCurrentEvent() << ENDL; return false; }
  //auto caloClusters = _evStore->GetObject<Herd::CaloClusters>("caloClusters");
  //if (!caloClusters) { COUT(DEBUG) << "caloClusters not present for event " << GetEventLoopProxy()->GetCurrentEvent() << ENDL; return false; }


  float calototedep = std::accumulate(caloHits->begin(), caloHits->end(), 0.,[](float sum, const Herd::Hit &hit) { return sum + hit.EDep(); });
  int calonhits =     std::accumulate(caloHits->begin(), caloHits->end(), 0.,[](int n, const Herd::Hit &hit) { if( hit.EDep()>0) return n+1; });
  _processstore->calonhits = calonhits;
  _processstore->calototedep = calototedep;
  //COUT(INFO)<<caloClusters->size()<<ENDL;
  //if( caloClusters ) _processstore->calonclusters = (int)caloClusters->size();

  if(calohitscutmc){
    auto mcTruth = _evStore->GetObject<Herd::MCTruth>("mcTruth");
    if (!mcTruth) {COUT(ERROR) << "mcTruth not present for event " << GetEventLoopProxy()->GetCurrentEvent() << ENDL;return false;}
    float mcmom = std::sqrt(mcTruth->primaries.at(0).initialMomentum * mcTruth->primaries.at(0).initialMomentum);
    if(calonhits < std::pow(10,(((1+log10(2))/3.))*std::log10(mcmom) + (2 - (1+std::log10(2))/3.)) ){ SetFilterResult(FilterResult::REJECT); }
  }

return true;
}

bool CaloGlob::Finalize() {
  const std::string routineName("CaloGlob::Finalize");

  auto globStore = GetDataStoreManager()->GetGlobalDataStore("globStore");
  if (!globStore) { COUT(ERROR) << "Global data store not found." << ENDL; return false; }

  return true;
}

//***************************

CaloGlobStore::CaloGlobStore(const std::string &name) :
  Algorithm{name}
   {
  }

  bool CaloGlobStore::Initialize() {
  const std::string routineName("CaloGlobStore::Initialize");
  Reset();
  return true;
}

  bool CaloGlobStore::Process() {
  const std::string routineName("CaloGlobStore::Process");
  return true;
}
  bool CaloGlobStore::Finalize() {
  const std::string routineName("CaloGlobStore::Finalize");
  return true;
}
bool CaloGlobStore::Reset() {
  const std::string routineName("CaloGlobStore::Finalize");

  calonhits=0;
  calototedep=0;
  calonclusters=0;

  return true;
}

