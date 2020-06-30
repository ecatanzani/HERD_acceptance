#ifndef MCTRUTHPROCESS_H_
#define MCTRUTHPROCESS_H_

#include "algorithm/Algorithm.h"

// HerdSoftware headers
#include "dataobjects/MCTruth.h"
#include "dataobjects/CaloHits.h"
#include "dataobjects/TrackInfoForCalo.h"
#include "dataobjects/StkIntersections.h"
#include "dataobjects/CaloGeoParams.h"

using namespace EA;

class TH1F;
class TH2F;
class TH3F;
class TGraph;
class TGraph2D;
class TVector3;
class MCtruthProcessStore;


class MCtruthProcess : public Algorithm {
public:

  MCtruthProcess(const std::string &name);
  bool Initialize();
  bool Process();
  bool Finalize();

private:

  //Store pointer
  std::shared_ptr<MCtruthProcessStore> _processstore;
  
  // Algorithm parameters
  bool filterenable;
  int minstkintersections;
  bool printcalocubemap;
  float mincalotrackx0;
  bool notfrombottom;

  // Created global objects
  // std::shared_ptr<TH1F> _histo; // Objects to be pushed on global store must be held by a shared_ptr
  std::shared_ptr<TGraph> _gdiscarded;
  std::shared_ptr<TH1F> _hstkintersections;
  std::shared_ptr<TH2F> _hgencthetaphi;
  std::shared_ptr<TH3F> _hgencoo;
  std::shared_ptr<TGraph2D> _ggencoo;
  std::shared_ptr<TGraph2D> _gcaloentry;
  std::shared_ptr<TGraph2D> _gcaloexit;
  std::shared_ptr<TH2F>_hcaloentryexitdir;
  std::shared_ptr<TH1F>_hshowerlength[Herd::RefFrame::NDirections][Herd::RefFrame::NDirections];
  std::shared_ptr<TH1F>_hshowerlengthall;

  //std::shared_ptr<TH2F> _hgencoo;

  // Utility variables
    void PrintCaloCubeMap();

  observer_ptr<EventDataStore> _evStore; // Pointer to the event data store

  TVector3 InterceptX(double, const TVector3 &, const TVector3 &) const;
};

class MCtruthProcessStore : public Algorithm {
public:
  MCtruthProcessStore(const std::string &name);
  bool Initialize();
  bool Process();
  bool Finalize();
  bool Reset();

  int mcNdiscarded;
  float mcDir[3];
  float mcCoo[3];
  float mcMom;
  float mcPhi;
  float mcCtheta;
  int mcStkintersections;
  float mcTracklengthcalox0;
  float mcTracklengthlysox0;
  float mcTrackcaloentry[3];
  float mcTrackcaloexit[3];
  int mcTrackcaloentryplane;
  int mcTrackcaloexitplane;

private:
  
};

#endif /* MCTRUTHPROCESS_H_ */
