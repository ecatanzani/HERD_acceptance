#include "MCtruthProcess.h"

// Root headers
#include "TMath.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TVector3.h"
#include "TGraph.h"
#include "TGraph2D.h"

// C/C++ standard headers
#include <numeric>

RegisterAlgorithm(MCtruthProcess);
RegisterAlgorithm(MCtruthProcessStore);


MCtruthProcess::MCtruthProcess(const std::string &name) :
  Algorithm{name},
  printcalocubemap{false},
  filterenable{true},
  minstkintersections{-1},
  mincalotrackx0{-999},
  notfrombottom{true}
   {
     DefineParameter("minstkintersections", minstkintersections);
     DefineParameter("printcalocubemap",    printcalocubemap);
     DefineParameter("filterenable",        filterenable);
     DefineParameter("mincalotrackx0",      mincalotrackx0);
     DefineParameter("notfrombottom",       notfrombottom);

  }

bool MCtruthProcess::Initialize() {
  const std::string routineName("MCtruthProcess::Initialize");

  // Setup the filter                                                                                                                                                                                                                       
  if (filterenable) SetFilterStatus(FilterStatus::ENABLED); else SetFilterStatus(FilterStatus::DISABLED);

  _evStore = GetDataStoreManager()->GetEventDataStore("evStore");
  if (!_evStore) {
    COUT(ERROR) << "Event data store not found." << ENDL;
    return false;
  }

  if(printcalocubemap) PrintCaloCubeMap();

  _processstore = std::make_shared<MCtruthProcessStore>("MCtruthProcessStore");
  COUT(INFO) << "InitializedProcessStore::" <<_processstore << ENDL;

  // Create the histogram
  _gdiscarded       = std::make_shared<TGraph>();
  _gdiscarded->SetNameTitle("gdiscarded","Discarded Events before simulated");
  _hgencthetaphi = std::make_shared<TH2F>("hgencthetaphi","MCtruth Generation;cos(#theta);Phi (rad)",1000,-1,1,100,-TMath::Pi(),+TMath::Pi());
  _hgencoo       = std::make_shared<TH3F>("hgencoo",      "MCtruth Generation;X(cm);Y(cm);Z(cm)",    100,-500,500,100,-500,500,100,-500,500);
  _hstkintersections = std::make_shared<TH1F>("hstkintersections", "Inyersection of track with STK;Occurrence", 101,-1.5,99.5);
  _ggencoo       = std::make_shared<TGraph2D>();
  _ggencoo->SetNameTitle("ggencoo","Generation Coordinates;X(cm);Y(cm);Z(cm)");
  _gcaloentry       = std::make_shared<TGraph2D>();
  _gcaloentry->SetNameTitle("gcaloentry","CALO Entry point;X(cm);Y(cm);Z(cm)");
  _gcaloexit       = std::make_shared<TGraph2D>();
  _gcaloexit->SetNameTitle("gcaloexit","CALO Exit point;X(cm);Y(cm);Z(cm)");

  _hcaloentryexitdir = std::make_shared<TH2F>("hcaloentryexitdir","CALO Entry (X) - Exit (Y)",Herd::RefFrame::NDirections+1, -0.5, Herd::RefFrame::NDirections+0.5, Herd::RefFrame::NDirections+1, -0.5, Herd::RefFrame::NDirections+0.5);
  _hcaloentryexitdir->GetXaxis()->SetBinLabel(_hcaloentryexitdir->GetNbinsX(),"NONE");
  for(int ibin=1; ibin<_hcaloentryexitdir->GetNbinsX(); ibin++) _hcaloentryexitdir->GetXaxis()->SetBinLabel(ibin, Herd::RefFrame::DirectionName[ibin-1].c_str());
  _hcaloentryexitdir->GetYaxis()->SetBinLabel(_hcaloentryexitdir->GetNbinsY(),"NONE");
  for(int ibin=1; ibin<_hcaloentryexitdir->GetNbinsY(); ibin++) _hcaloentryexitdir->GetYaxis()->SetBinLabel(ibin, Herd::RefFrame::DirectionName[ibin-1].c_str());
  
  _hshowerlengthall  = std::make_shared<TH1F>(Form("hshowerlengthall"), Form("Shower Lenght (X0) All"),200,0,100);
  for(int indir=0; indir < Herd::RefFrame::NDirections; indir++){
    for(int outdir=0; outdir < Herd::RefFrame::NDirections; outdir++){
      _hshowerlength[indir][outdir]  = std::make_shared<TH1F>(Form("hshowerlength_%d_%d",indir,outdir), Form("Shower Lenght (X0) [%s-%s]",Herd::RefFrame::DirectionName[indir].c_str(),Herd::RefFrame::DirectionName[outdir].c_str()),200,0,100);
      //COUT(DEBUG)<<Form("%s %s",_hshowerlength[indir][outdir]->GetName(),_hshowerlength[indir][outdir]->GetTitle())<<ENDL;
    }
  }
  return true;
}

bool MCtruthProcess::Process() {
  const std::string routineName("MCtruthProcess::Process");

  //Add the ProcessStore object for this event to the event data store
  //_processstore->Reset();
  _evStore->AddObject("MCtruthProcessStore",_processstore);

  //Set Filter Status
  SetFilterResult(FilterResult::ACCEPT);

  auto mctruth = _evStore->GetObject<Herd::MCTruth>("mcTruth");
  if (!mctruth) { COUT(DEBUG) << "MCTruth not present for event " << GetEventLoopProxy()->GetCurrentEvent() << ENDL; return false; }
  auto calotrack = _evStore->GetObject<Herd::TrackInfoForCalo>("trackInfoForCaloMC");
  if (!calotrack) { COUT(DEBUG) << "TrackInfoForCalo  not present for event " << GetEventLoopProxy()->GetCurrentEvent() << ENDL; return false; }
  auto stkintersections = _evStore->GetObject<Herd::StkIntersections>("stkIntersectionsMC");
  if (!stkintersections) { COUT(DEBUG) << "StkIntersections not present for event " << GetEventLoopProxy()->GetCurrentEvent() << ENDL; return false; }

  _gdiscarded->SetPoint(_gdiscarded->GetN(), _gdiscarded->GetN(), mctruth->nDiscarded);

  const auto &primary = mctruth->primaries.at(0);
  Herd::Point gencoo = primary.initialPosition;
	TVector3 genmom (primary.initialMomentum[Herd::RefFrame::Coo::X],primary.initialMomentum[Herd::RefFrame::Coo::Y],primary.initialMomentum[Herd::RefFrame::Coo::Z]);
  Double_t genctheta = genmom.CosTheta();
  Double_t genphi = genmom.Phi();
  Double_t mom = genmom.Mag();
  _hgencoo->Fill(gencoo[Herd::RefFrame::Coo::X],gencoo[Herd::RefFrame::Coo::Y],gencoo[Herd::RefFrame::Coo::Z]);
  _ggencoo->SetPoint(_ggencoo->GetN(),gencoo[Herd::RefFrame::Coo::X],gencoo[Herd::RefFrame::Coo::Y],gencoo[Herd::RefFrame::Coo::Z]);
  _hgencthetaphi->Fill(genctheta,genphi);

  int nstkintersections = static_cast<int>(stkintersections->intersections.size());
  _hstkintersections->Fill(nstkintersections);
  
  _processstore->mcDir[0] = primary.initialMomentum[Herd::RefFrame::Coo::X] / genmom.Mag();
  _processstore->mcDir[1] = primary.initialMomentum[Herd::RefFrame::Coo::Y] / genmom.Mag();
  _processstore->mcDir[2] = primary.initialMomentum[Herd::RefFrame::Coo::Z] / genmom.Mag();
  _processstore->mcNdiscarded = mctruth->nDiscarded;
  _processstore->mcCoo[0] = primary.initialPosition[Herd::RefFrame::Coo::X];
  _processstore->mcCoo[1] = primary.initialPosition[Herd::RefFrame::Coo::Y];
  _processstore->mcCoo[2] = primary.initialPosition[Herd::RefFrame::Coo::Z];
  _processstore->mcMom = genmom.Mag();
  _processstore->mcPhi = genmom.Phi();
  _processstore->mcCtheta = genmom.CosTheta();
  _processstore->mcStkintersections = nstkintersections;

  //Check number of intersections with STK
  if( nstkintersections < minstkintersections)  { SetFilterResult(FilterResult::REJECT); }
    
	Herd::RefFrame::Direction entrydir = calotrack->entrancePlane;
	Herd::RefFrame::Direction exitdir = calotrack->exitPlane;
  _hcaloentryexitdir->Fill( entrydir==Herd::RefFrame::Direction::NONE ? _hcaloentryexitdir->GetNbinsX()-0.5 : static_cast<int>(entrydir), exitdir==Herd::RefFrame::Direction::NONE ? _hcaloentryexitdir->GetNbinsY()-0.5 : static_cast<int>(exitdir));
	
  //Check MC track entrance plane
  if( notfrombottom) {
    if(calotrack->entrancePlane == Herd::RefFrame::Direction::Zneg) SetFilterResult(FilterResult::REJECT);
    } 

  float calotracklengthx0=-1;
  if( !(entrydir==Herd::RefFrame::Direction::NONE && exitdir==Herd::RefFrame::Direction::NONE) ){
	  _gcaloentry->SetPoint(_gcaloentry->GetN(), calotrack->entrance[Herd::RefFrame::Coo::X],calotrack->entrance[Herd::RefFrame::Coo::Y],calotrack->entrance[Herd::RefFrame::Coo::Z]);
	  _gcaloexit->SetPoint(_gcaloexit->GetN(), calotrack->exit[Herd::RefFrame::Coo::X],calotrack->exit[Herd::RefFrame::Coo::Y],calotrack->exit[Herd::RefFrame::Coo::Z]);
	  _hshowerlength[static_cast<int>(entrydir)][static_cast<int>(exitdir)]->Fill(calotrack->trackLengthCaloX0);
    _hshowerlengthall->Fill(calotrack->trackLengthCaloX0);
    calotracklengthx0 = calotrack->trackLengthCaloX0;

    _processstore->mcTracklengthcalox0 = calotrack->trackLengthCaloX0;
    _processstore->mcTracklengthlysox0 = calotrack->trackLengthLYSOX0;

    _processstore->mcTrackcaloentry[0] = calotrack->entrance[Herd::RefFrame::Coo::X];
    _processstore->mcTrackcaloentry[1] = calotrack->entrance[Herd::RefFrame::Coo::Y];
    _processstore->mcTrackcaloentry[2] = calotrack->entrance[Herd::RefFrame::Coo::Z];
    _processstore->mcTrackcaloexit[0] = calotrack->exit[Herd::RefFrame::Coo::X];
    _processstore->mcTrackcaloexit[1] = calotrack->exit[Herd::RefFrame::Coo::Y];
    _processstore->mcTrackcaloexit[2] = calotrack->exit[Herd::RefFrame::Coo::Z];
    _processstore->mcTrackcaloentryplane = static_cast<int>(entrydir);
    _processstore->mcTrackcaloexitplane = static_cast<int>(exitdir);
	  }

  //Check MC track length
  if(calotrack->trackLengthCaloX0<mincalotrackx0) SetFilterResult(FilterResult::REJECT);

  return true;
}

bool MCtruthProcess::Finalize() {
  const std::string routineName("MCtruthProcess::Finalize");

  auto globStore = GetDataStoreManager()->GetGlobalDataStore("globStore");
  if (!globStore) {COUT(ERROR) << "Global data store not found." << ENDL;return false;}

  globStore->AddObject(_hgencoo->GetName(), _hgencoo);
  globStore->AddObject(_hgencthetaphi->GetName(),_hgencthetaphi);
  globStore->AddObject(_ggencoo->GetName(),_ggencoo);
  globStore->AddObject(_gcaloentry->GetName(),_gcaloentry);
  globStore->AddObject(_gcaloexit->GetName(),_gcaloexit);
  globStore->AddObject(_gdiscarded->GetName(),_gdiscarded);
  globStore->AddObject(_hstkintersections->GetName(),_hstkintersections);
  globStore->AddObject(_hshowerlengthall->GetName(), _hshowerlengthall);
  globStore->AddObject(_hcaloentryexitdir->GetName(),_hcaloentryexitdir);
  for(int indir=0; indir < Herd::RefFrame::NDirections; indir++){
    for(int outdir=0; outdir < Herd::RefFrame::NDirections; outdir++){
      globStore->AddObject(_hshowerlength[indir][outdir]->GetName(), _hshowerlength[indir][outdir]);
    }
  } 


  return true;
}

void MCtruthProcess::PrintCaloCubeMap(){
  const std::string routineName("MCtruthProcess::PrintCaloCubeMap");
  auto globStore = GetDataStoreManager()->GetGlobalDataStore("globStore");
  if (!globStore) {COUT(ERROR) << "Global data store not found." << ENDL;}
  auto caloGeoParams = globStore->GetObject<Herd::CaloGeoParams>("caloGeoParams");
  if (!caloGeoParams) {COUT(ERROR) << "Event data store not found." << ENDL;}
  else
  {
    for(unsigned int icube=0; icube<caloGeoParams->NCubes(); icube++){
    printf("[%04u]\t%.2f\t%.2f\t%.2f\t%.2f\n", icube, caloGeoParams->CubeSize(), caloGeoParams->Position(icube)[Herd::RefFrame::Coo::X],caloGeoParams->Position(icube)[Herd::RefFrame::Coo::Y],caloGeoParams->Position(icube)[Herd::RefFrame::Coo::Z]);
  }
  }

}


//***************************

MCtruthProcessStore::MCtruthProcessStore(const std::string &name) :
  Algorithm{name}
   {
  }

  bool MCtruthProcessStore::Initialize() {
  const std::string routineName("MCtruthProcessStore::Initialize");
  Reset();
  return true;
}

  bool MCtruthProcessStore::Process() {
  const std::string routineName("MCtruthProcessStore::Process");
  return true;
}
  bool MCtruthProcessStore::Finalize() {
  const std::string routineName("MCtruthProcessStore::Finalize");
  return true;
}
bool MCtruthProcessStore::Reset() {
  const std::string routineName("MCtruthProcessStore::Finalize");

  mcNdiscarded = -1;
  for(int idir=0; idir<3; idir++) mcDir[idir] = -999.;
  for(int idir=0; idir<3; idir++) mcCoo[idir] = -999.;
  mcMom = -999.;
  mcPhi = -999.;
  mcCtheta = -999.;
  mcStkintersections = -999.;
  mcTracklengthcalox0 = -999.;
  mcTracklengthlysox0 = -999.;

  for(int idir=0; idir<3; idir++) mcTrackcaloentry[idir] = -999.;
  for(int idir=0; idir<3; idir++) mcTrackcaloexit[idir] = -999.;
  mcTrackcaloentryplane = -1;
  mcTrackcaloexitplane = -1;

  return true;
}

