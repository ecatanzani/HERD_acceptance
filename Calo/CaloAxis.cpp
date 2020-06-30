// Example headers
#include "CaloAxis.h"
#include "dataobjects/CaloGeoParams.h"

// Root headers
#include "TMatrixD.h"
#include "TMatrixDSym.h"
#include "TMatrixDSymEigen.h"

// C/C++ standard headers
#include <numeric>
#include <cmath>

namespace Herd{
RegisterAlgorithm(CaloAxis);
RegisterAlgorithm(CaloAxisStore);


CaloAxis::CaloAxis(const std::string &name) :
  Algorithm{name},
  filterenable{true},
  process_clusters{true}
   {
    DefineParameter("filterenable",  filterenable); 
    DefineParameter("process_clusters", process_clusters);
    DefineParameter("edepthreshold", edepthreshold);
  }

bool CaloAxis::Initialize() {
  const std::string routineName("CaloAxis::Initialize");

  _evStore = GetDataStoreManager()->GetEventDataStore("evStore");      if (!_evStore)   { COUT(ERROR) << "Event data store not found." << ENDL; return false; }
  _globStore = GetDataStoreManager()->GetGlobalDataStore("globStore"); if (!_globStore) { COUT(ERROR) << "Global data store not found." << ENDL; return false; }

  _processstore = std::make_shared<CaloAxisStore>("CaloAxisStore");
  
  // Setup the filter                                                                                                                                                                                                                       
  if (filterenable) SetFilterStatus(FilterStatus::ENABLED); else SetFilterStatus(FilterStatus::DISABLED);

  hhitedep = std::make_shared<TH1F>("hhitedep", "Hit.Edep()", 1500, log10(1e-10),log10(1e+5));

  return true;
}

bool CaloAxis::Process() {
  const std::string routineName("CaloAxis::Process");

  //Add the ProcessStore object for this event to the event data store
  //_processstore->Reset();
  _evStore->AddObject("CaloAxisStore",_processstore);

  //Set Filter Status
  SetFilterResult(FilterResult::ACCEPT);

  //Set vector sizes
  caloaxisinfos.clear();
  //ShowerEigenvalues.clear();
  //ShowerEigenvectors.clear();

  if( process_clusters ){
    auto caloclusters = _evStore->GetObject<CaloClusters>("caloClusters");
    if (!_globStore) {COUT(ERROR) << "Global data store not found." << ENDL; return false;}
    BuildAxis( *caloclusters );
  }
  else{
    auto calohits = _evStore->GetObject<CaloHits>("caloHitsMC");
    BuildAxis( *calohits );
  }

  _processstore->caloaxishits   = (unsigned int)caloaxisinfos.at(0).ShowerHits;
  _processstore->caloaxiscog[0] = (float)caloaxisinfos.at(0).ShowerCOG[RefFrame::Coo::X];
  _processstore->caloaxiscog[1] = (float)caloaxisinfos.at(0).ShowerCOG[RefFrame::Coo::Y];
  _processstore->caloaxiscog[2] = (float)caloaxisinfos.at(0).ShowerCOG[RefFrame::Coo::Z];
  _processstore->caloaxisdir[0] = (float)caloaxisinfos.at(0).ShowerDir[RefFrame::Coo::X];
  _processstore->caloaxisdir[1] = (float)caloaxisinfos.at(0).ShowerDir[RefFrame::Coo::Y];
  _processstore->caloaxisdir[2] = (float)caloaxisinfos.at(0).ShowerDir[RefFrame::Coo::Z];
  for(int i=0; i<3; i++)
    {
     _processstore->caloaxiseigval[i]    = (float)caloaxisinfos.at(0).ShowerEigenvalues[i];
     _processstore->caloaxiseigvec[i][0] = (float)caloaxisinfos.at(0).ShowerEigenvectors[i][RefFrame::Coo::X];
     _processstore->caloaxiseigvec[i][1] = (float)caloaxisinfos.at(0).ShowerEigenvectors[i][RefFrame::Coo::Y];
     _processstore->caloaxiseigvec[i][2] = (float)caloaxisinfos.at(0).ShowerEigenvectors[i][RefFrame::Coo::Z];
    }
  return true;
}

bool CaloAxis::DummyCaloCluster(){

  return true;     
}

bool CaloAxis::BuildAxis(CaloClusters caloclusters){

  for( auto const& calohits: caloclusters){
    BuildAxis(calohits);
  }

  return true;
}

bool CaloAxis::BuildAxis(CaloHits calohits){
  const std::string routineName("CaloAxis::BuildAxis");

  auto caloGeoParams = _globStore->GetObject<CaloGeoParams>("caloGeoParams");

  std::vector<std::array<double,4>> x;
  CaloAxisInfo caloaxisinfo;


  //Loop on hits, select hits with edep>threshold and create a vector containing coordinates [0,1,2] and weigths [3]
  //x[0]: Coo::X
  //x[1]: Coo::Y
  //x[2]: Coo::Z
  //x[3]: Hit::EDep
  int n=0;
  double sumw=0; 
  for (auto &hit : calohits) {
    
    if(hit.EDep() > edepthreshold){
      //std::array<double,4> x_i;
      Point pos = caloGeoParams->Position(hit.VolumeID());
      //x_i[0] = (pos[RefFrame::Coo::X]);
      //x_i[1] = (pos[RefFrame::Coo::Y]);
      //x_i[2] = (pos[RefFrame::Coo::Z]);
      //x_i[3] = (hit.EDep());
      hhitedep->Fill( log10(hit.EDep()));

      //x.push_back(x_i);
      x.push_back(std::array<double,4>{pos[RefFrame::Coo::X],pos[RefFrame::Coo::Y],pos[RefFrame::Coo::Z],hit.EDep()});
      sumw += hit.EDep();
      n++;
    }
  }  
  //Here we assume that calo hits are uncorrelated, so the covariance matrix is diagonal

  //Calculate the centroid of the energy deposit
  std::array<double,3> cog = {0,0,0}; 
  for(int i=0; i<n; i++) { for(int j=0; j<3; j++) { cog[j] += (1/sumw) * x[i][3] * x[i][j]; } }// printf("%f %f\n", x[i][j], cog[j]);} 
    
  //Reposition the hit point to the energy centroid
  for(int i=0; i<n; i++) { for(int j=0; j<3; j++) x[i][j] -= cog[j]; } 

 /*for(int i=0; i<n; i=i+10){
    printf("-- [%d][%d]  %f\n",i,0,X[i][0]);
    printf("-- [%d][%d]  %f\n",i,1,X[i][1]);
    printf("-- [%d][%d]  %f\n",i,2,X[i][2]);
    }
  printf("\n");
*/
  
  /*
  TMatrixD X(n,3);    //Data matrix (with respect to centroid)
  for(int i=0; i<n; i++){ for(int j=0; j<3; j++){ X[i][j] = x[i][j]; }}
  TMatrixD Xt(3,n);   //Data matrix transp3osed (with respect to centroid)
  for(int j=0; j<3; j++){ for(int i=0; i<n; i++){ Xt[j][i] = x[i][j]; }}
  TMatrixD XtX(Xt,TMatrixD::kMult,X);  //Xt * X
  */

/*printf("%f %f\n", Xt[0][1],  X[1][0]);
printf("%f %f\n", Xt[1][20], X[20][1]);
printf("%f %f\n", Xt[2][40], X[40][2]);
*/
 
/* for(int i=0; i<3; i++){
    printf("-- %f %f %f\n",XtX[i][0],XtX[i][1],XtX[i][2]);
    }

  printf("\n");
*/
  //Calculate the covariance matrix as   C_i,j = (1/(sumw-1)) * sum_k [ (w_k) Xt_i,k * X_k,j ]
  //https://en.wikipedia.org/wiki/Sample_mean_and_covariance#Weighted_samples
 // https://en.wikipedia.org/wiki/Weighted_least_squares
 
 /* TMatrixDSym C(3);
  for(int i=0; i<3; i++){
    for(int j=0; j<3; j++){
      C[i][j] = 0; for(int k=0; k<n; k++) C[i][j] += (1/(sumw)) * x[k][3] * XtX[i][j];
    }
  }*/

  TMatrixDSym C(3);
  for(int i=0; i<3; i++){
    for(int j=0; j<3; j++){
      C[i][j] = 0; for(int k=0; k<n; k++) C[i][j] += (1/(sumw)) * x[k][3] * x[k][i] * x[k][j];
    }
  }
/*
  for(int i=0; i<3; i++){
    printf("%f\t%f\t%f\n",C[i][0],C[i][1],C[i][2]);
    }
  printf("\n");
  */

  //Get eigenvalues of the covariance matrix
  TMatrixDSymEigen E(C);
  TVectorD eigvaluesvec = E.GetEigenValues();
  TMatrixD eigvecmatrix = E.GetEigenVectors();

  //for(int i=0; i<3; i++) printf("[%d] %f\t {%f,%f,%f}\n", i, eigvaluesvec[i], eigvecmatrix[i][0], eigvecmatrix[i][1], eigvecmatrix[i][2]);
  //printf("\n");
  //Sort Eigenvector by decreasing eigenvalues
  TVectorD eigvec0( TVectorD(TMatrixTColumn_const<double>(eigvecmatrix, 0) ) ); //could be done more elegantly....
  TVectorD eigvec1( TVectorD(TMatrixTColumn_const<double>(eigvecmatrix, 1) ) );
  TVectorD eigvec2( TVectorD(TMatrixTColumn_const<double>(eigvecmatrix, 2) ) );
  double eigval[3];   for(int i=0; i<3; i++){ eigval[i] = eigvaluesvec[i]; }
  std::vector<std::pair<double,TVectorD>> eigvec = { {eigval[0],eigvec0}, {eigval[1],eigvec1}, {eigval[2],eigvec2} };
  std::sort(eigvec.begin(),eigvec.end(), []( const std::pair<double,TVectorD>x, std::pair<double,TVectorD>y) {return x.first>y.first;});

 /*  
  if( eigval[0]<1e-10){
    printf("n::%d\n",n);
    C.Print();
  }
  */

  //Store values to containers
  caloaxisinfo.ShowerHits = (unsigned int)n;

  caloaxisinfo.ShowerCOG[RefFrame::Coo::X] = cog[0];
  caloaxisinfo.ShowerCOG[RefFrame::Coo::Y] = cog[1];
  caloaxisinfo.ShowerCOG[RefFrame::Coo::Z] = cog[2];

  double mag=0; for(int i=0; i<3; i++) mag+= pow(eigvec.at(0).second[i],2); mag=sqrt(mag); 
  caloaxisinfo.ShowerDir[RefFrame::Coo::X] = eigvec.at(0).second[0]/mag;
  caloaxisinfo.ShowerDir[RefFrame::Coo::Y] = eigvec.at(0).second[1]/mag;
  caloaxisinfo.ShowerDir[RefFrame::Coo::Z] = eigvec.at(0).second[2]/mag;

  for(int i=0; i<3; i++) caloaxisinfo.ShowerEigenvalues.push_back( eigvec.at(i).first );
  for(int i=0; i<3; i++) { 
    Vec3D v(eigvec.at(i).second[0],eigvec.at(i).second[1],eigvec.at(i).second[2]);
    caloaxisinfo.ShowerEigenvectors.push_back(v); }
  
  caloaxisinfos.push_back(std::move(caloaxisinfo));
  x.clear();

  return true;
}


bool CaloAxis::Finalize() {
  const std::string routineName("CaloAxis::Finalize");

  auto globStore = GetDataStoreManager()->GetGlobalDataStore("globStore");
  if (!globStore) {COUT(ERROR) << "Global data store not found." << ENDL; return false;}

  globStore->AddObject(hhitedep->GetName(), hhitedep);
  return true;
}

//***************************

CaloAxisStore::CaloAxisStore(const std::string &name) :
  Algorithm{name}
   {
  }

  bool CaloAxisStore::Initialize() {
  const std::string routineName("CaloAxisStore::Initialize");
  Reset();
  return true;
}

  bool CaloAxisStore::Process() {
  const std::string routineName("CaloAxisStore::Process");
  return true;
}
  bool CaloAxisStore::Finalize() {
  const std::string routineName("CaloAxisStore::Finalize");
  return true;
}
bool CaloAxisStore::Reset() {
  const std::string routineName("CaloAxisStore::Finalize");

  caloaxishits=0;
  for(int i=0; i<3; i++){
    caloaxiscog[i] = -999.;
    caloaxisdir[i] = -999.;
    caloaxiseigval[i] = -999.;

    for(int j=0; j<3; j++)
      {
      caloaxiseigvec[i][j] = -999.;
      }
    }

  return true;
}

} //namespace Herd