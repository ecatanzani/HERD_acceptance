/*
 * CaloGeomFidVolumeAlgo.h
 *
 *  Created on: 21 Oct 2019
 *      Author: Lorenzo Pacini
 */

#ifndef HERD_CALOGEOMFIDVOLUMEALGO_H_
#define HERD_CALOGEOMFIDVOLUMEALGO_H_

#include "algorithm/Algorithm.h"
#include "common/DirectionsArray.h"
#include "dataobjects/Plane.h"
#include "dataobjects/CaloGeoParams.h"
#include "CaloGeomFidVolume.h"
#include <array>

using namespace EA;

class CaloGeomFidVolumeStore;


namespace Herd {

/*! @brief An algorithm which computes information about the track inside the Calo.
 * @class CaloGeomFidVolumeAlgo CaloGeomFidVolumeAlgo.h algorithms/geometry/CaloGeomFidVolumeAlgo.h
 *
 * <B>Needed event objects:</B>
 *
 *   name          |     type          |  store      | optional       | description
 * ----------------|-------------------|-------------|----------------|-------------------------
 * mcTruth         |    MCTruth        | evStore     |    no          | Info about MC truth
 *
 * <B>Produced event objects:</B>
 *
 *   name                       | type             |   store   | description
 * -----------------------------|------------------|-----------|----------------------------------------
 * trackInfoForCaloMC           | TrackInfoForCalo | evStore   | Container of information about the track for the Calo.
 *
 */
class CaloGeomFidVolumeAlgo : public Algorithm {
public:
  /*! @brief Constructor.
   *
   * @param name The name of the algorithm object.
   */
  CaloGeomFidVolumeAlgo(const std::string &name);

  /*! @brief Initializes the planes used to compute the track calo information.
   *
   * @return true if initialization is done without errors, false otherwise.
   */
  bool Initialize();

  /*! @brief Computes the track information, so far it use the MC truth.
   *
   * @return true if no error occurs during processing, false otherwise.
   */
  bool Process();

  /*! @brief Do nothing.
   *
   * @return true if no error occurs during finalization, false otherwise.
   */
  bool Finalize();

private:
  
  //Store pointer
  std::shared_ptr<CaloGeomFidVolumeStore> _processstore;

  observer_ptr<EventDataStore> _evStore; ///< Pointer to the event data store.
  //TrackInfoForCalo *_trackInfoCalo; ///< The TrackInfoForCalo object to fill with the computed information.
  DirectionsArray<Plane> plane; ///< The calo surface for each directions.
  DirectionsArray<float> m;         
  DirectionsArray<float> q;

  
  double mXNEGYNEG;
  double qXNEGYNEG;
  double mXPOSYNEG;
  double qXPOSYNEG;
  double mXNEGYPOS;
  double qXNEGYPOS;
  double mXPOSYPOS;
  double qXPOSYPOS;
  

  bool filterenable;
  bool checkext;
  bool checkint;

  const float _XSideBig;    // cm
  const float _XSideSmall;  // cm
  const float _YSideBig;    // cm
  const float _YSideSmall;  // cm
  const float _ZCaloCenter; // cm
  const float _ZCaloHeight; // cm
  const float _phiXY;
  float cubeside; //cm
  float alpha;    //fraction of cube size to be contained in the fiducuial volume
  float shrink;   //cubeside*alpha

  std::array<Point,8> pxy;
  std::array<Point,2> pz;

  bool CheckExt();
  bool CheckInt();
  void FillCoo(const Herd::Point p, float coo[2][3], int index);

 
};                                    

} // namespace Herd

class CaloGeomFidVolumeStore : public Algorithm {
public:
  CaloGeomFidVolumeStore(const std::string &name);
  bool Initialize();
  bool Process();
  bool Finalize();
  bool Reset();

  float calofidvolalpha;
  bool calofidvolpass;
  short calofidvolxpos;
  short calofidvolxneg;
  short calofidvolypos;
  short calofidvolyneg;
  short calofidvolzpos;
  short calofidvolzneg;
  short calofidvolxnegyneg;
  short calofidvolxposyneg;
  short calofidvolxnegypos;
  short calofidvolxposypos;
  float calofidvolxposEntry[2][3];
  float calofidvolxnegEntry[2][3];
  float calofidvolyposEntry[2][3];
  float calofidvolynegEntry[2][3];
  float calofidvolzposEntry[2][3];
  float calofidvolznegEntry[2][3];
  float calofidvolxnegynegEntry[2][3];
  float calofidvolxposynegEntry[2][3];
  float calofidvolxnegyposEntry[2][3];
  float calofidvolxposyposEntry[2][3];

  
private:
};

#endif /* HERD_CALOGEOMFIDVOLUMEALGO_H_ */
