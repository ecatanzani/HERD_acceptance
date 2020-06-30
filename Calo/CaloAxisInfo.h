/*
 * CaloAxisInfo.h
 *
 *  Created on: 27 April 2020
 *      Author: Valerio Vagelli
 */

/*! @file CaloAxisInfo.h CaloAxisInfo class declaration, */

#ifndef HERD_CALOAXISINFO_H_
#define HERD_CALOAXISINFO_H_

// HerdSoftware headersCalo
#include "common/Vec3D.h"
#include "dataobjects/Point.h"

// C/C++ standard headers
#include <numeric>
#include <vector>

namespace Herd {

/*! @brief Container of information about the track inside the Calo
 * @struct CaloAxisInfo CaloAxisInfo.h dataobjects/CaloAxisInfo.h
 *
 */
struct CaloAxisInfo {
  
  unsigned short ShowerHits;                           ///Number of hits used for shower axis reconstruction
  Point ShowerCOG;                          ///Shower Center of Gravity
  Point ShowerDir;                          ///Shower axis Directions
  std::vector<double> ShowerEigenvalues;    ///Eigenvalues of Shower Covariance matrix
  std::vector<Vec3D> ShowerEigenvectors;    ///Eignevectors of Shower Covariance matrix

  /*! @brief Default constructor.
   *
   * Creates a CaloAxisInfo with standard values
   *
   */
  CaloAxisInfo() { Reset(); };

  /*!
   * @brief Set the members to default values
   *
   */
  void Reset();
};

using CaloAxisInfos = std::vector<CaloAxisInfo>;

} // namespace Herd

#endif /* HERD_CALOAXISINFO_H_ */
