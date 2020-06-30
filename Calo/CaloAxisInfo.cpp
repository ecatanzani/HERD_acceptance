/*
 * CaloAxisInfo.cpp
 *
 *  Created on: 27 April 2020
 *      Author: Valerio Vagelli
 */

#include "CaloAxisInfo.h"

namespace Herd {

void CaloAxisInfo::Reset() {
  
  ShowerHits = 0;
  ShowerCOG = Point();        
  ShowerDir = Point();            
  ShowerEigenvalues.clear();  
  ShowerEigenvectors.clear(); 

}
    
} // namespace Herd
