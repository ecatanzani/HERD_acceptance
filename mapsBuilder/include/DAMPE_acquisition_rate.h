#ifndef DAMPE_ACQUISITION_RATE_H
#define DAMPE_ACQUISITION_RATE_H

#include <iostream>
#include <memory>

#include "TFile.h"
#include "TH2D.h"

extern std::shared_ptr<TH2D> GetDAMPEAcquisitionRate(const std::string acqPath);

#endif