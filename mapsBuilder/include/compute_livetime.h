#ifndef COMPUTE_LIVETIME_H
#define COMPUTE_LIVETIME_H

#include <memory>
#include <iostream>
#include <math.h>

#include "TH2D.h"
#include "TH1D.h"
#include "TFile.h"

#define deadtime 3e-3

std::shared_ptr<TH1D> GetLiveTime(const std::shared_ptr<TH2D> dampe_acquition_rate);

#endif