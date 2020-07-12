#ifndef FLUXFIT_H
#define FLUXFIT_H

#include "TFile.h"
#include "TF1.h"
#include "TKey.h"
#include "TMultiGraph.h"
#include "TGraphAsymmErrors.h"
#include "TMath.h"

TF1 GetFluxFittedFunc(const std::string fluxPath);

#endif