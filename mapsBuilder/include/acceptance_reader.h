#ifndef ACCEPTANCE_READER_H
#define ACCEPTANCE_READER_H

#include "TFile.h"
#include "TH1D.h"

#include <memory>
#include <iostream>

std::shared_ptr<TH1D> GetAcceptanceHisto(const std::string accPath);


#endif