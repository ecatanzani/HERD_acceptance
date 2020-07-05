#ifndef EXTRACTOR_H
#define EXTRACTOR_H

#include <memory>
#include <vector>

#include "TH2D.h"
#include "TRandom3.h"

#include "healpix_base.h"
#include "chealpix.h"

extern void extract_from_distribution(
    const int nside,
    const std::vector<std::shared_ptr<TH2D>> h_event_distribution,
    const std::vector<double> old_pointing,
    const std::vector<double> pointing,
    const TRandom3 &rgen,
    std::vector<std::vector<float>> &pixel_dataMap);

#endif