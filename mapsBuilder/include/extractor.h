#ifndef EXTRACTOR_H
#define EXTRACTOR_H

#include <memory>
#include <vector>

#include "TH1D.h"
#include "TF1.h"
#include "TH2D.h"
#include "TRandom3.h"

#include "healpix_base.h"
#include "chealpix.h"

extern void extract_from_distribution(
    const int nside,
    const std::vector<std::shared_ptr<TH2D>> h_event_distribution,
    const std::vector<double> old_pointing,
    const std::vector<double> pointing,
    const double geo_lat,
    TRandom3 &rgen,
    std::vector<std::vector<float>> &pixel_dataMap,
    const TF1 &fitFunc,
    const std::shared_ptr<TH1D> acceptance,
    const std::shared_ptr<TH1D> live_time);

extern std::vector<std::vector<double>> get_events_coordinate(
    const unsigned int number_to_extract, 
    const double base_angle, 
    const std::shared_ptr<TH2D> evDist,
    const std::vector<double> pointing);

extern unsigned int GetNumberToExtract(
    TF1 fitFunc,
    const std::shared_ptr<TH1D> acceptance,
    const std::shared_ptr<TH1D> live_time,
    const int binIdx,
    const double geo_lat,
    TRandom3 &rgen);

#endif