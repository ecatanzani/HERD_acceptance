#include "extractor.h"

#include "TMath.h"

#include <iostream>

std::vector<std::vector<double>> get_events_coordinate(
    const int number_to_extract, 
    const double base_angle, 
    const std::shared_ptr<TH2D> evDist,
    const std::vector<double> pointing)
{
    std::vector<std::vector<double>> coord(number_to_extract, std::vector<double>(2, 0));
    double costheta, phi;
    for (auto evt=0; evt<number_to_extract; ++evt)
    {
        evDist->GetRandom2(costheta, phi);
        phi += base_angle;
        coord[evt][0] = pointing[0] + TMath::ACos(costheta)*TMath::Cos(phi);
        coord[evt][1] = pointing[1] + TMath::ACos(costheta)*TMath::Sin(phi);
    }

    return coord;
}

void extract_from_distribution(
    const std::vector<std::shared_ptr<TH2D>> h_event_distribution,
    const std::vector<double> old_pointing,
    const std::vector<double> pointing,
    const TRandom3 &rgen)
{
    auto base_angle = TMath::ATan((pointing[1]-old_pointing[1])/(pointing[0]-old_pointing[0]))*TMath::RadToDeg();
    
    // Extract from energy bin map
    for (auto it=h_event_distribution.begin(); it!=h_event_distribution.end(); ++it)
    {
        auto number_to_extract = 10;    // To be calculated including the flux and the livetime simulation
        auto coord = get_events_coordinate(number_to_extract, base_angle, *it, pointing);
    }
}