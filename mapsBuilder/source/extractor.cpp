#include "extractor.h"

#include "TMath.h"
#include "TAxis.h"

#include <iostream>
#include <stdlib.h>

std::vector<std::vector<double>> get_events_coordinate(
    const unsigned int number_to_extract, 
    const double base_angle, 
    const std::shared_ptr<TH2D> evDist,
    const std::vector<double> pointing)
{
    std::vector<std::vector<double>> coord(number_to_extract, std::vector<double>(2, 0));
    double costheta, phi;
    for (unsigned int evt=0; evt<number_to_extract; ++evt)
    {
        evDist->GetRandom2(costheta, phi);
        //phi += TMath::Pi(); // unify respect to the base angle extracted from pointing information
        phi += base_angle;

        auto local_lat = TMath::ACos(costheta)*TMath::Sin(phi);
        auto local_lon = TMath::ACos(costheta)*TMath::Cos(phi);

        coord[evt][1] = pointing[1] + local_lat;
        coord[evt][0] = pointing[0] + local_lon;
        
        while (coord[evt][1]>TMath::Pi())
            coord[evt][1] -= TMath::Pi();

        while (coord[evt][1]<0)
            coord[evt][1] += TMath::Pi();

        while (coord[evt][0]>2*TMath::Pi())
            coord[evt][0] -= 2*TMath::Pi();

        while (coord[evt][0]<0)
            coord[evt][0] += 2*TMath::Pi();

        /*

        if (coord[evt][1]<0 || coord[evt][0]<0)
        {
            std::cout << "\n";
            std::cout << "\nTheta: " << TMath::ACos(costheta);
            std::cout << "\nPhi: " << myphi;
            std::cout << "\nBase angle: " << base_angle;
            std::cout << "\nLocal LAT: " << local_lat << "\tLocal LON: " << local_lon;
            std::cout << "\nLAT: " << pointing[1] << "\tLON: " << pointing[0];
            std::cout << "\nActual LAT: " << coord[evt][1] << "\tActual LON: " << coord[evt][0];
            std::cout << "\nLAT differeence: " << lat_difference;
            std::cout << "\nLON difference: " << lon_difference;
            std::cout << "\n";
            exit(123);
        }
        */
    }

    return coord;
}

void extract_from_distribution(
    const int nside,
    const std::vector<std::shared_ptr<TH2D>> h_event_distribution,
    const std::vector<double> old_pointing,
    const std::vector<double> pointing,
    const double geo_lat,
    TRandom3 &rgen,
    std::vector<std::vector<float>> &pixel_dataMap,
    const TF1 &fitFunc,
    const std::shared_ptr<TH1D> acceptance,
    const std::shared_ptr<TH1D> live_time)
{
    auto base_angle = TMath::ATan2((pointing[1]-old_pointing[1]),(pointing[0]-old_pointing[0]));
    //std::cout << "\nBase angle: " << base_angle;

    // Extract from energy bin map
    for (auto it=h_event_distribution.begin(); it!=h_event_distribution.end(); ++it)
    {
        // Compute the number of points to extract
        auto number_to_extract = GetNumberToExtract(
            fitFunc, 
            acceptance, 
            live_time,
            (int)std::distance(h_event_distribution.begin(), it),
            geo_lat,
            rgen);

        auto coord = get_events_coordinate(number_to_extract, base_angle, *it, pointing);
        for (unsigned int idx=0; idx<coord.size(); ++idx)
        {
            long hpix = 0;
            //std::cout << "\nTO HEALPIX: Actual LAT: " << coord[idx][1] << "\tActual LON: " << coord[idx][0];
            ang2pix_ring(nside,coord[idx][1],coord[idx][0],&hpix);   
            ++pixel_dataMap[std::distance(h_event_distribution.begin(), it)][hpix];
        }
    }
}

unsigned int GetNumberToExtract(
    TF1 fitFunc,
    const std::shared_ptr<TH1D> acceptance,
    const std::shared_ptr<TH1D> live_time,
    const int binIdx,
    const double geo_lat,
    TRandom3 &rgen)
{   
    std::unique_ptr<TH1D> h_acqRate (static_cast<TH1D*>(acceptance->Clone("h_acqRate")));
    auto multiplyStatus = h_acqRate->Multiply(&fitFunc);
    if (!multiplyStatus)
        std::cout << "\nERROR during TF1 multiplication" << std::endl;
    double energyWidth = acceptance->GetBinWidth(binIdx);
    double bin_rate = h_acqRate->GetBinContent(binIdx)*energyWidth;
    double live_time_value = live_time->GetBinContent(live_time->GetXaxis()->FindBin(geo_lat));  
    double mean_number_of_events = bin_rate * live_time_value;
    
    unsigned int events = rgen.Poisson(mean_number_of_events);

    /*
    std::cout << "\nBin rate: " << bin_rate;
    std::cout << "\nLivetime: " << live_time_value;
    std::cout << "\nNumber of points: " << events << std::endl;
    */
   
    return events;
}