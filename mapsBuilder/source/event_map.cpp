#include "main.h"
#include "acceptance_simu.h"
#include "extractor.h"
#include "unique.h"

#include "healpix_base.h"
#include "chealpix.h"

#include <memory>

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TRandom3.h"
#include "TMath.h"

#define nside 32

void evaluateEventMap(
    const std::string evtPath,
    const std::string telPath,
    const std::string outPath,
    const bool verbose,
    AnyOption &opt,
    const TF1 &fitFunc,
    const std::shared_ptr<TH1D> acceptance,
    const std::shared_ptr<TH1D> live_time)
{   
    UInt_t rand_seed = 2;
    TRandom3 rgen(rand_seed);
    auto npix = nside2npix(nside); 

    TFile ang_event_file(evtPath.c_str(), "READ");
    if (ang_event_file.IsZombie())
    {
        std::cerr << "\n\nError reading input acceptance file: " << evtPath << std::endl;
        exit(100);
    }
    
    std::vector<std::shared_ptr<TH2D>> h_event_distribution (acceptance->GetNbinsX());
    
    std::vector<std::vector<float>> pixel_dataMap (acceptance->GetNbinsX());
    float iss_pointing[npix];
    for (auto idx=0; idx < npix; ++idx)
        iss_pointing[idx] = 0;

    init_data_maps(pixel_dataMap, npix);

    for (auto it=h_event_distribution.begin(); it!=h_event_distribution.end(); ++it)
    {   
        std::string histo_tmp_name = "h_angularDistribution_energyBin_" + std::to_string(std::distance(h_event_distribution.begin(), it));
        *it = std::shared_ptr<TH2D>(static_cast<TH2D*>(ang_event_file.Get(histo_tmp_name.c_str())));
        (*it)->SetDirectory(0);
    }
    
    ang_event_file.Close();

    TFile tree_file(telPath.c_str(), "READ");
    if (tree_file.IsZombie())
    {
        std::cerr << "\n\nError reading input telemetry file: " << telPath << std::endl;
        exit(100);
    }

    std::shared_ptr<TTree> RTItree(static_cast<TTree*>(tree_file.Get("RTI_tree")));
    std::vector<double> pointing(2, 0);
    double geo_lat = 0;
    RTItree->SetBranchAddress("glat", &pointing[1]);
    RTItree->SetBranchAddress("glon", &pointing[0]);
    RTItree->SetBranchAddress("geo_lat", &geo_lat);
    
    RTItree->GetEntry(0);

    std::vector<double> old_pointing(2);
    
    for (auto idx=0; idx<RTItree->GetEntries(); ++idx)
    {
        if ( (idx+1)%1000==0)
            std::cout << "\n[Event Loop]: " << idx+1;
        
        RTItree->GetEntry(idx);
         
        // Use healpix convention for pointing
        pointing[1] += 90;
        pointing[1] *= TMath::DegToRad();

        pointing[0] += 180;
        pointing[0] *= TMath::DegToRad();

        // Fill pointing map
        long hpix;
        ang2pix_ring(nside,pointing[1],pointing[0],&hpix);
        ++iss_pointing[hpix];

        if (!idx)
        {
            old_pointing = pointing;
            continue;
        }
        
        extract_from_distribution(
            nside,
            h_event_distribution,
            old_pointing,
            pointing,
            geo_lat,
            rgen,
            pixel_dataMap,
            fitFunc,
            acceptance,
            live_time);
    }

    auto mapsPath = uniqueOutFile(
        outPath, 
        opt, 
        0,
        true);

    write_healpix_map(
        iss_pointing,
        nside,
        mapsPath.c_str(),
        0,
        "G");
    
    write_final_maps(
        pixel_dataMap, 
        outPath, 
        opt);
    
}

void init_data_maps(std::vector<std::vector<float>> &pixel_dataMap, const long npix)
{
    for (unsigned int idx=0; idx<pixel_dataMap.size(); ++idx)
    {
        pixel_dataMap[idx].resize(npix);
        for (unsigned int sIdx=0; sIdx<pixel_dataMap[idx].size(); ++sIdx)
            pixel_dataMap[idx][sIdx] = 0;
    }
}

void write_final_maps(
    std::vector<std::vector<float>> &pixel_dataMap, 
    const std::string outputPath, 
    AnyOption &opt)
{   
    std::vector<std::vector<float>> pixel_dataMap_sum;

    // Sum maps
    for (unsigned int idx_map=0; idx_map<pixel_dataMap.size(); ++idx_map)
    {
        std::vector<float> tmp_sum_map (pixel_dataMap[idx_map].size(), 0);
        for (unsigned int idx_map_to_sum = idx_map; idx_map_to_sum<pixel_dataMap.size(); ++idx_map_to_sum)
            for(unsigned int single_map_idx=0; single_map_idx<pixel_dataMap[idx_map_to_sum].size(); ++single_map_idx)
                tmp_sum_map[single_map_idx] += pixel_dataMap[idx_map_to_sum][single_map_idx];
        pixel_dataMap_sum.push_back(tmp_sum_map);
    }
    
    float fmaps[pixel_dataMap_sum.size()][nside2npix(nside)];
    for(auto it=pixel_dataMap_sum.begin(); it!=pixel_dataMap_sum.end(); ++it)
    {
        auto mapsPath = uniqueOutFile(
            outputPath, 
            opt, 
            std::distance(pixel_dataMap_sum.begin(), it));

        for (unsigned int idx=0; idx<(*it).size(); ++idx)
            fmaps[std::distance(pixel_dataMap_sum.begin(), it)][idx] = (*it)[idx];

        write_healpix_map(
            fmaps[std::distance(pixel_dataMap_sum.begin(), it)],
            nside,
            mapsPath.c_str(),
            0,
            "G");
    }
}