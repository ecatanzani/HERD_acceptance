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
    const std::string accPath,
    const std::string telPath,
    const std::string outPath,
    const bool verbose,
    AnyOption &opt)
{   
    UInt_t rand_seed = 2;
    TRandom3 rgen(rand_seed);
    auto npix = nside2npix(nside); 

    TFile acceptance_file(accPath.c_str(), "READ");
    if (acceptance_file.IsZombie())
    {
        std::cerr << "\n\nError reading input acceptance file: " << accPath << std::endl;
        exit(100);
    }

    std::shared_ptr<TH1D> h_calo_filtered_fidvolume(static_cast<TH1D*>(acceptance_file.Get("h_calo_filtered_fidvolume")));
    std::shared_ptr<TH1D> h_mcgenspectrum(static_cast<TH1D*>(acceptance_file.Get("h_mcgenspectrum")));
    h_calo_filtered_fidvolume->SetDirectory(0);
    h_mcgenspectrum->SetDirectory(0);
    std::vector<std::shared_ptr<TH2D>> h_event_distribution (h_calo_filtered_fidvolume->GetNbinsX());
    
    std::vector<std::vector<float>> pixel_dataMap (h_calo_filtered_fidvolume->GetNbinsX());
    float iss_pointing[npix];
    for (auto idx=0; idx < npix; ++idx)
        iss_pointing[idx] = 0;

    init_data_maps(pixel_dataMap, npix);

    for (auto it=h_event_distribution.begin(); it!=h_event_distribution.end(); ++it)
    {   
        std::string histo_tmp_name = "h_angularDistribution_energyBin_" + std::to_string(std::distance(h_event_distribution.begin(), it));
        *it = std::shared_ptr<TH2D>(static_cast<TH2D*>(acceptance_file.Get(histo_tmp_name.c_str())));
        (*it)->SetDirectory(0);
    }
    
    acceptance_file.Close();

    TFile tree_file(telPath.c_str(), "READ");
    if (tree_file.IsZombie())
    {
        std::cerr << "\n\nError reading input telemetry file: " << telPath << std::endl;
        exit(100);
    }

    std::shared_ptr<TTree> RTItree(static_cast<TTree*>(tree_file.Get("RTI_tree")));
    std::vector<double> pointing(2, 0);
    RTItree->SetBranchAddress("glat", &pointing[1]);
    RTItree->SetBranchAddress("glon", &pointing[0]);
    
    RTItree->GetEntry(0);

    std::vector<double> old_pointing(2);
    
    for (auto idx=0; idx<RTItree->GetEntries(); ++idx)
    {
        if ( (idx+1)%1000==0)
            std::cout << "\n[Event Loop]: " << idx+1;
        
        RTItree->GetEntry(idx);
        
        //std::cout <<"\nLat before: " << pointing[1] << std::endl;

         // Use healpix convention for pointing
        pointing[1] += 90;
        pointing[1] *= TMath::DegToRad();

        pointing[0] += 180;
        pointing[0] *= TMath::DegToRad();
        
        //std::cout << "\nTheta: " << pointing[1] << "\tPhi: " << pointing[0] << std::endl;

        // Fill pointing map
        long hpix;
        ang2pix_ring(nside,pointing[1],pointing[0],&hpix);
        ++iss_pointing[hpix];

        if (!idx)
        {
            old_pointing = pointing;
            continue;
        }

        /*
        extract_from_distribution(
            nside,
            h_event_distribution,
            old_pointing,
            pointing,
            rgen,
            pixel_dataMap);
        */
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

    /*
    write_final_maps(
        pixel_dataMap, 
        outPath, 
        opt);
    */
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
    float fmaps[pixel_dataMap.size()][nside2npix(nside)];
    for(auto it=pixel_dataMap.begin(); it!=pixel_dataMap.end(); ++it)
    {
        auto mapsPath = uniqueOutFile(
            outputPath, 
            opt, 
            std::distance(pixel_dataMap.begin(), it), false);

        for (unsigned int idx=0; idx<(*it).size(); ++idx)
            fmaps[std::distance(pixel_dataMap.begin(), it)][idx] = (*it)[idx];

        write_healpix_map(
            fmaps[std::distance(pixel_dataMap.begin(), it)],
            nside,
            mapsPath.c_str(),
            0,
            "G");
    }
}