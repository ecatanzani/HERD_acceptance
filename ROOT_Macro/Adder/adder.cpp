#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"

#include <vector>
#include <memory>
#include <iostream>

void addAngularDistributions(const char* LE_dataFile, const char* HE_dataFile)
{
    TFile myInFile_LE(LE_dataFile, "READ");
	if (myInFile_LE.IsZombie())
	{
		std::cerr << "\n\nError reading input ROOT file: " << LE_dataFile << std::endl;
		exit(123);
	}

    auto numberOfBins = static_cast<TH1D*>(myInFile_LE.Get("h_calo_filtered_fidvolume"))->GetNbinsX();

    std::vector<std::shared_ptr<TH2D>> h_LE_event_distribution (numberOfBins);

    for (auto it=h_LE_event_distribution.begin(); it!=h_LE_event_distribution.end(); ++it)
    {   
        std::string histo_tmp_name = "h_angularDistribution_energyBin_" + std::to_string(std::distance(h_LE_event_distribution.begin(), it));
        *it = std::shared_ptr<TH2D>(static_cast<TH2D*>(myInFile_LE.Get(histo_tmp_name.c_str())));
        (*it)->SetDirectory(0);
        (*it)->Sumw2();
    }

    myInFile_LE.Close();

    TFile myInFile_HE(HE_dataFile, "READ");
	if (myInFile_HE.IsZombie())
	{
		std::cerr << "\n\nError reading input ROOT file: " << HE_dataFile << std::endl;
		exit(123);
	}

    std::vector<std::shared_ptr<TH2D>> h_HE_event_distribution (numberOfBins);

    for (auto it=h_HE_event_distribution.begin(); it!=h_HE_event_distribution.end(); ++it)
    {   
        std::string histo_tmp_name = "h_angularDistribution_energyBin_" + std::to_string(std::distance(h_HE_event_distribution.begin(), it));
        *it = std::shared_ptr<TH2D>(static_cast<TH2D*>(myInFile_HE.Get(histo_tmp_name.c_str())));
        (*it)->SetDirectory(0);
        (*it)->Sumw2();
    }

    myInFile_HE.Close();

    for (auto it=h_LE_event_distribution.begin(); it!=h_LE_event_distribution.end(); ++it)
        (*it)->Add(h_HE_event_distribution[std::distance(h_LE_event_distribution.begin(), it)].get());
    
    TFile outFile("eventAngularDistribution.root", "RECREATE");
    if (outFile.IsZombie())
    {
        std::cerr << "\n\nError writing output TFile" << std::endl;
        exit(123);
    }

    for (auto it=h_LE_event_distribution.begin(); it!=h_LE_event_distribution.end(); ++it)
        (*it)->Write();    

    outFile.Close();
}