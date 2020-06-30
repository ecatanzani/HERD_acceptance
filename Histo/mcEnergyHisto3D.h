#ifndef MCENERGYHISTO3D_H
#define MCENERGYHISTO3D_H

#include "algorithm/Algorithm.h"

// HerdSoftware headers
#include "dataobjects/MCTruth.h"

using namespace EA;

class TH3D;

class mcEnergyHisto3D : public Algorithm
{
public:
    mcEnergyHisto3D(const std::string &name);
    bool Initialize();
    bool Process();
    bool Finalize();

private:
    std::vector<double> energy_axispar;
    std::vector<double> polar_axispar;
    std::vector<double> azimut_axispar;

    std::vector<double> energy_binning;
    std::vector<double> polar_binning;
    std::vector<double> azimut_binning;
    bool logaxis;
    std::string title;

    void GenerateEnergyLogBinning();
    void GenerateEnergyBinning();
    void GeneratePolarBinning();
    void GenerateAzimutBinning();

    // Created global objects
    std::shared_ptr<TH3D> histo; // Objects to be pushed on global store must be held by a shared_ptr

    // Utility variables
    observer_ptr<EventDataStore> _evStore; // Pointer to the event data store
};

#endif