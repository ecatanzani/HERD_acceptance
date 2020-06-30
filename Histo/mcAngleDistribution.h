#ifndef mcAngleDistribution_H_
#define mcAngleDistribution_H_

#include "algorithm/Algorithm.h"

// HerdSoftware headers
#include "dataobjects/MCTruth.h"

using namespace EA;

class TH2D;

namespace Herd
{

    class mcAngleDistribution : public Algorithm
    {
    public:
        mcAngleDistribution(const std::string &name);
        bool Initialize();
        bool Process();
        bool Finalize();

    private:
        std::vector<double> energy_axispar;
        std::vector<double> polar_axispar;
        std::vector<double> azimuth_axispar;
        std::vector<double> energy_binning;
        bool logaxis;
        std::string title;

        void GenerateLogEnergyBinning();
        void GenerateEnergyBinning();
        int getCurrentEnergyBin(double energy);

        // Created global objects
        std::vector<std::shared_ptr<TH2D>> histo;

        // Utility variables
        observer_ptr<EventDataStore> _evStore; // Pointer to the event data store
    };

} // namespace Herd

#endif