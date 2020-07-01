#include "mcAngleDistribution.h"

// Root headers
#include "TH2D.h"

// HerdSoftware headers
#include "dataobjects/Line.h"
#include "dataobjects/MCTruth.h"
#include "dataobjects/Momentum.h"
#include "dataobjects/Point.h"

// C/C++ standard headers
#include <numeric>
#include <cmath>
#include <cstdlib>
#include <string>

namespace Herd
{

	RegisterAlgorithm(mcAngleDistribution);

	mcAngleDistribution::mcAngleDistribution(const std::string &name) : Algorithm{name},
																		energy_axispar{100, 1e+1, 1e+4},
																		polar_axispar{100, -1, 1},
																		azimuth_axispar{180, -M_PI, M_PI},
																		logaxis{false},
																		title("title")
	{

		DeclareConsumedObject("mcTruth", ObjectCategory::EVENT, "evStore");

		DefineParameter("energy_axispar", energy_axispar);
		DefineParameter("polar_axispar", polar_axispar);
		DefineParameter("azimuth_axispar", azimuth_axispar);
		DefineParameter("logaxis", logaxis);
		DefineParameter("title", title);
	}

	bool mcAngleDistribution::Initialize()
	{
		const std::string routineName("mcAngleDistribution::Initialize");

		_evStore = GetDataStoreManager()->GetEventDataStore("evStore");
		if (!_evStore)
		{
			COUT(ERROR) << "Event data store not found." << ENDL;
			return false;
		}

		// Check the user setting for the axis
		if (energy_axispar.size() != 3)
		{
			COUT(ERROR) << "The energy axis must be specified by exactly 3 parameters" << ENDL;
			return false;
		}
		if (polar_axispar.size() != 3)
		{
			COUT(ERROR) << "The polar axis must be specified by exactly 3 parameters" << ENDL;
			return false;
		}
		if (azimuth_axispar.size() != 3)
		{
			COUT(ERROR) << "The azimut axis must be specified by exactly 3 parameters" << ENDL;
			return false;
		}

		if ((float)((int)(energy_axispar[0])) != energy_axispar[0])
		{
			COUT(ERROR) << "The number of energy bins is not an integer." << ENDL;
			return false;
		}
		if ((float)((int)(polar_axispar[0])) != polar_axispar[0])
		{
			COUT(ERROR) << "The number of polar bins is not an integer." << ENDL;
			return false;
		}
		if ((float)((int)(azimuth_axispar[0])) != azimuth_axispar[0])
		{
			COUT(ERROR) << "The number of azimut bins is not an integer." << ENDL;
			return false;
		}

		if (energy_axispar[0] < 0)
		{
			COUT(ERROR) << "The number of energy bins is a negative value." << ENDL;
			return false;
		}
		if (polar_axispar[0] < 0)
		{
			COUT(ERROR) << "The number of polar bins is a negative value." << ENDL;
			return false;
		}
		if (azimuth_axispar[0] < 0)
		{
			COUT(ERROR) << "The number of azimut bins is a negative value." << ENDL;
			return false;
		}

		if (energy_axispar[1] >= energy_axispar[2])
		{
			COUT(ERROR) << "The lower axis limit is greater or equal to the upper limit." << ENDL;
			return false;
		}
		if (polar_axispar[1] >= polar_axispar[2])
		{
			COUT(ERROR) << "The lower axis limit is greater or equal to the upper limit." << ENDL;
			return false;
		}
		if (azimuth_axispar[1] >= azimuth_axispar[2])
		{
			COUT(ERROR) << "The lower axis limit is greater or equal to the upper limit." << ENDL;
			return false;
		}

		if (logaxis)
			GenerateLogEnergyBinning();
		else
			GenerateEnergyBinning();

		// Define the number of TH2 histos as the number of energy bins
		histo.resize(energy_axispar[0]);
		
		// Create the histos
		for (auto ebIdx=0; ebIdx< histo.size(); ++ebIdx)
		{
			std::string histo_name = "h_" + GetName() + "_energyBin_" + std::to_string(ebIdx);
			histo[ebIdx] = std::make_shared<TH2D> (
				histo_name.c_str(), 
				title.c_str(), 
				polar_axispar[0], polar_axispar[1], polar_axispar[2],
				azimuth_axispar[0], azimuth_axispar[1], azimuth_axispar[2]);
		}

		return true;
	}

	bool mcAngleDistribution::Process()
	{
		const std::string routineName("mcAngleDistribution::Process");

		auto mcTruth = _evStore->GetObject<Herd::MCTruth>("mcTruth");
		if (!mcTruth)
		{
			COUT(DEBUG) << "MCTruth not present for event " << GetEventLoopProxy()->GetCurrentEvent() << ENDL;
			return false;
		}

		Momentum &Mom = mcTruth->primaries[0].initialMomentum;
		Point &Pos = mcTruth->primaries[0].initialPosition;
		Line MCtrack(Pos, Mom);

		auto theta_deg = MCtrack.Polar() * 180 / M_PI;
		auto costheta = cos(MCtrack.Polar());
		auto phi_deg = MCtrack.Azimuth() * 180 / M_PI;
		auto phi = MCtrack.Azimuth();
		auto mcmom = std::sqrt(mcTruth->primaries[0].initialMomentum * mcTruth->primaries[0].initialMomentum);
		auto ebIdx = getCurrentEnergyBin(mcmom);
		
		histo[ebIdx]->Fill(costheta, phi);

		return true;
	}

	bool mcAngleDistribution::Finalize()
	{
		const std::string routineName("mcAngleDistribution::Finalize");
		auto globStore = GetDataStoreManager()->GetGlobalDataStore("globStore");
		if (!globStore)
		{
			COUT(ERROR) << "Global data store not found." << ENDL;
			return false;
		}
		
		// Book histos to the global store
		for (auto hIdx=0; hIdx<histo.size(); ++hIdx)
			globStore->AddObject(histo[hIdx]->GetName(), histo[hIdx]);

		return true;
	}

	void mcAngleDistribution::GenerateLogEnergyBinning()
	{
		energy_binning.resize((int)energy_axispar[0] + 1);
		double log_interval = (log10(energy_axispar[2]) - log10(energy_axispar[1])) / (int)energy_axispar[0];
		for (auto bIdx = 0; bIdx <= (int)energy_axispar[0]; ++bIdx)
			energy_binning[bIdx] = pow(10, log10(energy_axispar[1]) + bIdx * log_interval);
	}

	void mcAngleDistribution::GenerateEnergyBinning()
	{
		energy_binning.resize((int)energy_axispar[0] + 1);
		double interval = (energy_axispar[2] - energy_axispar[1]) / (int)energy_axispar[0];
		for (auto bIdx = 0; bIdx <= (int)energy_axispar[0]; ++bIdx)
			energy_binning[bIdx] = energy_axispar[1] + bIdx * interval;
	}
	
	int mcAngleDistribution::getCurrentEnergyBin(double energy)
	{
		for (auto ebIdx=0; ebIdx<energy_binning.size()-1; ++ebIdx)
			if (energy_binning[ebIdx] <= energy && energy_binning[ebIdx+1] >= energy)
				return ebIdx;
	}

} // namespace Herd