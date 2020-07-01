// Example headers
#include "mcEnergyHisto.h"

// Root headers
#include "TH1D.h"

// C/C++ standard headers
#include <numeric>
#include <cmath>
#include <cstdlib>
#include <string>

RegisterAlgorithm(mcEnergyHisto);

mcEnergyHisto::mcEnergyHisto(const std::string &name) : Algorithm{name},
												   axispar{100., 0., 100.},
												   logaxis{false},
												   title("title")
{
	DefineParameter("axispar", axispar);
	DefineParameter("logaxis", logaxis);
	DefineParameter("title", title);
}

bool mcEnergyHisto::Initialize()
{
	const std::string routineName("mcEnergyHisto::Initialize");

	_evStore = GetDataStoreManager()->GetEventDataStore("evStore");
	if (!_evStore)
	{
		COUT(ERROR) << "Event data store not found." << ENDL;
		return false;
	}

	// Check the user setting for the axis
	if (axispar.size() != 3)
	{
		COUT(ERROR) << "The axis must be specified by exactly 3 parameters" << ENDL;
		return false;
	}
	if ((float)((int)(axispar[0])) != axispar[0])
	{
		COUT(ERROR) << "The number of bins is not an integer." << ENDL;
		return false;
	}
	if (axispar[0] < 0)
	{
		COUT(ERROR) << "The number of bins is a negative value." << ENDL;
		return false;
	}
	if (axispar[1] >= axispar[2])
	{
		COUT(ERROR) << "The lower axis limit is greater or equal to the upper limit." << ENDL;
		return false;
	}
	
	if (logaxis)
		GenerateLogBinning();
	else
		GenerateBinning();

	// Create the histogram
	std::string histo_name = "h_" + GetName();
	histo = std::make_shared<TH1D>(histo_name.c_str(), title.c_str(), binning.size()-1, &(binning[0]));
	histo->GetXaxis()->SetTitle("MC Momentum (GV()");

	return true;
}

bool mcEnergyHisto::Process()
{
	const std::string routineName("mcEnergyHisto::Process");
	
	auto mctruth = _evStore->GetObject<Herd::MCTruth>("mcTruth");
	if (!mctruth)
	{
		COUT(DEBUG) << "MCTruth not present for event " << GetEventLoopProxy()->GetCurrentEvent() << ENDL;
		return false;
	}
	auto mcmom = std::sqrt(mctruth->primaries.at(0).initialMomentum * mctruth->primaries.at(0).initialMomentum);
	histo->Fill(mcmom);

	return true;
}

bool mcEnergyHisto::Finalize()
{
	const std::string routineName("mcEnergyHisto::Finalize");
	auto globStore = GetDataStoreManager()->GetGlobalDataStore("globStore");
	if (!globStore)
	{
		COUT(ERROR) << "Global data store not found." << ENDL;
		return false;
	}

	globStore->AddObject(histo->GetName(), histo);
	
	return true;
}

void mcEnergyHisto::GenerateLogBinning()
{
	binning.resize((int)axispar[0]+1);
	double log_interval = (log10(axispar[2]) - log10(axispar[1])) / (int)axispar[0];
	for (auto bIdx = 0; bIdx <= (int)axispar[0]; ++bIdx)
		binning[bIdx] = pow(10, log10(axispar[1]) + bIdx * log_interval);
}

void mcEnergyHisto::GenerateBinning()
{
	binning.resize((int)axispar[0]+1);
	double interval = (axispar[2] - axispar[1]) / (int)axispar[0];
	for (auto bIdx = 0; bIdx <= (int)axispar[0]; ++bIdx)
		binning[bIdx] = axispar[1] + bIdx * interval;
}
