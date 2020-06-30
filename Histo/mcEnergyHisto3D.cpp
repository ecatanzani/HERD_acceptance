#include "mcEnergyHisto.h"

// Root headers
#include "TH3D.h"

// HerdSoftware headers
#include "dataobjects/Line.h"

// C/C++ standard headers
#include <numeric>
#include <cmath>
#include <cstdlib>
#include <string>

RegisterAlgorithm(mcEnergyHisto3D);

mcEnergyHisto3D::mcEnergyHisto3D(const std::string &name) : Algorithm{name},
                                                            energy_axispar{100., 0., 100.},
                                                            polar_axispar{100., 0., 180.},
                                                            aximut_axispar{180., 0., 360.},
                                                            logaxis{false},
                                                            title("title")
{
    DefineParameter("energy_axispar", energy_axispar);
    DefineParameter("polar_axispar", polar_axispar);
    DefineParameter("azimut_axispar", azimut_axispar);
    DefineParameter("logaxis", logaxis);
    DefineParameter("title", title);
}

bool mcEnergyHisto3D::Initialize()
{
	const std::string routineName("mcEnergyHisto3D::Initialize");

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
    if (azimut_axispar.size() != 3)
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
    if ((float)((int)(azimut_axispar[0])) != azimut_axispar[0])
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
    if (azimut_axispar[0] < 0)
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
    if (azimut_axispar[1] >= azimut_axispar[2])
	{
		COUT(ERROR) << "The lower axis limit is greater or equal to the upper limit." << ENDL;
		return false;
	}
	
	if (logaxis)
		GenerateLogEnergyBinning();
	else
		GenerateEnergyBinning();

    GeneratePolarBinning();
    GenerateAzimutBinning();

	// Create the histogram
	std::string histo_name = "h_" + GetName();
	histo = std::make_shared<TH3D>(
        histo_name.c_str(), 
        title.c_str(), 
        energy_binning.size()-1, &(energy_binning[0]),
        polar_binning.size()-1, &(polar_binning[0]),
        azimut_binning.size()-1, &(azimut_binning[0]));

	histo->GetXaxis()->SetTitle("MC Momentum (GV()");
    histo->GetYaxis()->SetTitle("#theta");
    histo->GetZaxis()->SetTitle("#phi");

	return true;
}

bool mcEnergyHisto3D::Process()
{
	const std::string routineName("mcEnergyHisto3D::Process");
	
	auto mctruth = _evStore->GetObject<Herd::MCTruth>("mcTruth");
	if (!mctruth)
	{
		COUT(DEBUG) << "MCTruth not present for event " << GetEventLoopProxy()->GetCurrentEvent() << ENDL;
		return false;
	}
	
    Momentum &Mom = mcTruth->primaries[0].initialMomentum;
    Point Pos = mcTruth->primaries[0].initialPosition;
    Line MCtrack(Pos, Mom);

    auto theta = MCtrack.Polar();
    auto phi = MCtrack.Azimut();
    auto mcmom = std::sqrt(mctruth->primaries.at(0).initialMomentum * mctruth->primaries.at(0).initialMomentum);

    histo->Fill(mcmom, theta, phi);

	return true;
}

bool mcEnergyHisto3D::Finalize()
{
	const std::string routineName("mcEnergyHisto3D::Finalize");
	auto globStore = GetDataStoreManager()->GetGlobalDataStore("globStore");
	if (!globStore)
	{
		COUT(ERROR) << "Global data store not found." << ENDL;
		return false;
	}

	globStore->AddObject(GetName(), histo);
	
	return true;
}

void mcEnergyHisto3D::GenerateEnergyLogBinning()
{
	energy_binning.resize((int)energy_axispar[0]+1);
	double log_interval = (log10(energy_axispar[2]) - log10(energy_axispar[1])) / (int)energy_axispar[0];
	for (auto bIdx = 0; bIdx <= (int)energy_axispar[0]; ++bIdx)
		energy_binning[bIdx] = pow(10, log10(energy_axispar[1]) + bIdx * log_interval);
}

void mcEnergyHisto3D::GenerateEnergyBinning()
{
	energy_binning.resize((int)energy_axispar[0]+1);
	double interval = (energy_axispar[2] - energy_axispar[1]) / (int)energy_axispar[0];
	for (auto bIdx = 0; bIdx <= (int)energy_axispar[0]; ++bIdx)
		energy_binning[bIdx] = energy_axispar[1] + bIdx * interval;
}

void mcEnergyHisto3D::GeneratePolarBinning()
{
	polar_binning.resize((int)polar_axispar[0]+1);
	double interval = (polar_axispar[2] - polar_axispar[1]) / (int)polar_axispar[0];
	for (auto bIdx = 0; bIdx <= (int)polar_axispar[0]; ++bIdx)
		polar_binning[bIdx] = polar_axispar[1] + bIdx * interval;
}

void mcEnergyHisto3D::GenerateAzimutBinning()
{
	azimut_binning.resize((int)azimut_axispar[0]+1);
	double interval = (azimut_axispar[2] - azimut_axispar[1]) / (int)azimut_axispar[0];
	for (auto bIdx = 0; bIdx <= (int)azimut_axispar[0]; ++bIdx)
		azimut_binning[bIdx] = azimut_axispar[1] + bIdx * interval;
}