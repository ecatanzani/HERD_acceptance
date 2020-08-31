// Example headers
#include "mcGenSpectrum.h"

// Root headers
#include "TH1D.h"

// C/C++ standard headers
#include <numeric>
#include <cmath>

RegisterAlgorithm(mcGenSpectrum);

mcGenSpectrum::mcGenSpectrum(const std::string &name) : Algorithm{name},
														axispar{100., 1., 100000.},
														logaxis{true},
														title("title"),
														momrange{-1., -1.},
														index{-1}
{
	DefineParameter("axispar", axispar);
	DefineParameter("logaxis", logaxis);
	DefineParameter("title", title);
	DefineParameter("momrange", momrange);
	DefineParameter("index", index);
}

bool mcGenSpectrum::Initialize()
{
	const std::string routineName("mcGenSpectrum::Initialize");

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
	if (axispar[0] < 0.)
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
	histo = std::make_shared<TH1D>(histo_name.c_str(), title.c_str(), binning.size() -1, &(binning[0]));
	histo->GetXaxis()->SetTitle("MC Momentum (GV()");

	ngen = 0;

	return true;
}

bool mcGenSpectrum::Process()
{
	const std::string routineName("mcGenSpectrum::Process");

	auto mctruth = _evStore->GetObject<Herd::MCTruth>("mcTruth");
	if (!mctruth)
	{
		COUT(DEBUG) << "MCTruth not present for event " << GetEventLoopProxy()->GetCurrentEvent() << ENDL;
		return false;
	}
	ngen += mctruth->nDiscarded + 1;

	return true;
}

bool mcGenSpectrum::Finalize()
{
	const std::string routineName("mcGenSpectrum::Finalize");

	for (int bIdx = 1; bIdx < histo->GetNbinsX(); ++bIdx)
	{
		if (histo->GetBinLowEdge(bIdx + 1) < momrange[0]) 
			histo->SetBinContent(bIdx, 0);
		else if (histo->GetBinLowEdge(bIdx) > momrange[1]) 
			histo->SetBinContent(bIdx, 0);
		else
		{
			double lowedge = std::max(histo->GetBinLowEdge(bIdx), momrange[0]);
			double highedge = std::min(histo->GetBinLowEdge(bIdx + 1), momrange[1]);
			double w = 0;
			if (index == -1)
				w = (log10(highedge) - log10(lowedge)) / (log10(momrange[1]) - log10(momrange[0]));
			else
			{
				w = (1./(1-index))/(pow(highedge, -index + 1) - pow(lowedge, -index +1));
				w /= (1./(1-index))/(pow(momrange[1], -index + 1) - pow(momrange[0], -index +1))
			}
			histo->SetBinContent(bIdx, ngen * w);
		}
	}

	auto globStore = GetDataStoreManager()->GetGlobalDataStore("globStore");
	if (!globStore)
	{
		COUT(ERROR) << "Global data store not found." << ENDL;
		return false;
	}
	globStore->AddObject(histo->GetName(), histo);

	return true;
}

void mcGenSpectrum::GenerateLogBinning()
{
	binning.resize((int)axispar[0]+1);
	double log_interval = (log10(axispar[2]) - log10(axispar[1])) / (int)axispar[0];
	for (auto bIdx = 0; bIdx <= (int)axispar[0]; ++bIdx)
		binning[bIdx] = pow(10, log10(axispar[1]) + bIdx * log_interval);
}
void mcGenSpectrum::GenerateBinning()
{
	binning.resize((int)axispar[0]+1);
	double interval = (axispar[2] - axispar[1]) / (int)axispar[0];
	for (auto bIdx = 0; bIdx <= (int)axispar[0]; ++bIdx)
		binning[bIdx] = axispar[1] + bIdx * interval;
}
