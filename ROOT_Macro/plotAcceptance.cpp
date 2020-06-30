#include <iostream>
#include <memory>

#include "TFile.h"
#include "TH1D.h"
#include "TMath.h"

// Simulation parameters
#define genRadius 3		//Generation radius in meters

void buildAcceptance(const char* inFilePath, const char* outFilePath)
{
	TFile myInFile(inFilePath, "READ");
	if (myInFile.IsZombie())
	{
		std::cerr << "\n\nError reading input ROOT file: " << inFilePath << std::endl;
		exit(123);
	}

	TH1D* selected_events = static_cast<TH1D*>(myInFile.Get("hcalofidvolume"));
	TH1D* generated_events = static_cast<TH1D*>(myInFile.Get("h_mcgenspectrum"));
	
	selected_events->Sumw2();
	generated_events->Sumw2();

	selected_events->SetDirectory(0);
	generated_events->SetDirectory(0);

	myInFile.Close();

	auto genSurface = 4*TMath::Pi()*genRadius;
	auto scaleFactor = TMath::Pi()*genSurface;

	selected_events->Divide(generated_events);
	selected_events->Scale(scaleFactor);

	selected_events->SetName("acceptance");
	selected_events->SetTitle("Acceptance");

	TFile myOutFile(outFilePath, "RECREATE");
	if (myInFile.IsZombie())
	{
		std::cerr << "\n\nError writing output ROOT file: " << outFilePath << std::endl;
		exit(123);
	}

	selected_events->Write();

	myOutFile.Close();

}