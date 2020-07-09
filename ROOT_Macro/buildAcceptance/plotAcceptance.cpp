#include <iostream>
#include <memory>

#include "TFile.h"
#include "TH1D.h"
#include "TMath.h"

// Simulation parameters
#define genRadius 3		//Generation radius in meters

void buildAcceptance(const char* LE_dataFile, const char* HE_dataFile, const char* outFilePath)
{
	TFile myInFile_LE(LE_dataFile, "READ");
	if (myInFile_LE.IsZombie())
	{
		std::cerr << "\n\nError reading input ROOT file: " << LE_dataFile << std::endl;
		exit(123);
	}

	TH1D* selected_events_LE = static_cast<TH1D*>(myInFile_LE.Get("h_calo_filtered_fidvolume"));
	TH1D* generated_events_LE = static_cast<TH1D*>(myInFile_LE.Get("h_mcgenspectrum"));
	
	selected_events_LE->Sumw2();
	generated_events_LE->Sumw2();

	selected_events_LE->SetDirectory(0);
	generated_events_LE->SetDirectory(0);

	myInFile_LE.Close();

	TFile myInFile_HE(HE_dataFile, "READ");
	if (myInFile_HE.IsZombie())
	{
		std::cerr << "\n\nError reading input ROOT file: " << HE_dataFile << std::endl;
		exit(123);
	}

	TH1D* selected_events_HE = static_cast<TH1D*>(myInFile_HE.Get("h_calo_filtered_fidvolume"));
	TH1D* generated_events_HE = static_cast<TH1D*>(myInFile_HE.Get("h_mcgenspectrum"));
	
	selected_events_HE->Sumw2();
	generated_events_HE->Sumw2();

	selected_events_HE->SetDirectory(0);
	generated_events_HE->SetDirectory(0);

	myInFile_HE.Close();

	// Add histos

	TH1D* h_selected = static_cast<TH1D*>(selected_events_LE->Clone("h_selected"));
	h_selected->Add(selected_events_HE);

	TH1D* h_generated = static_cast<TH1D*>(generated_events_LE->Clone("h_generated"));
	h_generated->Add(generated_events_HE);

	auto genSurface = 4*TMath::Pi()*pow(genRadius,2);
	auto scaleFactor = TMath::Pi()*genSurface;

	h_selected->Divide(h_generated);
	h_selected->Scale(scaleFactor);

	h_selected->SetName("acceptance");
	h_selected->SetTitle("Acceptance");

	TFile myOutFile(outFilePath, "RECREATE");
	if (myOutFile.IsZombie())
	{
		std::cerr << "\n\nError writing output ROOT file: " << outFilePath << std::endl;
		exit(123);
	}

	h_selected->Write();

	myOutFile.Close();

}