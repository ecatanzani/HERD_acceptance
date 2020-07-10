#include "TFile.h"
#include "TF1.h"
#include "TH1D.h"

#include <memory>
#include <math.h>

#define deadtime 3e-3

void simuLiveTime(const char* fitFluxPath, const char* acceptancePath)
{
    TFile* inFile = TFile::Open(fitFluxPath, "READ");
    if (inFile->IsZombie())
    {
        std::cerr << "\n\nError reading input file: " << fitFluxPath << std::endl;
        exit(100);
    }

    std::unique_ptr<TF1> fitFunc (static_cast<TF1*>(inFile->Get("fitFunc")));

    inFile = TFile::Open(acceptancePath, "READ");
    if (inFile->IsZombie())
    {
        std::cerr << "\n\nError reading input file: " << acceptancePath << std::endl;
        exit(100);
    }

    std::unique_ptr<TH1D> h_acceptance (static_cast<TH1D*>(inFile->Get("acceptance")));
    h_acceptance->Sumw2();

    // Get acquisition rate bin per bin
    std::unique_ptr<TH1D> h_acqRate (static_cast<TH1D*>(h_acceptance->Clone("h_acqRate")));
    auto multiplyStatus = h_acqRate->Multiply(fitFunc.get());
    if (!multiplyStatus)
        std::cout << "\nERROR during TF1 multiplication" << std::endl;
    
    // Integrate acquisition rate in energy
    auto totalRate = h_acqRate->Integral("width");
    std::cout << "\n\nTotal rate (Hz): " << totalRate << std::endl;
    std::cout << "Total livetime (s): " << exp(-deadtime * totalRate) << std::endl;

    TFile outFile("lvtSimu.root", "RECREATE");
    if (outFile.IsZombie())
    {
        std::cerr << "\n\nError writing output file " << std::endl;
        exit(100);
    }

    h_acqRate->Write();
    fitFunc->Write();
    h_acceptance->Write();

    outFile.Close();
}