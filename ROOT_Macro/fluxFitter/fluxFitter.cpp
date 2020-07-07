#include "TFile.h"
#include "TF1.h"
#include "TKey.h"
#include "TMultiGraph.h"
#include "TGraphAsymmErrors.h"
#include "TMath.h"
#include "TApplication.h"
#include "TCanvas.h"

#include <memory>
#include <string>

void fluxFitter(const char* inFluxPath)
{   
    //TApplication app("app", 0, 0);
    
    TFile* inFile = TFile::Open(inFluxPath, "READ");
    if (inFile->IsZombie())
    {
        std::cerr << "\n\nError reading input file: " << inFluxPath << std::endl;
        exit(100);
    }

    std::unique_ptr<TMultiGraph> mg = std::make_unique<TMultiGraph>();
    mg->SetTitle("AMS electron flux; Energy(GeV); Flux");
    mg->SetName("flux_mg");
    
    for (auto&& keyAsObj : *inFile->GetListOfKeys())
    {
        auto key = static_cast<TKey*>(keyAsObj);
        if (!strcmp(key->GetClassName(), "TGraphAsymmErrors"))
            mg->Add(static_cast<TGraphAsymmErrors *>(key->ReadObj()));
    }

    TF1 fitFunc("fitFunc", "[0] + [1]*log10(x+1) + [2]*pow(log10(x+1),2) +[3]*pow(log10(x+1),3) +[4]*pow(log10(x+1),4) ", 1e-2, 1e+4);
    fitFunc.SetNpx(10000);
    mg->Fit("fitFunc", "IR");

    TFile outFile("flux_fit.root", "RECREATE");
    if (outFile.IsZombie())
    {
        std::cerr << "\n\nError writing output file" << std::endl;
        exit(100);
    }
    
    mg->Write();
    fitFunc.Write();

    outFile.Close();

    /*
    TCanvas c1("c1", "AMS02 Flux");
    c1.cd();
    mg->Draw();

    TCanvas c2("c2", "AMS02 Flux - Fitting function");
    c2.cd();
    fitFunc.Draw();

    app.Run();
    */
}