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

void mgBuilder(const char* inFluxPath, const char* inFluxPathE3)
{
    TFile* inFile = TFile::Open(inFluxPathE3, "READ");
    if (inFile->IsZombie())
    {
        std::cerr << "\n\nError reading input file: " << inFluxPathE3 << std::endl;
        exit(100);
    }

    std::unique_ptr<TMultiGraph> mgE3 = std::make_unique<TMultiGraph>();
    mgE3->SetTitle("All electron flux (AMS, DAMPE, CALET); Energy(GeV); Flux (#E^3)");
    mgE3->SetName("flux_mg_E3");

    for (auto&& keyAsObj : *inFile->GetListOfKeys())
    {
        auto key = static_cast<TKey*>(keyAsObj);
        if (!strcmp(key->GetClassName(), "TGraphAsymmErrors"))
            mgE3->Add(static_cast<TGraphAsymmErrors *>(key->ReadObj()));
    }

    inFile = TFile::Open(inFluxPath, "READ");
    if (inFile->IsZombie())
    {
        std::cerr << "\n\nError reading input file: " << inFluxPath << std::endl;
        exit(100);
    }

    std::unique_ptr<TMultiGraph> mg = std::make_unique<TMultiGraph>();
    mg->SetTitle("All electron flux (AMS, DAMPE, CALET); Energy(GeV); Flux");
    mg->SetName("flux_mg");

    for (auto&& keyAsObj : *inFile->GetListOfKeys())
    {
        auto key = static_cast<TKey*>(keyAsObj);
        if (!strcmp(key->GetClassName(), "TGraphAsymmErrors"))
            mg->Add(static_cast<TGraphAsymmErrors *>(key->ReadObj()));
    }

    TFile outFile("flux_mg.root", "RECREATE");
    if (outFile.IsZombie())
    {
        std::cerr << "\n\nError writing output file" << std::endl;
        exit(100);
    }
    
    mg->Write();
    mgE3->Write();

    outFile.Close();

}

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
    mg->SetTitle("All electron flux (AMS, DAMPE, CALET); Energy(GeV); Flux");
    mg->SetName("flux_mg");
    
    for (auto&& keyAsObj : *inFile->GetListOfKeys())
    {
        auto key = static_cast<TKey*>(keyAsObj);
        if (!strcmp(key->GetClassName(), "TGraphAsymmErrors"))
            mg->Add(static_cast<TGraphAsymmErrors *>(key->ReadObj()));
    }

    TF1 fitFuncE3("fitFuncE3", "pow(log10(x),3) * pow(10, [0] + [1]*log10(x) + [2]*pow(log10(x),2) +[3]*pow(log10(x),3) + [4]*pow(log10(x),4)) ", 5, 1e+4);
    fitFuncE3.SetNpx(10000);
    mg->Fit("fitFuncE3", "IR");

    TF1 fitFunc("fitFunc", "(pow(log10(x),3) * pow(10, [0] + [1]*log10(x) + [2]*pow(log10(x),2) +[3]*pow(log10(x),3) + [4]*pow(log10(x),4))) / pow(x,3)", 5, 1e+4);
    for (int nPar=0; nPar<fitFuncE3.GetNpar(); ++nPar)
        fitFunc.SetParameter(nPar, fitFuncE3.GetParameter(nPar));

    TFile outFile("flux_fit.root", "RECREATE");
    if (outFile.IsZombie())
    {
        std::cerr << "\n\nError writing output file" << std::endl;
        exit(100);
    }
    
    mg->Write();
    fitFuncE3.Write();
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